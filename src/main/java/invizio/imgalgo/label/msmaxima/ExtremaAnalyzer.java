package invizio.imgalgo.label.msmaxima;



import Jama.Matrix;
import invizio.imgalgo.label.WindowMaxima;
import net.imglib2.type.numeric.RealType;
import net.imglib2.util.Intervals;
import net.imglib2.util.Util;
import net.imglib2.view.Views;

//import net.imglib2.img.Img;
import net.imglib2.Cursor;
import net.imglib2.Interval;
import net.imglib2.Point;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.algorithm.localextrema.LocalExtrema;
import net.imglib2.algorithm.localextrema.LocalExtrema.LocalNeighborhoodCheck;
import net.imglib2.algorithm.neighborhood.Neighborhood;
import net.imglib2.algorithm.neighborhood.RectangleShape;
import net.imglib2.algorithm.neighborhood.RectangleShape.NeighborhoodsAccessible;


/**
 * 
 * @author Benoit Lombardot
 * 
 */




public class ExtremaAnalyzer < T extends RealType<T> > {
	
	RandomAccessibleInterval<T> img1;
	double[] res1;
	WindowMaxima.ExtremaType extremaType;
	
	int nDim;
	int size;
	RandomAccess<T> imgRA1;
	RandomAccess<Neighborhood<T>> imgNeighRA1;
	LocalNeighborhoodCheck< Point, T > localNeighborhoodCheck;
	
	long[] position;
	boolean initialized = false;
	
	double[] neighData1;
	Matrix Mhessian;
	Matrix Mgradient;
	double[] deltaPosition;
	Double eigValRatioCriteria;
	
	/**
	 * 
	 * @param img
	 * @param res
	 * @param extremaType wether minima or maxima (on the isExtrema method depend on it)
	 */
	public ExtremaAnalyzer( RandomAccessibleInterval<T> img1, double[] res1, WindowMaxima.ExtremaType extremaType)
	{
		//this.neighData1 = getNeighborhoodInImage( convertToLongPos(pt.getPosition()), neighs_RA1);
		this.img1 = img1;
		this.res1 = res1;
		this.nDim = img1.numDimensions();
		this.size = nDim;
		// get the neighborhoodRandomAccessible on img;
		this.extremaType = extremaType; 
		this.imgRA1 = img1.randomAccess();
		int R1 = 1;
		this.imgNeighRA1 = getNeighborhoodRandomAcessible( img1, R1 );
		
		T val = Util.getTypeFromInterval(img1).createVariable();
		val.setZero();
		switch(extremaType)
		{
			case MINIMA:	
				localNeighborhoodCheck = new LocalExtrema.MinimumCheck< T >( val );
				break;
			default: // case MAXIMA:
				localNeighborhoodCheck  = new LocalExtrema.MaximumCheck< T >( val );
				break;
		}
	}
	
	
	/**
	 * Set the position to analyze in the image. the doubl[] interface is provided for convenience, value will be cast to long integer 
	 * @param position
	 * 
	 */
	public void setPosition(double[] position)
	{
		if (position.length != nDim)
			return;
		
		long[] longPosition = new long[nDim];
		for(int d=0; d<nDim; d++ )
			longPosition[d] = (long) position[d];
		
		this.setPosition(longPosition);
	}
	
	
	
	/**
	 * Set the position to analyze in the image
	 * @param position: coordinates of a point the same dimensionality as the input image
	 */
	public void setPosition(long[] position)
	{
		if (position.length != nDim)
			return;
		
		this.position = position;
		this.neighData1 = getNeighborhoodInImage(position, imgNeighRA1);
		this.Mgradient = null;
		this.Mhessian = null;
		this.deltaPosition = null;
		this.eigValRatioCriteria = null;
	}
	
	
	
	/**
	 * @return the subpixel position of the extrema in pixel unit
	 */
	public double[] getOptimizedPosition()
	{
		if (position==null)
			return null;
		
		deltaPosition = getDeltaPosition();
		double[] optimizedPosition = new double[nDim];
		for(int d=0; d<nDim; d++)
		{ 
			optimizedPosition[d] = position[d] + deltaPosition[d];
		}
		
		return optimizedPosition;
	}
	
	
	
	/**
	 * Calculate the eigenvalue ratio criteria in 2D or 3D according to point dimensionality 
	 * lambda1/lambda3<r <=in 2D=> tr(H)^2/Det(H) < (r+1)^2/r  from Lowe, Distinctive Image Features from Scale-Invariant Keypoints, 2004
	 * lambda1/lambda3<r <=in 3D=> tr(H)^2/Det(H) < (2r+1)^3/r^2 from Allaire, Full Orientation Invariance and Improved Feature Selectivity of 3D SIFT with Application to Medical Image Analysis, 2008
	 * @return return the eigenvalue ratio criteria
	 */
	public double getEigValRatioCriteria()
	{
		computeEigValRatioCriteria();
		return eigValRatioCriteria;
	}
	
	protected void computeEigValRatioCriteria()
	{
		if ( eigValRatioCriteria != null)
			return;
			
		if(Mhessian == null)
			computeHessian();
		double det = Mhessian.det();
		if(det==0)
			eigValRatioCriteria = Double.NaN;
		else
			eigValRatioCriteria = Math.pow(Mhessian.trace(),nDim) / det;
		
	}
	
	
	/**
	 * 
	 * @param ratio
	 * @return
	 */
	public boolean isEigValRatioSuperiorTo( double ratioThreshold )
	{
		computeEigValRatioCriteria();
		double criteriaThreshold = Math.pow((nDim-1)*ratioThreshold+1,nDim)/Math.pow(ratioThreshold,nDim-1);		
		return eigValRatioCriteria > criteriaThreshold;
	}
	
	
	
	public boolean isExtrema()
	{
		imgRA1.setPosition(position); 
		Point p = localNeighborhoodCheck.check( imgRA1, imgNeighRA1.get() );
		if ( p != null )
			return true;
		return false;
	}
	
	
	/**
	 * @return the translation from the point to its subpixel position in pixel unit
	 * remark: dpos = inv(H(x))*G
	 */
	protected double[] getDeltaPosition()
	{	
		if (deltaPosition!=null)
			return deltaPosition;
		if (Mhessian==null)
			this.computeHessian();
		if (Mgradient==null)
			this.computeGradient();
		
		Matrix Mdpos = new Matrix(size,1);
		if (Mhessian.rank()==size) 	// test if Mhessian is invertible
			Mdpos = Mhessian.inverse().times(Mgradient).times(-1);
		else 
			return null;
					
		double[] dpos = Mdpos.getColumnPackedCopy();
		for(int d=0; d<nDim; d++)
		{ 
			dpos[d] /= res1[d];
		}

		
		return dpos;
	}
	
	
	
	
	// the delta x at one scale should be adapted to the resolution ratios
	protected void computeGradient()
	{
		int center = (int)(Math.pow(3,nDim)/2);
		Mgradient = new Matrix(size,1);
		for(int d = 0; d<nDim; d++ )
		{
			double val = ( neighData1[center+move(1,d)] - neighData1[center+move(-1,d)] ) / (2*res1[d]);
			Mgradient.set(d, 0, val);
		}
		
	}
	
	
	
	// the delta x at one scale should be adapted to the resolution ratios
	protected void computeHessian()
	{
		double[] hessian = new double[size*size];
		int center = (int)(Math.pow(3,nDim)/2);
		for (int d=0; d<nDim; d++)
		{
			hessian[d+size*d] = ( -2*neighData1[center] +neighData1[center+move(1,d)] +neighData1[center+move(-1,d)] )/(res1[d]*res1[d]);
			for (int f=d+1; f<nDim; f++)
			{
				hessian[d+f*size] = ( ( neighData1[center+move( 1,d)+move(1,f)] - neighData1[center+move( 1,d)+move(-1,f)] ) / (2*res1[f]) 
						             -( neighData1[center+move(-1,d)+move(1,f)] - neighData1[center+move(-1,d)+move(-1,f)] ) / (2*res1[f]) ) / (2*res1[d]); 
				hessian[f+d*size] = hessian[d+f*size];
			}
		}
		Mhessian = new Matrix(hessian, size);
	}
	
	
	
	
	// get the neighborhood of radius r of a point pt 
	protected double[] getNeighborhoodInImage(long[] position, RandomAccess<Neighborhood<T>> neighs_RA)
	{
		int R=1;
		int nDim = position.length;
		double[] data = new double[(int)Math.pow( 2*R+1, nDim)];
		
		neighs_RA.setPosition(position);
		Cursor<T> c0 = neighs_RA.get().cursor();
		int count=0;
		while( c0.hasNext())
			data[count++] =  c0.next().getRealFloat();
		
		return data;
	}
	
	
	
	protected int move(int mov, int dim )
	{	
		return mov*(int)Math.pow(3, dim);
	}
	
	
	
	// create a neighborhood random accessible for a neighborhood of radius R and and a
	/**
	 * create a neighborhood random accessible on image img
	 * @param img, image on which to set the neighborhood
	 * @return
	 */
	protected   RandomAccess<Neighborhood<T>>   getNeighborhoodRandomAcessible( RandomAccessibleInterval<T> img)
	{
		int R=1;
		return getNeighborhoodRandomAcessible( img, R );
	}
	
	/**
	 * 
	 * @param img
	 * @param R
	 * @return
	 */
	protected   RandomAccess<Neighborhood<T>>   getNeighborhoodRandomAcessible( RandomAccessibleInterval<T> img, int R )
	{
		boolean skipCenter = false;
		final RectangleShape shape = new RectangleShape( R, skipCenter );
		NeighborhoodsAccessible<T> neighs = shape.neighborhoodsRandomAccessible(Views.extendBorder(img));
		
		Interval interval = Intervals.expand(img, 0);
		RandomAccess<Neighborhood<T>> neighs_RA = neighs.randomAccess(interval);
		
		return neighs_RA;
	}
	
	
	
	
	
	// helper for debug
//	protected static Img<FloatType> getGaussAtPos(int r,int nDim,double[] c0, double ds)
//	{
//		ImagePlus imp;
//		if(nDim==2)
//			imp = IJ.createImage("Untitled", "32-bit black", 2*r+1, 2*r+1, 1);
//		else
//			imp = IJ.createImage("Untitled", "32-bit black", 2*r+1, 2*r+1, 2*r+1);
//		Img<FloatType> img1 = ImagePlusAdapter.convertFloat(imp);
//		// go through the pixel and draw a gaussian centered on c0
//		Cursor<FloatType> cursor = img1.cursor();
//		long[] pos = new long[nDim];
//		while( cursor.hasNext()){
//			FloatType val = cursor.next();
//			cursor.localize(pos);
//			float D2 = 0;
//			for(int d=0; d<nDim; d++){D2 += (pos[d]-c0[d])*(pos[d]-c0[d]); }
//			D2 += ds*ds;
//			val.setReal( Math.exp( -D2/4 ) );
//		}
//		return img1;
//	}
	
	public static void main(final String... args)
	{
//		new ij.ImageJ();
//		
//		int nDim=3;
//		double[] res = new double[] {1,1,1};
//
//		int r=5;
//		double[] cref = new double[] {r,r,r};
//		double[] c0 = new double[] {r+0.3,r+0.4,r-0.25};
//		double ds=0.2;
//		
//		Img<FloatType> img1 = getGaussAtPos( r,nDim, c0, ds);
//		
//		ImgLib2Utils.convertImgToImagePlus( img1 , "test", "Grays", new int[] {2*r+1,2*r+1,2*r+1,1,1}, null).show();
//
//		// instantiate pointN_NeighCalculation at cref
//		ExtremaAnalyzer<FloatType> extremaAnalyzer = new ExtremaAnalyzer<FloatType>( img1, res, WindowedMaxima.ExtremaType.MAXIMA);
//		extremaAnalyzer.setPosition( cref );
//		if( extremaAnalyzer.isExtrema() )
//		{
//			double[] cOptimized = extremaAnalyzer.getOptimizedPosition();
//			double eigValRatio = extremaAnalyzer.getEigValRatioCriteria();
//			System.out.println(Arrays.toString(cOptimized));
//			System.out.println(eigValRatio);
//		}
//		
		
		
	}
	
	
	
	
}
