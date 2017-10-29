package invizio.imgalgo.label.msmaxima;


import java.util.ArrayList;
import java.util.List;

//import mpicbg_scicomp.imgTools.core.util.ImgLib2Utils;
import net.imglib2.Cursor;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
//import net.imglib2.img.Img;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.NumericType;
import net.imglib2.view.Views;


/**
 * 
 * @author Benoit Lombardot
 * 
 */


public class DOGPyramidHandler < T extends  NumericType< T > & NativeType< T > > {

	// parameters
	protected final RandomAccessibleInterval<T> img0;				// input image
	private double[] resolution0;
	private int nScalePerOctave;				// number of scale per Octave
	private int nSuppScalePerOctave;			// number of scale stored without downsampling (beyond the number of scale per octave ) 
		
		
	// deducted from parameters
	private int nDim;
	private int nOctave;						// number of octave in the pyramid (deducted from image size )
	private int maxLevel; 						// total max number of level in the pyramid (deducted from the number of octave)  
	private List< List<RandomAccessibleInterval<T>> > imgPyramid;	// image at each level of an octave, each octave of the pyramid
	
	GaussianPyramidHandler<T> gaussianPyramid;
	private boolean initialized=false;
	private int nThread;
		
		
	/**
	 * 
	 * @param img0 
	 * rk: resolution is set to default value 1 for every dimension of img0
	 * rk: nScalePerOctave is set to 3   
	 * rk: nSuppScalePerOctave is set to 0 
	 */	
	public DOGPyramidHandler(RandomAccessibleInterval<T> img0)
	{
		this.img0 = img0;
		this.nScalePerOctave = 3;
		this.nSuppScalePerOctave = 0;

		this.nDim = img0.numDimensions();
		this.resolution0 = new double[nDim];
		for(int d=0; d<nDim; d++){ this.resolution0[d]=1; }
		
		this.nThread =  Runtime.getRuntime().availableProcessors();	 
	}

	
	
	/**
	 * 
	 * @param img0 
	 * @param resolution : resolution of the voxel in the image
	 * rk: nScalePerOctave is set to 2   
	 * rk: nSuppScalePerOctave is set to 0 
	 */
	public DOGPyramidHandler(RandomAccessibleInterval<T> img0, double[] resolution )
	{
		this.img0 = img0;
		this.nScalePerOctave = 2;
		this.nSuppScalePerOctave = 0;

		this.nDim = img0.numDimensions();
		this.resolution0 = resolution;
		
		this.nThread =  Runtime.getRuntime().availableProcessors();	
	}
	
	
	
	/**
	 * 
	 * @param img0 
	 * @param resolution : resolution of the voxel in the image
	 * @param nScalePerOctave :   
	 * @param nSuppScalePerOctave : 
	 */
	public DOGPyramidHandler(RandomAccessibleInterval<T> img0, double[] resolution, int nScalePerOctave, int nSuppScalePerOctave )
	{
		this.img0 = img0;
		this.nScalePerOctave = nScalePerOctave;
		this.nSuppScalePerOctave = nSuppScalePerOctave;

		this.nDim = img0.numDimensions();
		this.resolution0 = resolution;
		
		this.nThread =  Runtime.getRuntime().availableProcessors();	
	}
	
	
	
	/**
	 * 
	 * @param i level in the pyramid
	 * @return The size of voxel at the ith level of the pyramid
	 * rk: images in the same octave will have the same voxel size 
	 */	
	public double[] getResolution(int level)
	{
		initializePyramid();
		if (level<0 | level>=maxLevel ){ return null; }
		return gaussianPyramid.getResolution( level );
	}
	
	
	/**
	 * 
	 * @param i level in the pyramid
	 * @return The downsampling factors of the level i image wrt. the input image in each image dimensions 
	 * rk: images in the same octave will have the same downsampling factors
	 */
	public double[] getDownsampling(int level)
	{	
		initializePyramid();
		if (level<0 | level>=maxLevel ){ return null; }
		return gaussianPyramid.getDownsampling(level);
	}
	
	
	/**
	 * 
	 * @param i level in the pyramid
	 * @return The scale (detail size) of the level i image (in pixel) 
	 */
	public double getPixelScale(int level)
	{
		initializePyramid();
		if (level<0 | level>=maxLevel ){ return 0; }
		return gaussianPyramid.getPixelScale( level );
	}
	
	
	/**
	 * 
	 * @param i level in the pyramid
	 * @return The scale (detail size) of the level i image (in image resolution unit) 
	 */
	public double getPhysicalScale(int level)
	{
		initializePyramid();
		if (level<0 | level>=maxLevel ){ return 0; }
		return gaussianPyramid.getPhysicalScale( level );
	}
	
	/**
	 * 
	 * @return the maximum level of the pyramid (does not take into account the redundant level present at multiple scale)
	 */
	public int getMaxLevel()
	{
		initializePyramid();
		return maxLevel;
	}
	
	
	/**
	 * 
	 * @param i level of the pyramid
	 * @return return the pyramid image at level i
	 * rk: if nSuppScalePerOctave>0 the function returns the image with lowest resolution (i.e. largest octave)
	 */
	public RandomAccessibleInterval<T> getLevel(int i)
	{
		// todo check that the request makes sense
		initializePyramid();
		if (i<0 | i>=maxLevel ){ return null; }
		int oct = i/nScalePerOctave;
		int octLevel = i%nScalePerOctave;
		generateLevel(oct,octLevel);
		
		return imgPyramid.get(oct).get(octLevel);
	}
	
	
	/**
	 * 
	 * @param oct : octave of the pyramid
	 * @param octLevel : level in the octave
	 * @return return the pyramid image in octave oct and octave level octLevel
	 * rk: this is the only way to generate and retrieve levels when octLevel>=nScalePerOctave
	 */
	public RandomAccessibleInterval<T> getLevel(int oct, int octLevel)
	{
		initializePyramid();
		if ( oct<0 | oct>=nOctave | octLevel<0 | octLevel>=(nScalePerOctave + nSuppScalePerOctave) ){ return null; }
		generateLevel(oct,octLevel);
		return imgPyramid.get(oct).get(octLevel);
	}
	
	/**
	 * adjust the number of thread used by the gaussian blurring (Imglib2)
	 * @param nThread
	 */
	public void setNThread(int nThread)
	{
		this.nThread = Math.min(  Math.max(1, nThread) ,  Runtime.getRuntime().availableProcessors()  );
		gaussianPyramid.setNThread(this.nThread);
	}
	
	
	
	///////////////////////////////////////////////////////////////////////////////////
	// private functions //////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////
	
	private void initializePyramid(){
		
		if ( initialized ){ return; }
		gaussianPyramid = new GaussianPyramidHandler<T>(this.img0, this.resolution0, this.nScalePerOctave, this.nSuppScalePerOctave + 1 );
		maxLevel = gaussianPyramid.getMaxLevel() - 1;
		nOctave = gaussianPyramid.getMaxLevel()/nScalePerOctave;
		
		
		//initialize imgPyramid
		imgPyramid = new ArrayList< List<RandomAccessibleInterval<T>> >(nOctave);
		for(int oct=0; oct<nOctave; oct++)
		{	
			imgPyramid.add( new ArrayList<RandomAccessibleInterval<T>>(nScalePerOctave+nSuppScalePerOctave) );
			for(int octLevel=0; octLevel<nScalePerOctave+nSuppScalePerOctave; octLevel++)
			{
				imgPyramid.get(oct).add(octLevel,null);
			}
		}
		
		initialized= true;	
	}
	
	
	
		
	private void generateLevel(int oct, int octLevel){
		
		//IJ.log("oct,octLevel: "+oct+","+octLevel);
		
		if ( imgPyramid.get(oct).get(octLevel)!=null ){ return; }
		
		
		// make sure that all prior levels of the dog pyramid are already calculated
		if (octLevel>0)
		{
			if (  octLevel==nScalePerOctave  &&  oct<(maxLevel/nScalePerOctave-1)  ) 
			{	
				gaussianPyramid.getLevel(oct+1,0);
			}
			generateLevel(oct, octLevel-1);
		}
		else if(oct>0)
		{
			generateLevel(oct-1, nScalePerOctave - 1 );
		}
		
		RandomAccessibleInterval<T> g1 = gaussianPyramid.getLevel(oct,octLevel);
		RandomAccessibleInterval<T> g2 = gaussianPyramid.getLevel (oct, octLevel+1);
		RandomAccessibleInterval<T> dog = subtractImg(g1, g2); 
		
		//ImgLib2Utils.convertImgToImagePlus(g1, "gauss_"+oct+"_"+octLevel, "Grays", null, null ).show();
		//ImgLib2Utils.convertImgToImagePlus(g2, "gauss_"+oct+"_"+(octLevel+1), "Grays", null, null ).show();
		//ImgLib2Utils.convertImgToImagePlus(dog, "dog_"+oct+"_"+octLevel, "Grays", null, null ).show();
		
		// does subtractImg change the actual content of g1 in gaussianPyramid ?
		// does imgPyramid.get(oct).get(octLevel) is just a reference to an image in gaussian Pyramid
		// if yes to both then
		// When calculating a level all prior level should already be calculated (including the nSuppScale level).
		// also this would avoid having 2 full pyramid in memory, only dog pyramid refering to the modified image in the gaussian pyramid would remain
		
		imgPyramid.get(oct).set( octLevel , dog  );
		
		
		
	}
	
	
	/**
	 * subtract image I2 from I1 in I1
	 * @param I1 an image
	 * @param I2 an image
	 * @return the updated image I1
	 */
	private static < T extends NumericType< T > & NativeType< T > > RandomAccessibleInterval<T> subtractImg(RandomAccessibleInterval<T> I1, RandomAccessibleInterval<T> I2)
	{
		// assumes image are the same size
		 
		Cursor<T> c1 = Views.iterable( I1 ).cursor();
		RandomAccess<T> c2 = I2.randomAccess();
		int[] pos = new int[I1.numDimensions()];
		while( c1.hasNext() )
		{	
			T t1 = c1.next();
			c1.localize(pos);
			c2.setPosition(pos);
			t1.mul(-1);
			t1.add( c2.get() );
		}
		
		return I1;
	}
	
	
	
/////////////////////////////////////////////////////////////
// for testing  /////////////////////////////////////////////

	public static void main(final String... args)
	{
//		//start imageJ
//		new ij.ImageJ();
//
//		// load an image
//		IJ.open("F:\\project_data\\test data\\blobs.tif");
//		ImagePlus imp = IJ.getImage();
//		Img<FloatType> input = ImagePlusAdapter.convertFloat(imp);
//
//		// display the value of sigma and resolution for different resolution
//		double[] res0 = new double[] {1,1};
//		int nScalePerOctave = 2;
//		int nSuppScalePerOctave = 1;
//		int maxOctaveLevel = nScalePerOctave + nSuppScalePerOctave;
//
//		DOGPyramidHandler<FloatType> dogPyramid = new DOGPyramidHandler<FloatType>(input, res0, nScalePerOctave, nSuppScalePerOctave);
//
//		System.out.println( "level\t"+"oct\t" + "octLevel\t" + "pixScale\t" + "phyScale\t" + "downSampling\t" + "Resolution\t" );
//
//		int maxLev = dogPyramid.getMaxLevel();
//		System.out.println("max level: "+ maxLev );
//
//		for( int level = 0; level<4; level++) 
//		{	
//			double[] ds = dogPyramid.getDownsampling(level);
//			double[] res = dogPyramid.getResolution(level);
//			double phyScale = dogPyramid.getPhysicalScale(level);
//			double pixScale = dogPyramid.getPhysicalScale(level);
//
//			int oct = level/nScalePerOctave;
//			int octLevel = level%nScalePerOctave; 
//			System.out.println( level+"\t"+oct+"\t"+octLevel+"\t"+
//					pixScale+"\t"+
//					phyScale+"\t"+
//					Arrays.toString(ds)+"\t"+
//					Arrays.toString(res)+"\t"
//					);
//
//			Img<FloatType> imgLevel = dogPyramid.getLevel(level); 
//			ImgLib2Utils.convertImgToImagePlus(imgLevel, imp.getTitle()+"_"+level+"_"+oct+"_"+octLevel, "Grays", imp.getDimensions(), imp.getCalibration() ).show();
//
//			
//			if( octLevel == nScalePerOctave-1)
//			{
//				int equivLevel = level;
//				octLevel++;
//				equivLevel++;
//				while (octLevel<maxOctaveLevel && equivLevel<maxLev)
//				{
//					ds = dogPyramid.getDownsampling(equivLevel);
//					res = dogPyramid.getResolution(equivLevel);
//					phyScale = dogPyramid.getPhysicalScale(equivLevel);
//					pixScale = dogPyramid.getPhysicalScale(equivLevel);
//
//					System.out.println( equivLevel+"\t"+oct+"\t"+octLevel+"\t"+
//							pixScale+"\t"+
//							phyScale+"\t"+
//							Arrays.toString(ds)+"\t"+
//							Arrays.toString(res)+"\t"
//							);
//
//					Img<FloatType> imgSuppLevel = dogPyramid.getLevel(oct, octLevel); 
//					ImgLib2Utils.convertImgToImagePlus(imgSuppLevel, imp.getTitle()+"_"+equivLevel+"_"+oct+"_"+octLevel, "Grays", imp.getDimensions(), imp.getCalibration() ).show();
//
//					octLevel++;
//					equivLevel++;
//				}
//			}
//		}


	}
	
	
}
