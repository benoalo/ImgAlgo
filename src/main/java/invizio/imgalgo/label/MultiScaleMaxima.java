package invizio.imgalgo.label;

//for debug function
//import ij.IJ;
//import ij.ImagePlus;
//import ij.gui.OvalRoi;
//import ij.gui.Overlay;
//import ij.plugin.Duplicator;
//import ij.plugin.frame.RoiManager;
//import java.awt.Color;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import ij.IJ;
import ij.ImagePlus;
import invizio.imgalgo.Attributes;
import invizio.imgalgo.Point_A;
import invizio.imgalgo.RealPoint_A;
import invizio.imgalgo.label.WindowMaxima.ExtremaType;
import invizio.imgalgo.label.msmaxima.DOGPyramidHandler;
import invizio.imgalgo.label.msmaxima.ScaleSpaceExtremaAnalyzer;
import invizio.imgalgo.util.Label;
import net.imagej.ImageJ;
import net.imglib2.Point;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.integer.IntType;
import net.imglib2.type.numeric.real.FloatType;


/**
 * 
 * @author Benoit Lombardot
 * 
 */

/*
 * TODO:
 * 	[x] make the class derive from default label Algorithm
 * 	[x] replace Img by RAI
 *  [x] provide constructor giving access to all the parameter of the preprocessing
 *  [x] provide setter for parameter allowing to select parameters
 *  [-] check if the code is factory agnostic (this prob depends wheter the pyramide construction is factory agnostic)
 *  
 *  [-] check if the min/max scale passed as argument should be physical or pixel based
 *  
 */



/**
* 
* This class allow to detect extrema of intensity in an image scale space.
* by default the maxima of the original image are detected. It is possible
* to detect minima by calling @link{setExtremaType()}. Image can be 2D or 3D and have anistropic voxel 
* 
* nscalePerOctave: integer with minimum value 1 and default is 3 
* ExtremaType: minima or maxima, default maxima 
* @author Benoit lombardot, scientific computing facility, MPI-CBG
* @param <T> type of the input image
*/
public class MultiScaleMaxima< T extends RealType<T> & NativeType< T > > extends DefaultLabelAlgorithm<T>
{
	
	
	////////////////////////////////////////////
	// global variable /////////////////////////
	
	// parameter for the initial detection. these can be set via the constructor
	double[] resolution;							// the resolution of the input image voxel 
	int nScalePerOctave=4;							// number of image in an octave from one image to the next the scale is multiplied by 2^(1/nScalePerOctave)
	WindowMaxima.ExtremaType extremaInOrigImage = WindowMaxima.ExtremaType.MAXIMA;
													// type of extrema in to detect in the image
	WindowMaxima.ExtremaType extremaDOG = WindowMaxima.ExtremaType.MINIMA;
													// type of extrema to detect in the dog (always inverted wrt extremaInOrigImage)
	
	
	// parameter to select the maxima
	double minScale=1;								// minimum scale at which to detect extrema
	double maxScale=Double.MAX_VALUE;				// maximum scale at which to detect extrema
	double threshold=0; 							// absolute threshold in the original image 
	double hMin=0;									// dynamics of the peaks as measured in the dog pyramid (after intensity rescaling)
	double anistropyMax=2; 							// eratio
	
	
	/** @param minScale minScale and maxScale indicate the range of scale at which maxima can be detected*/
	public void setMinScale(double minScale) {
		this.minScale = minScale;
	}

	/** @param maxScale minScale and maxScale indicate the range of scale at which maxima can be detected*/
	public void setMaxScale(double maxScale) {
		this.maxScale = maxScale;
	}

	/** @param threshold on the intensity in the input image */
	public void setThreshold(double threshold) {
		this.threshold = threshold;
	}

	/** @param hMin a threshold on maxima dynamics (i.e. absolute value in scale space. higher value are selected) */
	public void sethMin(double hMin) {
		this.hMin = Math.max( 0 , hMin );
	}
	
	/**
	 * this parameter helps getting rid of detection on ridge and edge
	 * @param anistropyMax a threshold on hessian eigenvalue ratio (lower value are selected)
	 */
	public void setMaxAnisotropy(double anistropyMax) {
		this.anistropyMax = Math.max( 1, anistropyMax );
	}
	
	
	
	
	// data structure and variable the processing
	List<Point_A> maxList;							// list with all the detected extrema
	List<RealPoint_A> realMaxList;					// list with all the detected extrema
	int nDim;
	
	
	
	
	
	//////////////////////////////////////////////
	// constructors //////////////////////////////
	
	public MultiScaleMaxima(RandomAccessibleInterval<T> input)
	{
		super(input);
		this.nDim = input.numDimensions();
		
		double[] resolution = new double[nDim];
		for( int d=0;d<nDim; d++)
			resolution[d]=1;
		
		this.resolution = resolution;
		this.minScale = getMinPhysicalScale();
		this.maxScale = getMaxPhysicalScale();
		this.setExtremaType( WindowMaxima.ExtremaType.MAXIMA );
		
		this.maxList = new ArrayList<Point_A>();
		this.realMaxList = new ArrayList<RealPoint_A>();
		
	}
	
	
	public MultiScaleMaxima(RandomAccessibleInterval<T> input, double[] resolution)
	{
		super(input);
		this.nDim = input.numDimensions();
		this.resolution = resolution;
		this.minScale = getMinPhysicalScale();
		this.maxScale = getMaxPhysicalScale();
		this.setExtremaType( WindowMaxima.ExtremaType.MAXIMA );
		
	}
	
	
	public MultiScaleMaxima(RandomAccessibleInterval<T> input, double[] resolution, double minScale, double maxScale, int nScalePerOctave)
	{
		super(input);
		
		this.resolution = resolution;
		this.minScale = Math.max( minScale , getMinPhysicalScale() );
		this.maxScale = Math.min( maxScale , getMaxPhysicalScale() ) ;
		this.nScalePerOctave = Math.max(1, nScalePerOctave);
		this.setExtremaType( WindowMaxima.ExtremaType.MAXIMA );
		this.nDim = input.numDimensions();
	}
	
	
	
	
	
	////////////////////////////////////////////
	// getter //////////////////////////
	
	// getter
	public double[] getResolution() 	{ return resolution; }
	public double[] getscaleRange() 	{ return new double[] {minScale, maxScale}; }
	
	
	
	
	/**
	 * @return a the list of extrema detected in the scale-space for user set parameters (minScale, maxScale, threshold, hMin, anistropyMax) 
	 */
	public List<Point_A> getExtrema()
	{ 
		if( maxList == null ) {
			run();
		}
		
		return filterMaxima( maxList ); 
	}

	
	public List<RealPoint_A> getRealExtrema()
	{ 
		if( realMaxList == null )
			run();
		return filterMaxima( realMaxList ); 
	}
	

	

	private <L extends Attributes> List<L> filterMaxima(List<L> points){
		
		
		List<L> points_filtered = new ArrayList<L>();
		
		// lambda1/lambda2<r <= in 2D => tr(H)^2/Det(H) < (r+1)^2/r  from Lowe, Distinctive Image Features from Scale-Invariant Keypoints, 2004
		// lambda1/lambda3<r <= in 3D => tr(H)^3/Det(H) < (2r+1)^3/r^2 from Allaire, Full Orientation Invariance and Improved Feature Selectivity of 3D SIFT with Application to Medical Image Analysis, 2008
		// the following formula would work for dim 2 and 3 it is not checked for superior dimension
		// lambda1/lambda_nDim<r <= in 2D or 3D => tr(H)^nDim/Det(H) < ((nDim-1)r+1)^nDim/r^(nDim-1)
		double anisoCrit =  Math.pow((nDim-1)*anistropyMax+1,nDim)/Math.pow(anistropyMax,nDim-1);
		
		for( L pt : points)
		{	
			//System.out.println( pt.toString( ) );
			
			if( Math.abs(pt.getAttribute("value")) 		>= 	hMin 		&&
				pt.getAttribute("intensity") 			>= 	threshold 	&&
				pt.getAttribute("scale")				>= 	minScale	&&
				pt.getAttribute("scale")				<= 	maxScale	&&
				pt.getAttribute("anisoCrit")			<=	anisoCrit		) 
			{	
				points_filtered.add(  pt );
			}
		}
		return points_filtered;
	}

		


	
	
	
	/**
	 * 
	 * @return
	 */
	@Override
	public RandomAccessibleInterval<IntType> getLabelMap()
	{ 
		if( labelMap==null ) {
			long[] dimensions = new long[nDim];
			input.dimensions( dimensions );
			labelMap = Label.toLabelMap( dimensions , this.getExtrema() );
			
			System.out.println("filtered extrema: " + getExtrema().size() );
			System.out.println("unfiltered extrema: " + maxList.size() );
		}
		return labelMap;
	}
	
	
	
	
	
	public double	getMinPhysicalScale() {
		
		double min_minScale= Double.MAX_VALUE;
		for(int d=0; d<nDim; d++) {
			min_minScale = Math.min( min_minScale , resolution[d] );
		}
		return min_minScale;
	}
	
	
	
	public double	getMaxPhysicalScale() {
		
		double max_maxScale= Double.MIN_VALUE;
		for(int d=0; d<nDim; d++) {
			max_maxScale = Math.max( max_maxScale , input.dimension(d)*resolution[d]/2.0 );
		}
		return max_maxScale;
	}
		
	
	//////////////////////////////////////////////////////
	// helper ////////////////////////////////////////////

	
	/**
	 *  @param extremaInOrigImage the type of extrema to detect in the input image
	 */
	private void setExtremaType(WindowMaxima.ExtremaType extremaInOrigImage) {	
		
		this.extremaInOrigImage = extremaInOrigImage; 
		if(extremaInOrigImage.equals(WindowMaxima.ExtremaType.MINIMA) )
		{ 
			this.extremaDOG = WindowMaxima.ExtremaType.MAXIMA;
		}
		else
		{
			this.extremaDOG = WindowMaxima.ExtremaType.MINIMA;
		}
		
	}
	
	
	
	
	/////////////////////////////////////////////////////
	// public methods  //////////////////////////////////
	
	@Override
	protected void process()
	{
		this.maxList = new ArrayList<Point_A>();
		this.realMaxList = new ArrayList<RealPoint_A>();
		
		// calculate the dog pyramid
		int nSuppScalePerOctave = 2; // will save us some image resampling when detecting maxima in the scale space dimension
		DOGPyramidHandler<T> dogPyramid = new DOGPyramidHandler<T>(input, resolution, nScalePerOctave, nSuppScalePerOctave);
		
		RandomAccess<T> inputRA = input.randomAccess();
		
		double IntensityNormFactor = 1 / (  Math.pow( 2.0 , 1.0/(double)nScalePerOctave )  -  1  ); // cf lowe sift, difference of gaussian
		
		// find the start level corresponding to minScale
		int level=0;
		while(  dogPyramid.getPhysicalScale(level)<minScale  &&  level<dogPyramid.getMaxLevel()  )
		{
			//System.out.println("Skipped scale "+ dogPyramid.getPhysicalScale(level) +" (2^"+level/(float)nScalePerOctave+")" );
			level++;
		}

		
		while( dogPyramid.getPhysicalScale( level ) <= maxScale  && level<(dogPyramid.getMaxLevel()-1) )
		{
			double scale = dogPyramid.getPhysicalScale( level );
			//System.out.println("process detection in scale "+scale+" (2^"+ (level/(float)nScalePerOctave) +")" );
			
			int octave = level/nScalePerOctave;
			int octLevel = level%nScalePerOctave;
			if (octLevel<1 & octave>0)
			{
				octave--;
				octLevel = octLevel+nScalePerOctave;
			}
			RandomAccessibleInterval<T> dogImg0 = dogPyramid.getLevel(octave, octLevel-1); // return null if the level does not exist 
			RandomAccessibleInterval<T> dogImg1 = dogPyramid.getLevel(octave, octLevel);
			RandomAccessibleInterval<T> dogImg2 = dogPyramid.getLevel(octave, octLevel+1);
			
			// calculate extrema candidate at that level if necessary
			float extrema_threshold = 0f;
			if( extremaDOG.equals( WindowMaxima.ExtremaType.MAXIMA ) )
				extrema_threshold = (float)hMin;
			else	
				extrema_threshold = (float)(-hMin);
			
			int neighborhood_radius = 1;
			WindowMaxima<T> maxDetector = new WindowMaxima<T>( dogImg1, extrema_threshold, neighborhood_radius, extremaDOG, WindowMaxima.NeighborhoodType.SQUARE );
			List<Point> peaks = maxDetector.getExtrema();
			//List<Point> peaks = WindowMaxima.getExtrema(dogImg1, extrema_threshold, neighborhood_radius, extremaDOG);
				
			// test if a candidate are maxima in the scale space.
			double[] resolution1 = dogPyramid.getResolution(octave*nScalePerOctave);
			ScaleSpaceExtremaAnalyzer<T> scaleExtremaAnalyzer = new ScaleSpaceExtremaAnalyzer<T>( dogImg0, dogImg1, dogImg2, resolution1, scale, nScalePerOctave, extremaDOG );
			
			double[] currentLevelDownSampling = dogPyramid.getDownsampling(octave*nScalePerOctave); // level of downSampling in each dimension WRT the input image
			RandomAccess<T> dogImg1RA = dogImg1.randomAccess(); 
			//int count = 0;
			long[] inputdims = new long[nDim];
			input.dimensions(inputdims);
			for(Point pt : peaks)
			{
				long[] positionDogImg1 = new long[nDim];
				pt.localize(positionDogImg1);
				scaleExtremaAnalyzer.setPosition( positionDogImg1 );
				if ( scaleExtremaAnalyzer.isExtrema() )
				{
					//count++;
					// rescale the point to match with original image coordinate
					double[] position = new double[nDim];
					long[] positionLong = new long[nDim];
					for(int d=0; d<nDim; d++){
						position[d] = positionDogImg1[d]*currentLevelDownSampling[d];	
					}
					
					double anisoCrit = scaleExtremaAnalyzer.getEigValRatioCriteria();
					double optimizedScale = scaleExtremaAnalyzer.getOptimizedScale();
					double[] optimizedPosition = scaleExtremaAnalyzer.getOptimizedPosition();
					for(int d=0; d<nDim; d++){
						double pos = optimizedPosition[d] * currentLevelDownSampling[d];
						if( pos<0 )
							pos=0;
						if( pos>=(double) inputdims[d] )
							pos=(double)inputdims[d]-1;
						optimizedPosition[d] = pos;
						positionLong[d] = (long) pos;
					}
					dogImg1RA.setPosition(positionDogImg1);
					double value = dogImg1RA.get().getRealDouble();
					value = value*IntensityNormFactor;
					
					inputRA.setPosition(positionLong);
					double intensity = inputRA.get().getRealDouble();
					
					
					Map<String, Double> attributes = new HashMap<String, Double>();
					attributes.put("intensity", intensity	);
					attributes.put("scale", 	scale		);
					attributes.put("value", 	value		);
					attributes.put("anisoCrit", anisoCrit	);
					Point_A pt2 = new Point_A(positionLong, attributes);
					maxList.add( pt2 );
					
					
					Map<String, Double> attributes2 = new HashMap<String, Double>();
					attributes2.put("intensity",intensity		);
					attributes2.put("scale", 	optimizedScale	);
					attributes2.put("value", 	value			);
					attributes2.put("anisoCrit",anisoCrit		);
					RealPoint_A pt3 = new RealPoint_A(optimizedPosition, attributes2);
					realMaxList.add( pt3 );
					
				}
			}
			
			// debug /////////////////////////
			//System.out.println("octave,octaveLevel, level,numpeaks  "+octave+" , "+octLevel+" , "+level+" , "+ count );// debug
			//System.out.println("I norm factor: "+ IntensityNormFactor );// debug
			//System.out.println("down sampling: "+ Arrays.toString(currentLevelDownSampling) );// debug
			
			////////////////////////////////////////////////
			level++;
		}		
		//System.out.println("numpeaks: "+ maxList.size() );// debug
		
	}
	
		
	
	
	
	
	
	
	////////////////////////////////////////////
	// utils for debug  ////////////////////////
	
	/*
	private void displayPeaks(Img<T> img, List<Point> peaks, String title)
	{
		long[] dims = new long[img.numDimensions()];
		img.dimensions(dims);
		ImageJimgLib2Connector.floatImageToImagePlus( img , title, "Grays", new int[] {(int)dims[0],(int)dims[1],1,1,1}).show();
		
		ImagePlus imp = IJ.getImage();
		long[] pos = new long[img.numDimensions()];
		float[] px = new float[peaks.size()];
		float[] py = new float[peaks.size()];
		int count = 0;
		for(Point pt : peaks)
		{	pt.localize(pos);
			px[count]=pos[0];
			py[count]=pos[1];
			count++;
		}
		
		PointRoi roi = new PointRoi(px, py);
		imp.setRoi(roi);
	}*/

	/*
	private void displayImg(Img<T> img, String title)
	{
		long[] dims = new long[img.numDimensions()];
		img.dimensions(dims);
		ImageJimgLib2Connector.floatImageToImagePlus( img , title, "Grays", new int[] {(int)dims[0],(int)dims[1],1,1,1}).show();
	}
	*/
	
	/* private void displayScalePeaks(Img<T> img, List<PointNS> peaks, String title)
	{
		long[] dims = new long[img.numDimensions()];
		img.dimensions(dims);
		ImageJimgLib2Connector.floatImageToImagePlus( img , title, "Grays", new int[] {(int)dims[0],(int)dims[1],1,1,1}).show();
		
		ImagePlus imp = IJ.getImage();
		float[] px = new float[peaks.size()];
		float[] py = new float[peaks.size()];
		int count = 0;
		for(PointNS pt : peaks)
		{	;
			px[count]= (float)pt.getPosition(0);
			py[count]= (float)pt.getPosition(1);
			count++;
		}
		
		PointRoi roi = new PointRoi(px, py);
		imp.setRoi(roi);
	}
	*/
	
	
	
//	public static void drawScaleCircles(ImagePlus imp , List<PointN> extrema)
//	{
//		drawScaleCircles( imp , extrema, Color.blue);
//	}
//
//	public static void drawScaleCircles(ImagePlus imp , List<PointN> extrema, Color color)
//	{
//		if (extrema.size()==0){ return; }
//		
//		Overlay ov = imp.getOverlay();
//		if( ov == null )
//			ov = new Overlay();
//		
//		int nDim = extrema.get(0).getnDim();
//		if(nDim == 2) 
//			for( PointN pt : extrema)
//			{
//				double[] pos = pt.getPosition();
//				double r = pt.getFeature("scale");
//				OvalRoi roi = new OvalRoi(pos[0]-r-0.5, pos[1]-r-0.5, 2*r+1, 2*r+1);
//				roi.setStrokeColor(color);
//				ov.add(roi);
//				//IJ.log("(Det,Tr)  ("+pt.getFeature("hDet")+","+pt.getFeature("hTr")+")");
//				//IJ.log("ratioCrit "+pt.getFeature("eigvalRatioCriterion") );
//				IJ.log("scale "+pt.getFeature("scale") );
//				IJ.log("value "+pt.getFeature("value") );
//			}
//		else if (nDim==3)
//		{
//			for( PointN pt : extrema)
//			{
//				double[] pos = pt.getPosition();
//				double r = pt.getFeature("scale");
//				OvalRoi roi = new OvalRoi(pos[0]-r-0.5, pos[1]-r-0.5, 2*r+1, 2*r+1);
//				roi.setPosition( (int)pos[2]+1 );
//				roi.setStrokeColor(color);
//				ov.add(roi);
//			
//			}
//		}
//		
//		imp.setOverlay(ov);
//		imp.show();
//	}
//	
//	public static void AddCirclesToRoiManager(List<PointN> extrema, java.awt.Color color)
//	{
//		if (extrema.size()==0){ return; }
//		
//		RoiManager roiManager = RoiManager.getInstance();
//		if (roiManager==null)
//			roiManager = new RoiManager();
//			
//		int nDim = extrema.get(0).getnDim();
//		if(nDim == 2) 
//			for( PointN pt : extrema)
//			{
//				double[] pos = pt.getPosition();
//				double r = pt.getFeature("scale");
//				OvalRoi roi = new OvalRoi(pos[0]-r-0.5, pos[1]-r-0.5, 2*r+1, 2*r+1);
//				roi.setStrokeColor(color);
//				roiManager.addRoi(roi);
//			}
//		else if (nDim==3)
//		{
//			for( PointN pt : extrema)
//			{
//				double[] pos = pt.getPosition();
//				double r = pt.getFeature("scale");
//				OvalRoi roi = new OvalRoi(pos[0]-r-0.5, pos[1]-r-0.5, 2*r+1, 2*r+1);
//				roi.setPosition( (int)pos[2]+1 );
//				roi.setStrokeColor(color);
//				roiManager.addRoi(roi);
//			}
//		}
//		
//		roiManager.setVisible(true);
//	}
//	
//	public static void drawCircles(ImagePlus imp , List<double[]> extrema, int r, java.awt.Color color)
//	{
//		Overlay ov = new Overlay();
//		int nDim = extrema.get(0).length;
//		if(nDim == 2) 
//			for( double[] pos : extrema)
//			{
//				OvalRoi roi = new OvalRoi(pos[0]-r-0.5, pos[1]-r-0.5, 2*r+1, 2*r+1);
//				roi.setStrokeColor(color);
//				ov.add(roi);
//			}
//		else if (nDim==3)
//		{
//			for( double[] pos : extrema)
//			{
//				OvalRoi roi = new OvalRoi(pos[0]-r-0.5, pos[1]-r-0.5, 2*r+1, 2*r+1);
//				roi.setPosition( (int)pos[2]+1 );
//				roi.setStrokeColor(color);
//				ov.add(roi);
//			
//			}
//		}
//		
//		imp.setOverlay(ov);
//		imp.show();
//	}
//	
//
//	
//	public static void main(final String... args)
//	{
//		new ij.ImageJ();
//		
//		double[] res0 = new double[] {1,1};
//		double minScale = 4;
//		double maxScale = 64;
//		int nScalePerOctave = 3;
//		// maxima selection
//		double noise = 10;
//		//double eRatio = 100;
//		
//		//IJ.run("Blobs (25K)");
//		IJ.open("F:\\project_data\\blobs.tif");
//		
//		ImagePlus imp = IJ.getImage();
//		IJ.run(imp, "Gaussian Blur...", "sigma=1 stack");
//		Img<FloatType> input = ImagePlusAdapter.convertFloat(imp);
//		
//		MultiScaleMaxima<FloatType> MSmax = new MultiScaleMaxima<FloatType>(input, res0, minScale, maxScale);
//		MSmax.setNScalePerOctave(nScalePerOctave);
//		MSmax.setExtremaType(WindowedMaxima.ExtremaType.MAXIMA);
//		List<PointN> maxList = MSmax.getExtrema(noise,false); 
//		//List<PointN> realMaxList = MSmax.getExtrema(noise,true); 
//		
//		IJ.log("max:" + maxList.size());
//		//drawScaleCircles(imp , maxList );
//		AddCirclesToRoiManager(maxList, Color.red);
//		//drawScaleCircles(imp , realMaxList, Color.red );
//		
//		
//	}
	
	
	public static void main(final String... args)
	{
		
		ImageJ ij = new ImageJ();
		ij.ui().showUI();
		
		ImagePlus imp = IJ.openImage("F:\\projects\\blobs32.tif");
		//ImagePlus imp = IJ.openImage("C:/Users/Ben/workspace/testImages/blobs32.tif");
		ij.ui().show(imp);
		
		
		Img<FloatType> img = ImageJFunctions.wrap(imp);
		img = (Img<FloatType>) ij.op().filter().gauss(img, 1.0);

		float threshold = 100;
		float hMin = 100;
		float sMin = 8;
		float sMax = 64;
		double[] pixelSize = new double[] { 1d , 1d };
		MultiScaleMaxima<FloatType> labeler = new MultiScaleMaxima<FloatType>( img , pixelSize );
		
		labeler.setThreshold(threshold);
		labeler.sethMin(hMin);
		labeler.setMinScale(Math.max( sMin, labeler.getMinPhysicalScale() )  );
		labeler.setMaxScale(Math.min( sMax, labeler.getMaxPhysicalScale() )  );
		//labeler.setMaxAnisotropy( 2.0 );
		
		RandomAccessibleInterval<IntType> labelMap = labeler.getLabelMap();
		
		String str = labelMap==null ? "null" : labelMap.toString();
		
		System.out.println("hello Ms-Maxima:" + str );
		System.out.println("min possible scale:" +  labeler.getMinPhysicalScale() );
		System.out.println("max possible scale:" +  labeler.getMaxPhysicalScale() );
		
		ij.ui().show(labelMap);
		
		System.out.println("done!");
	}
	
	
	
}
