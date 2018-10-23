package invizio.imgalgo.label;


import net.imagej.ImageJ;
import net.imglib2.Cursor;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.algorithm.neighborhood.DiamondTipsShape;
import net.imglib2.algorithm.neighborhood.Neighborhood;
import net.imglib2.algorithm.neighborhood.RectangleShape;
import net.imglib2.algorithm.stats.ComputeMinMax;
import net.imglib2.img.Img;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.IntegerType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.integer.IntType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;
import ij.IJ;
import ij.ImagePlus;
import invizio.imgalgo.HierarchicalFIFO;
import invizio.imgalgo.label.hwatershed.HWatershedHierarchy;
import invizio.imgalgo.label.hwatershed.HWatershedHierarchy.ConnectivityHWS;
import invizio.imgalgo.util.Pixel;
import invizio.imgalgo.util.RAI;


/**
 * 
 * @author Benoit Lombardot
 * 
 */

	
	
public class SeededWatershed < T extends RealType<T> & NativeType<T>, U extends RealType<U>> extends DefaultLabelAlgorithm< T > {

	// IFT seeded watershed
	// with quantized priority queue	
	
	
	public enum WatershedConnectivity
	{
		FACE(Pixel.Connectivity.FACE),
		FULL(Pixel.Connectivity.FULL);
		
		Pixel.Connectivity conn;
		
		WatershedConnectivity(Pixel.Connectivity conn)
		{			
			this.conn = conn;
		}
		
		Pixel.Connectivity getConn()
		{
			return conn;
		}
	}
	
	
	
	RandomAccessibleInterval<U> seed;
	float threshold = 0;
	WatershedConnectivity connectivity = WatershedConnectivity.FACE;
	
	
	public SeededWatershed( RandomAccessibleInterval<T> input, RandomAccessibleInterval<U> seed, float threshold, WatershedConnectivity connectivity) {
		
		super( input );
		
		this.getScaleFactor();
		labelMap = RAI.convertToInteger( input , scaleFactor, offset );
		
		this.seed = seed;
		this.threshold = (threshold+offset)*scaleFactor;
		this.connectivity = connectivity;
		
	}
	
	
	
	
	
	// 2016-06-02 (thresh=100, seed=Hmaximax  with h=5)
	// algorithm				| t1-head	| blobs		|
	// full, current 			| 13.1s		| 42ms		|
	// face, in place			| 1.0s		| 25ms		|
	// full, in place			| 2.5s		| 33ms		|
	// face, in place, tresh	| 0.8s		| 19ms		|
	// ==> can't reproduce these results for t1-head in place 
	
	// 2017-06-28: priority queue modified to handle float in any range (no proper sorting of float, all pix in one level are considered equivalent)
	// algorithm				| t1-head		| blobs		|
	// face, in place			| 	4.5s 		| 18ms		| 
	// full, in place			| 	9.0s		| 25ms		| 
	
	@Override
	protected void process()
	{
		
		int ndim = labelMap.numDimensions();
		long[] dimensions = new long[ndim]; labelMap.dimensions(dimensions);
		
		
		float min = Float.MAX_VALUE;
		float max = Float.MIN_VALUE;
		Cursor<IntType> Cursor = Views.iterable( labelMap ).cursor() ;
		while ( Cursor.hasNext() )
		{
			float val = Cursor.next().getRealFloat();
			if(val<min)
				min=val;
			else if( val>max )
				max = val;
		}	
		min = Math.max(min, threshold);
		
		// create a priority queue
		HierarchicalFIFO Q = new HierarchicalFIFO( min, max, (int)Math.max(this.minNumberOfLevel, max-min));
		
		
		
		// fill the queue
		final RandomAccess< IntType > input_RA = labelMap.randomAccess();
		final Cursor< U > seed_cursor = Views.iterable(seed).localizingCursor();
		long idx=-1;
		long[] pos = new long[ndim];
		while( seed_cursor.hasNext() )
		{
			float valSeed = seed_cursor.next().getRealFloat(); 
			input_RA.setPosition( seed_cursor );
			IntType pInput = input_RA.get();
			float pVal = pInput.getRealFloat();

			if ( pVal>=min)
			{	if ( valSeed>0)
				{
					seed_cursor.localize(pos);
					idx = Pixel.getIdxFromPos(pos, dimensions);
					Q.add( idx, (int)pVal );
					pInput.setReal(min-1-valSeed);
				}
			}
			else
			{
				pInput.setReal(min-1);
			}
		}
		
		// extend input and seeds to to deal with out of bound
		IntType outOfBound = labelMap.randomAccess().get().createVariable(); 
		outOfBound.setReal(min-1);
		RandomAccess< IntType > input_XRA = Views.extendValue(labelMap, outOfBound ).randomAccess();
		
		RandomAccess<Neighborhood<IntType>> inputNeigh_RAccess= null;
		if( connectivity == WatershedConnectivity.FACE )
		{
			DiamondTipsShape shape = new DiamondTipsShape( 1 );
			DiamondTipsShape.NeighborhoodsAccessible<IntType> inputNeigh_RA =  shape.neighborhoodsRandomAccessible( Views.extendValue(labelMap, outOfBound ) );
			inputNeigh_RAccess= inputNeigh_RA.randomAccess();
		}
		else{ //( connectivity == WatershedConnectivity.FULL ){
			RectangleShape shape = new RectangleShape( 1, true );
			RectangleShape.NeighborhoodsAccessible<IntType> inputNeigh_RA =  shape.neighborhoodsRandomAccessible( Views.extendValue(labelMap, outOfBound ) );
			inputNeigh_RAccess= inputNeigh_RA.randomAccess();
		}
		long[] nPos = new long[ndim];
		long nIdx=0;
        while( Q.HasNext() )
		{ 	
        	final HierarchicalFIFO.Item pixel  = Q.Next();
        	final long pIdx = pixel.getIndex();
			final float pVal = Q.getCurrent_value();
			
			final long[] posCurrent = new long[ndim];
			Pixel.getPosFromIdx( pIdx, posCurrent, dimensions);
			
			input_XRA.setPosition(posCurrent);
			float pLabel = input_XRA.get().getRealFloat();
			
			inputNeigh_RAccess.setPosition(posCurrent);
			Neighborhood<IntType> neighborhood_U = inputNeigh_RAccess.get();
			Cursor<IntType> neighCursor = neighborhood_U.localizingCursor();
			while( neighCursor.hasNext() ){
				final IntType n = neighCursor.next();
				
				final float nVal = n.getRealFloat();
				if ( nVal>=min ) // is not queued yet and is in bound?
				{
					neighCursor.localize(nPos);
					nIdx = Pixel.getIdxFromPos(nPos, dimensions);
					Q.add( nIdx, Math.min(pVal, nVal) );
					n.setReal(pLabel);
				}
			}
			
		}
        
        final IntType minU = labelMap.randomAccess().get().createVariable();
        minU.setReal(min-1);
        final IntType minusOneU = labelMap.randomAccess().get().createVariable();
        minusOneU.setReal(-1);
        Cursor<IntType> input_cursor2 = Views.iterable( labelMap ).cursor();
        while( input_cursor2.hasNext() )
		{
        	IntType p = input_cursor2.next();
        	if (p.getRealFloat()>=(min-1) )
        	{
        		p.setReal(0);
        	}
        	else
        	{
        		p.sub(minU);
            	p.mul(minusOneU);
        	}
		}
        
		return;
	}
	
	
	
	
//	@Deprecated
//	public static <T extends RealType<T>, U extends RealType<U>> RandomAccessibleInterval<T> watershedInPlace(RandomAccessibleInterval<T> input, RandomAccessibleInterval<U> seed)
//	{
//		float threshold = Float.NEGATIVE_INFINITY; 
//		return watershedInPlace(input, seed, threshold, WatershedConnectivity.FULL);
//	}
//	
//	@Deprecated
//	public static <T extends RealType<T>, U extends RealType<U>> RandomAccessibleInterval<T> watershedInPlace(RandomAccessibleInterval<T> input, RandomAccessibleInterval<U> seed, float threshold)
//	{
//		return watershedInPlace(input, seed, threshold, WatershedConnectivity.FULL);
//	}
//	
//	
//	public static <T extends RealType<T>, U extends RealType<U>, V extends RealType<V>> 
//		RandomAccessibleInterval<V> watershed(RandomAccessibleInterval<T> input, RandomAccessibleInterval<U> seed)
//	{
//		float threshold = Float.NEGATIVE_INFINITY; 
//		boolean inPlace = false;
//		return watershed(input, seed, threshold, inPlace);
//	}
//	
//	public static <T extends RealType<T>, U extends RealType<U>, V extends RealType<V>> 
//		RandomAccessibleInterval<V> watershed(RandomAccessibleInterval<T> input, RandomAccessibleInterval<U> seed, float threshold)
//	{
//		boolean inPlace = false;
//		return watershed(input, seed, threshold, inPlace);
//	}
//	
//	public static <T extends RealType<T>, U extends RealType<U>, V extends RealType<V>> 
//		RandomAccessibleInterval<V> watershed(RandomAccessibleInterval<T> input, RandomAccessibleInterval<U> seed, boolean inPlace)
//	{
//		float threshold = Float.NEGATIVE_INFINITY; 
//		return watershed(input, seed, threshold, inPlace);
//	}
	
	
//	public static <T extends RealType<T>, U extends RealType<U>, V extends RealType<V>> 
//		RandomAccessibleInterval<V> watershed(RandomAccessibleInterval<T> input0, RandomAccessibleInterval<U> seed, float threshold, boolean inPlace)
//	{
//		// if type is RGB --> return the input image with a warning: only gray level images are processed
//		// this should happen before that function
//		//RandomAccessibleInterval<V> input = null;
//		
//		if( ! inPlace ){ // create a container copying the input image but in 32 bit int type
//			
//			//Object type = Util.getTypeFromInterval( input0 );
//			
//			int nDim = input0.numDimensions();
//			long[] dims = new long[nDim];
//			input0.dimensions(dims);
//			long nPixel = 1;
//			for (int d=0; d<nDim ; d++){
//				nPixel *= dims[d];
//			}
//			
//			Img<IntType> input;
//			if( nPixel < Math.pow(2, 31) ){
//				input = new ArrayImgFactory<IntType>().create( input0, new IntType());
//			}
//			else{
//				input = new CellImgFactory<IntType>().create( input0, new IntType());
//			}
//			
//			final Cursor< IntType > out = input.localizingCursor();
//			final RandomAccess< T > in = input0.randomAccess();
//
//			while( out.hasNext() )
//			{
//				out.fwd();
//				in.setPosition( out );
//				out.get().setReal( in.get().getRealFloat() );
//			}
//			return (RandomAccessibleInterval<V>) watershedInPlace(input, seed, threshold, WatershedConnectivity.FACE);
//
//		}
//		else // if inplace
//		{ 
//			RandomAccessibleInterval<T> input =  input0;
//			return (RandomAccessibleInterval<V>) watershedInPlace(input, seed, threshold, WatershedConnectivity.FACE);
//		}
//		
//	}
//	
	
	public static void main(final String... args)
	{
		
		ImageJ ij = new ImageJ();
		ij.ui().showUI();
		
		//Dataset data = (Dataset) ij.io().open("F:\\projects\\blobs32.tif");
		//Img<FloatType> img = (Img<FloatType>) data.getImgPlus();
		
		ImagePlus imp0 = IJ.openImage("C:/Users/Ben/workspace/testImages/sampleNoise_std50_blur10.tif"); //blobs32.tif");
		Img<FloatType> img = ImageJFunctions.wrap(imp0);
		float threshold = 0.5f;
		float hMin = 0.01f;
		
		HMaxima<FloatType> maxima = new HMaxima<FloatType>( img, threshold, hMin );
		RandomAccessibleInterval<IntType> seeds = maxima.getLabelMap();
		
		ij.ui().show(seeds);
		
		SeededWatershed<FloatType,IntType> watershed = new SeededWatershed<FloatType,IntType>( img , seeds, threshold, WatershedConnectivity.FACE );
		
		RandomAccessibleInterval<IntType> labelMap = watershed.getLabelMap();
		
		
		ij.ui().show(labelMap);
		
		
//		new ij.ImageJ();
//		//ImageJ ij = new ImageJ();
//		//ij.ui().showUI();
//		
//		IJ.open("F:\\projects\\t1-head.tif");
//		//IJ.open("/Users/lombardo/workspace/test_images/t1-head.tif");
//		//IJ.open("/Users/lombardo/workspace/test_images/blobs.tif");
//		//IJ.open("/Users/lombardo/workspace/test_images/mitosis.tif");
//		ImagePlus imp = IJ.getImage();
//		IJ.run(imp, "Gaussian Blur...", "sigma=2 stack");
//		
//
//		Img<FloatType> input = ImagePlusAdapter.convertFloat(imp);
//		
//		Img<FloatType> inputSeed = ImagePlusAdapter.convertFloat(imp);
//		HMaximaLabeling maxDetector = new HMaximaLabeling();
//		Img<IntType> seed = maxDetector.HMaxima(inputSeed,5);
//		
//		
//		int nIter=10;	
//		String lut = "3-3-2 RGB";
//		RandomAccessibleInterval<FloatType> ws4 = null;
//		boolean inPlace = false;
//		
//		DebugHelper.trackDeltaTime(null);
//		long startTime = System.nanoTime();
//		for (int i=0; i<nIter; i++)
//		{
//			//ws4 = watershedInPlace( input.copy(), seed, thresh/255f, WatershedConnectivity.FACE);
//			ws4 = watershed( input, seed, inPlace);
//			DebugHelper.trackDeltaTime(null);
//		}
//		long stopTime = System.nanoTime();
//
//		ImagePlus imp4 = ImageCreationUtilities.convertImgToImagePlus(ws4, "c4 from normalized float", lut, imp.getDimensions(), imp.getCalibration());
//		imp4.show() ;
//		ImagePlus impInput = ImageCreationUtilities.convertImgToImagePlus(input, "input after watershed with inplace="+inPlace, lut, imp.getDimensions(), imp.getCalibration());
//		impInput.show() ;
//		
//		long mean_dt = (stopTime-startTime)/(nIter*1000000);
//		System.out.println("Average time per iteration: " + mean_dt );
	}
	

}
