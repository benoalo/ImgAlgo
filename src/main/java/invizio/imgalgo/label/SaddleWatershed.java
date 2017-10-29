package invizio.imgalgo.label;



import java.io.IOException;

import ij.IJ;
import net.imagej.ImageJ;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.algorithm.morphology.distance.DistanceTransform;
import net.imglib2.algorithm.morphology.distance.DistanceTransform.DISTANCE_TYPE;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgs;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.integer.IntType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.util.Intervals;


/**
 * 
 * @author Benoit Lombardot
 * 
 */



public class SaddleWatershed < T extends RealType<T> & NativeType<T> > extends DefaultLabelAlgorithm<T> {

	float threshold;
	public RandomAccessibleInterval<FloatType> distMap;
	
	
	public SaddleWatershed(RandomAccessibleInterval<T> input , float threshold ) {
		
		super(input);
		this.threshold = threshold;
		
	}
	
	
	@Override
	protected void process()
	{
		
		distMap = ArrayImgs.floats( Intervals.dimensionsAsLongArray( input ) );
		final double lambda = 40.0;
		DistanceTransform.transform( input, distMap, DISTANCE_TYPE.L1, lambda );
		
	}
	
	
	
	public static void main(final String... args) throws IOException
	{
		ImageJ ij = new ImageJ();
		ij.ui().showUI();
		
		//Dataset dataset = (Dataset) ij.io().open("C:/Users/Ben/workspace/testImages/blobs32.tif"); //rleLabelTest/toLabel2.tif");
		//@SuppressWarnings("unchecked")
		//Img<FloatType> img = (Img<FloatType>) dataset.getImgPlus();
		//Img<FloatType> img2 = ij.op().convert().float32( img );
		Img<FloatType> img =ImageJFunctions.wrap( IJ.openImage("C:/Users/Ben/workspace/testImages/blobs32.tif") );
		//Img<FloatType> img =ImageJFunctions.wrap( IJ.openImage("F:/projects/noise_std50_blur10_thresh40.tif") );
		//img = (Img<FloatType>) ij.op().math().multiply( img, new FloatType(-1) );
		
		
		float threshold = 100;
		int nIter=1;
		int warmup=0;
		
		long dt1 = 0;
		long min_dt1 = Long.MAX_VALUE;
		SaddleWatershed<FloatType> labeler1 = null;
		RandomAccessibleInterval<FloatType> distMap = null;
		for(int i=0 ; i<nIter ; i++) 
		{
			labeler1 = new SaddleWatershed<FloatType>( img , threshold);
			//labelMap1 = labeler1.getLabelMap();
			labeler1.run();
			distMap = labeler1.distMap; 
			long dt = labeler1.getProcessTime();
			if( i>=warmup )
				dt1 += dt;
			min_dt1 = min_dt1 > dt ? dt : min_dt1;
			//System.out.println("iter "+ i + ": " + dt );
		}
		System.out.println("dt1 " + (dt1/(nIter-warmup)));
		System.out.println("min1 " + (min_dt1));
		
		ij.ui().show(img);
		ij.ui().show(distMap);
	}
}
