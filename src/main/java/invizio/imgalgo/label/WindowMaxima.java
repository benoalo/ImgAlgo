package invizio.imgalgo.label;

import java.util.ArrayList;
import java.util.List;

import ij.IJ;
import ij.ImagePlus;
import invizio.imgalgo.util.Label;
import net.imagej.ImageJ;
import net.imglib2.Cursor;
import net.imglib2.Interval;
import net.imglib2.Point;
import net.imglib2.RandomAccessible;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.algorithm.localextrema.LocalExtrema;
import net.imglib2.algorithm.localextrema.LocalExtrema.LocalNeighborhoodCheck;
import net.imglib2.algorithm.neighborhood.CenteredRectangleShape;
import net.imglib2.algorithm.neighborhood.HyperSphereShape;
import net.imglib2.algorithm.neighborhood.Neighborhood;
import net.imglib2.algorithm.neighborhood.RectangleShape;
import net.imglib2.algorithm.neighborhood.Shape;
import net.imglib2.img.Img;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.integer.IntType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.util.Intervals;
import net.imglib2.util.Util;
import net.imglib2.view.Views;


/**
 * 
 * @author Benoit Lombardot
 * 
 */

// Ideas:
//  [-] use histogram technic to speed up the process with rectangular neighborhood (similar to IJ1 rank filter)
//  [-] 


public class WindowMaxima < T extends RealType<T> & NativeType<T> > extends DefaultLabelAlgorithm<T>  {

	// parameter
	float threshold=0;
	int neighborhoodRadius = 1;
	int[] span;
	ExtremaType extremaType = ExtremaType.MAXIMA;
	NeighborhoodType neighborhoodType= NeighborhoodType.SQUARE;
	
	// data structure
	List<Point> extrema;
	
	
	public enum ExtremaType{
		MINIMA,
		MAXIMA;
	}
	
	public enum NeighborhoodType{
		RECTANGLE,
		SQUARE,
		SPHERE;
	}
	
	public WindowMaxima(RandomAccessibleInterval<T> input, float threshold, int neighborhoodRadius, ExtremaType extremaType, NeighborhoodType neighborhoodType) {
		super(input);
		this.threshold = threshold;
		this.neighborhoodRadius = neighborhoodRadius;
		this.extremaType = extremaType;
		this.neighborhoodType = neighborhoodType;
	}


	public WindowMaxima(RandomAccessibleInterval<T> input, float threshold) {
		super(input);
		this.threshold = threshold;
		this.neighborhoodRadius = 1;
		this.extremaType = ExtremaType.MAXIMA;
		this.neighborhoodType = NeighborhoodType.SQUARE;
	}

	
	public WindowMaxima(RandomAccessibleInterval<T> input, float threshold, int[] span,  ExtremaType extremaType) {
		super(input);
		this.threshold = threshold;
		this.span = span;
		this.neighborhoodRadius = 1;
		this.extremaType = extremaType;
		this.neighborhoodType = NeighborhoodType.RECTANGLE;
	}
	
	@Override
	protected void process() {
		
		extrema = new ArrayList<Point>();
		
		RandomAccessible< T > imgX = Views.extendBorder(input);
		Interval interval = Intervals.expand(input, 0);
		
		
		LocalNeighborhoodCheck< Point, T > localNeighborhoodCheck;
		T val = Util.getTypeFromInterval(input).createVariable();
		val.setReal(threshold);
		switch(extremaType)
		{
			case MINIMA:	
				localNeighborhoodCheck = new LocalExtrema.MinimumCheck< T >( val );
				break;
			default: // case MAXIMA:
				localNeighborhoodCheck  = new LocalExtrema.MaximumCheck< T >( val );
				break;
		}
		
		int nDim = input.numDimensions();
		boolean skipCenter = true; //implementation of the extrema checked is such that we don't care of equal values
		Shape shape;
		switch(neighborhoodType)
		{
		case RECTANGLE:
			if( span == null) {
				span = new int[nDim];
				for( int d=0; d<nDim; d++) {
					span[d] = Math.max(1,(int)neighborhoodRadius );
				}
			}
			shape = new CenteredRectangleShape( span, skipCenter );
			break;
			
		case SQUARE:	
			shape = new RectangleShape( neighborhoodRadius, skipCenter );
			break;

		default: // case SPHERE:
			shape = new HyperSphereShape( neighborhoodRadius);
			break;
		
		}
		
		final Cursor< T > center = Views.flatIterable( input ).cursor();
		
		
		
		for ( final Neighborhood< T > neighborhood : shape.neighborhoods( Views.interval( imgX, interval)) ) 
		{
			center.fwd();
			Point p = localNeighborhoodCheck.check( center, neighborhood );
			if ( p != null )
				extrema.add( p );
		}
		
		return ;
	}
	
	
	@Override
	public RandomAccessibleInterval<IntType> getLabelMap()
	{
		if( labelMap == null ) {
			long[] dims = new long[input.numDimensions()];
			input.dimensions(dims);
			
			labelMap =  Label.toLabelMap(dims, this.getExtrema() );
		}
		
		return labelMap;		
	}
	
	public List<Point> getExtrema() {
		
		if (extrema == null)
			run();
		
		return extrema;
		
	}
	
	
	
	public static void main(final String... args)
	{
		
		ImageJ ij = new ImageJ();
		ij.ui().showUI();
		
		//ImagePlus imp = IJ.openImage("F:\\projects\\blobs32.tif");
		ImagePlus imp = IJ.openImage("C:/Users/Ben/workspace/testImages/blobs32.tif");
		ij.ui().show(imp);
		
		
		Img<FloatType> img = ImageJFunctions.wrap(imp);
		float threshold = 100;
		int[] span = new int[] { 1 , 1 };
		WindowMaxima<FloatType> labeler = new WindowMaxima<FloatType>( img, threshold, span, ExtremaType.MAXIMA);
		
		RandomAccessibleInterval<IntType> labelMap = labeler.getLabelMap();
		
		String str = labelMap==null ? "null" : labelMap.toString();
		
		System.out.println("hello D-maxima:" + str );
		ij.ui().show(labelMap);
		
		System.out.println("done!");
	}


	
	
	
}
