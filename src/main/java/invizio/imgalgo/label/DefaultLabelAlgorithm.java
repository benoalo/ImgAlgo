package invizio.imgalgo.label;


import net.imglib2.Cursor;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.RandomAccess;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.cell.CellImgFactory;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.integer.IntType;
import net.imglib2.view.Views;

/**
 * 
 * @author Benoit Lombardot
 * 
 */

/*
 * 
 * TODO:
 * 	[x] handle in place and modification of the input in the defaultLabelAlgorithm
 * 		[x] add a constructor
 * 		[x] implement convertToInteger(...)
 * 		[x] implement duplicate(...)
 *  [x] Add the missing MultiScaleMaxima 
 *	[x] make all label algorithm extend that class
 *		[x] HMaxima
 *		[x] AreaMaxima
 *		[x] SeededWatershed
 *		[x] HWatershed
 *		[x] CCL
 *		[x] Maxima
 *		[x] WindowMaxima
 *		[x] MultiScaleMaxima
 *	[-] make a minimal test for each label algorithm
 *		[-] HMaxima
 *		[-] AreaMaxima
 *		[-] SeededWatershed
 *		[-] HWatershed
 *		[-] CCL
 *		[-] Maxima
 *		[-] Window Maxima
 *		[-] MultiScale Maxima
 * 	[-] Create getExtrema methods for all maxima detection function
 * 
 * 
 * label algorithm:
 * 	- always return an integer type
 *  - never modify the input
 *  	- if process() modifies the input it is replaced by a copy of itself
 *  	- if the process can be done in place the input is copied and discretized in labelMap (i.e. the output)
 *  		- in case the input is a real type to avoid loosing all the information in the process image with range 
 *  		inferior to 255 are rescaled to fit that range.
 * 
 * 
 */

public class DefaultLabelAlgorithm < T extends RealType<T> & NativeType<T> >  implements Runnable{
	
	// 
	public enum InputProcessing {
		NONE ,
		MODIFIED ;
	}
	
	protected RandomAccessibleInterval<T> input=null;
	protected boolean isProcessed = false;
	long processTime=-1;
	int progress = 0;
	protected long numberOfLabels=0; 
	int minNumberOfLevel = 255;
	
	protected RandomAccessibleInterval<IntType> labelMap = null;
	
	boolean duplicate = false;
	boolean inPlace = false;
	
	public DefaultLabelAlgorithm( RandomAccessibleInterval<T> input, InputProcessing inputProcessing ) {
		
		switch( inputProcessing )
		{
		case NONE:
			this.input = input;
			break;
		
//		case INPLACE:
//			this.labelMap = convertToInteger(input);
//			scaleFactorProcessed= false;
//			break;
			
		case MODIFIED:
			this.input = duplicate(input);
			break;
		}
		
		
	}
	
	public DefaultLabelAlgorithm( RandomAccessibleInterval<T> input ) {
		
		this.input = input;
		
	}
	
	@Override
	public void run() {
		long StartTime = System.nanoTime();
		process();
		processTime = System.nanoTime() - StartTime;
		isProcessed=true;
	}
	
	
	protected void process() {
		// to be implemented in daughter class
	}
	
	
	public RandomAccessibleInterval<IntType> getLabelMap() {
		if( ! isProcessed )
			run();
	
		return labelMap;
	}

	
	
	public int getProgress() {
		return progress;
	}

	
	public long getProcessTime() {
		return processTime;
	}
	
	
	public long getNumberOfLabels() {
		return numberOfLabels;
	}
	
	
	// TODO: guarantee that the factory will be the same
	private RandomAccessibleInterval<T> duplicate( RandomAccessibleInterval<T> input ){
		
		
		// determine an appropriate factory for the output
		long nPixel = Views.iterable( input ).size();
		T valT = input.randomAccess().get().createVariable();
		RandomAccessibleInterval<T> output = null;
		if( nPixel < Math.pow(2, 31) ){
			output = new ArrayImgFactory<T>().create( input, valT);
		}
		else{
			output = new CellImgFactory<T>().create( input, valT );
		}

		// copy input into output
		final Cursor< T > out = Views.iterable( output ).cursor();
		final RandomAccess< T > in = input.randomAccess();

		while( out.hasNext() )
		{
			out.fwd();
			in.setPosition( out );
			out.get().set( in.get() );
		}

		return output;
	}

	
	
	
	
	
	
	

}
