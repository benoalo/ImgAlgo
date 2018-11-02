package invizio.imgalgo.util;

import net.imglib2.Cursor;
import net.imglib2.FinalInterval;
import net.imglib2.IterableInterval;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.ImgFactory;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.integer.IntType;
import net.imglib2.util.Util;
import net.imglib2.view.Views;

/**
 * 
 * @author Benoit Lombardot
 * 
 */

public class RAI {

	/**
	 * Duplicate a RandomAccessibleInterval<T>
	 * @param input 
	 * @return a copy of the input with the same type
	 */
	// TODO: guarantee that the factory will be the same
	public static <T extends NativeType<T>> RandomAccessibleInterval<T> duplicate( RandomAccessibleInterval<T> input ){
		
		// create an appropriate factory for the output
		T sampleT = input.randomAccess().get().createVariable();
		ImgFactory<T> imgFactory = Util.getArrayOrCellImgFactory( input, sampleT );
		RandomAccessibleInterval<T> output = imgFactory.create( input );
		
		
		// copy input into output (and adapt the process depending whether the two RAI have the same iteration order)
		
		final IterableInterval<T> input_I = Views.iterable( input );
		final IterableInterval<T> output_I = Views.iterable( output );
		
		if(  input_I.iterationOrder().equals( output_I.iterationOrder() )  )
		{
			final Cursor< T > in = input_I.cursor();
			final Cursor< T > out = output_I.cursor();
			while( out.hasNext() )
			{
				out.next().set( in.next() );
			}
		}
		else
		{
			final Cursor< T > out = output_I.cursor();
			final RandomAccess< T > in = input.randomAccess();
	
			while( out.hasNext() )
			{
				out.fwd();
				in.setPosition( out );
				out.get().set( in.get() );
			}
		}
		
		return output;
	}
	
	
	
	/**
	 * this function down-sample the input image and return an image which size
	 * is input.dimension(i)*ratio[i])
	 * 
	 * @param input
	 *            an image to down-sample
	 * @param DownSampleFactor
	 *            the down-sampling factor for each dimension (value should be
	 *            superior to 1 for down-sampling)
	 * @return a down-sampled version of input
	 */
	public static <T extends NativeType<T> > RandomAccessibleInterval<T> decimate(RandomAccessibleInterval<T> input, double[] DownSampleFactor) {
		// create a new image
		int nDim = input.numDimensions();
		long[] out_dims = new long[nDim];
		for (int i = 0; i < nDim; i++) {
			out_dims[i] = (int) (input.dimension(i) / DownSampleFactor[i]);
		}
		
		// create an empty RAI for the output
		T sampleT = input.randomAccess().get().createVariable();
		ImgFactory<T> imgFactory = Util.getArrayOrCellImgFactory( new FinalInterval(out_dims), sampleT );
		RandomAccessibleInterval<T> output = imgFactory.create( out_dims );
		
		
		
		// iterate through the image
		RandomAccess<T> input_RA = input.randomAccess();
		Cursor<T> out_cursor = Views.iterable( output ).localizingCursor();
		int[] out_pos = new int[nDim];
		int[] in_pos = new int[nDim];
		while (out_cursor.hasNext()) {
			out_cursor.fwd();
			out_cursor.localize(out_pos);
			for (int d = 0; d < nDim; d++) {
				in_pos[d] = (int) (out_pos[d] * DownSampleFactor[d]);
			}
			input_RA.setPosition(in_pos);
			out_cursor.get().set( input_RA.get() );
		}

		return output;
	}
	
	
	
	/**
	 * Converts an image to an integer image with a certain scaling and offset 
	 * @param input
	 * @param scaleFactor
	 * @param offset
	 * @return a rescaled image
	 */
	public static <T extends RealType<T> > RandomAccessibleInterval<IntType> convertToInteger( RandomAccessibleInterval<T> input, float scaleFactor, float offset ){
		
		
		// determine an appropriate factory for the output
		ImgFactory<IntType> imgFactory = Util.getArrayOrCellImgFactory( input, new IntType(0) );
		RandomAccessibleInterval<IntType> output =  imgFactory.create( input);

		// copy input into output
		// could such conversion benefit more infrastructure if using Convert api in imglib2  
		
		if(  Views.iterable(input).iterationOrder().equals( Views.iterable(output).iterationOrder() )  )
		{
			final Cursor< T > in = Views.iterable(input).cursor();
			final Cursor< IntType > out = Views.iterable(output).cursor();
			while( out.hasNext() )
			{
				out.next().setReal( ( in.next().getRealFloat() + offset ) * scaleFactor  );
			}
		}
		else
		{
			final Cursor< IntType > out = Views.iterable(output).cursor();
			final RandomAccess< T > in = input.randomAccess();
	
			while( out.hasNext() )
			{
				out.fwd();
				in.setPosition( out );
				out.get().setReal( ( in.get().getRealFloat() + offset ) * scaleFactor );
			}
		}
		
		return output;
	}
	
	/**
	 * Creates a new image with a certain scaling and offset of pixel values in input 
	 * @param input
	 * @param scaleFactor
	 * @param offset
	 * @return a rescaled image
	 */
	public static <T extends RealType<T> > RandomAccessibleInterval<T> scale( RandomAccessibleInterval<T> input, float scaleFactor, float offset ){
		
		
		// determine an appropriate factory for the output
		ImgFactory<T> imgFactory = Util.getSuitableImgFactory(input, input.randomAccess().get().createVariable() );
		RandomAccessibleInterval<T> output =  imgFactory.create( input );

		// copy input into output
		// could such conversion benefit more infrastructure if using Convert api in imglib2  
		
		if(  Views.iterable(input).iterationOrder().equals( Views.iterable(output).iterationOrder() )  )
		{
			final Cursor< T > in = Views.iterable(input).cursor();
			final Cursor< T > out = Views.iterable(output).cursor();
			while( out.hasNext() )
			{
				out.next().setReal( ( in.next().getRealFloat() + offset ) * scaleFactor  );
			}
		}
		else
		{
			final Cursor< T > out = Views.iterable(output).cursor();
			final RandomAccess< T > in = input.randomAccess();
	
			while( out.hasNext() )
			{
				out.fwd();
				in.setPosition( out );
				out.get().setReal( ( in.get().getRealFloat() + offset ) * scaleFactor );
			}
		}
		
		return output;
	}
	
	
}
