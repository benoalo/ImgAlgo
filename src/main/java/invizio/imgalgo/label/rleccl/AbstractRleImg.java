package invizio.imgalgo.label.rleccl;

import java.util.ArrayList;
import java.util.List;

import net.imglib2.Cursor;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.type.numeric.RealType;
import net.imglib2.view.Views;

/**
 * 
 * @author Benoit Lombardot
 * 
 */

public class AbstractRleImg implements RleImg {

	List< List<PixelRun> > rleLineList = null;
	int nDim;
	long[] dim;
	long[] rleDim;
	int nLines=1;
	
	
	
	@Override
	public int getNumberOfLines() {
		return nLines;
	}

	@Override
	public long[] getRleDim() {
		return rleDim;
	}

	@Override
	public long[] getDim() {
		return dim;
	}

	@Override
	public List<PixelRun> getLine(int index)
	{
		return rleLineList.get(index);
	}
	
	
	@Override
	public boolean isInBound( long[] pos )
	{
		for(int d = 0 ; d < nDim-1 ; d++ )
		{
			int coordinate = (int)pos[d];
			if( coordinate < 0 || coordinate >= rleDim[d] )
			{
				return false;
			}	
		}	
		return true;
	}

	@Override
	public void setLine(int index, List<PixelRun> runLine)
	{
		rleLineList.set( index, runLine );
	}

	
	protected <T extends RealType<T>> void initialize(RandomAccessibleInterval<T> rai, Float threshold)
	{
		
		// create the list of PixelRun for each line
		
		
		Cursor<T> cursor = Views.flatIterable( rai ).cursor();
		long lineLength = dim[0];
		
		for(  int lineIndex = 0  ;  lineIndex < nLines  ;  lineIndex++  )
		{
			List<PixelRun> rleLine = new ArrayList<PixelRun>();
			int start=-1;
			boolean prev = false;
			for( int pixIndex = 0 ; pixIndex<lineLength ; pixIndex++)
			{
				final T t = cursor.next();
				if( t.getRealFloat() > threshold )
				{
					if( !prev ) { 
						start=pixIndex;
					}
					prev = true;
				}
				else
				{
					if( prev ) {
						rleLine.add( new PixelRun(start, pixIndex) );
						start=-1;
					}
					prev = false;
				}
			}
			if( start>=0) // modify from > to >= (2018-06-22)
				rleLine.add( new PixelRun(start, (int)lineLength) );
			
			this.setLine(lineIndex, rleLine);
			//rleLineList.add( rleLine );
			
		}
	}
	
}
