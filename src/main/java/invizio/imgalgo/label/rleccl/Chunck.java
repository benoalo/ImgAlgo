package invizio.imgalgo.label.rleccl;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Future;

import invizio.imgalgo.label.RleCCL_multithreaded;
//import invizio.imgalgo.label.RleCCL_multithreaded2;
import net.imglib2.AbstractInterval;
import net.imglib2.Interval;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.RealType;

/**
 * 
 * @author Benoit Lombardot
 * 
 */

public class Chunck < T extends RealType<T> & NativeType<T> > extends AbstractInterval
{
	public Chunck<T> parent = null;
	public List<Chunck<T>> children = new ArrayList<Chunck<T>>();
	public boolean isProcessed = false;
	public boolean isActive = false;
	public Future<?> future;
	
	public int id = -1;
	public RleCCL_multithreaded<T>.RleWorker worker;
	//public RleCCL_multithreaded2<T>.RleWorker worker2;
	
	public int val;
	public int cutDim = -1;
	
	
	Chunck( Interval interval )
	{
		super( interval );
	} 
	
	Chunck( long[] min, long[] max)
	{
		super( min , max);
	}
	
	@Override
	public String toString(){
		String str = "";
		for(int d=0; d<this.numDimensions(); d++ )
			str += "[" + this.min(d) + "," + this.max(d) + "]" ;
		
		return str;
	}
	
}
