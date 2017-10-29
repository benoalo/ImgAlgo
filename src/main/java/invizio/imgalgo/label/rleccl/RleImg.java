package invizio.imgalgo.label.rleccl;

import java.util.List;

/**
 * 
 * @author Benoit Lombardot
 * 
 */


public interface RleImg {

	public int getNumberOfLines();
	
	public long[] getRleDim();
	
	public long[] getDim();
	
	public List<PixelRun> getLine(int index);
	
	public void setLine(int index , List<PixelRun> runLine );
	
	public boolean isInBound( long[] pos );
	
}
