package invizio.imgalgo.label.rleccl;

import java.util.List;

import invizio.imgalgo.util.Pixel;
import net.imglib2.Interval;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.type.numeric.RealType;
import net.imglib2.view.Views;

/**
 * 
 * @author Benoit Lombardot
 * 
 */

public class RleImgView extends AbstractRleImg
{
	
	RleImg sourceRleImg;
	long[] origin;
	Interval interval;
	
	public RleImgView( RleImg sourceRleImg, Interval interval)
	{
		this.sourceRleImg = sourceRleImg;
		this.interval = interval;
		
		nDim = sourceRleImg.getDim().length;
		
		dim = new long[nDim];
		origin = new long[nDim];
		for(int d=0; d<nDim ; d++) { 
			dim[d] = interval.max(d) - interval.min(d) + 1;
			origin[d] = interval.min(d);
		}
		
		
		rleDim = new long[nDim-1];
		nLines=1;
		for(int d=1; d<nDim ; d++) {
			rleDim[d-1] = dim[d];
			nLines *= dim[d];
		}
		
		
	}

	public <T extends RealType<T>> RleImgView( RleImg sourceRleImg, RandomAccessibleInterval<T> raiView, float threshold)
	{
		this( sourceRleImg, (Interval) raiView );		
		this.initialize(raiView, threshold);
		
	}

	@Override
	public List<PixelRun> getLine(int index) {
		
		// convert i to a position in the source rleImg;
		long[] posInSource = Pixel.getPosFromIdx( index, rleDim );
		for ( int d=0 ; d<nDim-1 ; d++)
			posInSource[d] += origin[d+1];
		
		int indexInSource = (int) Pixel.getIdxFromPos(posInSource, sourceRleImg.getRleDim() );
		
		return sourceRleImg.getLine( indexInSource );
	}

	
	@Override
	public void setLine(int index, List<PixelRun> runLine )
	{
		// convert i to a position in the source rleImg;
		long[] posInSource = Pixel.getPosFromIdx( index, rleDim );
		for ( int d=0 ; d<nDim-1 ; d++)
			posInSource[d] += origin[d+1];
		
		int indexInSource = (int) Pixel.getIdxFromPos(posInSource, sourceRleImg.getRleDim() );
		
		sourceRleImg.setLine( indexInSource, runLine );
		
	}
	
	
	public RleImg getSourceRleImg() {
		return sourceRleImg;
	}

}
