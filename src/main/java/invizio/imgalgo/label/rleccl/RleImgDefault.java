package invizio.imgalgo.label.rleccl;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import invizio.imgalgo.util.Pixel;
import net.imglib2.Cursor;
import net.imglib2.Interval;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.type.numeric.RealType;
import net.imglib2.view.Views;

/**
 * 
 * @author Benoit Lombardot
 * 
 */


public class RleImgDefault extends AbstractRleImg{
	
	
	public RleImgDefault( RleImg  rleImg0, RleImg  rleImg1, int mergeDim)
	{
		nDim = rleImg0.getDim().length;
		
		dim = new long[nDim];
		for( int d=0; d<nDim ; d++ ) {
			dim[d] = rleImg0.getDim()[d] ; 
		}
		dim[mergeDim] += rleImg1.getDim()[mergeDim];
		
		nLines = rleImg0.getNumberOfLines() +  rleImg1.getNumberOfLines() ;
		
		rleDim = new long[nDim-1];
		for( int d=1; d<nDim ; d++ ) {
			rleDim[d-1] = dim[d] ; 
		}

		rleLineList = new ArrayList< List<PixelRun> >(nLines);
		for( int i=0 ; i<nLines; i++ ) {
			rleLineList.add( null );
		}
		
		final long[] origin0 = new long[nDim-1];
		setRunLines( rleImg0, origin0 );
		
		long[] origin1 = new long[nDim-1];
		origin1[mergeDim-1] = rleImg0.getRleDim()[mergeDim-1];
		setRunLines( rleImg1, origin1 );
		
		//System.out.println( "origin0=" + Arrays.toString(origin0) );
		//System.out.println( "origin1=" + Arrays.toString(origin1) );
			
		
	}

	
	private void setRunLines( RleImg rleImg0 , long[] origin0 )
	{
		long[] rleDim0 = rleImg0.getRleDim();
		long[] pos = new long[nDim-1];
		for( int i0=0 ; i0<rleImg0.getNumberOfLines() ; i0++ ) {
			
			// convert i0 to a position in the merged image
			final long[] pos0 = Pixel.getPosFromIdx(i0, rleDim0);
			for(int d=0 ; d<nDim-1 ; d++) {
				pos[d] = origin0[d] + pos0[d];
			}
			
			// convert pos to an index in the merged image
			int idx = (int) Pixel.getIdxFromPos( pos , rleDim );
			//System.out.println("idx="+idx+ "  ; i0="+i0);
			// set run line in the new image
			rleLineList.set(idx, rleImg0.getLine( i0 ) );
			//if( i0==126 )
			//{
			//	List<PixelRun> lastLine = rleImg0.getLine( i0 );
				//System.out.println( "number of run : " + lastLine.size() ); 
			//}
		}
	}
	
	
	
	public RleImgDefault( Interval interval )
	{
		nDim = interval.numDimensions();
		dim = new long[nDim];
		interval.dimensions( dim );
		
		
		nLines = 1;
		rleDim = new long[nDim-1];
		for( int d=1 ; d<nDim ; d++)
		{
			rleDim[d-1] = dim[d]; 
			nLines = nLines * (int)dim[d] ;
		}
		
		rleLineList = new ArrayList< List<PixelRun> >(nLines);
		for(  int lineIndex = 0  ;  lineIndex < nLines  ;  lineIndex++  ){
			rleLineList.add(null);
		}
	}
	
	
	
	public <T extends RealType<T>> RleImgDefault( RandomAccessibleInterval<T> rai , float threshold){
		
		this(rai);
		this.initialize( rai, threshold );
		
	}
	
}

