package invizio.imgalgo.label;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import ij.IJ;
import invizio.imgalgo.label.rleccl.PixelRun;
import invizio.imgalgo.label.rleccl.RleImg;
import invizio.imgalgo.label.rleccl.RleImgDefault;
import invizio.imgalgo.util.Pixel;
import net.imagej.ImageJ;
import net.imglib2.Cursor;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.img.ImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.integer.IntType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.util.Util;
import net.imglib2.view.Views;

/**
 * 
 * @author Benoit Lombardot
 * 
 */

public class RleCCL < T extends RealType<T> & NativeType<T> > extends DefaultLabelAlgorithm<T>  {
	
	//private float threshold;
	private List<Integer> parent; 
	private RleImg runImg = null;
	private float threshold;
	
	private final int nDim;	
	private long[] rleDim;
	private  long[][] neighDeltaPos;
	private  int[] neighLinearOffset;
	
	
	
	public List<Integer> getLabelEquivalence() {
		return parent;
	}

	public RleImg getRleImg() {
		return runImg;
	}

	
	
	
	public RleCCL(RandomAccessibleInterval<T> input , float threshold ) {
		
		super(input);
		
		this.threshold = threshold;
		parent = new ArrayList<Integer>();
		nDim = input.numDimensions();
		
		
	}
	
	
	public RleCCL( RleImg runImg )
	{
		super( null );
		this.runImg = runImg;
		parent = new ArrayList<Integer>();
		nDim = runImg.getDim().length;
	}
		
	
	
	
	@Override
	protected void process()
	{
		// create a run representation of the image
		if ( runImg == null )
			runImg = new RleImgDefault( input , threshold );
		
		rleDim = runImg.getRleDim();
		
		neighDeltaPos = Pixel.getConnectivityPos(nDim-1, Pixel.Connectivity.LEXICO_FULL);
		neighLinearOffset = Pixel.getIdxOffsetToCenterPix(neighDeltaPos, rleDim );
        
		        
		
		// label each run line according the neighboring run lines already processed
		int label=1;
		parent.add(0);
		parent.add(1);
		
		for( int lineIndex=0 ; lineIndex < runImg.getNumberOfLines() ; lineIndex++ )
		{
			final List<PixelRun> currentLine = runImg.getLine( lineIndex );
			List<List<PixelRun>> neighborLines = getNeighborLines( lineIndex );
			label = updateLabelEquivalence( currentLine , neighborLines, label );			
		}
		parent.remove( parent.size()-1);
		
		
		// flatten the label look up
		numberOfLabels=0;
		for( int i=1 ; i<parent.size() ; i++ )
		{
			final int parent_i = findRoot(i);//parent.get(i); 
			
			if( parent_i != i)
			{
				parent.set( i , parent.get( parent_i ) );
			}
			else {
				++numberOfLabels;
			}
		}
		
	}
	
	
	
	
	private int updateLabelEquivalence( final List<PixelRun> currentLine , List<List<PixelRun>> neighborLines, int label )
	{
		int nNeigh = neighLinearOffset.length;
		int[] nRunIdx = new int[nNeigh];

		for( PixelRun cRun : currentLine )
		{
			cRun.label = label;
			for( int neighLineIdx=0 ; neighLineIdx<neighborLines.size() ; neighLineIdx++ )
			{
				List<PixelRun> neighborLine = neighborLines.get( neighLineIdx );
				int numberOfNeighborRun = neighborLine.size();
				if (nRunIdx[neighLineIdx] >= numberOfNeighborRun ) {
					continue ;
				}
				PixelRun nRun = neighborLine.get( nRunIdx[neighLineIdx]);
				while( cRun.end >= nRun.start )
				{
					if( cRun.start <= nRun.end)
					{
						if( cRun.label<label)
							union( nRun , cRun );
						else
							cRun.label = findRoot(nRun.label);
					}
					if( nRun.end <= cRun.end  && (++nRunIdx[neighLineIdx]) < numberOfNeighborRun) {
						// then go to the next neighbor run 
						nRun = neighborLine.get( nRunIdx[neighLineIdx] );
					}
					else {
						// there are no other neighbor run in contact with current run
						// go to the next neighbor line
						break;
					}
					
				}
					
			}
			if( cRun.label >= label )
				parent.add( ++label );
		}
		
		
		return label;
	}
	
	
	@Override
	public RandomAccessibleInterval<IntType> getLabelMap() {
		if( ! isProcessed )
			run();
		
		if ( labelMap==null ) {
			long StartTime = System.nanoTime();
			createLabelMap();
			processTime += System.nanoTime() - StartTime;
		}
		return labelMap;
	}
	
	
	
	
	private void createLabelMap()
	{
		// go through the runs and create a label image
		
		int newLabel = 0;
		for( int i=1 ; i<parent.size() ; i++  ) {
			int pari = parent.get(i);
			if( pari == i )
				parent.set(i, ++newLabel );
			else
				parent.set(i, parent.get(pari) );
		}

		
		
		ImgFactory<IntType> imgFactory = Util.getArrayOrCellImgFactory( input, new IntType(0) );
		labelMap = imgFactory.create( input, new IntType(0) );
					
		final Cursor<IntType> cursor = Views.flatIterable( labelMap ).cursor();
		final int lineLength = (int) input.dimension(0);
		int prevRunEnd = 0;
		for( int lineIndex=0 ; lineIndex < runImg.getNumberOfLines() ; lineIndex++ )
		{
			
			final List<PixelRun> currentLine = runImg.getLine( lineIndex );
			for( PixelRun cRun : currentLine )
			{
				// position cursor at the beginning of the run
				final long jump = cRun.start - prevRunEnd;
				if( jump>0)
					cursor.jumpFwd( jump );
				
				// iterate till the end of the run
				final int runLength =  cRun.end - cRun.start;
				final int value = parent.get( cRun.label );
				for( int i=0 ; i<runLength ; i++ )
					cursor.next().set( value );
				prevRunEnd = cRun.end;
			}
			
			// go to the beginning of next line
			cursor.jumpFwd( lineLength - prevRunEnd);
			prevRunEnd = 0;
		}
	}
	

	
	
	private List<List<PixelRun>> getNeighborLines(int lineIndex)
	{
		final int nNeigh = neighLinearOffset.length;
		List<List<PixelRun>> neighborLines = new ArrayList<List<PixelRun>>();
		long[] currentPos = new long[nDim-1];
		Pixel.getPosFromIdx(lineIndex, currentPos, rleDim );
		for( int n=0 ; n<nNeigh ; n++)
		{
			final long[] neighPos = new long[nDim-1];
			for(int d=0; d<nDim-1 ; d++)
				neighPos[d] = currentPos[d] + neighDeltaPos[n][d]; 
			if( runImg.isInBound( neighPos ) )	
			{
				final List<PixelRun> neighLine = runImg.getLine( lineIndex + neighLinearOffset[n] );
				neighborLines.add( neighLine );
			}
		}
		return neighborLines;
	}

	
	
	private int findRoot(int p)
	{
		int parent_p = parent.get(p);
		if( parent_p == p )
			return p;
		int r = findRoot(parent_p );
		parent.set( p , r );
		return r;
	}
	
	
	
	private void union(PixelRun nRun, PixelRun cRun)
	{
		final int r = findRoot(nRun.label);
		if( cRun.label < r )
			parent.set(r, cRun.label);
		else if( cRun.label > r ) 
		{
			parent.set(cRun.label, r);
			cRun.label = r;
		}
		
	}
	
	
	
		
	
	
	
	

	
	
	
	
	public static void main(final String... args) throws IOException
	{
		ImageJ ij = new ImageJ();
		ij.ui().showUI();
		
		//Dataset dataset = (Dataset) ij.io().open("C:/Users/Ben/workspace/testImages/blobs32.tif"); //rleLabelTest/toLabel2.tif");
		//@SuppressWarnings("unchecked")
		//Img<FloatType> img = (Img<FloatType>) dataset.getImgPlus();
		//Img<FloatType> img2 = ij.op().convert().float32( img );
		//Img<FloatType> img =ImageJFunctions.wrap( IJ.openImage("C:/Users/Ben/workspace/testImages/blobs32.tif") );
		//Img<FloatType> img =ImageJFunctions.wrap( IJ.openImage("F:/projects/2DEmbryoSection_Mette_contourMaskInv.tif") );
		Img<FloatType> img = ImageJFunctions.wrap( IJ.openImage("C:/Users/Ben/workspace/testImages/noise2000_std50_blur10.tif") );
		
		//img = (Img<FloatType>) ij.op().math().multiply( img, new FloatType(-1) );
		
		
		float threshold = 0.5f;
		int nIter=20;
		int nWarmup=10;
		
		long dt1 = 0;
		RleCCL<FloatType> labeler1 = null;
		RandomAccessibleInterval<IntType> labelMap1 = null;
		for(int i=0 ; i<nIter ; i++) 
		{
			labeler1 = new RleCCL<FloatType>( img , threshold);
			labelMap1 = labeler1.getLabelMap();
			if(i>=nWarmup) {
				long dt = labeler1.getProcessTime();
				dt1 += dt;
			}
		}
		System.out.println("dt1 " + (dt1/(nIter-nWarmup)));

		
		long dt2 = 0;
		CCL<FloatType> labeler2 = null;
		RandomAccessibleInterval<IntType> labelMap2 = null;
		for(int i=0 ; i<nIter ; i++) 
		{
			labeler2 = new CCL<FloatType>( img , threshold);
			labelMap2 = labeler2.getLabelMap();
			if(i>=nWarmup) {
				long dt = labeler2.getProcessTime();
				dt2 += dt;
			}
		}
		System.out.println("dt2 " + (dt2/(nIter-nWarmup)));

		
		//ImageJFunctions.wrap(img, "image").show();
		//ImageJFunctions.wrap(labelMap, "label map, n="+labeler.getNumberOfLabels() ).show();
		
		ij.ui().show(img);
		ij.ui().show(labelMap1);
		ij.ui().show(labelMap2);
		
	}
	
	
	
}
