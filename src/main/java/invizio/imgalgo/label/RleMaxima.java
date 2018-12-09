package invizio.imgalgo.label;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.python.antlr.runtime.misc.Stats;

import ij.IJ;
import ij.ImagePlus;
import ij.measure.Measurements;
import ij.plugin.Duplicator;
import ij.plugin.ImageCalculator;
import ij.process.ImageStatistics;
import ij.process.StackStatistics;
import invizio.imgalgo.label.rleccl.PixelRun;
import invizio.imgalgo.label.rleccl.RleImg;
import invizio.imgalgo.label.rleccl.RleMaxImg;
import invizio.imgalgo.label.rleccl.ValuedPixelRun;
import invizio.imgalgo.util.Pixel;
import net.imagej.ImageJ;
import net.imglib2.Cursor;
import net.imglib2.RandomAccess;
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


// update the run creation to select only potential max run
// update the run labeling to label only the real maxima run

public class RleMaxima < T extends RealType<T> & NativeType<T> > extends DefaultLabelAlgorithm<T>  {
	
	//private float threshold;
	private List<Integer> parent; 
	private RleImg runImg = null;
	private float threshold;
	
	private final int nDim;	
	private long[] rleDim;
	private  long[][] neighDeltaPos;
	private  int[] neighLinearOffset;
	
	private boolean labelSetAtParentRoot=false;
	
	
	public List<Integer> getLabelEquivalence() {
		return parent;
	}

	public RleImg getRleImg() {
		
		if( ! isProcessed  )
			run();
		
		setLabelAtParentRoots();
		
		// update runImg label to their parent root
		for( int lineIndex=0 ; lineIndex < runImg.getNumberOfLines() ; lineIndex++ )
		{
			final List<PixelRun> currentLine = runImg.getLine( lineIndex );
			final List<PixelRun> newRunLine = new ArrayList<PixelRun>();
			for( PixelRun cRun0 : currentLine )
			{	
				final ValuedPixelRun cRun = (ValuedPixelRun) cRun0;
				final int labelValue = parent.get( cRun.label );
				if( labelValue>0 ) {
					// position cursor at the beginning of the run
					cRun0.label = labelValue;
					newRunLine.add(cRun0);
				}
			}
			// set the new runline
			runImg.setLine(lineIndex, newRunLine);
		}
		
		return runImg;
	}

	
	
	
	public RleMaxima(RandomAccessibleInterval<T> input , float threshold ) {
		
		super(input);
		
		this.threshold = threshold;
		parent = new ArrayList<Integer>();
		nDim = input.numDimensions();
		
		
	}
	
	
	public RleMaxima( RleImg runImg )
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
			runImg = new RleMaxImg( input , threshold );
		
		rleDim = runImg.getRleDim();
		neighDeltaPos = Pixel.getConnectivityPos(nDim-1, Pixel.Connectivity.INVERT_LEXICO_FULL);
		neighLinearOffset = Pixel.getIdxOffsetToCenterPix(neighDeltaPos, rleDim );
		long[][] neighDeltaPos2 = Pixel.getConnectivityPos(nDim-1, Pixel.Connectivity.LEXICO_FULL);
		int[] neighLinearOffset2 = Pixel.getIdxOffsetToCenterPix(neighDeltaPos2, rleDim );
		
		// test if run is possible maxima
		// strategy 1: look at the complete neighborhood of a run for validation in the 1st pass
		// strategy 2: look one neighbor line at each pass (n neighbor line --> n real raster pass )
		// strategy 3: for each run line do: for each neighbor line do: check run valididity 
		//             --> 1 pass, trying to limit the number of jump between neighbor lines
		// strategy 4: do one raster pass through the pixels, checking all neighbor lines runs for validity
		//			   --> 1 real raster pass, tool to get neigh line already set up
		        
		// strategy 4
		
		
		final RandomAccess<T> input_RA = input.randomAccess();
		long[] pos = new long[nDim];
		final int lineLength = (int)input.dimension(0);
		for( int lineIndex=0 ; lineIndex < runImg.getNumberOfLines() ; lineIndex++ )
		{
			List<List<PixelRun>> neighborLines = getNeighborLines( lineIndex );
			List<List<PixelRun>> neighborLines2 = getNeighborLines( lineIndex, neighLinearOffset2, neighDeltaPos2);
			
			long[] linePos = new long[nDim-1];
			Pixel.getPosFromIdx(lineIndex, linePos, rleDim );
			for( int d=1; d<nDim; d++) {
				pos[d] = linePos[d-1];
			}
			
			// make a list of neighbor runs sorted by start pixel
			List<ValuedPixelRun> neighRunList = new ArrayList<ValuedPixelRun>();
			for(List<PixelRun> nruns : neighborLines) {
				for( PixelRun nrun0 : nruns ) {
					final ValuedPixelRun nrun = (ValuedPixelRun) nrun0;
					if( nrun.valid ) {
						nrun.label=-1;
						neighRunList.add( nrun );
					}
				}
			}
			for(List<PixelRun> nruns : neighborLines2) {
				for( PixelRun nrun0 : nruns ) {
					final ValuedPixelRun nrun = (ValuedPixelRun) nrun0;
					if( nrun.valid )
						neighRunList.add( nrun );
				}
			}
			neighRunList.sort( (PixelRun r1, PixelRun r2)->r1.start-r2.start );
			
			// for each pix in map test validity of mapped runs
			//final Iterator<ValuedPixelRun> runIterator = neighRunList.iterator(); 
			//while( runIterator.hasNext() ) {
			for( ValuedPixelRun run : neighRunList ) {
						
				//final ValuedPixelRun run = runIterator.next();
				int pos0 = Math.max(run.start-1, 0);
				pos[0] = pos0;
				input_RA.setPosition(pos);
				final int iMax = Math.min(run.end+1, lineLength);
				final float runValue = run.value;
				int countPlateau=0;
				if( run.label==-1)
					countPlateau = run.plato;
				else
					countPlateau = run.postplato;
				int iPrev=-10;
				for(int i=pos0; i<iMax; i++) {
					final float pixValue = input_RA.get().getRealFloat();
					
					if( pixValue > runValue ) {
						run.valid=false;
						break;
					}
					// count the number of plateau
					if ( pixValue == runValue ) {
						if( i>iPrev+1 )
							countPlateau++;
						iPrev=i;
					}
						
					input_RA.fwd(0); // move one pixel along dimension 0
				}
				if( run.label==-1)
					run.plato = countPlateau;
				else
					run.postplato = countPlateau;
				run.label=0;
			}
		}
		
				
		// label each run line according the neighboring run lines already processed
		
		neighDeltaPos = Pixel.getConnectivityPos(nDim-1, Pixel.Connectivity.LEXICO_FULL);
		neighLinearOffset = Pixel.getIdxOffsetToCenterPix(neighDeltaPos, rleDim );
		
		
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
		
		// final handling of post plato issue
		for( int lineIndex=0 ; lineIndex < runImg.getNumberOfLines() ; lineIndex++ )
		{
			final List<PixelRun> currentLine = runImg.getLine( lineIndex );
			for(PixelRun run0: currentLine) {
				ValuedPixelRun run = (ValuedPixelRun) run0;
				if ( run.postplato>=1 )
					parent.set( findRoot(run.label) , 0);
			}
		}
		
		
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
		final int nNeigh = neighLinearOffset.length;
		final int[] nRunIdx = new int[nNeigh];

		for( PixelRun cRun0 : currentLine )
		{
			final ValuedPixelRun cRun = (ValuedPixelRun) cRun0;  //
			if ( ! cRun.valid )
				cRun.label = 0;
			else
				cRun.label = label;
			
			// int cRunNeighCount = 0;
			for( int neighLineIdx=0 ; neighLineIdx<neighborLines.size() ; neighLineIdx++ )
			{
				final List<PixelRun> neighborLine = neighborLines.get( neighLineIdx );
				int numberOfNeighborRun = neighborLine.size();
				if (nRunIdx[neighLineIdx] >= numberOfNeighborRun ) {
					continue ;
				}
				ValuedPixelRun nRun = (ValuedPixelRun) neighborLine.get( nRunIdx[neighLineIdx]);
				while( cRun.end >= nRun.start )
				{
					if( cRun.start <= nRun.end)
					{
						if( cRun.value == nRun.value)
						{	
							cRun.plato--;
							nRun.postplato--;
							
							if( cRun.label<label)
								union( nRun , cRun );
							else
								cRun.label = findRoot(nRun.label);
						}
						else if( cRun.value > nRun.value ) {
							final int r = findRoot(nRun.label);
							parent.set(r, 0);
						}
						
					}
					if( nRun.end <= cRun.end  && (++nRunIdx[neighLineIdx]) < numberOfNeighborRun) {
						// then go to the next neighbor run 
						nRun = (ValuedPixelRun) neighborLine.get( nRunIdx[neighLineIdx] );
					}
					else {
						// there are no other neighbor run in contact with current run
						// go to the next neighbor line
						break;
					}
					
				}
					
			}
			//if ( cRunNeighCount<=0 && cRun.plato ) {
			if ( cRun.plato>=1 ) {
				//cRun.label=0;
				parent.set( findRoot(cRun.label) , 0);
				cRun.valid=false;
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
		setLabelAtParentRoots();	
		
		ImgFactory<IntType> imgFactory = Util.getArrayOrCellImgFactory( input, new IntType(0) );
		labelMap = imgFactory.create( input );
					
		final Cursor<IntType> cursor = Views.flatIterable( labelMap ).cursor();
		final int lineLength = (int) input.dimension(0);
		int prevRunEnd = 0;
		for( int lineIndex=0 ; lineIndex < runImg.getNumberOfLines() ; lineIndex++ )
		{
			final List<PixelRun> currentLine = runImg.getLine( lineIndex );
			for( PixelRun cRun0 : currentLine )
			{	
				final ValuedPixelRun cRun = (ValuedPixelRun) cRun0;
				final int labelValue = parent.get( cRun.label );
				if( labelValue>0 ) {
					// position cursor at the beginning of the run
					final long jump = cRun.start - prevRunEnd;
					if( jump>0)
						cursor.jumpFwd( jump );
					
					// iterate till the end of the run
					final int runLength =  cRun.end - cRun.start;
					for( int i=0 ; i<runLength ; i++ )
						cursor.next().set( labelValue );
					prevRunEnd = cRun.end;
				}
			}
			
			// go to the beginning of next line
			cursor.jumpFwd( lineLength - prevRunEnd);
			prevRunEnd = 0;
		}
	}
	
	private void setLabelAtParentRoots()
	{
		if( labelSetAtParentRoot )
			return;
		
		int newLabel = 0;
		for( int i=1 ; i<parent.size() ; i++  ) {
			int pari = parent.get(i);
			if( pari == i )
				parent.set(i, ++newLabel );
			else
				parent.set(i, parent.get(pari) );
		}
		labelSetAtParentRoot = true;
		
		return;
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

	private List<List<PixelRun>> getNeighborLines(int lineIndex, int[] neighLinearOffset, long[][] neighDeltaPos)
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
	
	
	
	private void union(ValuedPixelRun nRun, ValuedPixelRun cRun)
	{
		final int r = findRoot(nRun.label);
		if( cRun.label < r ) {
			parent.set(r, cRun.label);
		}
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
		Img<FloatType> img =ImageJFunctions.wrap( IJ.openImage("C:/Users/Ben/workspace/testImages/blobs32.tif") );
		//Img<FloatType> img =ImageJFunctions.wrap( IJ.openImage("F:/projects/2DEmbryoSection_Mette_contourMaskInv.tif") );
		//Img<FloatType> img = ImageJFunctions.wrap( IJ.openImage("C:/Users/Ben/workspace/testImages/noise2000_std50_blur10.tif") );
		//Img<FloatType> img = ImageJFunctions.wrap( IJ.openImage("C:/Users/Ben/workspace/testImages/t1-head.tif") );
		//Img<FloatType> img = ImageJFunctions.wrap( IJ.openImage("C:/Users/Ben/workspace/testImages/rleMaxTest/sample_test_rle_max7.tif") );
		
		
		//img = (Img<FloatType>) ij.op().math().multiply( img, new FloatType(-1) );
		
		
		float threshold = -100f;
		int nIter=10;
		int nWarmup = 1;

		long dt1 = 0;
		RleMaxima<FloatType> labeler1 = null;
		RandomAccessibleInterval<IntType> labelMap1 = null;
		for(int i=0 ; i<nIter ; i++) 
		{
			labeler1 = new RleMaxima<FloatType>( img , threshold);
			labelMap1 = labeler1.getLabelMap();
			if( i>=nWarmup ) {
				long dt = labeler1.getProcessTime();
				dt1 += dt;
			}
		}
		System.out.println(  "dt1 " + ( dt1/(nIter-nWarmup) )  );

		
		long dt2 = 0;
		Maxima<FloatType> labeler2 = null;
		RandomAccessibleInterval<IntType> labelMap2 = null;
		for(int i=0 ; i<nIter ; i++) 
		{
			labeler2 = new Maxima<FloatType>( img , threshold);
			labelMap2 = labeler2.getLabelMap();
			if( i>=nWarmup ) {
				long dt = labeler2.getProcessTime();
				dt2 += dt;
			}
		}
		System.out.println(  "dt2 " + ( dt2/(nIter-nWarmup) )  );
		
		
		// test that label map select the same pixels
		final int nDim= labelMap1.numDimensions();
		double deltaValue = Double.MIN_VALUE;
		ImageStatistics stats = null;
		ImagePlus imp3 = null;
		
		ImagePlus imp1 = ImageJFunctions.wrap(labelMap1,"label map 1");
		imp1 = (new Duplicator()).run(imp1);
		IJ.setThreshold(imp1 ,0.5, 65535);
		
		ImagePlus imp2 = ImageJFunctions.wrap(labelMap2,"label map 2");
		imp2 = (new Duplicator()).run(imp2);
		IJ.setThreshold(imp2 ,0.5, 65535);

		if ( nDim == 2 ) { 
			IJ.run(imp1, "Convert to Mask", "");
			IJ.run(imp2, "Convert to Mask", "");
			imp3 = (new ImageCalculator()).run("Subtract create 32-bit", imp1, imp2);
			stats = imp3.getStatistics();		
		}
		else if( nDim>2 ) {
			IJ.run(imp1, "Convert to Mask", "background=Dark black");
			IJ.run(imp2, "Convert to Mask", "background=Dark black");
			imp3 = (new ImageCalculator()).run("Subtract create 32-bit stack", imp1, imp2);
			stats = new StackStatistics(imp3);
		}
		deltaValue = stats.max - stats.min;
		ij.log().log(0, "image difference is: " + deltaValue );
		
		ij.ui().show(img);
		ij.ui().show(labelMap1);
		ij.ui().show(labelMap2);
		ij.ui().show(imp3);
		
		
	}
	
	
	
}
