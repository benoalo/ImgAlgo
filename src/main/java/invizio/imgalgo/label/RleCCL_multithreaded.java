package invizio.imgalgo.label;


import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;
import java.util.Queue;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorCompletionService;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;

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
import ij.IJ;
import invizio.imgalgo.label.rleccl.Chunck;
import invizio.imgalgo.label.rleccl.MinSurfaceChuncker;
import invizio.imgalgo.label.rleccl.PixelRun;
import invizio.imgalgo.label.rleccl.RleImg;
import invizio.imgalgo.label.rleccl.RleImgDefault;
import invizio.imgalgo.label.rleccl.RleImgView;
import invizio.imgalgo.util.Pixel;

/**
 * 
 * @author Benoit Lombardot
 * 
 */


/*
 * TODO:
 *  [x] v bug in 2D
 *  [x] neighborhood issue in nD
 *  [-] test various shape and cut in 3D
 *  
 *  The trade off for the the threads seem to be a few ms. to be worthe multithreading we need tasks in 10s of ms
 *  most of the time only the initial construction of runs and runs labeling is worth that effort
 *  remark: the main speed gain is on
 *  
 */

public class RleCCL_multithreaded  < T extends RealType<T> & NativeType<T> > extends DefaultLabelAlgorithm<T>
{

	private int numberOfThreads;
	private final float threshold;
	private final int nDim;
	
	
	public RleCCL_multithreaded(RandomAccessibleInterval<T> input, float threshold )
	{
		
		super(input);
		
		this.numberOfThreads = Runtime.getRuntime().availableProcessors();
		this.threshold = threshold;
		this.nDim = input.numDimensions();	
				
	}
	
	
	public RleCCL_multithreaded(RandomAccessibleInterval<T> input, float threshold, int numberOfThreads )
	{
		
		super(input);
		
		this.numberOfThreads = numberOfThreads;
		this.threshold = threshold;
		this.nDim = input.numDimensions();	
				
	}
	
	
	
	@Override
	protected void process() 
	{
		long StartTime  = System.nanoTime();
		
		// define chuncks in the data
		int numberOfChuncks = numberOfThreads;
		int[] dimsToChunck = new int[nDim-1];
		for(int i=0; i<nDim-1; i++){ 
			dimsToChunck[i] = i+1 ; 
		}
		MinSurfaceChuncker<T> chuncker = new MinSurfaceChuncker<T>( input, dimsToChunck, numberOfChuncks);
		List< Chunck<T> > tree = chuncker.getChuncks();
		
		long dt = System.nanoTime() - StartTime;
		System.out.println("0 - Create chunck, dt = " + dt);
		
		StartTime = System.nanoTime();
		
		// label individual chuncks
		
		//toMerge = new LinkedList<Chunck<T>>();
		ExecutorCompletionService<Chunck<T>> executorCS = new ExecutorCompletionService<Chunck<T>>( Executors.newFixedThreadPool( numberOfThreads ) );
		
		dt = System.nanoTime() - StartTime;
		System.out.println("1 - instantiate executor service, dt = " + dt);
		
		StartTime = System.nanoTime();
		
		
		for( Chunck<T> chunck : tree )
		{
			if ( chunck.children.size() == 0  )
			{
				RleLabeler worker = new RleLabeler( chunck , input );
				chunck.worker = worker;
				executorCS.submit( worker );
			}	
		}
		
		
		for(int i=0; i<tree.size(); i++)
		{
			Chunck<T> chunck=null;
			try {
				chunck = executorCS.take().get();
			} catch (InterruptedException | ExecutionException e) {
				e.printStackTrace();
			}
			if( chunck == null ) {
				continue;
			}
			chunck.isProcessed = true;
			if( chunck.parent != null) 
			{
				Chunck<T> chunck0 = chunck.parent.children.get(0); 
				Chunck<T> chunck1 = chunck.parent.children.get(1); 
				if( chunck0.isProcessed && chunck1.isProcessed )
				{
					RleMerger merger = new RleMerger(chunck.parent);
					chunck.parent.worker = merger;
					executorCS.submit( merger );				
				}
			}
			else {
				break;
			}
		}
		
		System.out.println("3 - Chunck labelling and merging, dt = " + (System.nanoTime()-StartTime) );
		
	
		StartTime = System.nanoTime();
		// Create label image from chuncks and unified parent
		//	- create labelMap container
		ImgFactory<IntType> imgFactory = Util.getArrayOrCellImgFactory( input, new IntType(0) );
		labelMap = imgFactory.create( input, new IntType(0) );

		Chunck<T> rootChunck = null;
		for( Chunck<T> chunck : tree )
			if ( chunck.parent == null  )
				rootChunck = chunck;
		RleImg runImg = rootChunck.worker.runImg;
		List<Integer> parent = rootChunck.worker.parent;
		
		//
		int newLabel = 0;
		for( int i=1 ; i<parent.size() ; i++  ) {
			int pari = parent.get(i);
			if( pari == i ) {
				//System.out.println( "pre/final:" + i + "/" + parent.get(i) + "  -> root");
				parent.set(i, ++newLabel );
			}
			else {
				//System.out.println( "pre/final:" + i + "/" + parent.get(pari) );
				parent.set(i, parent.get(pari) );
			}
		}
		
		//  - write the image 
		//List<Future<?>> futures2 = new ArrayList<Future<?>>();
		int count=0;
		for( Chunck<T> chunck : tree )
		{
			if ( chunck.children.size() == 0  )
			{
				RleImg partialRleImg = new RleImgView( runImg , chunck );
				RandomAccessibleInterval<IntType> partialLabelMap = Views.interval(labelMap, chunck);
				LabelWriter worker = new LabelWriter( partialLabelMap, partialRleImg , parent ); 
				executorCS.submit( worker );
				count++;
			}
			
		}
		
		// wait for the threads to finish
		for( int i=0; i<count; i++)
		{
			try {
				executorCS.take();
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}
		
		//LabelWriter worker = new LabelWriter( labelMap, runImg , parent );
		//worker.call();
		System.out.println("4 - final labeling, dt = " + (System.nanoTime()-StartTime) );
		
	}
	
	
	public abstract class RleWorker implements Callable<Chunck<T>>{

		Chunck<T> chunck;
		RleImg runImg;
		List<Integer> parent;
		RandomAccessibleInterval<T> chunckedRai;
		
		public void run() {
			// implement in daugther class
		}
		
		private RleImg getRleImg() {
			return runImg;
		}

		private List<Integer> getLabelEquivalence() {
			return parent;
		}
		
	}
	
	
	class RleLabeler extends RleWorker {
		
		RleLabeler( Chunck<T> chunck, RandomAccessibleInterval<T> rai )
		{
			this.chunck = chunck;
			this.chunckedRai = Views.interval( rai , chunck );
		}
		
		@Override
		public Chunck<T> call() throws Exception {
			//long startTime = System.nanoTime();
			RleCCL<T> labeler = new RleCCL<T>( chunckedRai , threshold );
			labeler.run();
			//long dt = labeler.getProcessTime();
			runImg = labeler.getRleImg();
			parent = labeler.getLabelEquivalence();
			
			//Integer parent = chunck.parent==null ? chunck.id : chunck.parent.id;
			//System.out.println("labeling: "+chunck.id+"  ;  par="+parent + "  ;  interval: "+chunck.toString() + "  ;  t="+dt );
			//System.out.println("labeling: "+chunck.id+"  ;  par="+parent + "  ;  interval: "+chunck.toString() );
			//System.out.println("runImg: dim="+Arrays.toString(runImg.getRleDim())+"  ;  nLines="+runImg.getNumberOfLines() );
				

			return chunck;
			/*
			chunck.isProcessed = true;
			
			if( chunck.parent == null )
				return chunck;
			
			Chunck<T> chunck0 = chunck.parent.children.get(0); 
			Chunck<T> chunck1 = chunck.parent.children.get(1); 
			if( chunck0.isProcessed && chunck1.isProcessed && !chunck.parent.isActive)
			{
				chunck.parent.isActive=true; //
				RleMerger merger = new RleMerger(chunck.parent);
				chunck.parent.worker = merger;
				toMerge.add(chunck.parent);
				System.out.println("labeling: "+chunck.id+"  ;  par="+parent + "  ; merging started  ;  t="+(System.nanoTime()-startTime) );				
			}
			else {
				System.out.println("labeling: "+chunck.id+"  ;  par="+parent + "  ; merging not started  ;  t="+(System.nanoTime()-startTime) );
			}
		
			*/
		}
		
		
	}
	
	
	
	
	
	class RleMerger extends RleWorker
	{
		
		//long[] rleDim;
		int mergeDim;
		private  long[][] neighDeltaPos;
		private  int[] neighLinearOffset;
		private  int nNeigh;
		
		public RleMerger(Chunck<T> chunck)
		{
			this.chunck = chunck;
			this.mergeDim = chunck.cutDim;
			// define the neighborhoods
		}

		@Override
		public Chunck<T> call() {
			//long startTime = System.nanoTime();
			
			Chunck<T> chunck0 = chunck.children.get(0); 
			Chunck<T> chunck1 = chunck.children.get(1); 
			RleImg rleImg0 = chunck0.worker.getRleImg();
			RleImg rleImg1 = chunck1.worker.getRleImg();
			List<Integer> labelEq0 = chunck0.worker.getLabelEquivalence();
			List<Integer> labelEq1 = chunck1.worker.getLabelEquivalence();
			
						
			// merge chuncks together
	        
			////////////////////////////////////////////////////////
			//	1 - merge parent
	        ////////////////////////////////////////////////////////
			//System.out.println( Arrays.toString(labelEq0.toArray()) );
			//System.out.println( Arrays.toString(labelEq1.toArray()));
			
			
			int nLabel0 = labelEq0.size()-1;
			labelEq1.remove(0);
			for( int i=0; i<labelEq1.size() ; i++ ) {
				labelEq1.set(i, labelEq1.get(i)+nLabel0 );
			}
			parent = new ArrayList<Integer>();
			for(Integer val : labelEq0)
				parent.add( val );
			for(Integer val : labelEq1)
				parent.add( val );
			
			//System.out.println(Arrays.toString(parent.toArray()));
			
	        ////////////////////////////////////////////////////////
			//	2 - offset the runs in volume 1
	        ////////////////////////////////////////////////////////

			int nRunLines1 = rleImg1.getNumberOfLines();
			for( int i=0; i<nRunLines1; i++ ) {
				for(PixelRun run : rleImg1.getLine(i) ) {
					run.label += nLabel0; 
				}
			}
			
			
	        ////////////////////////////////////////////////////////
			//	3 - rework equivalence at the boundary between chuncks
	        ////////////////////////////////////////////////////////
			
			//		- define neighborhood
			long[] rleDim0 = rleImg0.getRleDim();
			long[][] neighDeltaPos0 = Pixel.getConnectivityPos(nDim-1, Pixel.Connectivity.FULL);
			neighDeltaPos = new long[(int)Math.pow(3, nDim-2)][];
			int count=0;
			for(int i=0; i<neighDeltaPos0.length; i++) {
				if( neighDeltaPos0[i][mergeDim-1]==-1 ) {
					neighDeltaPos[count] = neighDeltaPos0[i];
					count++;
				}
			}
			neighLinearOffset = Pixel.getIdxOffsetToCenterPix(neighDeltaPos, rleDim0 );
	        nNeigh = neighLinearOffset.length;

	        // 		- analyze the first hyperslice in vol1 along the merging dimension
	        //		- for each line of runs in the hyperslice, perform union find
	        long[] rleDim1 = rleImg1.getRleDim();
			
	        long[] hDim = new long[ rleDim0.length];
	        long nPos = 1;
	        for ( int d=0 ; d<rleDim0.length ; d++ ) {
	        	if ( d == mergeDim-1 )
	        		hDim[d]=1;
	        	else
	        		hDim[d] = rleDim0[d];
	        	nPos *= hDim[d];
	        }
	        
	        //System.out.println("hDim="+Arrays.toString(hDim) +"  ;  nPos="+nPos);
	        
	        long[] pos0 = new long[rleDim0.length];
	        for( int ih=0 ; ih<nPos ; ih++) {
	        	// convert i1 to a position in rleImg1, no translation needed since hyperslice and rleImg1 have the same origin
				final long[] pos1 = Pixel.getPosFromIdx(ih, hDim);
				// convert pos1 to index in rleImg1 
				int i1 = (int) Pixel.getIdxFromPos( pos1 , rleDim1 );
				// get the current line of runs
				final List<PixelRun> currentLine1 = rleImg1.getLine( i1 );
				
				//debug
				//System.out.println("ih="+ih+"  ;  posh="+Arrays.toString(pos1));
				//System.out.println("i1="+i1+"  ;  pos1="+Arrays.toString(pos1));
				//System.out.println("currentLine.size="+ currentLine1.size());
				
				
				//get the equivalent of pos1 in rleImg0
				pos0 = pos1;
				pos0[mergeDim-1] += rleImg0.getDim()[mergeDim];
				
				// get the neighbor lines of current line in rleImg1
				final List<List<PixelRun>> neighborLines0 = getNeighborLines(pos0, rleDim0, rleImg0);

				//System.out.println("pos0="+Arrays.toString(pos0));
				//System.out.println("nNeighbor="+ neighborLines0.size());

				// for each neighbor union find run in the lines (this will update the parent array)
				updateLabelEquivalence( currentLine1 , neighborLines0);
				
	        }
	        
	        
	        ////////////////////////////////////////////////////////
	        //		- flatten the merged parent
	        ////////////////////////////////////////////////////////
	        
	        for( int i=1 ; i<parent.size() ; i++ )
			{
				final int root_i = findRoot(i);
				parent.set( i , root_i );
				//if( parent_i != i)
				//{
				//	parent.set( i , parent.get( root_i ) );
				//}
			}
	        
	        ////////////////////////////////////////////////////////
	        //	- create a unified rleImg
	        ////////////////////////////////////////////////////////
	        
	        runImg = new RleImgDefault( rleImg0 , rleImg1, chunck.cutDim );
			
			return chunck;
	        /*
			chunck.isProcessed = true;
			
			long dt = System.nanoTime() - startTime ;
			
			// for debug
			Integer parent = chunck.parent==null ? chunck.id : chunck.parent.id;
			//System.out.println("merging: "+chunck.id+"  ;  par="+parent + "  ;  interval: "+chunck.toString() );
			
			
			if( chunck.parent == null )
			{	
				System.out.println("merging: "+chunck.id+"  ;  par="+parent + "  ;  all mergin done!  ;  dt: "+(System.nanoTime()-startTime) );
				return;
			}
			
			if( chunck.parent.children.get(0).isProcessed && chunck.parent.children.get(1).isProcessed)
			{
				RleMerger merger = new RleMerger(chunck.parent);
				chunck.parent.worker = merger;
				//futures.add( executor.submit( merger ) );
				System.out.println("merging: "+chunck.id+"  ;  par="+parent + "  ;  merging started  ;  dt: "+(System.nanoTime()-startTime) );
			}
			else {
				System.out.println("merging: "+chunck.id+"  ;  par="+parent + "  ;  merging not started  ;  dt: "+(System.nanoTime()-startTime) );
			}
			*/
		}
		
		
		
		private List<List<PixelRun>> getNeighborLines(long[] currentPos, long[] rleDim0, RleImg rleImg0)
		{
			
			int nRleDim = rleDim0.length;
			List<List<PixelRun>> neighborLines = new ArrayList<List<PixelRun>>();
			int i0 = (int)Pixel.getIdxFromPos( currentPos, rleDim0);
			for( int n=0 ; n<nNeigh ; n++)
			{
				final long[] neighPos = new long[nRleDim];
				for(int d=0; d<nRleDim ; d++)
					neighPos[d] = currentPos[d] + neighDeltaPos[n][d]; 
				if( rleImg0.isInBound( neighPos ) )	
				{
					final List<PixelRun> neighLine = rleImg0.getLine( i0 + neighLinearOffset[n] );
					neighborLines.add( neighLine );
				}
			}
			return neighborLines;
		}

		
		private void updateLabelEquivalence( final List<PixelRun> currentLine , List<List<PixelRun>> neighborLines )
		{
			int[] nRunIdx = new int[nNeigh];

			for( PixelRun cRun : currentLine )
			{
				//int label = cRun.label;
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
							//if( cRun.label<label)
								union( nRun , cRun );
							//else
							//	cRun.label = findRoot(nRun.label);
						}
						if( nRun.end <= cRun.end  && (++nRunIdx[neighLineIdx]) < numberOfNeighborRun) {
							// then go to the next neighbor run 
							nRun = neighborLine.get( nRunIdx[neighLineIdx] );
						}
						else {
							// there are no other neighbor run in contact with current run
							// go to the next neighbor line
							break;
						}	// end if neighbor run finishes before current run
					}	// end while nrun in contact with crun in current neighborline
				}	// end for neighbor line
				
			}	// end for run in currentline
			
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
			final int nRoot = findRoot(nRun.label);
			int cRoot = findRoot(cRun.label);
			
			if( cRoot < nRoot ) {
				parent.set(nRoot, cRoot);
				//cRun.label = cRoot;
			}
			else if( cRoot > nRoot ) 
			{
				parent.set(cRoot, nRoot);
				cRun.label = nRoot;
			}
			
		}

		
		
	}
	
	
	class LabelWriter implements Callable<Chunck<T>>{
		
		List<Integer> parent;
		RleImg runImg;
		 RandomAccessibleInterval<IntType> partialLabelMap;
		 
		public LabelWriter( RandomAccessibleInterval<IntType> partialLabelMap, RleImg runImg, List<Integer> parent) {
			
			this.parent = parent;
			this.runImg = runImg;
			this.partialLabelMap = partialLabelMap;
			
		}

		@Override
		public Chunck<T> call() {
			long startTime = System.nanoTime();
			final Cursor<IntType> cursor = Views.flatIterable( partialLabelMap ).cursor();
			final int lineLength = (int) partialLabelMap.dimension(0);
			int prevRunEnd = 0;
			
			//System.out.println("labelWriter: nLines="+runImg.getNumberOfLines());
			//System.out.println("labelWriter: n run="+runImg.getLine(126).size() );
			
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
			System.out.println("chunck label writing : "+(System.nanoTime()-startTime));
			
			return null;
		}
		
	}
	
	
	public static void main(final String... args) throws IOException
	{
		ImageJ ij = new ImageJ();
		ij.ui().showUI();
		
		//Dataset dataset = (Dataset) ij.io().open("C:/Users/Ben/workspace/testImages/blobs32.tif");
		//@SuppressWarnings("unchecked")
		//Img<FloatType> img = (Img<FloatType>) dataset.getImgPlus();
		
		Img<FloatType> img =ImageJFunctions.wrap( IJ.openImage("C:/Users/Ben/workspace/testImages/noise5000_std50_blur10.tif") );
		//Img<FloatType> img =ImageJFunctions.wrap( IJ.openImage("F:/projects/blobs32.tif") );
		float threshold = 0.5f;
		int nThread = 2;
		int nIter=100;
		int warmup = 10;
		
		long dt1 = 0;
		long min_dt1 = Long.MAX_VALUE;
		RandomAccessibleInterval<IntType> labelMap1=null;
		long dt2 = 0;
		long min_dt2 = Long.MAX_VALUE;
		RandomAccessibleInterval<IntType> labelMap2=null;
		//RleCCL<FloatType> labeler1 = null;
		//RandomAccessibleInterval<IntType> labelMap1 = null;
		for(int i=0 ; i<nIter ; i++) 
		{
			RleCCL_multithreaded<FloatType> labelParallel = new RleCCL_multithreaded<FloatType>(img, threshold, nThread);
			labelMap1 = labelParallel.getLabelMap();
			long dt = labelParallel.getProcessTime();
			if( i>=warmup)
				dt1 += dt;
			min_dt1 = min_dt1 > dt ? dt : min_dt1;
			
			RleCCL<FloatType> label = new RleCCL<FloatType>(img, threshold);
			labelMap2 = label.getLabelMap();
			dt = label.getProcessTime();
			if( i>=warmup)
				dt2 += dt;
			min_dt2 = min_dt2 > dt ? dt : min_dt2;
			//System.out.println("iter "+ i + ": " + dt );
		}
		System.out.println("dt1 " + (dt1/(nIter-warmup)));
		System.out.println("min1 " + (min_dt1));
		System.out.println("dt2 " + (dt2/(nIter-warmup)));
		System.out.println("min2 " + (min_dt2));

		
		ImageJFunctions.show( labelMap1 );
		ImageJFunctions.show( labelMap2 );
		
		System.out.println("done!");
		
		
		
		
	}
	
	
	

	
	
	
	
	
	
	
	
	
}
