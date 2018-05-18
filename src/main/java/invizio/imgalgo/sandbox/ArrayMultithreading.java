package invizio.imgalgo.sandbox;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;



public class ArrayMultithreading {

	
	// use case 1:
	//		- copy a source array to mutiple destination arrays
	//		- copy Img to multiple Img
	
	// use case 2:
	//		- write multiple source arrays to a destination array
	//		- write multiple source arrayImgs to a destination arrayImg
	
	int sourceSize;
	int nThreads;
	int[] source;
	List<int[]> dests;
	
	
	public  void initCopyOnetoMany(int sourceSize, int nThreads)
	{
		this.nThreads = nThreads;		
		this.source = new int[sourceSize];
		for( int i=0; i<sourceSize; i++)
		{
			source[i] = i;
		}
		
		dests = new ArrayList<int[]>();
		for( int i=0; i<nThreads; i++)
		{
			int destSize = sourceSize / nThreads ; 
			if ( i == nThreads-1 ) {
				destSize = sourceSize - sourceSize / nThreads*(nThreads-1);
			}
			dests.add( new int[destSize]);
		}
		
	}
	
	public  int runCopyOnetoMany()
	{

		ExecutorService executor = Executors.newFixedThreadPool(nThreads);
		int start=0;
		List<Future<Integer>> futures = new ArrayList<Future<Integer>>();
		for( int i=0; i<nThreads; i++)
		{
			int[] dest = dests.get(i);
			int len  = dest.length;
			Future<Integer> future = executor.submit( new copyWorker(source, dest, start, len ) );
			futures.add(future);
			start += len;
		}
		
		
		Integer sum=0;
		for( Future<Integer> future : futures )
		{
			try {
				sum += future.get();
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
		
		return sum;
		
	}
	
	
	class copyWorker implements Callable<Integer> {
		int[] source; 
		int[] dest;
		int start;
		int length;
		
		public copyWorker(int[] source, int[] dest, int start, int length) {
			this.source = source;
			this.dest = dest;
			this.start = start;
			this.length = length;
		}
		
		@Override
		public Integer call() throws Exception {
			
			int sum=0;
			for( int i=0; i<length; i++)
			{
				int val = source[start+i];
				dest[i] = val;
				for( int j=1; j<1000; j++)
					sum += val + (start+i) / j ;
			}
			//System.out.println( sum );
			return sum;
		}
		
	}

	
	public static void main(final String... args)
	{
		int warmupRuns = 5;
		int runs = 10;
		int size = 1*1000000;
		int nThreads = Runtime.getRuntime().availableProcessors();
		
		int [] nThreadsList = new int[] {1, 2, 4, 5, 6, 9, 12};
		for(int j=0; j<nThreadsList.length; j++)
		{
			nThreads = nThreadsList[j]; // (int) Math.pow(2, j);
			
			long dt0 = 0;
			for(int i=0; i<warmupRuns; i++)
			{
				ArrayMultithreading testCopyOneToMany = new ArrayMultithreading();
				testCopyOneToMany.initCopyOnetoMany(size, nThreads);
				long start =  System.nanoTime();
				testCopyOneToMany.runCopyOnetoMany();
				dt0 += System.nanoTime() - start;
			}
			
			
			int sum = 0;
			long dt = 0;
			for(int i=0; i<runs; i++)		
			{
				ArrayMultithreading testCopyOneToMany = new ArrayMultithreading();
				testCopyOneToMany.initCopyOnetoMany(size, nThreads);
				
				long start =  System.nanoTime();
				sum += testCopyOneToMany.runCopyOnetoMany();
				dt  += System.nanoTime() - start ;
			}
			
			System.out.println( "\nnThreads: "+ nThreads );
			System.out.println( "warmup, time per run: "+ (dt0/warmupRuns) );
			System.out.println( "time per run: "+ (dt/runs) );
			System.out.println( "sum: "+ sum );
		}

		
	}
	
}
