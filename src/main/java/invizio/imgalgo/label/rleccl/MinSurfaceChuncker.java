package invizio.imgalgo.label.rleccl;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;
import java.util.Queue;

import net.imglib2.Interval;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.RealType;



/**
 * 
 * @author Benoit Lombardot
 * 
 */


//TODO: 
//	[-] check that the dimension does not get smaller than one
	


public class MinSurfaceChuncker < T extends RealType<T> & NativeType<T> >
{
	
	Interval interval0;
	long[] dims;
	int nDim;
	
	int[] dimsToSplit;
	int numberOfChuncks;
	
	List<Chunck<T>> tree;
	boolean treeProcessed = false;
	
	
	
	public List<Chunck<T>> getChuncks(){
		
		if( ! treeProcessed )
			process();
		
		return tree;
	}
	
	
	
	public MinSurfaceChuncker( Interval interval0 ,int[] dimsToSplit, int numberOfChuncks)
	{
		this.numberOfChuncks = numberOfChuncks;
		this.dimsToSplit = dimsToSplit; 
		this.interval0 = interval0;
		
		nDim = interval0.numDimensions();
		dims = new long[nDim];
		interval0.dimensions(dims);	
	}
	
	
	
	protected void process(){
		
		int[] splitPerDim = defineSplitPerDim( dimsToSplit );
		
		tree = createMergingTree( dimsToSplit, splitPerDim);
		
	}

	
	
	private int[] defineSplitPerDim( int[] dimsToSplit )
	{
		long vol = 1;
		for(int d=0; d<nDim; d++) {
			vol *= dims[d];
		}
		
		int numberOfSplitableDim = dimsToSplit.length;
		int[] splitPerDim = new int[numberOfSplitableDim];
		long[] surfPerCut = new long[numberOfSplitableDim];
		for(int i=0; i<numberOfSplitableDim ; i++) {
			splitPerDim[i] = 1;
			surfPerCut[i] = vol/dims[dimsToSplit[i]];
			//System.out.println("cuDim:"+dimsToSplit[i]+"  ;  surf:"+surfPerCut[i]);

		}
		
		
		if ( numberOfSplitableDim == 1 ){
			splitPerDim[0] = numberOfChuncks;
		}
		else if ( numberOfSplitableDim == 2)
		{
			List<Integer> factors = getFactors( numberOfChuncks );
			long minSurface = numberOfChuncks * surfPerCut[0];
		
			// brute force minimisation
			int nFactor = factors.size();
			int nCase = (int)Math.pow(2, nFactor);
			for( int i=0; i<nCase ; i++){
				int mul1=1;
				int mul2=1;
				int val = i;
				for(int n=nFactor-1 ; n>=0 ; n--) {
					int xx = val / (int)Math.pow(2, n);
					val -= xx*(int)Math.pow(2, n);
					final int factor = factors.get(n); 
					if( xx==1 )
						mul1 *= factor;
					else
						mul2 *= factor;
				}
				long surface =  (mul1-1) * surfPerCut[0] + (mul2-1) * surfPerCut[1];
				//System.out.println("surface:"+surface);
				if (surface <= minSurface ) {
					minSurface = surface;
					splitPerDim[0] = mul1;
					splitPerDim[1] = mul2;
				}	
			}
		}
		
		//System.out.println("split per dim :"+ Arrays.toString(splitPerDim));
		// brute for force minimization is possible for larger dimensionality
		// but might become impractical when the number of factors grows
		// i.e. the process should remain 2 order of magnitude faster than the image processing
		
		return splitPerDim;
	}
	
	
	
	
	
	private List<Chunck<T>> createMergingTree(int[] dimsToSplit, int[] splitPerDim)
	{ 
		
		int id = 0;
		List<Chunck<T>> tree = new ArrayList<Chunck<T>>();
		Chunck<T> initChunck = new Chunck<T>(interval0 );
		initChunck.id = id;
		tree.add( initChunck );
		
		for( int i=0 ; i<dimsToSplit.length ; i++  )
		{
			int dim = dimsToSplit[i];
			Queue<Chunck<T>> Q = new LinkedList<Chunck<T>>();
			
			// find the leafs of the tree, add them to the queue and initialize their value
			for(Chunck<T> chunck : tree) {
				if( chunck.children.size() == 0 ) {
					if( splitPerDim[ i ]>1 ) {
						chunck.val = splitPerDim[ i ];
						Q.add(chunck);
					}
				}
			}
			
			// build the tree
			while( ! Q.isEmpty() )
			{
				Chunck<T> chunck = Q.poll();
				int val = chunck.val;
				int val0 = val / 2;
				int val1 = val - val0;
				
				long[] min0 = new long[nDim];
				long[] max0 = new long[nDim];
				long[] min1 = new long[nDim];
				long[] max1 = new long[nDim];
				
				chunck.min(min0);
				chunck.min(min1);
				chunck.max(max0);
				chunck.max(max1);
				
				long width = (max0[dim]-min0[dim])/val*val0;
				max0[dim] = min0[dim]+width;
				min1[dim] = max0[dim]+1;
				
				Chunck<T> chunck0 = new Chunck<T>( min0, max0);
				chunck0.id = ++id;
				chunck0.val = val0;
				chunck0.parent = chunck;
				Chunck<T> chunck1 = new Chunck<T>( min1, max1);
				chunck1.id = ++id;
				chunck1.val = val1;
				chunck1.parent = chunck;
				
				tree.add(chunck0);
				tree.add(chunck1);
				
				if ( val0>1 )
					Q.add(chunck0);
				if ( val1>1 )
					Q.add(chunck1);
				
				chunck.children.add(chunck0);
				chunck.children.add(chunck1);	
				chunck.cutDim = dim;
			} // end while

		} // end for splitDim
					
		return tree;
	}
	
	
	
	private int[] prime = new int[] {2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97};
	List<Integer> divisor;
	
	// decompose number in a product of prime numbers
	// helper function to estimate the min surface chuncks splitting in 2 dimensions 
	private List<Integer> getFactors(int number)
	{
		// decompose number in a product of prime numbers
		divisor = new ArrayList<Integer>();
		int idPrime = 0;
		while( number > 1 )
		{
			if( (number % prime[idPrime]) != 0 ){
				idPrime++;
			}
			else {
				number = number/prime[idPrime];
				divisor.add(prime[idPrime]);
			}
		}
		return divisor;
	}
	
	
	
}
