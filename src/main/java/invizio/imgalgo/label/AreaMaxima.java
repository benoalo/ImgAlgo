package invizio.imgalgo.label;

import invizio.imgalgo.util.Pixel;
import invizio.imgalgo.util.RAI;
import net.imglib2.Cursor;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.array.ArrayImg;
import net.imglib2.img.basictypeaccess.IntAccess;
import net.imglib2.img.basictypeaccess.array.IntArray;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.outofbounds.OutOfBounds;
import net.imglib2.outofbounds.OutOfBoundsConstantValueFactory;
import net.imglib2.outofbounds.OutOfBoundsFactory;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.integer.IntType;
import net.imglib2.type.numeric.real.AbstractRealType;
import net.imglib2.util.Fraction;
import net.imglib2.view.Views;


/**
 * 
 * @author Benoit Lombardot
 * 
 */

/*
 * TODO:
 *	[-] handle when HMaxima is handled (since the code are twins,  only the attribute change) 
 * 	[-] make the code img factory agnostic
 * 		[-] replace the use of parent isActivePeak and crriteria (if next point is possible 
 * 			then both would be resolved at the same time, plus a good speed boost)
 * 		[-] queue could be replaced by hierarchical FIFO + a LIFO queue could also be setup for the last step of the algorithm 
 *  [-] run the algorithm in place
 * 	[x] extend a LabelAlgorithm class
 *  [x] use RAI rather than Img		

 */



public class AreaMaxima < T extends RealType<T> & NativeType<T> > extends DefaultLabelAlgorithm<T>  {

	
	// parameter
	float AreaThresh=0;
	float threshold=0;
	
	// data structure for the algorithm
	int[] parent;
	boolean[] is_ActivePeak;
	int[] criteria;
	
	
	
	public AreaMaxima(RandomAccessibleInterval<T> input, float threshold, float AreaMin) {
		
		super( input );
		
		T sample = input.randomAccess().get().createVariable();
		this.getScaleFactor();
		if( scaleFactor>1) //&& sample instanceof AbstractRealType)
			this.input = RAI.scale( input , scaleFactor, offset );
		
		//ImageJFunctions.wrap(this.input, "test").show();
		
		this.threshold = (threshold+offset)*scaleFactor;
		this.AreaThresh = AreaMin;
		//System.out.println("Threshold:"+this.threshold);
		//System.out.println("AreaThreshold:"+this.AreaThresh);
		//System.out.println("offset:"+this.offset);
		//System.out.println("scalefactor:"+this.scaleFactor);
		
		
		
	}
	
	@Override
	protected void process()
	{
				
		long numPixel = Views.iterable( input ).size();
		int ndim = input.numDimensions();
 		long[] dimensions = new long[ndim];
 		input.dimensions(dimensions);
 		
 		// get image min and max
		Cursor<T> in_cursor = Views.iterable( input ).cursor();
		float aux = in_cursor.next().getRealFloat();
		float min = aux;
		float max = aux;

		while ( in_cursor.hasNext() )
		{
			float val = in_cursor.next().getRealFloat();
			if(val<min)       {  min = val;  }
			else if( val>max ){  max = val;  }
		}
		
		
		//if( max-min < this.minNumberOfLevel ) {
		//	input = this.duplicate(input );
		//}
		
		
		
		min = Math.max((int)min, (int)threshold+1); // the threshold is excluded from the pixel that will be processed
		max = Math.max(min, max);
		
		
		//////////////////////////////////////////////////////////////////////////////////////
		// build an ordered list of the pixel ////////////////////////////////////////////////

		 // get image histogram
		int nlevel = (int)max - (int)min + 1;
        int[] histo_Int = new int[ nlevel ];
        
        // create a flat iterable cursor
// 		long[] minindim = new long[ ndim ], maxindim = new long[ ndim ];
//        for ( int d = 0; d < ndim; ++d ){   minindim[ d ] = 0 ;    maxindim[ d ] = dimensions[d] - 1 ;  }
//        FinalInterval interval = new FinalInterval( minindim, maxindim );
        
        final Cursor< T > in_flatcursor = Views.flatIterable( Views.interval( input, input)).cursor();
        

        while ( in_flatcursor.hasNext() )
		{
        	int level = ((int)in_flatcursor.next().getRealFloat() - (int)min); 
        	if(level>=0)
        	{
        		histo_Int[ level ]++;
        	}
		}

        // get each level start to build the sorted list
        int[] posInLevel = new int[ nlevel ];
        posInLevel[nlevel-1] = 0;
        for (int i = nlevel - 1; i > 0; i--)
            posInLevel[ i-1 ] = posInLevel[ i ] + histo_Int[ i ];

        // build a sorted list // higher level are put first // for each level pix are in lexicographic order
        int[] Sorted_Pix = new int[(int) numPixel];
        in_flatcursor.reset(); 
        int idx = 0; 
        

        while ( in_flatcursor.hasNext() )
        {
        	int level = ((int)in_flatcursor.next().getRealFloat() - (int)min);
        	if(level>=0)
        	{
        		Sorted_Pix[posInLevel[level]] = idx;
        		posInLevel[level]++;
        	}
        	idx++;
        }
            
		//////////////////////////////////////////////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////

		
		// define image dealing with out of bound
        T outOfBoundValue = input.randomAccess().get().createVariable();
 		final OutOfBoundsFactory< T, RandomAccessibleInterval< T >> oobImageFactory =  new OutOfBoundsConstantValueFactory< T, RandomAccessibleInterval< T >>( outOfBoundValue );
 		final OutOfBounds< T > input_X = oobImageFactory.create( input );
 		RandomAccess<T> input_RA = input.randomAccess();
 				
 		// define the connectivity
 		long[][] neigh = Pixel.getConnectivityPos(ndim, Pixel.Connectivity.FULL);
 		int[] n_offset = Pixel.getIdxOffsetToCenterPix(neigh, dimensions);
 		long[][] neighMov = Pixel.getSuccessiveMove(neigh);
 		int nNeigh = neigh.length;
 		
 		// first pass in sorted order of the pixel
 		parent = new int[(int) numPixel ];
 		int pval, nval, rval;
 		int nidx, ridx;
 		long[] position = new long[ndim];
 		
 		is_ActivePeak = new boolean[(int) numPixel ];
 		for( int i = 0; i < is_ActivePeak.length; i++){  is_ActivePeak[i] = false;}
 		criteria = new int[(int) numPixel ];
 		
 		for(int pidx: Sorted_Pix)
		{
			Pixel.getPosFromIdx(pidx, position, dimensions);//getPosFromIdx( pidx, dimensions, position);
			input_X.setPosition(position);
			pval = (int)input_X.get().getRealFloat();
			
			if(pval<=threshold){ break; }
			
			parent[pidx]=pidx;
			is_ActivePeak[pidx]=true;
			criteria[pidx] = 1;
			
			//boolean is_ActiveAux = true;
			
			// for each neighbor
			for(int i=0; i<nNeigh; i++)
			{
				input_X.move(neighMov[i]);
				if( input_X.isOutOfBounds() ){  continue;  }
				nval = (int)input_X.get().getRealFloat();
				nidx = pidx + n_offset[i];
				
				if ( (pval < nval)  |  ((pval == nval) & (nidx < pidx)) ) // test if n was already processed 
				{	
					ridx = FindRoot(nidx);
					if( ridx == pidx) { continue; }
					if ( ! is_ActivePeak[ridx])
					{
						is_ActivePeak[pidx] = false;
						continue;
					}
					
					Pixel.getPosFromIdx(ridx, position, dimensions);//getPosFromIdx( ridx, dimensions, position);
					input_RA.setPosition(position);
					rval = (int)input_RA.get().getRealFloat();
					
					union(ridx, pidx, rval, pval, AreaThresh);
									
				}
			}	
		}
		
		// label the image
		int current_label=1;
		for (int i = Sorted_Pix.length - 1; i >= 0; i--)
        {
			idx = Sorted_Pix[i];
            if ( parent[idx] != idx ) {
                parent[idx] = parent[parent[idx]];
            }
            else{ 
            	if( criteria[idx]>=AreaThresh){
            		//parent[idx] = criteria[idx]; // to color with the peak volume
            		parent[idx] = current_label; // to color with label
            		current_label++;
            	}
            	else{
            		parent[idx] = 0;
            	}           	
            }
        }
		numberOfLabels = current_label-1;
		
		// create an output image from the label array
		final IntAccess access = new IntArray( parent );
		final Fraction frac = new Fraction(1,1);
		final ArrayImg<IntType, IntAccess> array =	new ArrayImg<IntType, IntAccess>( access, dimensions, frac );// create a Type that is linked to the container
		final IntType linkedType = new IntType( array );
		// pass it to the DirectAccessContainer
		array.setLinkedType( linkedType );	
		
		labelMap = array;

		return;

	
	}
	
	
	private int FindRoot(int n)
    {
        if (parent[n] != n)
        {
            return FindRoot(parent[n]);
        }
        else { return n; }
    }
	
	
	private void union(int r, int p, int valr, int valp, float AreaThresh)
	{
		int cr = criteria[r];
		if (  ( valr == valp && cr<=AreaThresh)  |  ( cr<AreaThresh )  )
		{
			criteria[p] = criteria[p] + criteria[r];
			criteria[r]=0;
			parent[r] = p;
			is_ActivePeak[r] = false;
			//System.out.println("criteria: "+ criteria[p] + "  ;  valp: "+valp);
		}
		else // without this else. a new region are restarted indefinitely when previous one reaches the max Area 
		  is_ActivePeak[p]=false;  
		
		return;
	}

	
//	protected static void getPosFromIdx(int idx, long[] dimensions, long[] position)
//	{
//		int ndim = dimensions.length;
//		for( int d=0; d<ndim; d++)
//		{
//			position[d] = idx % dimensions[d];
//			idx /=  dimensions[d] ;
//		}
//	}
	
	
	
//	public <T extends RealType<T> > Img<IntType> AreaMaxima(Img<T> input, final int AreaThresh)
//	{
//		float threshold = (float)input.firstElement().createVariable().getMinValue();
//		return AreaMaxima(input, AreaThresh, threshold);
//	}



	
	
	
}
