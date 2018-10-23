package invizio.imgalgo.label.hwatershed;


import java.io.IOException;

import ij.IJ;
import ij.ImagePlus;
import invizio.imgalgo.HierarchicalFIFO;
import invizio.imgalgo.label.DefaultLabelAlgorithm;
import invizio.imgalgo.label.Maxima;
import invizio.imgalgo.label.hwatershed.Tree;
import invizio.imgalgo.util.Pixel;
import invizio.imgalgo.util.RAI;
import net.imagej.Dataset;
import net.imagej.ImageJ;
import net.imglib2.Cursor;
import net.imglib2.FinalInterval;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.algorithm.stats.ComputeMinMax;
import net.imglib2.img.Img;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.integer.IntType;
import net.imglib2.type.numeric.real.AbstractRealType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;

/**
 * 
 * @author Benoit Lombardot
 * 
 */


/*
 * TODO:
 *  [x] make algorithm class derive from an abstract algorithm class that extends runnable and has a progress getter
 *
 * 	[-] make the HWatershedHierarchy img factory agnostic
 * 		[-] make isDequeued a randomAcssessible interval booleanType
 *  [x] make input random Accessible
 *  [x] check that output are always IntType (should be true for all label algorithm)
 *  [x] bring getPosFromIndex to a Util class
 *  [x] bring input image conversion to a Util class (conversion/copy done in mother class)
 *  [-] study the need for Connectivity (all imglib2 algo does not need to depend on Connectivty)
 *  [x] make the scaling of the input seemless from outside the class
 *  	[x] adapt input parameter
 *  	[x] adapt the tree features (Imax and HCriteria)
 *  [-] possible bug when threshold is non zeros 
 *  
 */


public class HWatershedHierarchy <T extends RealType<T> & NativeType<T>> extends DefaultLabelAlgorithm<T> {
	
	// Connectivity enum that limits the settable connectivity to Face and Full while still pointing to the standard connectivity enum
	public enum ConnectivityHWS
	{
		FACE(Pixel.Connectivity.FACE),
		FULL(Pixel.Connectivity.FULL);
		
		Pixel.Connectivity conn;
		
		ConnectivityHWS(Pixel.Connectivity conn)
		{			
			this.conn = conn;
		}
		
		Pixel.Connectivity getConn()
		{
			return conn;
		}
	}
	
	private float threshold;
	private ConnectivityHWS connectivity;
	private Tree segmentTree;
	
	
	public HWatershedHierarchy(RandomAccessibleInterval<T> input, float threshold, ConnectivityHWS connectivity)
	{
		super( input );
		
		this.getScaleFactor();
		labelMap = RAI.convertToInteger( input , scaleFactor, offset );
				
		this.threshold = (float) Math.floor((threshold + offset ) * scaleFactor);
		this.connectivity = connectivity;
	}
	
	
	
	public Tree getSegmentTree() {
		if( ! isProcessed )
			run();
		return segmentTree;
	}

	
	
	
	
	protected void process()
	{
		
		//////////////////////////////////////////////////////////////////////
		// initialisation ////////////////////////////////////////////////////
		
		IntType Tmin = labelMap.randomAccess().get().createVariable();
		IntType Tmax = Tmin.createVariable();
		ComputeMinMax.computeMinMax(labelMap, Tmin, Tmax);
				
		float min = Math.max(threshold, Tmin.getRealFloat());
		float max = Tmax.getRealFloat();
		
		// get local maxima (8/26 connected by default)
		Maxima<IntType> maxLabeler = new Maxima<IntType>(labelMap, min);
		RandomAccessibleInterval<IntType> seed = maxLabeler.getLabelMap();	
		IntType TnSeeds = new IntType(0);
		IntType Tdummy = new IntType(0);
		ComputeMinMax.computeMinMax(seed, Tdummy, TnSeeds);
		int nLeaves = (int) TnSeeds.getRealFloat();
		
		
		HierarchicalFIFO Q = new HierarchicalFIFO( min, max);
		
		int ndim = labelMap.numDimensions();
		long[] dimensions = new long[ndim]; labelMap.dimensions(dimensions);
		
		// create a flat iterable cursor
		long[] minInt = new long[ ndim ], maxInt = new long[ ndim ];
		for ( int d = 0; d < ndim; ++d ){   minInt[ d ] = 0 ;    maxInt[ d ] = dimensions[d] - 1 ;  }
		FinalInterval interval = new FinalInterval( minInt, maxInt );
		final Cursor< IntType > input_cursor = Views.flatIterable( Views.interval( labelMap, interval)).cursor();
		final Cursor< IntType > seed_cursor = Views.flatIterable( Views.interval( seed, interval)).cursor();
		
		// initialize tree and node features arrays
		double[] hCriteria = new double[2*nLeaves];
		double[] Imax = new double[2*nLeaves];
		int[] parent = new int[2*nLeaves];
		int[][] children = new int[2*nLeaves][];
		for(int i=0; i<hCriteria.length; i++)
		{
			children[i] = new int[] {-1,-1};
			parent[i]=i;
			hCriteria[i]=0;
			Imax[i]=min;
		}
		
		// fill the queue
		long idx=-1;
		while( input_cursor.hasNext() )
		{
			++idx;
			IntType pInput = input_cursor.next();
			int pVal = pInput.getInteger();
			int valSeed = seed_cursor.next().getInteger();
			if ( pVal>=min)
			{
				if ( valSeed>0)
				{
					Q.add( idx, pVal );
					pInput.setReal(min-1-valSeed);
					Imax[valSeed]= pVal;
				}
			}
			else
			{
				pInput.setReal(min-1);
			}
		}
		
		
		// extend input and seeds to to deal with out of bound
		IntType outOfBoundT = labelMap.randomAccess().get().createVariable(); 
		outOfBoundT.setReal(min-1);
		RandomAccess< IntType > input_XRA = Views.extendValue(labelMap, outOfBoundT ).randomAccess();
		RandomAccess< IntType > input_XRA2 = input_XRA.copyRandomAccess();
		
		// define the connectivity
		long[][] neigh = Pixel.getConnectivityPos(ndim, connectivity.getConn() );
		int[] n_offset = Pixel.getIdxOffsetToCenterPix(neigh, dimensions);
		long[][] dPosList = Pixel.getSuccessiveMove(neigh);
		int nNeigh = n_offset.length;
		
		
		boolean[] isDequeued = new boolean[ (int) Views.iterable( labelMap ).size() ];
		for (int i=0; i<isDequeued.length; i++)
			isDequeued[i]=false;
		
		/////////////////////////////////////////////////////////////////////////////////////
		// building the watershed and the tree //////////////////////////////////////////////
		
		int newNode = nLeaves;
		while( Q.HasNext() )
		{ 	
						
			final HierarchicalFIFO.Item pixel = Q.Next();
			final long pIdx = pixel.getIndex();
			final float pVal = Q.getCurrent_value();
			
			final long[] posCurrent = new long[ndim];
			Pixel.getPosFromIdx(pIdx, posCurrent, dimensions);
			input_XRA.setPosition(posCurrent);
			IntType p = input_XRA.get();
			int pLeaf = (int)(min - 1 - p.getRealFloat());
			int pNode = findRoot(pLeaf, parent);
			isDequeued[(int)pIdx]=true;
			
			// loop on neighbors			
			input_XRA2.setPosition(posCurrent);
			for( int i=0; i<nNeigh; i++)
			{
				final long nIdx = pIdx + n_offset[i];
				
				input_XRA2.move(dPosList[i]);
				final IntType n = input_XRA2.get();
				final float nVal = n.getRealFloat();
				
				if ( nVal != (min-1) ) // if n is in-bound
				{	
					if( isDequeued[(int)nIdx] ) // p is the lowest point 
					{	
						int nLeaf = (int)(min - 1 - nVal);
						int nNode = findRoot(nLeaf, parent);
						
						if( nNode != pNode ) // 2 distincts nodes are meeting and p is the saddle : merge Nodes
						{							
							newNode++;
							double Hn = Imax[nNode]-pVal;
							hCriteria[nNode] = Hn;
							double Hp = Imax[pNode]-pVal;
							hCriteria[pNode] = Hp;

							//merge the node with smallest H with first neighbor node that has higher dynamics
							int node1, node2;
							double H1, H2;
							if (Hp == Hn){
								node1 = pNode;
								node2 = nNode;
							}
							else if ( Hp < Hn ){
								node1 = pNode;
								node2 = nLeaf;
								H1 = Hp;
								H2 = hCriteria[node2];
								while( H2 <= H1 )
								{	
									node2 = parent[node2];
									H2 = hCriteria[node2];
								}
								
								double HMerge;
								if( parent[node2]==node2 )
									HMerge = Double.POSITIVE_INFINITY;
								else
									HMerge = Math.min( hCriteria[children[parent[node2]][0]], hCriteria[children[parent[node2]][1]] );
								
								while( H1 > HMerge )
								{
									node2 = parent[node2];
									H2 = hCriteria[node2];
									if( parent[node2]==node2 )
										HMerge = Double.POSITIVE_INFINITY;
									else
										HMerge = Math.min( hCriteria[children[parent[node2]][0]], hCriteria[children[parent[node2]][1]] );
								}
							}
							else{ // if( Hn <= Hp )
								node1 = nNode;
								node2 = pLeaf;
								H1 = Hn;
								
								H2 = hCriteria[node2];
								while( H2 <= H1 )
								{	
									node2 = parent[node2];
									H2 = hCriteria[node2];
								}
								double HMerge;
								if( parent[node2]==node2 )
									HMerge = Double.POSITIVE_INFINITY;
								else
									HMerge = Math.min( hCriteria[children[parent[node2]][0]], hCriteria[children[parent[node2]][1]] );
								
								while( H1 > HMerge )
								{
									node2 = parent[node2];
									H2 = hCriteria[node2];
									if( parent[node2]==node2 )
										HMerge = Double.POSITIVE_INFINITY;
									else
										HMerge = Math.min( hCriteria[children[parent[node2]][0]], hCriteria[children[parent[node2]][1]] );
								}
							}
							mergeNodes(node1, node2, newNode, parent, children);
							pNode=findRoot(newNode, parent);
						
							Imax[newNode]= Math.max(Imax[node1], Imax[node2]);
							hCriteria[newNode] =  Math.max(hCriteria[node1], hCriteria[node2]); //Imax[newNode]-pVal;
						}
					}
					
					if ( nVal>=min ) // is not queued yet and is in bound?
					{
						//Q.add( nIdx, (int)nVal ); // if that policy is chosen, one should make sure that Q is able to increase current level (i.e. nVal could be superior to pVal)
						Q.add( nIdx, Math.min(pVal, nVal) );
						n.setReal(min -1 - pLeaf);
					}
				}
				
				
			} // end loop on neighbor
			
		} // end while
		
		// for root nodes adjust there height to Imax(rootLabel)-min. 
		for(int i=0 ; i<parent.length; i++)
		{
			if( hCriteria[i]>0 & parent[i]==i)
			{
				hCriteria[i] = Imax[i]-min;
			}
		}
		
		//////////////////////////////////////////////////////////////////////////////////
		// final pass on the label image /////////////////////////////////////////////////
		
		// convert the input to label image (label L is stored in input with value min-1-L all other value should be equal to min-1 )
		final IntType minT = labelMap.randomAccess().get().createVariable();
        minT.setReal(min-1);
        final IntType minusOneT = labelMap.randomAccess().get().createVariable();
        minusOneT.setReal(-1);
        Cursor<IntType> input_cursor2 = Views.iterable( labelMap ).cursor();
        while( input_cursor2.hasNext() )
		{
        	IntType p = input_cursor2.next();
        	if (p.getRealFloat()>=(min-1) )
        	{
        		p.setReal(0);
        	}
        	else
        	{
        		p.sub(minT);
            	p.mul(minusOneT);
        	}
		}
		
        // hCriteria and Imax should be rescaled so everything is seemless from outside the class
        for(int i=0 ; i<Imax.length; i++)
		{
        	Imax[i] =  Imax[i] / scaleFactor - offset ;
        	hCriteria[i] = hCriteria[i] / scaleFactor;  
		}
        segmentTree = new Tree(parent, children);
        segmentTree.setFeature("dynamics", hCriteria );
        segmentTree.setFeature("Imax", Imax );
        
        
        this.numberOfLabels = nLeaves;
        
        return;
        // at the end, input was transformed to a label image
        // hCriteria contains the dynamics of each peak
        // parent link nodes to their parent node, if label L has no parent, parent[L]=0 
        // these can be used to build any hMap on the fly.
	}
	
	
	protected static void mergeNodes(int node1, int node2, int newNode, int[] parent, int[][] children)
	{
		
		// update children of parent(node1) if needed
		int par1 = parent[node1];
		if( par1 != node1){ // node1 is not a root
			for(int i=0; i<2; i++){
				if ( children[par1][i] == node1  ){
					children[par1][i] = newNode;
					parent[newNode] = par1;
				}
			}
		}
							
		// update children of parent(node2) if needed
		int par2 = parent[node2];
		if( par2 != node2){ // node2 is not a root
			for(int i=0; i<2; i++){
				if ( children[par2][i] == node2  ){	
					children[par2][i] = newNode;
					parent[newNode] = par2;
				}
			}
		}
		parent[node1] = newNode;
		parent[node2] = newNode;
		children[newNode][0] = node1;
		children[newNode][1] = node2;
		
		return;
	}
	
	
	
	protected static int findRoot(int label, int[] parent)
	{
		// find the root of label 
		int labelParent = parent[label];
		while (labelParent!=label)  
		{
			label = labelParent;
			labelParent=parent[label];
		}
		return label;
	}
	
	


		
	// implementation pseudo code
	//
	// Input: I
	// Output: watershed label, Parent, Hcriteria, Imax,
	// algo will discretize value
	//
	// I_localMin = local_minima(I)
	// Imin = min(I); Imax = max(I)
	// Imin = max(Imin, Thresh)
	// Q = Queue labeled pixel in a priority queue with FIFO policy
	// define out of bound to have value Imin-1
	// initialize parent to -I+Imin if labeled, 0 otherwise // size wise parent will 2*max(I_localMin)-1
	// initialize Hcriteria to 0
	// initialize isDequeued to false
	// while Q.hasNext()
	//
	//	p=Q.next()
	//  isdeQueued[p]=true;
	//
	//	pLeaf= p.getlabel()
	//	pRoot= findRoot(pLeaf)
	//
	//	pVal = p.getVal();
	//	for n in N(p)
	//     if n was already dequeued // meaning that p is the lowest point()
	//      nLeaf = n.getLabel() 
	//		nNode = finRoot(nLabel)
	//		if nNode != pNode // 2 nodes are meeting at that point, create a new node 
	//			maxNode++ // replace maxNode by newNode
	//			Hcriteria(nNode) = Imax(nNode)-pVal  
	//			Hcriteria(pNode) = Imax(pNode)-pVal
	//			
	//			//merge the smaller peak with first neighbor node that has higher dynamics
	//			if H(nNode)>=H(pNode)
	//              node1 = pLeaf
	//				node = nLeaf
	//				while(H(par(node))<=H(node1))
	//					node=parent(node)
	//				merge(node, node1, maxNode)
	//			else 
	//				node1 = nLeaf
	//				node = pLeaf
	//				while(H(parent(node))<=H(node1))
	//					node=parent(node)
	//				merge(node, node1, maxNode)
	//			Imax(maxNode) = max( Imax(node), Imax(node1) )
	//			H[maxNode] = Imax[maxNode]-pVal;
	//
	//		if nVal>Imin // not queued and in-bound
	//			Q.add(n, nVal)
	//			I(n)=pLeaf
	//
	// merge(node1, node2, newNode) // updates parent, children
	// {
	//	parent(newNode) = newNode
	//	
	//	// update children of parent(node1) if needed
	//	p1 = parent[node1]
	//	if p1!=node1 // node1 is note root
	//		for (i=0, i<2; i++)
	//			if children[p1][i] == node1
	//				children[p1][i] = newNode
	//
	//	// update the children of parent(node2) if needed 
	//	p2 = parent[node2]
	//	if ...
	//
	//	parent(node1) = newNode
	//	parent(node2) = newNode
	// }
	//	

	
	public static void main(final String... args) throws IOException
	{
		ImageJ ij = new ImageJ();
		ij.ui().showUI();
		
		//Dataset data = (Dataset) ij.io().open("F:\\projects\\blobs32.tif");
		//Img<FloatType> img = (Img<FloatType>) data.getImgPlus();
		
		ImagePlus imp = IJ.openImage("F:\\projects\\blobs32.tif");
		
		ij.ui().show(imp);
		
		Img<FloatType> img = ImageJFunctions.wrap(imp);
		float threshold = 8f;
		HWatershedHierarchy<FloatType> watershed = new HWatershedHierarchy<FloatType>( img , threshold, ConnectivityHWS.FACE );
		
		RandomAccessibleInterval<IntType> labelMap = watershed.getLabelMap();
		
		
		ij.ui().show(labelMap);
		
		/*
		// image to flood
		new ij.ImageJ();
		//IJ.open("F:\\projects\\blobs.tif");
		IJ.open("F:\\projects\\2D_8peaks.tif");
		ImagePlus imp = IJ.getImage();
		//IJ.run(imp, "Gaussian Blur...", "sigma=2");
		ImagePlusImgConverter impImgConverter = new ImagePlusImgConverter(imp);
		Img<FloatType> input = impImgConverter.getImgFloatType();
		
		float threshold = Float.NEGATIVE_INFINITY;
		HWatershedLabeling<FloatType> maxTreeConstructor = new HWatershedLabeling<FloatType>(input, threshold, Connectivity.FULL);
		
		Img<IntType> output = maxTreeConstructor.getLabelMapMaxTree();
		ImagePlus imp_out = impImgConverter.getImagePlus(output);
		imp_out.show();
		
		//impImgConverter.getImagePlus( maxTreeConstructor.labelMapDebug ).show();
		
		int[] parents = maxTreeConstructor.getTree().getParentsAsArray();
		Map<Integer,Node> treeNodes = maxTreeConstructor.getTree().getNodes();
		double[] dynamics =  maxTreeConstructor.getTree().getFeature("dynamics");
		for( Node node : treeNodes.values())
		{
			int id = node.getId();
			int pId= id;
			if ( node.getParent()!= null )
				pId = node.getParent().getId();
			
			String str = "Id:"+id+"  ;  parent:"+pId+"  ;  children:";
			for(Node nodeC : node.getChildren() )
				str = str+nodeC.getId()+", ";
			str = str+"  ;  dyn:"+dynamics[id];
			System.out.println(str);
		}
		
		//System.out.println("Parents: " + Arrays.toString(parents));
		
		double[] attributes = maxTreeConstructor.getTree().getFeature("dynamics");
		System.out.println("Attributes: " + Arrays.toString(attributes));
		
		
		//maxTreeConstructor.getImgLabeling();
		
		//float[] attributes2 = new float[] {1, 2, 4, 8, 16, 32};
		float[] attributes2 = new float[] {0};
		//Arrays.copyOf(attributes, attributes.length);
		//Arrays.sort(attributes2);
		for( float val : attributes2 )
		{
			
			float h= val;
			Img<IntType> hMax = maxTreeConstructor.getHFlooding2(h);
			ImagePlus imp_hMax = impImgConverter.getImagePlus(hMax);
			imp_hMax.setTitle("hMax (h="+h+")");
			imp_hMax.show();
			IJ.run(imp_hMax, "3-3-2 RGB", "");
			IJ.setMinAndMax(imp_hMax, 0, parents.length);
			
		}
		
		*/
		
		
	}
	

}
