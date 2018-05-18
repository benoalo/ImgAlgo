package invizio.imgalgo.label;



import net.imagej.ImageJ;

/*
Author: Benoit Lombardot, Scientific Computing Facility, MPI-CBG, Dresden  

Copyright 2017 Max Planck Institute of Molecular Cell Biology and Genetics, Dresden, Germany

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following 
conditions are met:

1 - Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2 - Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer 
in the documentation and/or other materials provided with the distribution.

3 - Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived 
from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE 
COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES 
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) 
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/


import net.imglib2.Cursor;
import net.imglib2.IterableInterval;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.integer.IntType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;

import ij.IJ;
import ij.ImagePlus;

import invizio.imgalgo.label.hwatershed.HWatershedHierarchy;
import invizio.imgalgo.label.hwatershed.Tree;
//import invizio.imgalgo.label.hwatershed.TreeUtils;
import invizio.imgalgo.util.RAI;
import invizio.imgalgo.label.hwatershed.HTreeLabeling;

/**
 * 
 * @author Benoit Lombardot
 * 
 */


/*
 * TODO: 
 *  [x] make algorithm class derive from an abstract algorithm class that extends runnable and has a progress getter
 * 	[-] check HWatershed is img factory agnostic
 *  [x] make input random Accessible
 *  [x] check that output are always IntType (should be true for all label algorithm)
 *  [x] bring input image conversion to a Util class (conversion/copy done in mother class)
 *  [-] two separated region can have the same label ==> for a peak i, if Imax(i)-T < hMin the peak should not be in the label map 
 *  
 *  Rk: do hMin and threshold need to be adapted to the rescaling done in segment hierarchy (they should not have to)
 *  
 */



// if pos is updated tree labeling does not change
// if hMin is updated the segmentMap slice is constant but still need to be relabeled. currently we don't keep a copy and have to redo the clicking

// ideally multithread the copy of the segmentMap hyperslice (look how to optimize that)
// as well as the filling of the labelMap 



public class HWatershed <T extends RealType<T> & NativeType<T>> extends DefaultLabelAlgorithm<T>  {

	Tree segmentTree0;
	RandomAccessibleInterval<IntType> segmentMap0;
	
	
	// parameters
	float threshold = 0;
	float hMin = 0;
	float peakFlooding = 100;
	boolean allowSpliting = false;
	
	boolean createSegmentHierarchy = true;
	boolean updateTreeLabeling = true;
	boolean updateLabelMap = true;
	
	// data structure
	int[] nodeIdToLabel;  	// current tree labeling
	int[] nodeIdToLabelRoot;
	double[] peakThresholds;
	HTreeLabeling treeLabeler;
	
	
	 /**
	  *  @param threshold , all pixel below threshold are set to 0 
	  */
	public void setThreshold(float threshold) {
		this.threshold = threshold;
		updateLabelMap = true;
		isProcessed = false;
	}

	/**
	 * @param hMin the minimum  dynamics of a peak not to be merged
	 */
	public void sethMin(float hMin) {
		this.hMin = hMin;
		updateTreeLabeling = true;
		updateLabelMap = true;
		isProcessed = false;
	}

	/**
	 *  @param peakFlooding , percent of the peak that will be flooded  (percent between label maximum and the threshold) 
	 */
	public void setPeakFlooding(float peakFlooding) {
		
		peakFlooding = Math.max( 0 		, peakFlooding );
		peakFlooding = Math.min( 100 	, peakFlooding );
		this.peakFlooding = peakFlooding;
		updateLabelMap = true;
		isProcessed = false;
	}
	
	public void allowSpliting(boolean allowSpliting) {
		
		this.allowSpliting = allowSpliting;
		
		updateTreeLabeling = true;
		updateLabelMap = true;
		isProcessed = false;
	}
	
	
	
	/**
	 * 
	 * @param segmentTree a tree where each leaf correspond to a label in the segment map and each other not to a merging of these region
	 * @param segmentMap0 a labelmap defining non overlapping region
	 * @param intensity0 a graylevel image the same size as segmentMap0
	 */
	public HWatershed(RandomAccessibleInterval<T> input, RandomAccessibleInterval<IntType> segmentMap0, Tree segmentTree ){
		
		super(input);
		
		this.segmentTree0 = segmentTree;
		this.segmentMap0 = segmentMap0;
		createSegmentHierarchy = false;
	}
	
	
	
	/**
	 * build a hierarchy from the input image and makes it available to build HWatershed based on the parameter set ()
	 * @param intensity0 a graylevel image 
	 */
	public HWatershed(RandomAccessibleInterval<T> input ){
		
		super(input);
		createSegmentHierarchy = true;
	}
	
	
	@Override
	protected void process() {
		
		createSegmentHierarchy();
		updateTreeLabeling();
		updateLabelMap();
	}
	
	
	/**
	 *  create the segment Hierarchy (full watershed label map and segment tree )
	 */
	private void createSegmentHierarchy(){
		
		if( ! createSegmentHierarchy ) {	return;		}
		
		T threshold0T = input.randomAccess().get().createVariable();
		threshold0T.setReal( threshold0T.getMinValue() );
		float threshold0 = threshold0T.getRealFloat();
		HWatershedHierarchy<T> SegmentHierarchyBuilder = new HWatershedHierarchy<T>(input , threshold0,  HWatershedHierarchy.ConnectivityHWS.FACE );
		
		segmentMap0 = SegmentHierarchyBuilder.getLabelMap();
		segmentTree0 = SegmentHierarchyBuilder.getSegmentTree();
		
		createSegmentHierarchy = false;
	}
	
	
	
	/**
	 *  update the tree labeling
	 */
	private void updateTreeLabeling(){
		
		if( ! updateTreeLabeling ) {	return;		}
		
		treeLabeler = new HTreeLabeling(segmentTree0);
		
		int nNodes = segmentTree0.getNumNodes();
		nodeIdToLabel = new int[nNodes];
		nodeIdToLabelRoot = new int[nNodes];
		peakThresholds = new double[nNodes];
		this.numberOfLabels = treeLabeler.getLabeling(hMin, threshold, peakFlooding, allowSpliting, nodeIdToLabel, nodeIdToLabelRoot, peakThresholds);
		
		//boolean makeNewLabels = false;
		//nodeIdToLabel =  TreeUtils.getTreeLabeling(segmentTree0, "dynamics", hMin, makeNewLabels );
		
		updateTreeLabeling = false;
	}

	
	
	/**
	 * Relabel the label map according to the tree labeling.
	* @return a label image corresponding to the current tree labeling, threshold, percentFlooding parameters
	 */
	private void updateLabelMap(){
		if( ! updateLabelMap ) {	return;		}
		
		IterableInterval<T> input_iterable = Views.iterable( input );
		labelMap = RAI.duplicate(segmentMap0); //segmentMap0.factory().create(dims, segmentMap0.randomAccess().get().createVariable() );
		fillLabelMap2(labelMap, input_iterable);
		
		updateLabelMap = false;
	}
	
	
	
	@Deprecated
	protected void fillLabelMap( final RandomAccessibleInterval<IntType> segmentMap, final IterableInterval<T> intensity ){
		
		
		double[] Imax = segmentTree0.getFeature("Imax");
		double minI = Imax[0]; // Imax is set to min of the image for the node 0  (it represents the background)
		// create a new (continuous) labeling where only the region visible in the current labeling (and threshold) are numbered
		int nNodes = segmentTree0.getNumNodes();
		int[] nodeIdToLabel2 = new int[nNodes];
		int[] labelled = new int[nNodes];
		int currentLabel = 0;
		for(int i=0; i<nNodes ; i++)
		{
			int parent = nodeIdToLabel[i];
			double ImaxParent = Imax[parent]; 
			if( ImaxParent >= threshold && ImaxParent>minI )
			{
				int parentLabel = labelled[parent];
				if( parentLabel > 0 ) // parent is already labelled, label the node according to its parent
				{
					nodeIdToLabel2[i] = parentLabel;
				}
				else // parent is not labeled yet, create a new label
				{
					currentLabel++;
					labelled[parent] = currentLabel;
					nodeIdToLabel2[i] = currentLabel;
				}
			}
			// else do nothing since the nodeIdToLabel2 is initialized to 0 
		}
		this.numberOfLabels = currentLabel;
		
		
		int nNode = Imax.length;
		float[] peakThresholds = new float[nNode];
		for(int i=0;i<nNode; i++)
			peakThresholds[i] =  threshold + ((float)Imax[i]-threshold)*(1f-peakFlooding/100f);	
		
		Cursor<IntType> cursor = Views.iterable( segmentMap ).cursor();
		Cursor<T> cursorImg = intensity.cursor();
		while( cursor.hasNext() )
		{
			T imgPixel = cursorImg.next();
			float val =imgPixel.getRealFloat();
			
			IntType pixel = cursor.next();
			if(  val >= threshold )
			{
				int node = (int)pixel.getRealFloat();
				int labelRoot = segmentTree0.getNodes().get(node).labelRoot;
				if(  val >= peakThresholds[labelRoot]  )
				{	
					int label = nodeIdToLabel2[node];
					float finalVal = (float)label;
					pixel.setReal( finalVal );
				}
				else
					pixel.setReal( 0.0 );
			}
			else
				pixel.setReal( 0.0 );
		}
		
		return;
	}
	
	
	protected void fillLabelMap2( final RandomAccessibleInterval<IntType> segmentMap, final IterableInterval<T> intensity  ){
		
		
		//int nNodes = segmentTree0.getNumNodes();
		//int[] nodeIdToLabel = new int[nNodes];
		//int[] nodeIdToLabelRoot = new int[nNodes];
		//double[] peakThresholds = new double[nNodes];
		
		//this.nLabels = treeLabeler.getLabeling(hMin, threshold, percentFlooding, keepOrphanPeak, nodeIdToLabel, nodeIdToLabelRoot, peakThresholds);
		
		Cursor<IntType> cursor =  Views.iterable( segmentMap ).cursor();
		Cursor<T> cursorImg = intensity.cursor();
		while( cursor.hasNext() )
		{
			T imgPixel = cursorImg.next();
			float val =imgPixel.getRealFloat();
			
			IntType pixel = cursor.next();
			if(  val >= threshold )
			{
				final int nodeId = (int)pixel.getRealFloat();
				final int labelRoot = nodeIdToLabelRoot[nodeId];
				if(  val >= peakThresholds[labelRoot]  )
				{	
					final int label = nodeIdToLabel[nodeId];
					pixel.setReal( (float)label );
				}
				else
					pixel.setReal( 0.0 );
			}
			else
				pixel.setReal( 0.0 );
		}
		
		return;
		
	}
	
	
	
	
	// code to lazily create single slice of the labelmap; 
	
	
	RandomAccessibleInterval<IntType> labelMapSlice; // current hyperslice
	IterableInterval<T> inputSlice; // current hyperslice
	
	public RandomAccessibleInterval<IntType> getLabelMap( int dim, long pos ) {
		
		process(dim, pos);
	
		return labelMapSlice;
	}
	
	protected void process(int dim, long pos) {
		
		createSegmentHierarchy();
		updateTreeLabeling();
		updateLabelMap( dim, pos );
	}
	
	
	/**
	 *  Relabel a particular slice of the label map according to the tree labeling.	
	 *  @param dim , dimension along which the slice will be cut
	 *  @param pos , position of the slice on long the dimension
	 */
	private void updateLabelMap( int dim, long pos){
		
		int nDims = segmentMap0.numDimensions();
		
		if (nDims>2)
		{	
			labelMapSlice = RAI.duplicate( Views.hyperSlice(segmentMap0, dim, pos) );
			inputSlice = (IterableInterval<T>) Views.hyperSlice(input, dim, pos);
		}
		else{
			labelMapSlice =  RAI.duplicate( segmentMap0 ); //segmentMap0.factory().create(dims, segmentMap0.randomAccess().get().createVariable() );
			inputSlice = Views.iterable( input );
		}
		
		fillLabelMap2( labelMapSlice, inputSlice );
		
	}
	
	
	
	
	public static void main(final String... args)
	{
		
		ImageJ ij = new ImageJ();
		ij.ui().showUI();
		
		ImagePlus imp = IJ.openImage("F:\\projects\\blobs32.tif");
		//ImagePlus imp = IJ.openImage("C:/Users/Ben/workspace/testImages/blobs32.tif");
		ij.ui().show(imp);
		
		
		Img<FloatType> img = ImageJFunctions.wrap(imp);
		float threshold = 100;
		Float hMin = 50f;
		Float peakFlooding = 100f;
	
		HWatershed<FloatType> hWatershed = new HWatershed<FloatType>( img );
		hWatershed.sethMin((float)(double)hMin);
		hWatershed.setThreshold((float)((double)threshold));
		hWatershed.setPeakFlooding((float)(double)peakFlooding);
		
		RandomAccessibleInterval<IntType> labelMap = (Img<IntType>) hWatershed.getLabelMap();
		ij.ui().show( labelMap );
	
	}
	
	
	
	
	
}





