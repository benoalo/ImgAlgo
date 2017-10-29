package invizio.imgalgo;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import net.imglib2.RealPoint;
import net.imglib2.Point;




/**
 * N-dimensional point to which arbitrary attributes can be attached by mean of a dictionary  
 * @author Benoit Lombardot
 *
 */

/* TODO:
 * 	[-] make that version depends on RealPoint
 * 	[-] make utils to cast to Point and vice et versa
 * 	[-] use this new point in MSmax
 */

public class Point_A extends Point implements Attributes {

		Map<String,Double> attributes; 
		
		/**
		 * 
		 * @param position an array with n values indicating the position of an n dimensionnal point
		 * @param attributes2 a dictionary associating a string attributes to double values
		 */
		public Point_A(long[] position, Map<String, Double> attributes2)
		{
			super( position );
			this.attributes = attributes2;
		}
		
		
		public Point_A(int[] position, Map<String, Double> attributes)
		{
			super( position );
			this.attributes = attributes;
		}
		
		@Override
		public Map<String,Double> getAttributes() { 
			
			return attributes; 
		}
		
		
		@Override		
		public double getAttribute(String attribute) {
			
			return attributes.get(attribute); 
		}
		
		
		@Override
		public void   setAttribute(String attribute, double value) {
			
			attributes.put( attribute, value); 
		}
		
		public double getDistanceTo(Point_A p)
		{
			double sum = 0;
			for (int d = 0; d < this.numDimensions() ; d ++)
			{
				sum += Math.pow(  ( this.getDoublePosition(d) - p.getDoublePosition(d) )  ,  2  );
			}
			return Math.sqrt(sum);
		}
	
		
		@Override
		public String toString() {
			String str="pos=";
			
			for( int d=0;d<this.numDimensions() ; d++)
				str += this.position[d] +" , ";
			
			str += "  ;  ";
			for(Entry<String,Double> item : attributes.entrySet() )
				str += item.getKey() +"="+ item.getValue() + " , ";
			
			return str;
		}
		
		
		
	public static List<Point> toPoint(List<RealPoint> realPoints) {
		
		List<Point> points = new ArrayList<Point>();
		int nDim = points.get(0).numDimensions();
		for( RealPoint realPoint : realPoints ) {
			Point point = new Point( nDim );
			for(int d=0; d<nDim; d++) {
				final long position = (int) realPoint.getDoublePosition(d);
				point.setPosition(position, d);
			}
			points.add( point );
		}
		
		return points;
	}
	
	/**
	 * TODO
	 * @param points a list of Point
	 * @return a list of RealPoint
	 */
	public static List<RealPoint> toRealPoint(List<Point> points) {
		List<RealPoint> realPoints = new ArrayList<RealPoint>();

		return realPoints;
	}


	
}
