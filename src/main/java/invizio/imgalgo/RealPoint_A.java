package invizio.imgalgo;



import java.util.Map;
import net.imglib2.RealPoint;


/**
 * N-dimensional point to which arbitrary attributes can be attached by mean of a dictionary  
 * @author Benoit Lombardot Benoit
 * 
 */

/* TODO:
 * 	[-] make utils to cast to Point and vice et versa
 */

public class RealPoint_A extends RealPoint  implements Attributes {

		Map<String,Double> attributes; 
		
		/**
		 * 
		 * @param position an array with n values indicating the position of an n dimensionnal point
		 * @param attributes a dictionary associating a string attributes to double values
		 */
		public RealPoint_A(double[] position, Map<String, Double> attributes)
		{
			super( position );
			this.attributes = attributes;
		}
		
		
		public RealPoint_A(float[] position, Map<String, Double> attributes)
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
			
			attributes.put(attribute, value);
		}
		
		
		public double getDistanceTo(RealPoint_A p)
		{
			double sum = 0;
			for (int d = 0; d < this.numDimensions() ; d ++)
			{
				sum += Math.pow(  ( this.getDoublePosition(d) - p.getDoublePosition(d) )  ,  2  );
			}
			return Math.sqrt(sum);
		}
	
		
		
	
}
