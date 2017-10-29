package invizio.imgalgo;

/**
 * @author Benoit Lombardot
 */

import java.util.Map;

public interface Attributes{

	public Map<String,Double> getAttributes();
	
	public double getAttribute(String attribute);
	
	public void setAttribute(String attribute, double value);
	
}
