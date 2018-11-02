package invizio.imgalgo.label.rleccl;

/**
 * 
 * @author Benoit Lombardot
 * 
 */


public class ValuedPixelRun extends PixelRun {

	public float value=0;
	public boolean valid=true;
	public boolean plato=false;
	public boolean postplato=false;
	
	public ValuedPixelRun(int start, int end, float value) {		
		super(start, end);
		this.value = value;
	}	
}
