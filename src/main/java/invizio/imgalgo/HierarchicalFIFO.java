package invizio.imgalgo;


import java.util.ArrayDeque;
import java.util.ArrayList;


/**
 * @author Benoit Lombardot
 */


public class HierarchicalFIFO {
	
	
	private int current_level;
	public int getCurrent_level() {
		return current_level;
	}

	float current_value;
	public float getCurrent_value() {
		return current_value;
	}
	
	private void update_currentValue(){
		current_value = min + (max-min)*current_level/(nbin-1);
	}
	
	public float getMin() {
		return min;
	}

	private final float min, max;
	private int max_level, nbin;
	
	private ArrayList<ArrayDeque<Item>> QueueList; // could use a treeMap to sort it with FloatType ( comput cost ? )
	
	
	public HierarchicalFIFO(float min, float max, int nbin)
	{
		QueueList = new ArrayList<ArrayDeque<Item>>(nbin);
		for(int i=0; i<nbin; i++)
			QueueList.add( new ArrayDeque<Item>() );
		this.min = min;
		this.max = max;
		this.max_level = nbin-1;
		this.nbin = nbin;
		current_level = max_level;
		update_currentValue();
	}
	
	public HierarchicalFIFO(float min, float max)
	{
		int nbin = (int)max - (int)min +1;
		QueueList = new ArrayList<ArrayDeque<Item>>(nbin);
		for(int i=0; i<nbin; i++)
			QueueList.add( new ArrayDeque<Item>() );
		this.min = min;
		this.max = max;
		this.max_level = nbin-1;
		this.nbin = nbin;
		current_level = max_level;
		update_currentValue();
	}
	
	
	public void add(long idx, float val)
	{
		int level = (int)( (val - min) / (max - min) * (nbin-1) ) ;
		QueueList.get(level).add(  new Item( val, idx )  );
		//current_level = Math.max(current_level,level);
	}
	
	
	public boolean HasNext()
	{
		while( QueueList.get(current_level).isEmpty() & current_level>0){
			current_level--;
			update_currentValue();
		}
		
		return !QueueList.get(current_level).isEmpty();
	}
	
	
	public Item Next()
	{	
		return QueueList.get(current_level).poll();	
	}


	public class Item
	{
		float value;
		long index;
		
		public Item(float value, long index){
			this.value = value;
			this.index = index;
		}
		
		public float getValue(){
			return value;
		}
		
		public long getIndex(){
			return index;
		}
		
	}
	
}