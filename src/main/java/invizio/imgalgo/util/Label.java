package invizio.imgalgo.util;


import java.util.List;


import net.imglib2.Localizable;

import net.imglib2.RandomAccess;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgs;
import net.imglib2.type.numeric.integer.IntType;

/**
 * 
 * @author Benoit Lombardot
 * 
 */


// Utility for labelMaps
public class Label {
	
	
	/*
	 * TODO:
	 * 	[-] if ? implements Point_A and has a scale attribute one can sort the point per scale and print disk corresponding to scale size 
	 */
	public static <L extends Localizable> Img<IntType> toLabelMap(long[] dims, List<L> points) {
		
		Img<IntType> resultImg = ArrayImgs.ints(dims);
		RandomAccess<IntType> ra = resultImg.randomAccess();
	
		int label = 1;
		for (Localizable point : points) {
			ra.setPosition( point );
			ra.get().set( label );
			label++;
		}
	
		return resultImg;
	}
	
		
}
