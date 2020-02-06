package mumfordShah;
/**
 * Class for the representation of a point with an left-sided gradient
 * 
 * @author Kilian Hohm
 * 
 */
public class Point {

	private int position;
	private double value;
	private double gradient;

	/**
	 * Constructor for creating a new Point
	 * 
	 * @param pos
	 *            position (x-axis)
	 * @param val
	 *            value (y-axis)
	 * @param grad
	 *            left-sided gradient
	 */
	public Point(int pos, double val, double grad) {
		position = pos;
		value = val;
		gradient = grad;
	}

	/**
	 * @return returns the position of the Point
	 */
	public int getPosition() {
		return position;
	}

	/**
	 * @return returns the Value of the Point
	 */
	public double getValue() {
		return value;
	}

	/**
	 * @return returns the left sided gradient of the Point
	 */
	public double getGradient() {
		return gradient;
	}

	/**
	 * changes the left-sided gradient to a new value
	 * 
	 * @param grad
	 *            new left-sided gradient
	 */
	public void setGradient(double grad) {
		gradient = grad;
	}
	
	public String toString(){
		return "(" + String.valueOf(position) + "," + String.valueOf(value) + "," + String.valueOf(gradient) + ")";
	}

}