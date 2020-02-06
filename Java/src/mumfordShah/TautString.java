package mumfordShah;
import java.util.LinkedList;
import java.util.ListIterator;

/**
 * Java class for calculating the tautstring
 * 
 * @author Kilian Hohm
 * 
 */

public class TautString {
	private double[] F; // cumsum
	private double[] s; // quadcumsum
	private double gamma;
	private int n;
	private int iter;
	private LinkedList<Point> sx;
	private LinkedList<Point> sv;
	private LinkedList<Point> compResult;
	private LinkedList<Point> tempResult;
	double tvSolution = 0; // isn't calcualated with variable gamma (always gamma[0] is used)
	double tempTvSolution = 0;
	

	/**
	 * Constructor
	 * 
	 * @param F
	 *            cumulative sum of the data which defines the shape of the tube
	 * @param gamma
	 *            parameter which defines the width of the tube
	 */
	public TautString(double[] F, double[] f, double gamma) {
		this.F = F;
		this.s = quadCumSum(f);
		this.gamma = gamma;
		n = F.length;

		sx = new LinkedList<Point>();
		sx.addLast(new Point(0, F[0], Double.NEGATIVE_INFINITY));

		sv = new LinkedList<Point>();
		sv.addLast(new Point(0, F[0], Double.POSITIVE_INFINITY));
		
		iter = 1;

		compResult = new LinkedList<Point>();
		tempResult = new LinkedList<Point>();
	}

	/**
	 * Calculate the tautstring only for the complete interval F[0],...,F[n-1]
	 */
	public LinkedList<Point> run() { // doesn't calculate for every interval only for the
						// complete

		for (int i = 1; i < n; i++) {
			processNext(i, n, 0);
		}

		// If there isn't any point in result we have to add (0,F[0])
		if (compResult.isEmpty())
			compResult.add(new Point(0, F[0], 0));

		// finally after all points are processed add the last point F(n) to the
		// linked list results
		Point p = compResult.getLast();
		double gradient = (F[n-1] - p.getValue())/(n-1 - p.getPosition());
		compResult.addLast(new Point(n - 1, F[n - 1], gradient));
		
		return compResult;
	}
	
	public LinkedList<Point> runNext(){
		if(iter < n){
			LinkedList<Point> tempSV = cloneList(sv);
			LinkedList<Point> tempSX = cloneList(sx);
			
			
			processNext(iter, iter, 1); // 
			
			// If there isn't any point in result we have to add (0,F[0])
			if (tempResult.isEmpty())
				tempResult.add(new Point(0, F[0], 0));

			// finally after all points are processed add the last point F(n) to
			// the linked list results
			Point resultLast = tempResult.getLast();
			double grad = (F[iter] - resultLast.getValue())
					/ (iter - resultLast.getPosition());
			tempResult.addLast(new Point(iter, F[iter], grad));

			sv = tempSV;
			sx = tempSX;
			
			LinkedList<Point> tempResult = cloneList(this.tempResult);
			this.tempResult = cloneList(compResult);

			processNext(iter, n, 1);
			iter++;

			return tempResult;
		}
		else if(iter == n){
			// If there isn't any point in result we have to add (0,F[0])
			if (tempResult.isEmpty())
				tempResult.add(new Point(0, F[0], 0));

			// finally after all points are processed add the last point F(n) to the
			// linked list results
			Point resultLast = tempResult.getLast();
			double grad = (F[n - 1] - resultLast.getValue())
					/ (n - 1 - resultLast.getPosition());
			tempResult.addLast(new Point(n - 1, F[n - 1], grad));

			
			
			if (compResult.isEmpty())
				compResult.add(new Point(0, F[0], 0));

			// finally after all points are processed add the last point F(n) to the
			// linked list results
			resultLast = compResult.getLast();
			grad = (F[n - 1] - resultLast.getValue())
					/ (n - 1 - resultLast.getPosition());
			compResult.addLast(new Point(n - 1, F[n - 1], grad));
			
			iter++;
			return tempResult;
			
		}
		else{
			throw new IndexOutOfBoundsException("There isn't any new point to process!");
		}
	}
	
	public double runNextGetValue(){
		if(iter < n){
			LinkedList<Point> tempSV = cloneList(sv);
			LinkedList<Point> tempSX = cloneList(sx);
			
			
			processNext(iter, iter, 1); // 
			
			// If there isn't any point in result we have to add (0,F[0])
			if (tempResult.isEmpty())
				tempResult.add(new Point(0, F[0], 0));

			// finally after all points are processed add the last point F(n) to
			// the linked list results
			
			
			Point resultLast = tempResult.getLast();
			double grad = (F[iter] - resultLast.getValue())
					/ (iter - resultLast.getPosition());
			tempResult.addLast(new Point(iter, F[iter], grad));
			
			
			
			// calculate tempTVSolution (if changes they have to be made at 4 positions in the code)
			double tempTVSolution = 0;
			if(resultLast.getPosition()>=0){
				int firstPosition = -1;
				if(tempResult.size() == 1){
					firstPosition = resultLast.getPosition()+1;
				}
				else if(tempResult.size() > 1){
					firstPosition = resultLast.getPosition()+1;
				}
					
//				for(int l = firstPosition; l <= iter; l++)
//					tempTVSolution = tempTVSolution + Math.pow(f[l-1] - grad, 2) ; // (f(l) - tautLösung(l))^2 with tautLösung(l) = -grad
				tempTVSolution += (iter-firstPosition+1)*Math.pow(grad, 2) - 2*grad*(F[iter]-F[firstPosition-1]) + (s[iter]-s[firstPosition-1]);
				if(tempResult.size() > 2)
					tempTVSolution = tempTVSolution + gamma*Math.abs(grad-resultLast.getGradient());
			}
			
			tempTVSolution += tempTvSolution;

			tempTvSolution = tvSolution;

			sv = tempSV;
			sx = tempSX;
			
			this.tempResult = cloneList(compResult);

			processNext(iter, n, 1);
			iter++;

			return tempTVSolution;
		}
		else if(iter == n){
			// If there isn't any point in result we have to add (0,F[0])
			if (tempResult.isEmpty())
				tempResult.add(new Point(0, F[0], 0));

			// finally after all points are processed add the last point F(n) to the
			// linked list results
			Point resultLast = tempResult.getLast();
			double grad = (F[n - 1] - resultLast.getValue())
					/ (n - 1 - resultLast.getPosition());
			tempResult.addLast(new Point(n - 1, F[n - 1], grad));




			if (compResult.isEmpty())
				compResult.add(new Point(0, F[0], 0));

			// finally after all points are processed add the last point F(n) to the
			// linked list results
			resultLast = compResult.getLast();
			grad = (F[n - 1] - resultLast.getValue())
					/ (n - 1 - resultLast.getPosition());
			compResult.addLast(new Point(n - 1, F[n - 1], grad));

			// calculate tempTVSolution
			double tempTVSolution = 0;
			if(resultLast.getPosition()>=0){
				int firstPosition = -1;
				if(compResult.size() == 1){
					firstPosition = resultLast.getPosition()+1;
				}
				else if(compResult.size() > 1){
					firstPosition = resultLast.getPosition()+1;
				}
//				for(int l = firstPosition; l <= iter; l++)
//					tempTVSolution = tempTVSolution + Math.pow(f[l-1] - grad, 2) ; // (f(l) - tautLösung(l))^2 with tautLösung(l) = -grad
				tempTVSolution += (iter-firstPosition+1)*Math.pow(grad, 2) - 2*grad*(F[iter]-F[firstPosition-1]) + (s[iter]-s[firstPosition-1]);
				if(compResult.size() > 2)
					tempTVSolution = tempTVSolution + gamma*Math.abs(grad-resultLast.getGradient());
			}

			iter++;
			return tvSolution + tempTVSolution;
			
		}
		else{
			throw new IndexOutOfBoundsException("There isn't any new point to process!");
		}
	}

	/*
	 * (non-Javadoc) Print what is left in sx and sv, and all results of the
	 * intervals
	 * 
	 * @see java.lang.Object#toString()
	 */
	public String toString() {
		String resultString = "SX: ";
	
		for (int i = 0; i < sx.size(); i++) {
			resultString = resultString + "(" + sx.get(i).getPosition() + ","
					+ sx.get(i).getValue() + ") ";
		}
	
		resultString = resultString + "\nSV: ";
	
		for (int i = 0; i < sv.size(); i++) {
			resultString = resultString + "(" + sv.get(i).getPosition() + ","
					+ sv.get(i).getValue() + ") ";
		}
	
		resultString = resultString + "\nresult: \n";
	
		for (int j = 0; j < compResult.size(); j++) {
				resultString = resultString + "("
						+ compResult.get(j).getPosition() + ","
						+ compResult.get(j).getValue() + ") ";
			resultString = resultString + "\n";
		}
	
		return resultString;
	}

	/**
	 * This function clones the passed LinkedList of points. There is a new
	 * instance for each point of the old LinkedList, so if you change a point
	 * in the new list it won't be changed in the old one
	 * 
	 * @param temp
	 *            LinkedList<Point> that will be cloned
	 * @return a completely new clone of temp
	 */
	private LinkedList<Point> cloneList(LinkedList<Point> temp) {
		LinkedList<Point> res = new LinkedList<Point>();
		ListIterator<Point> iter = temp.listIterator();
		while (iter.hasNext()) {
			Point curr = iter.next();
			Point newP = new Point(curr.getPosition(), curr.getValue(),
					curr.getGradient());
			res.addLast(newP);
		}

		return res;
	}

	/**
	 * This functions executes the tasks which has to be done in the main loop
	 * of {@link tautString.run()} or {@link tautString.runIntervals()}
	 * 
	 * @param i
	 *            current position
	 * @param j
	 *            total length of the intervall
	 * @param k
	 *            use result.get(k) to save the resultstring
	 */
	private void processNext(int i, int j, int k) { // add the next point to sx
													// and sv and
													// testForOverlaps
		calcSX(i, j); // Calculate convex minorant sx_i
		calcSV(i, j); // Calculate concave majorant sv_i
		testForOverlap(k, j); // Test and process overlaps (argument:
								// result.get(0) will be used to enter in
								// solution)
	}

	/**
	 * Calculates the convex hull of the upper bound by adding the next point
	 * F[i]
	 * 
	 * @param i
	 *            current position
	 * @param j
	 *            total length
	 */
	private void calcSX(int i, int j) {
		// create new point gradient isn't important yet (will be changed later)
		Point curr = null;
		if (i < j - 1)
			curr = new Point(i, F[i] + gamma, 0);
		else
			curr = new Point(i, F[i], 0);

		while (sx.size() >= 1) {
			// get last point of sx (previous point of curr)
			Point prev = sx.getLast();

			// calculate the current gradient by building the difference
			// quotient of the current and previous point
			double currGrad = (curr.getValue() - prev.getValue())
					/ (curr.getPosition() - prev.getPosition());

			// Test if current gradient is bigger than previous gradient (then
			// sx is convex)
			if (currGrad > prev.getGradient()) {
				curr.setGradient(currGrad); // Set gradient of current point to
											// currGrad
				break;
			}

			// If sx with the current point is no longer convex delete last
			// point in sx
			sx.removeLast();
		}

		sx.addLast(curr);

	}

	/**
	 * Calculates the concave hull of the lower bound by adding the next point
	 * F[i]
	 * 
	 * @param i
	 *            current position
	 * @param j
	 *            total length
	 */
	private void calcSV(int i, int j) {
		// create new point gradient isn't important yet (will be changed later)
		Point curr = null;
		if (i < j - 1)
			curr = new Point(i, F[i] - gamma, 0);
		else
			curr = new Point(i, F[i], 0);

		while (sv.size() >= 1) {
			// get last point of sx (previous point of curr)
			Point prev = sv.getLast();

			// calculate the current gradient by building the difference
			// quotient of the current and previous point
			double currGrad = (curr.getValue() - prev.getValue())
					/ (curr.getPosition() - prev.getPosition());

			// Test if current gradient is smaller than previous gradient (then
			// sv is concave)
			if (currGrad < prev.getGradient()) {
				curr.setGradient(currGrad); // Set gradient of current point to
											// currGrad
				break;
			}

			// If sx with the current point is no longer convex delete last
			// point in sx
			sv.removeLast();
		}
		sv.addLast(curr);

	}

	/**
	 * This functions test if after calculating sx and sv there is an overlap.
	 * If there is it will be processing from left to right sx/sv and saving the
	 * results in the specified list until ther won't be an overlap.
	 * 
	 * @param i
	 *            use result.get(k) to save the resultstring
	 * @param j
	 *            total length of the intervall
	 */
	private void testForOverlap(int i, int j) {
		// If one of the lists (or both) only contain one element it cannot be
		// tested for overlaps
		if (sx.size() <= 1 || sv.size() <= 1)
			return;

		// Point in sx at position 1 (important gradient information contained)
		Point up = sx.get(1); // there are at least two elements in sx
		double gradUp = up.getGradient(); // sx'(0++)
		double valUp = up.getValue();
		int posUp = up.getPosition();

		// Point in sv at position 1 (important gradient information contained)
		Point down = sv.get(1); // there are at least two elements in sv
		double gradDown = down.getGradient(); // sv'(0++)
		double valDown = down.getValue();
		int posDown = down.getPosition();

		// As long as there is an overlap concerning the first segment of sx and
		// sv fix it
		if (gradDown >= gradUp) { // sv'(0++) > sx'(0++)

			Point sxFirst = sx.getFirst();
			Point resultLast = null;
			
			double tempTVSolution = 0;

			if (compResult.isEmpty()) {
				sxFirst.setGradient(0);

			} else {
				resultLast = compResult.getLast();
				double grad = (sxFirst.getValue() - resultLast.getValue())
						/ (sxFirst.getPosition() - resultLast.getPosition());
				
				sxFirst.setGradient(grad);
			}
			
			
			
			tempTVSolution = 0;
			boolean tester = false;

			if(i == 1){
				if (tempResult.isEmpty()) {
					sxFirst.setGradient(0);

				} else {
					resultLast = tempResult.getLast();
					double grad = (sxFirst.getValue() - resultLast.getValue())
							/ (sxFirst.getPosition() - resultLast.getPosition());
					
					sxFirst.setGradient(grad);
				}
				if(tempResult.isEmpty())
					resultLast = sxFirst;
				else
					resultLast = tempResult.getLast();
				
				// to calculate TV solution use -grad because data is inserted backwards
				int firstPosition = resultLast.getPosition()+1;
				
				tempResult.addLast(sxFirst);
				
//				for(int l = firstPosition; l <=sxFirst.getPosition(); l++)
//					tempTVSolution = tempTVSolution + Math.pow(f[l-1] - sxFirst.getGradient(), 2) ; // (f(l) - tautLösung(l))^2 with tautLösung(l) = -grad
				tempTVSolution += (sxFirst.getPosition()-firstPosition+1)*Math.pow(sxFirst.getGradient(), 2) - 2*sxFirst.getGradient()*(F[sxFirst.getPosition()]-F[firstPosition-1]) + (s[sxFirst.getPosition()]-s[firstPosition-1]);
				if(tempResult.size() > 2)
					tempTVSolution = tempTVSolution + gamma*Math.abs(sxFirst.getGradient()-resultLast.getGradient());

				

				this.tempTvSolution += tempTVSolution;
				
				tester = true;
			}

			tempTVSolution = 0;
			if (j == n){
				if(compResult.isEmpty())
					resultLast = sxFirst;
				else
					resultLast = compResult.getLast();
				
				// to calculate TV solution use -grad because data is inserted backwards
				int firstPosition = resultLast.getPosition()+1;
				
				compResult.addLast(sxFirst);
				
//				for(int l = firstPosition; l <=sxFirst.getPosition(); l++)
//					tempTVSolution = tempTVSolution + Math.pow(f[l-1] - sxFirst.getGradient(), 2) ; // (f(l) - tautLösung(l))^2 with tautLösung(l) = -grad
				tempTVSolution += (sxFirst.getPosition()-firstPosition+1)*Math.pow(sxFirst.getGradient(), 2) - 2*sxFirst.getGradient()*(F[sxFirst.getPosition()]-F[firstPosition-1]) + (s[sxFirst.getPosition()]-s[firstPosition-1]);
				if(compResult.size() > 2)
					tempTVSolution = tempTVSolution + gamma*Math.abs(sxFirst.getGradient()-resultLast.getGradient());

				

				
				tvSolution = tvSolution + tempTVSolution;
				
				if(!tester){
					this.tempTvSolution += tempTVSolution;
				}
					
			}

			// remove first element in sx and sv because they are no longer
			// needed
			sx.removeFirst();
			sv.removeFirst();

			if (posUp < posDown) { // add Point up (now first in sx) to sv on
									// first position and calculate the new
									// gradient for down (second position in sv)
				Point newPoint = new Point(up.getPosition(), up.getValue(),
						Double.POSITIVE_INFINITY);
				sv.addFirst(newPoint);
				// change gradient sx'(0++) to -inf so it will be convex (for
				// the next calculation of sxi)
				sx.getFirst().setGradient(Double.NEGATIVE_INFINITY);

				// calcualate new Gradient (sv'(0++))
				double newGrad = (valDown - valUp) / (posDown - posUp);

				// update gradient of down
				down.setGradient(newGrad);

			} else if (posUp > posDown) {
				Point newPoint = new Point(down.getPosition(), down.getValue(),
						Double.NEGATIVE_INFINITY);
				sx.addFirst(newPoint);
				// change gradient sv'(0++) to +inf so it will be concave (for
				// the next calculation of svi)
				sv.getFirst().setGradient(Double.POSITIVE_INFINITY);
				// calculate new gradient (sx'(0++))
				double newGrad = (valUp - valDown) / (posUp - posDown);

				// update gradient of up
				up.setGradient(newGrad);
			} else { // is it possible that posUp == posDown? yes? only the last
						// point is in it

			}

			testForOverlap(i, j); // recursive Overlap testing
		}
	}

	/**
	 * Matlab connection to pass the result of the complete interval
	 * 
	 * @return returns the result for the complete interval F[0],...,F[n-1]
	 */
	public double[][] getCompResult() {
		double[][] res = new double[compResult.size()][3];
		ListIterator<Point> iter = compResult.listIterator();
		for (int j = 0; j < compResult.size(); j++) {
			Point curr = iter.next();
			res[j][1] = curr.getValue();
			res[j][0] = curr.getPosition();
			res[j][2] = curr.getGradient();
		}
		return res;
	}
	
	public double[][] getTempResult() {
		double[][] res = new double[tempResult.size()][3];
		ListIterator<Point> iter = tempResult.listIterator();
		for (int j = 0; j < tempResult.size(); j++) {
			Point curr = iter.next();
			res[j][1] = curr.getValue();
			res[j][0] = curr.getPosition();
			res[j][2] = curr.getGradient();
		}
		return res;
	}
	
	private double[] quadCumSum(double[] data){
		double[] res = new double[data.length+1];
		res[0] = 0;
		for(int i = 1; i < data.length+1; i++){
			res[i] = res[i-1] + data[i-1]*data[i-1];
		}
		return res;
	}
}


