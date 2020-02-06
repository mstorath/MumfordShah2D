package mumfordShah;
import java.util.LinkedList;

public class L1MumfordShah1D extends AbstractPotts1D {
	
	/*
	 * 		Computes the solution to the following 1D minimization problem
	 * 			
	 * 			min	
	 */
	
	double alpha;
	private double d1r, dlr;
	private int curL,curR;
	private TautString[] tautD1r, tautDlr;
	
	public L1MumfordShah1D(){
		double[][] f = new double[5][5];
		double gamma = 1;
		double alpha = 1;
		set(f,gamma,alpha);
	}
	
	public L1MumfordShah1D(double[] f, double gamma, double alpha){
		double[][] tempf = new double[1][1];
		tempf[0] = f;
		set(tempf,gamma,alpha);
	}
	
	public L1MumfordShah1D(double[][] f, double gamma, double alpha){
		set(f,gamma,alpha);
	}
	
	public void set(double[][] f,double gamma,double alpha){
		set(f,gamma);
		
		d1r = 0; dlr = 0; curL = 0; curR = 0;
		
		tautD1r = new TautString[channels];
		tautDlr = new TautString[channels];

		this.alpha = alpha;
		
	}

	@Override
	protected double[][] computeMu_LtoR(int l, int r) {
		double[][] tempf = new double[channels][r-l+1];
		double[][] tempCumsum = new double[channels][0];
		for(int dim = 0; dim < channels; dim++){
			for(int i = l-1; i < r; i++){
				tempf[dim][i-(l-1)] = f[dim][i];
			}
			tempCumsum[dim] = cumSum(tempf[dim],tempf[dim].length);
		}
		return getSolutionWithTaut(tempCumsum,tempf,alpha); //tempCumsum replace by Arrays.copyOfRange(m[row], l, r+1)
	}
	
	private double[][] getSolutionWithTaut(double[][] cumSum, double[][]tempf, double alpha){
		double[][] temp = new double[channels][cumSum[0].length-1];
		for(int dim = 0; dim < channels; dim++){
			TautString taut = new TautString(cumSum[dim],tempf[dim], alpha);
			LinkedList<Point> result = taut.run();

			result.removeFirst(); // Remove first element because it has no information

			int j = 0; // has to be outside of the for declaration or it will be reset to zero every time
					
			while(!result.isEmpty()){
				// one loop represents a constant piece of the result function
				Point p = result.pop(); // get and remove first element in list
				for(; j < p.getPosition() && j < temp[0].length; j++)
					temp[dim][j] = p.getGradient();
			}
		}
		return temp;
	}
	
	@Override
	protected double computeDlr(int l, int r) {
		// initialize tautD1r's and compute d11
		if(r == 1 && l == 1){
			d1r = 0;
			for(int dim = 0; dim < channels; dim++){
				double[] F = cumSum(f[dim], n);
				tautD1r[dim] = new TautString(F,f[dim],alpha);
				d1r += tautD1r[dim].runNextGetValue();
			}
			curR = r; curL = l;
		}
		// compute d1r
		else if (r == curR+1 && l == 1){
			d1r = 0;
			for(int dim = 0; dim < channels; dim++)
				d1r += tautD1r[dim].runNextGetValue();
			curR = r; curL = l;
		}
		// initialize tautDlr's and compute dr-1r
		else if(r == curR && l == curR){
			dlr = 0;
			
			// prepare data for "backwards" computation (remark that TV(data) == TV(dataBackwards))
			double[][] data = new double[channels][r-1];
			for(int dim = 0; dim < channels; dim++) for(int iter = 0; iter < r-1; iter++)
				data[dim][iter] = f[dim][r-1-iter];
			
			// initialize taut strings and compute the value
			for(int dim = 0; dim < channels; dim++){
				double[] Data = cumSum(data[dim],r-1);
				tautDlr[dim] = new TautString(Data,data[dim],alpha);
				dlr += tautDlr[dim].runNextGetValue();
			}
			curR = r; curL = l;
		}
		// compute dlr
		else if (r == curR && l == curL-1){
			dlr = 0;
			
			for(int dim = 0; dim < channels; dim++)
				dlr += tautDlr[dim].runNextGetValue();
			
			curR = r; curL = l;
		}
		// we don't have to compute something new
		else if(r == curR && l == curL-1);
		// else throw exception
		else
			throw(new IllegalArgumentException("l has to be curL-1"));
			

		// return the computed value
		if(l == 1)
			return d1r;
		else
			return dlr;
	}

}
