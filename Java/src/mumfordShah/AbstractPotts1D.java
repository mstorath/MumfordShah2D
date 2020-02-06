package mumfordShah;
import java.util.concurrent.Callable;


public abstract class AbstractPotts1D implements Callable<double[][]> {

	// MEMEBERS
	protected double gamma; // Potts constant
	protected int n; // n is the number of columns of the data field f
	protected int channels; // amount of color channels
	private int[] p; // stores the partition found by findBestPartition
	private double[][] result; // stores the result found by the algorithm (first dimension is the color channel, second dimension the data)
	double[][] f; // stores the original data
	private int jumpsAmount = -1;
	
	
	
	// CONSTRUCTORS
	public AbstractPotts1D(){
		this.gamma = 1;
	}
	
	public AbstractPotts1D(double[] f, double gamma){
		double[][] tempf = new double[1][1];;
		tempf[0] = f;
		set(tempf,gamma);
	}
	
	public AbstractPotts1D(double[][] f, double gamma){
		set(f,gamma);
	}
	
	// CONSTRUCTOR METHOD
	protected void set(double[][] f, double gamma){
		n = f[0].length;
		channels = f.length;
		
		// store data
		this.f = f;
		
		// create gamma vector
		this.gamma = gamma;
		
		p = new int[n+1];
		
		// create result (initially the 0-vector is the result)
		result = new double[channels][n];
	}
	
	
	
	// METHODS
	public void run(){
		// for loop over all rows of the data field f
			findBestPartition();
			segmentationFromPartition();
	}
	
	public double[][] getResult(){
		return result;
	}
	
	public double[][] call(){
		run();
		return result;
	}

	protected void findBestPartition(){
		double [] B = new double[n+1];
		double b, dlr;
		B[0] = -gamma;
		
		
		for(int r = 1; r < n+1; r++){
			B[r] = computeDlr(1,r);
			
			for(int l = r; l >= 2; l--){
				dlr = computeDlr(l,r);
				
				if(B[r] < gamma + dlr) // abortion criteria (acceleration)
					break;
				
				b = B[l-1] + gamma + dlr;
				if(b <= B[r]){
					B[r] = b;
					p[r] = l-1;
				}
			}
		}
		
	}
	
	protected void segmentationFromPartition(){
		int r = n;
		int l = p[r];
		
		while(r > 0){
			double[][] mu = computeMu_LtoR(l+1,r); // compute solution on data l+1,...,r
			for(int dim = 0; dim < channels; dim++){
				for(int t = 0; t < r-l; t++){
					result[dim][l+t] = mu[dim][t];
				}
			}
			//System.out.println(l);
			r = l;
			l = p[r];
			jumpsAmount++;
		}
	}
	
	protected double[] cumSum(double[] data, int n){
		double[] res = new double[n+1];
		res[0] = 0;
		for(int i = 1; i < n+1; i++){
			res[i] = res[i-1] + data[i-1];
		}
		return res;
	}
	
	public int getJumpsAmount(){
		return jumpsAmount;
	}
	
	public int[] getPartition(){
		return p;
	}

	
	protected abstract double computeDlr(int l, int r);

	protected abstract double[][] computeMu_LtoR(int l, int r);
	
	
}
