package mumfordShah;

import java.util.ArrayList;

public abstract class GaussElim {
	protected int n;
	protected double alpha;
	protected ArrayList<double[][]> diagonals;
	protected ArrayList<double[][]> factors;
	
	public GaussElim(int n, double alpha) {
		this.n = n;
		this.alpha = alpha;
		initialize();
	}
	
	protected abstract void initialize();
	
	public abstract double[] computeMu(double[] y);
	
	public double computeDlr(double[] y){
		double result = 0;
		
		double[] mu = computeMu(y);
		
		for(int idx = 0; idx < mu.length - 1; idx++)
			result += Math.pow((y[idx]-mu[idx]),2) + alpha * Math.pow((mu[idx+1]-mu[idx]),2);
		
		result += Math.pow((y[mu.length-1]-mu[mu.length-1]),2);
		
		return result;
	}
}
