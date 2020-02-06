package mumfordShah;

import java.util.ArrayList;

public class GaussL2Mum extends GaussElim {
	private double a,b,c;
	
	public GaussL2Mum(int n, double alpha){
		super(n,alpha);
		
		diagonals = new ArrayList<double[][]>(n);
		factors = new ArrayList<double[][]>(n);
		
		initialize();
	}
	
	public double[][] getDiagonals(int i){
		return diagonals.get(i);
	}
	
	protected void initialize(){
		a = alpha+1;
		b = -alpha;
		c = 2*alpha+1;
		
		diagonals = new ArrayList<double[][]>(n);
		factors = new ArrayList<double[][]>(n);
		
		double[][] D0 = {{1}};
		double[][] C0 = new double[0][0];
		diagonals.add(0,D0);
		factors.add(0,C0);
		
		double[][] D1 = {{a,a-b*b/a}};
		double[][] C1 = {{-b/a}};
		diagonals.add(1,D1);
		factors.add(1,C1);
		
		double[][] D,C;
		for(int idx = 2; idx<n; idx++){
			D = new double[1][idx+1];
			C = new double[1][idx];
			
			D[0][0] = a;
			
			for(int curIdx = 1; curIdx < idx; curIdx++){
				double multi = -b/D[0][curIdx-1];
				C[0][curIdx-1] = multi;
				D[0][curIdx] = c+b*multi;
			}

			double multi = -b/D[0][idx-1];
			C[0][idx-1] = multi;
			D[0][idx] = a+b*multi;
			
			diagonals.add(idx,D);
			factors.add(idx,C);
		}
	}
	
	
	public double[] computeMu(double[] y){
		double[][] factor,diagonal;
		double[] result,b;
		int len = y.length;
		if(len == 1) return new double[] {y[0]};
		diagonal = diagonals.get(len-1);
		factor = factors.get(len-1);
		
		result = new double[len];
		b = new double[len];
		
		b[0] = y[0];
		for(int idx = 1; idx < len; idx++)
			b[idx] = y[idx] + factor[0][idx-1]*b[idx-1];
		
		result[len-1] = b[len-1]/diagonal[0][len-1];
		for(int idx = len-2; idx >= 0; idx--)
			result[idx] = (b[idx]-this.b*result[idx+1])/diagonal[0][idx];
		
		
		return result;
	}
}
