package mumfordShah;

public class L2MumfordShah1DFaster extends AbstractPotts1D {
	
	/*
	 * 		Computes the solution to the following 1D minimization problem
	 * 			
	 * 			min	
	 */
	
	double alpha;
	private double[] prevD1r, prevDlr;
	private double[] d1r, dlr;
	private double[] prevU1r, prevUlr;
	private double[] u1r, ulr;
	private double[][] tempF;
	private int curL,curR;
	GaussElim gauss;
	double[] scriptF;
	double[][][] Dlr;

	public L2MumfordShah1DFaster(){
		double[][] f = new double[5][5];
		double gamma = 1;
		double alpha = 1;
		set(f,gamma,alpha);
	}

	public L2MumfordShah1DFaster(double[] f, double gamma, double alpha){
		double[][] tempf = new double[1][1];
		tempf[0] = f;
		set(tempf,gamma,alpha);
	}


	public L2MumfordShah1DFaster(double[][] f, double gamma, double alpha){
		set(f,gamma,alpha);
	}
	
	public L2MumfordShah1DFaster(double[][] f, double gamma, double alpha, GaussElim gauss){
		set(f,gamma,alpha,gauss);
	}

	public void set(double[][] f, double gamma, double alpha){
		GaussElim tempGauss = new GaussL2Mum(f[0].length,alpha);
		set(f,gamma,alpha,tempGauss);
	}
	
	public void set(double[][] f,double gamma,double alpha, GaussElim gauss){
		set(f,gamma);
		
		// initialize scriptF perhaps make it possible to hand over a pre computed scriptF since it is always the same
		scriptF = new double[n];
		scriptF[0] = 1;
		for(int idx = 1; idx < scriptF.length; idx++)
			scriptF[idx] = alpha * scriptF[idx-1] / (alpha + scriptF[idx-1]) + 1;
		
		curL = 0; curR = 0;
		
		this.alpha = alpha;
		
		this.gauss = gauss;
		
	}
	

	@Override
	protected double[][] computeMu_LtoR(int l, int r) {
		double[][] tempf = new double[channels][r-l+1];
		for(int dim = 0; dim < channels; dim++)
			for(int i = l-1; i < r; i++)
				tempf[dim][i-(l-1)] = f[dim][i];
		return computeSolution(tempf);
	}
	
	
	private double[][] computeSolution(double[][] data){

		double[][] result = new double[channels][data[0].length];

		for(int dim = 0; dim < channels; dim++){
			result[dim] = gauss.computeMu(data[dim]);
		}


		return result;
	}
	
	@Override
	protected double computeDlr(int l, int r) {
		if(r == 1 && l == 1){
			prevD1r = new double[channels];
			d1r = new double[channels];
			prevU1r = new double[channels];
			u1r = new double[channels];
			for(int dim = 0; dim < channels; dim++){
				u1r[dim] = f[dim][r-1];
			}
			curR = r; curL = l;
		}

		else if (r == curR+1 && l == 1){
			prevD1r = d1r;
			d1r = new double[channels];
			prevU1r = u1r;
			u1r = new double[channels];
			curR = r; curL = l;
		}

		else if(r == curR && l == curR){
			prevDlr = new double[channels];
			dlr = new double[channels];
			prevUlr = new double[channels];
			ulr = new double[channels];
			
			// prepare data for "backwards" computation
			tempF = new double[channels][r-1];
			for(int dim = 0; dim < channels; dim++) for(int iter = 0; iter < r-1; iter++)
				tempF[dim][iter] = f[dim][r-1-iter];
			
			for(int dim = 0; dim < channels; dim++){
				ulr[dim] = tempF[dim][0];
			}
			curR = r; curL = l;
		}

		else if (r == curR && l == curL-1){
			prevDlr = dlr;
			dlr = new double[channels];
			prevUlr = ulr;
			ulr = new double[channels];
			curR = r; curL = l;
		}
		// we don't have to compute something new
		else if(r == curR && l == curL-1);
		// else throw exception
		else
			throw(new IllegalArgumentException("l has to be curL-1"));
		
		// return the computed value
		if(l == 1){
			for(int dim = 0; dim < channels; dim++){
				u1r[dim] = ((scriptF[r-l] - 1) * prevU1r[dim] + f[dim][r-l])/scriptF[r-l];
				d1r[dim] = prevD1r[dim] + (scriptF[r-l] - 1)*Math.pow(u1r[dim] - prevU1r[dim], 2) + Math.pow(u1r[dim] - f[dim][r-l], 2);
			}
			return sum(d1r);
		}	
		else{
			for(int dim = 0; dim < channels; dim++){
				ulr[dim] = ((scriptF[r-l] - 1)*prevUlr[dim] + tempF[dim][r-l])/scriptF[r-l];
				dlr[dim] = prevDlr[dim] + (scriptF[r-l] - 1)*Math.pow(ulr[dim] - prevUlr[dim], 2) + Math.pow(ulr[dim] - tempF[dim][r-l], 2);
			}
			return sum(dlr);	
		}
	}
	
	private double sum(double[] data){
		double summedData = 0;
		for(int i = 0; i < data.length; i++)
			summedData += data[i];
		return summedData;
	}

}
