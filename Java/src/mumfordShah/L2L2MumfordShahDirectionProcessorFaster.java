package mumfordShah;

public class L2L2MumfordShahDirectionProcessorFaster extends
AbstractDirectionProcessor {
	
	private double alpha;

	// set and reset method
	public void set(double[][][] image, double gamma, double alpha, int[] direction, double directionWeight){
		set(image, gamma, direction, directionWeight);
		this.alpha = alpha;
	}
	
	public void setGauss(GaussElim gauss){
		if(gauss == null){
			this.L2L2 = false;
			this.gauss = null;
		}
		else{
			this.L2L2 = true;
			this.gauss = gauss;
		}
		
	}

	@Override
	protected AbstractPotts1D newAbstractPotts1D(double[][] rowCol) {
		return new L2MumfordShah1DFaster(rowCol,gamma,alpha);
	}

	@Override
	protected AbstractPotts1D newAbstractPotts1D(double[][] rowCol,
			GaussElim gauss) {
		return new L2MumfordShah1DFaster(rowCol,gamma,alpha,gauss);
	}

}
