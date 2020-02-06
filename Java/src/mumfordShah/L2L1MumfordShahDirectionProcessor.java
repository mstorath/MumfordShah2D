package mumfordShah;

public class L2L1MumfordShahDirectionProcessor extends
AbstractDirectionProcessor {
	
	private double alpha;

	// set and reset method
	public void set(double[][][] image, double gamma, double alpha, int[] direction, double directionWeight){
		set(image, gamma, direction, directionWeight);
		this.alpha = alpha;
	}

	@Override
	protected AbstractPotts1D newAbstractPotts1D(double[][] rowCol) {
		return new L1MumfordShah1D(rowCol,gamma,alpha);
	}

	@Override
	protected AbstractPotts1D newAbstractPotts1D(double[][] rowCol,
			GaussElim gauss) {
		return new L1MumfordShah1D(rowCol,gamma,alpha);
	}

}
