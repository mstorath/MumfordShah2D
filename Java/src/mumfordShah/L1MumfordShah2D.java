package mumfordShah;


public class L1MumfordShah2D extends AbstractPotts2DParallel{

	
	private double alpha;
	private L1MumfordShah1D mumShah;
	
	
	public L1MumfordShah2D(double[][] image, double gamma, double alpha, double mu0, double tau){
		super(image,gamma,mu0,tau);
		double[][][] tempImage = new double[1][1][1];
		tempImage[0] = image;
		set(tempImage, gamma, alpha, mu0, tau);
	}	
	
	public L1MumfordShah2D(double[][][] image, double gamma, double alpha, double mu0, double tau){
		super(image,gamma,mu0,tau);
		set(image, gamma, alpha, mu0, tau);
	}
	
	
	// set and reset method
	public void set(double[][][] image, double gamma, double alpha, double mu0, double tau){
		set(image,gamma,mu0,tau);
		this.alpha = alpha;

	}

	protected double[][] compute1DSolution(double[][] rowCol, double mu) {
		mumShah = new L1MumfordShah1D();
		mumShah.set(rowCol, 2*gamma/(1+mu), 2*alpha/(1+mu));
		mumShah.run();
		return mumShah.getResult();
	}


	@Override
	protected AbstractPotts1D newAbstractPotts1D(double[][] rowCol, double mu) {
		return new L1MumfordShah1D(rowCol, 2*gamma/(1+mu), 2*alpha/(1+mu));
	}

	@Override
	protected AbstractPotts1D newAbstractPotts1D(double[][] rowCol, double mu, GaussElim gauss) {
		// dummy function will never be called
		return null;
	}
}
