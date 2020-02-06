package mumfordShah;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

public abstract class AbstractPotts2DParallel {
	private Image img;
	protected double gamma;
	private double mu0;
	private double tau;
	protected int channels, nRow, nCol; 			// number of color channels (channels), height (nRow) and width (nCol) of the image
	private int iter;							// iteration counter
	private Image u,v,lambda;					// iterator variables
	private double eps;							// convergence parameter
	boolean verbose = true;
	boolean L2L2 = false;						// flag for L2L2 MumfordShah
	protected GaussElim gauss;
	private int numThreads = 1;
	

	public AbstractPotts2DParallel(double[][] image, double gamma, double mu0, double tau){
		double[][][] tempImage = new double[1][1][1];
		tempImage[0] = image;
		set(tempImage, gamma, mu0, tau);
	}
	
	public AbstractPotts2DParallel(double[][][] image, double gamma, double mu0, double tau){
		set(image, gamma, mu0, tau);
	}


	// set and reset method
	public void set(double[][][] image, double gamma, double mu0, double tau){
		// set parameters
		this.img = new Image(image);
		this.gamma = gamma;
		this.mu0 = mu0;
		this.tau = tau;
		nRow = img.getNRow();
		nCol = img.getNCol();
		channels = img.getNChannels();
		eps = 0.00001;
		iter = 0;
		


		// initialize iterators
		resetIterators();
	}
	
	public void setNumThreads(int num){
		this.numThreads = num;
	}

	public void resetIterators(){
		v = img.clone();
		u = new Image(channels,nRow,nCol);
		lambda = new Image(channels,nRow,nCol);
	}

	public boolean setEps(double eps){
		this.eps = eps;
		return true;
	}

	public double[][][] run(){

		double err = Image.computeError(u, v);
		double mu = mu0;

		while(err/(nRow*nCol*channels) > eps){
			// compute solutions horizontally
			ExecutorService executor = Executors.newFixedThreadPool(numThreads);
			List<Future<double[][]>> futureList = new ArrayList<Future<double[][]>>();
			for(int idx = 0; idx < nRow; idx++){
				double[][] row;

				row = Image.getRow(img, v, lambda, mu, idx);
				Callable<double[][]> toCall;
				if(!L2L2){
					toCall = newAbstractPotts1D(row, mu);
				}
				else{
					toCall = newAbstractPotts1D(row, mu, gauss);
				}
				Future<double[][]> future = executor.submit(toCall);
				futureList.add(future);
			}

			for(int idx = 0; idx < nRow; idx++){
				try{
					Future<double[][]> fut = futureList.get(idx);
					double[][] result = fut.get();
					u.setRow(result, idx);

					if(verbose){
						if(idx%128 == 0)
							System.out.print("\n");
						System.out.print("*");
					}
				} catch(Exception e) {
					e.printStackTrace();
				}
			}

			// compute solutions vertically
			futureList = new ArrayList<Future<double[][]>>();
			for(int idx = 0; idx < nCol; idx++){
				double[][] col;

				col = Image.getCol(img, u, lambda, mu, idx);
				Callable<double[][]> toCall;
				if(!L2L2){
					toCall = newAbstractPotts1D(col, mu);
				}
				else{
					toCall = newAbstractPotts1D(col, mu, gauss);
				}
				Future<double[][]> future = executor.submit(toCall);
				futureList.add(future);

			}

			for(int idx = 0; idx < nCol; idx++){
				try{
					Future<double[][]> fut = futureList.get(idx);
					double[][] result = fut.get();
					v.setCol(result, idx);

					if(verbose){
						if(idx%128 == 0)
							System.out.print("\n");
						System.out.print("*");
					}
				} catch(Exception e) {
					e.printStackTrace();
				}

			}

			executor.shutdown();

			// update dual parameter
			lambda.updateDualParam(mu,u,v);

			// update coupling parameter
			// mu = tau*mu;
			mu = mu0*Math.pow(iter*tau, 2);
			// mu = mu0*Math.pow(iter, tau);
			err = Image.computeError(u, v);
			iter++;

			if(verbose){
				System.out.print("\n");
				System.out.println(iter + ". Iteration:");
				System.out.println(" average quadratic error per pixel and channel =  " + err/(nRow*nCol*channels));
				System.out.println(" mu =  " + mu);
			}
			else{
				System.out.print("*");
			}


			if(iter > 100){
				System.out.println("Method didn't converge!");
				break;
			}

		}


		System.out.println();
		return u.getArray();
	}

	protected abstract AbstractPotts1D newAbstractPotts1D(double[][] rowCol, double mu);
	protected abstract AbstractPotts1D newAbstractPotts1D(double[][] rowCol, double mu, GaussElim gauss);
}
