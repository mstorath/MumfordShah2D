package mumfordShah;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

public abstract class AbstractDirectionProcessor {
	private Image img;
	protected double gamma;
	protected int channels, nRow, nCol; 			// number of color channels (channels), height (nRow) and width (nCol) of the image
	boolean verbose = true;
	boolean L2L2 = false;							// flag for L2L2 MumfordShah
	protected GaussElim gauss;
	int[] direction;
	double directionWeight;
	private int numThreads = 1;


	public AbstractDirectionProcessor(){

	}


	// set and reset method
	public void set(double[][][] image, double gamma, int[] direction, double directionWeight){
		// set parameters
		this.img = new Image(image);
		this.gamma = gamma;
		
		this.direction = direction;
		this.directionWeight = directionWeight;
		
		nRow = img.getNRow();
		nCol = img.getNCol();
		channels = img.getNChannels();
	}
	
	public void setNumThreads(int num){
		this.numThreads = num;
	}
	
	public double[][][] run(){

		// start executor service (perhaps configure to make it not multi threading comaptible)
		ExecutorService executor = Executors.newFixedThreadPool(numThreads);
		
		// compute the number of "direction rows" to be computed (depending on the direction it is a different amount) if a direction is not yet defined it just doesn't compute anything
		int amount = 0;
		if(direction[0] == 1 && direction[1] == 0)
			amount = nCol;
		else if(direction[0] == 1 && direction[1] == 1)
			amount = nRow + nCol - 1;
		else if(direction[0] == 2 && direction[1] == 1)
			amount = 2*nCol + nRow-2;
		else if(direction[0] == 1 && direction[1] == 2)
			amount = nCol + 2*(nRow - 1);
		

		// generate new future task list
		List<Future<double[][]>> futureList = new ArrayList<Future<double[][]>>();
		// iterate over specified direction
		for(int idx = 0; idx < amount; idx++){
			double[][] row;

			row = img.getDirection(direction, idx);

			// initialize new one dimensional potts
			Callable<double[][]> toCall;
			if(!L2L2){
				toCall = newAbstractPotts1D(row);
			}
			else{
				toCall = newAbstractPotts1D(row, gauss);
			}
			
			// add task to future list
			Future<double[][]> future = executor.submit(toCall);
			futureList.add(future);
		}

		// get the results from the future list
		for(int idx = 0; idx < amount; idx++){
			try{
				Future<double[][]> fut = futureList.get(idx);
				double[][] result = fut.get();
				img.setDirection(result, direction, idx);

			} catch(Exception e) {
				e.printStackTrace();
			}
		}

		// shutdown the executor
		executor.shutdown();

		// return the result
		return img.getArray();
	}

	protected abstract AbstractPotts1D newAbstractPotts1D(double[][] rowCol);
	protected abstract AbstractPotts1D newAbstractPotts1D(double[][] rowCol, GaussElim gauss);
}
