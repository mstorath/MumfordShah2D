package mumfordShah;

public class Image {

	double[][][] image;
	int channels, nRow, nCol;
	
	// Constructor for an image from an array
	Image(double[][][] image){
		setImage(image);
	}
	
	// Constructor for a zero image
	Image(int channels, int nRow, int nCol){
		this.image = new double[channels][nRow][nCol];
		this.channels = channels;
		this.nRow = nRow;
		this.nCol = nCol;
	}
	
	// get number of rows
	public int getNRow(){
		return nRow;
	}
	
	// get number of columns
	public int getNCol(){
		return nCol;
	}
	
	// get number of channels
	public int getNChannels(){
		return channels;
	}
	
	// get specified row of the image
	public double[][] getRow(int index){
		double[][] result = new double[channels][nCol];
		
		for(int dim = 0; dim < channels; dim++) for(int j = 0; j < nCol; j++)
			result[dim][j] = image[dim][index][j];
		
		return result;
	}
	
	// get specified column of the image
	public double[][] getCol(int index){
		double[][] result = new double[channels][nRow];
		
		for(int dim = 0; dim < channels; dim++) for(int j = 0; j < nRow; j++)
			result[dim][j] = image[dim][j][index];
		
		return result;
	}
	
	// set specified row of the image
	public boolean setRow(double[][] row, int index){
		// test if there is the right amount of elements in row
		boolean success = false;
		if(row[0].length == nCol && row.length == channels) success = true;
		
		// overwrite the elements of the specified row
		for(int dim = 0; dim < channels; dim++) for(int j = 0; j < nCol; j++)
			image[dim][index][j] = row[dim][j];
		
		return success;
	}
	
	// set specified column of the image
	public boolean setCol(double[][] col, int index){
		// test if there is the right amount of elements in col
		boolean success = false;
		if(col[0].length == nRow && col.length == channels) success = true;

		// overwrite the elements of the specified column
		for(int dim = 0; dim < channels; dim++) for(int j = 0; j < nRow; j++)
			image[dim][j][index] = col[dim][j];
		return success;
	}
	
	public boolean setImage(double[][][] image){
		this.image = image;
		channels = image.length;
		nRow = image[0].length;
		nCol = image[0][0].length;
		return true;
	}
	
	public Image clone(){
		double[][][] cloned = new double[channels][nRow][nCol];
		
		for(int dim = 0; dim < channels; dim++) for(int i = 0; i < nRow; i++) for(int j = 0; j < nCol; j++)
			cloned[dim][i][j] = image[dim][i][j];
		
		return new Image(cloned);
		
	}
	
	public double[][][] getArray(){
		return image;
	}
	
	public static double computeError(Image imgU, Image imgV) throws IllegalArgumentException{
		double[][][] u = imgU.getArray();
		double[][][] v = imgV.getArray();
		
		if(imgU.getNRow() != imgV.getNRow() || imgU.getNCol() != imgV.getNCol() || imgU.getNChannels() != imgV.getNChannels())
			throw(new IllegalArgumentException("Images have to be of the same dimension"));

		int channelsDim = imgU.getNChannels();
		int rows = imgU.getNRow();
		int cols = imgU.getNCol();
		
		double error = 0;		
		
		for(int dim = 0; dim < channelsDim; dim++) for(int i = 0; i < rows; i++) for(int j = 0; j < cols; j++)
			error += Math.pow(u[dim][i][j]-v[dim][i][j],2);
		
		return error;
	}

	public void updateDualParam(double mu, Image u, Image v) {
		double[][][] arrayU, arrayV;
		arrayU = u.getArray();
		arrayV = v.getArray();
		
		for(int dim = 0; dim < channels; dim++) for(int i = 0; i < nRow; i++) for(int j = 0; j < nCol; j++)
			image[dim][i][j] += mu*(arrayU[dim][i][j] - arrayV[dim][i][j]);
		
	}

	public static double[][] getRow(Image img, Image v, Image lambda, double mu, int idx) {
		double[][] arrayF, arrayV, arrayLambda, result;
		
		arrayF = img.getRow(idx);
		arrayV = v.getRow(idx);
		arrayLambda = lambda.getRow(idx);
		
		int n = arrayF[0].length;
		int channels = arrayF.length;
		
		result = new double[channels][n];
		
		for(int dim = 0; dim < channels; dim++) for(int i = 0; i < n; i++)
			result[dim][i] = (arrayF[dim][i] + mu*arrayV[dim][i] - arrayLambda[dim][i])/(1+mu);
		
		return result;
	}
	
	public static double[][] getCol(Image img, Image u, Image lambda, double mu, int idx) {
		double[][] arrayF, arrayU, arrayLambda, result;
		
		arrayF = img.getCol(idx);
		arrayU = u.getCol(idx);
		arrayLambda = lambda.getCol(idx);
		
		int n = arrayF[0].length;
		int channels = arrayF.length;
		
		result = new double[channels][n];
		
		for(int dim = 0; dim < channels; dim++) for(int i = 0; i < n; i++)
			result[dim][i] = (arrayF[dim][i] + mu*arrayU[dim][i] + arrayLambda[dim][i])/(1+mu);
		
		return result;
	}
	
	public double[][] getDirection(int[] direction, int idx){
		double[][] result = null;
		int x0, y0; // start coordinates
		int stepX, stepY; // step size in x/y direction
		int size;					// size of the result array
		
		int[] startCoord = computeStartCoordinates(direction, idx);
		x0 = startCoord[0];
		y0 = startCoord[1];
		
		size = computeSize(direction, x0, y0);
		
		stepX = direction[0];
		stepY = direction[1];
		
		result = new double[channels][size];
		
		int iter = 0;
		
		while(iter < size){
			for(int dim = 0; dim < channels; dim++){
				result[dim][iter] = image[dim][x0 + iter*stepX][y0 + iter*stepY];
			}
			iter++;
		}
		
		return result;
	}
	
	public void setDirection(double[][] data, int[] direction, int idx){
		int x0, y0; // start coordinates
		int stepX, stepY; // step size in x/y direction
		int size;					// size of the result array
		
		int[] startCoord = computeStartCoordinates(direction, idx);
		x0 = startCoord[0];
		y0 = startCoord[1];
		
		size = computeSize(direction, x0, y0);
		
		stepX = direction[0];
		stepY = direction[1];
				
		int iter = 0;
		
		while(iter < size){
			for(int dim = 0; dim < channels; dim++){
				image[dim][x0 + iter*stepX][y0 + iter*stepY] = data[dim][iter];
			}
			iter++;
		}
	}
	
	public int computeSize(int[] direction, int x0, int y0){
		int size = 0;
		int maxX, maxY;
		
		maxX = (int)Math.ceil((double)(nRow-x0)/direction[0]);
		
		if (direction[1] == 0)
			maxY = Integer.MAX_VALUE;
		else
			maxY = (int)Math.ceil((double)(nCol-y0)/direction[1]);
		
		size = Math.min(maxX, maxY);
		
		return size;
	}
	
	public int[] computeStartCoordinates(int[] direction, int idx){
		int result[] = new int[2];
		if (idx < nCol){
			result[0] = 0;
			result[1] = idx;
		}
		else if(direction[0] == 1 && direction[1] == 1){
			result[0] = idx-nCol;
			result[1] = 0;
		}
		else if(direction[0] == 2 && direction[1] == 1){
			if(idx < 2*nCol){
				result[0] = 1;
				result[1] = idx-nCol;
			}
			else{
				result[0] = 2 + (idx-2*nCol);
				result[1] = 0;
			}	
		}
		else if(direction[0] == 1 && direction[1] == 2){
			result[0] = 1 + (idx - nCol)/2;
			result[1] = (idx - nCol) % 2;
		}
		return result;
	}
}
