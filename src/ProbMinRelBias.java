import cern.colt.matrix.impl.SparseObjectMatrix2D;

/**
 * @author Ulfiani Primawati
 *
 * This is an implementation of the probabilistic thresholding technique
 * that minimizes the maximum normalized bias. This technique does not perform
 * randomized rounding. Instead, each coefficient is either retained or
 * discarded, according to the probabilities y(i), and y(i)s are selected
 * to minimize a desired error metric.
 *
 * This algorithm is based on a paper proposed by M. Garofalakis and
 * P. G. Gibbons in 2004, titled "Probabilistic Wavelet Synopses".
 */

public class ProbMinRelBias {

    public static double[] perturbAndCalcNorm(double[] wavelet, int[] nzArray, double[] data, double percentile){
        double nthPercentile = ProbMinRelVar.findPercentile(data, percentile);
        double[] minValues = ProbMinRelVar.applyPerturbationRule(wavelet, nzArray, data, nthPercentile);
        double[] norm = new double[2 * wavelet.length];

        for(int i = 0; i < wavelet.length; i++){
            norm[i] = Math.max(Math.abs(minValues[i]), nthPercentile);
        }

        // Added this
        int j = 0;
        for(int i = wavelet.length; i < norm.length; i++){
            norm[i] = Math.max(Math.abs(data[j]), nthPercentile);
            j++;
        }

        return norm;
    }

    public static void callMainFunction(double[] wavelet, double b, int q, double[] data, double percentile){
        int[] nzArray = ProbMinRelVar.constructNzArray(wavelet);
        int length_1 = wavelet.length;
        int length_2 = (int)Math.round(b * q) + 1;
        double[] norm = perturbAndCalcNorm(wavelet, nzArray, data, percentile);

        // Instantiate variables for calculating optimal y(i)
        double[][] mValues = new double[length_1][length_2];
        double[][] yValues = new double[length_1][length_2];
        double[][] leftAllot = new double[length_1][length_2];
        boolean[][] checked = new boolean[length_1][length_2];

        // Calculate optimal value for j = overall root of the tree
        double root;
        double rootSpace;
        int indexB = (int)Math.rint(b * q);

        mValues[0][(int)Math.rint(b * q)] = Double.POSITIVE_INFINITY;

        for(int l = 1; l <= q; l++){
            if(wavelet[0] == 0.0){
                root = 0;
                rootSpace = 0;
            }
            else{
                rootSpace = (double)l / (double)q;
                root = ((1.0 - rootSpace) * Math.abs(wavelet[0])) / norm[1];
            }

            for(int _b = 0; _b <= (int)Math.rint(b * q) - l; _b++){
                double __b = (double)_b / (double)q;
                double next = getOptimalNormBias(wavelet, __b, q, 1,
                        nzArray, norm, mValues, yValues, leftAllot, checked);
                if(root + next < mValues[0][indexB]){
                    mValues[0][indexB] = root + next;
                    yValues[0][indexB] = rootSpace;
                    leftAllot[0][indexB] = __b;
                }
            }

            if(wavelet[0] == 0.0)
                break;
        }

        checked[0][indexB] = true;

        // Put all optimal y values for each coefficient in an array
        double[] chosenY = new double[length_1];
        double[] bValue = new double[length_1];

        bValue[0] = b;

        chosenY[0] = yValues[0][(int)(b * q)];
        bValue[1] = leftAllot[0][(int)(b * q)];

        for(int i = 1; i < length_1; i++){
            chosenY[i] = yValues[i][(int)Math.round(bValue[i] * q)];
            if(i * 2 < length_1){
                bValue[2 * i] = leftAllot[i][(int)Math.round(bValue[i] * q)];
                bValue[2 * i + 1] = bValue[i] - bValue[2 * i] - chosenY[i];
            }
        }

        // Perform coin flips
        performCoinFlips(wavelet, chosenY);

        double sum = 0;
        for(double y : chosenY)
            sum += y;

        System.out.println(sum);
        System.out.println();

        for(double w : chosenY) {
            System.out.println(w);
        }
    }

    public static void performCoinFlips(double[] wavelet, double[] chosenY){
        for(int i = 0; i < wavelet.length; i++){
            if(!retainCoefficient(chosenY[i]))
                wavelet[i] = 0.0;
        }
    }

    public static boolean retainCoefficient(double yValue){
        double randomDouble = Math.random();

        if(yValue == 0)
            return false;

        return randomDouble <= yValue;
    }

    public static double getOptimalNormBias(double[] wavelet, double B, int q, int root,
                                            int[] nzArray, double[] norm,
                                            double[][] mValues, double[][]yValues, double[][] leftAllot,
                                            boolean[][] checked){

        int indexB = (int)Math.rint(B * q);

        if(root > wavelet.length - 1 || B - (double)nzArray[root] >= 1e-6)
            return 0;

        if(nzArray[root] > B * q)
            return Double.POSITIVE_INFINITY;

        if(checked[root][indexB])
            return mValues[root][indexB];

        mValues[root][indexB] = Double.POSITIVE_INFINITY;

        double rootLeft;
        double rootRight;
        double rootSpace;

        for(int l = 1; l <= q; l++){
            if(wavelet[root] == 0.0){
                rootLeft = 0.0;
                rootRight = 0.0;
                rootSpace = 0.0;
            }
            else if(2 * root > wavelet.length - 1){
                rootSpace = (double)l / (double)q;
                rootLeft = ((1.0 - rootSpace) * Math.abs(wavelet[root])) / norm[root];
                rootRight = rootLeft;
            }
            else{
                rootSpace = (double)l / (double)q;
                rootLeft = ((1.0 - rootSpace) * Math.abs(wavelet[root])) / norm[2 * root];
                rootRight = ((1.0 - rootSpace) * Math.abs(wavelet[root])) / norm[2 * root + 1];
            }

            for(int _b = 0; _b <= indexB - l; _b++){
                double __b = (double)_b / (double)q;
                double left = getOptimalNormBias(wavelet, __b, q, 2 * root,
                        nzArray, norm, mValues, yValues, leftAllot, checked);
                double right = getOptimalNormBias(wavelet, B - rootSpace - __b, q, 2 * root + 1,
                        nzArray, norm, mValues, yValues, leftAllot, checked);
                double max = Math.max(rootLeft + left, rootRight + right);

                if(max < mValues[root][indexB]){
                    mValues[root][indexB] = max;
                    yValues[root][indexB] = rootSpace;
                    leftAllot[root][indexB] = __b;
                }
            }

            if(wavelet[root] == 0)
                break;
        }

        checked[root][indexB] = true;
        return mValues[root][indexB];
    }

    public static double getOptimalNormBias(double[] wavelet, double b, int q, int root,
                                            int[] nzArray, double[] norm, SparseObjectMatrix2D objectMatrix2D){

        int indexB = (int)Math.rint(b * q);

        if(root > wavelet.length - 1 || b - (double)nzArray[root] >= 0)
            return 0;

        if(!(objectMatrix2D.getQuick(indexB, root) == null))
            return ((MinNSE)objectMatrix2D.getQuick(indexB, root)).getmValue();

        MinNSE obj = new MinNSE(Double.POSITIVE_INFINITY, 0, 0);
        objectMatrix2D.setQuick(indexB, root, obj);

        double rootLeft;
        double rootRight;
        double rootSpace;

        for(int l = 0; l <= q; l++){
            if(wavelet[root] == 0.0){
                rootLeft = 0.0;
                rootRight = 0.0;
                rootSpace = 0.0;
            }

            else{
                rootSpace = (double)l / (double)q;
                rootLeft = ((1.0 - rootSpace) * Math.abs(wavelet[root])) / norm[2 * root];
                rootRight = ((1.0 - rootSpace) * Math.abs(wavelet[root])) / norm[2 * root + 1];
            }

            for(int _b = 0; _b <= indexB - (int)Math.rint(rootSpace * q); _b++){
                double __b = (double)_b / (double)q;
                double left = getOptimalNormBias(wavelet, __b, q, 2 * root,
                        nzArray, norm, objectMatrix2D);
                double right = getOptimalNormBias(wavelet, b - rootSpace - __b, q, 2 * root + 1,
                        nzArray, norm, objectMatrix2D);
                double max = Math.max(rootLeft + left, rootRight + right);

                if(max < ((MinNSE)objectMatrix2D.getQuick(indexB, root)).getmValue()){
                    MinNSE newObj = new MinNSE(max, rootSpace, __b);
                    objectMatrix2D.setQuick(indexB, root, newObj);
                }
            }

            if(wavelet[root] == 0)
                break;
        }

        return ((MinNSE)objectMatrix2D.getQuick(indexB, root)).getmValue();
    }

    public static double[] callMainFunction1(double[] wavelet, double b, int q, double[] data, double percentile){

        long startingTime = System.currentTimeMillis();

        int[] nzArray = ProbMinRelVar.constructNzArray(wavelet);
        int length_1 = wavelet.length;
        int length_2 = (int)Math.round(b * q) + 1;
        double[] norm = perturbAndCalcNorm(wavelet, nzArray, data, percentile);

        SparseObjectMatrix2D objectMatrix2D = new SparseObjectMatrix2D(length_2, length_1);

        // Calculate optimal value for j = overall root of the tree
        double root;
        double rootSpace;
        int row = (int)Math.rint(b * q);
        int column = 0;

        MinNSE obj = new MinNSE(Double.POSITIVE_INFINITY, 0, 0);
        objectMatrix2D.setQuick(row, column, obj);

        for(int l = 1; l <= q; l++){
            if(wavelet[0] == 0.0){
                root = 0;
                rootSpace = 0;
            }
            else{
                rootSpace = (double)l / (double)q;
                root = ((1.0 - rootSpace) * Math.abs(wavelet[0])) / norm[1];
            }

            for(int _b = 0; _b <= (int)Math.rint(b * q) - (int)Math.rint(rootSpace * q); _b++){
                double __b = (double)_b / (double)q;
                double next = getOptimalNormBias(wavelet, __b, q, 1,
                        nzArray, norm, objectMatrix2D);
                if(root + next <= ((MinNSE)objectMatrix2D.getQuick(row, column)).getmValue()){
                    MinNSE newObj = new MinNSE(root + next, rootSpace, __b);
                    objectMatrix2D.setQuick(row, column, newObj);
                }
            }

            if(wavelet[0] == 0.0)
                break;
        }

        // Put all optimal y values for each coefficient in an array
        double[] chosenY = new double[length_1];
        double[] bValue = new double[length_1];

        bValue[0] = b;

        chosenY[0] = ((MinNSE)objectMatrix2D.getQuick(row, column)).getyValue();
        bValue[1] = ((MinNSE)objectMatrix2D.getQuick(row, column)).getLeftAllot();

        // Print mValue, yValue, leftAllot
//        for(int i = 0; i < length_2; i++){
//            //System.out.print(i + " ");
//            for(int j = 0; j < length_1; j++){
//                if(!(objectMatrix2D.getQuick(i, j) == null)) {
//                    MinNSE store = (MinNSE) objectMatrix2D.getQuick(i, j);
//                    System.out.print(store.getmValue() + "\t" + store.getyValue() + "\t" + store.getLeftAllot() + "\t");
//                }
//                else
//                    System.out.print(0 + "\t" + 0 + "\t" + 0 + "\t");
//            }
//            System.out.println();
//        }


        for(int i = 1; i < length_1; i++){
            int ithRow = (int)(Math.rint(bValue[i] * q));
            if(objectMatrix2D.getQuick(ithRow, i) == null) {
                chosenY[i] = 1;
                System.out.println("coefficient " + wavelet[i]);
            }
            else
                chosenY[i] = ((MinNSE)objectMatrix2D.getQuick(ithRow, i)).getyValue();

            if(i * 2 < length_1) {
                if(objectMatrix2D.getQuick(ithRow, i) == null)
                    bValue[i * 2] = 1;
                else
                    bValue[i * 2] = ((MinNSE)objectMatrix2D.getQuick(ithRow, i)).getLeftAllot();
                bValue[i * 2 + 1] = bValue[i] - bValue[i * 2] - chosenY[i];
            }
        }

        double sum = 0;
        for(int i = 0; i < length_1; i++){
            //System.out.println(chosenY[i]);
            sum += chosenY[i];
        }

        performCoinFlips(wavelet, chosenY);
        System.out.println("sum " + sum);

        System.out.println("Time for executing: " +
                (System.currentTimeMillis() - startingTime) + "milliseconds.");

        return chosenY;
    }
}
