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
        double[] norm = new double[wavelet.length];

        for(int i = 0; i < wavelet.length; i++){
            norm[i] = Math.max(Math.abs(minValues[i]), nthPercentile);
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

        for(double w : wavelet) {
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
}
