import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * @author Ulfiani Primawati
 *
 * This is an implementation of the probabilistic thresholding technique
 * that minimizes the maximum relative error. This algorithm
 * deterministically keeps the most important wavelet coefficients while
 * randomly rounding the other coefficients either down to zero or up to a
 * larger value (called a "rounding value").
 *
 * This algorithm is based on a paper proposed by M. Garofalakis and
 * P. G. Gibbons in 2004, titled "Probabilistic Wavelet Synopses".
 */

public class ProbMinRelVar {
    /**
     * This method performs perturbation rule. For each subtree Tj such that
     * (1) one of its child subtrees, say T(2j) has all zero coefficients,
     * (2) its other child subtree T(2j+1) has at least one nonzero coefficient, and
     * (3) the minimum data value in T2j is less than the minimum data value in T(2j+1),
     * we perturb c(2j) in wavelet. "wavelet" will contain the coefficients after
     * pertubation.
     *
     * @param wavelet the wavelet
     * @param nzArray an array that consists of the number of nonzero coefficients rooted at each subtree
     * @param data the original data
     * @param nthPercentile Nth percentile value
     */
    public static double[] applyPerturbationRule(double[] wavelet, int[] nzArray, double[] data, double nthPercentile){
        double delta = Math.min(0.01, nthPercentile / 100.0);
        double[] minValues = minEachSubtree(data);

        for(int i = (wavelet.length / 2) - 1; i >= 1; i--){
            // If the number of non-zero coefficients is greater than 0, then check its children's
            if(nzArray[i] != 0){
                boolean perturbed = false;
                double minLeftSubTree = minValues[2 * i];
                double minRightSubTree = minValues[2 * i + 1];
                int nzLeft = nzArray[2 * i];
                int nzRight = nzArray[2 * i + 1];
                if(minLeftSubTree < minRightSubTree && nzLeft == 0 && nzRight > 0){
                    //Apply perturbation to the left child
                    wavelet[2 * i] = delta * (double)flipCoin();
                    nzArray[2 * i] = nzArray[2 * i] + 1;
                    perturbed = true;
                }
                else if(minLeftSubTree > minRightSubTree && nzRight == 0 && nzLeft > 0){
                    //Apply perturbation to the right child
                    wavelet[2 * i + 1] = delta * (double)flipCoin();
                    nzArray[2 * i + 1] = nzArray[2 * i + 1] + 1;
                    perturbed = true;
                }

                // Update the rest of nzArray if perturbation performed
                if(perturbed){
                    int index = i;
                    int nextIndex;
                    while(index >= 1){
                        if(index % 2 == 0)
                            nextIndex = index / 2;
                        else
                            nextIndex = (index - 1) / 2;
                        nzArray[index] = nzArray[index] + 1;
                        index = nextIndex;
                    }
                    nzArray[0] = nzArray[0] + 1;
                }

            }
        }
        return minValues;
    }

    private static int flipCoin(){
        double randomDouble = Math.random();
        if(randomDouble < 0.5)
            return 1;
        else
            return -1;
    }

    public static double[] minEachSubtree(double[] data){
        double[] minimumValues = new double[data.length];
        double leftVal;
        double rightVal;

        // Find minimum data for each subtree rooted at the last level
        int j = data.length - 1;
        for(int i = data.length - 1; i >= data.length / 2; i--){
            rightVal = data[j];
            j--;
            leftVal = data[j];
            j--;
            minimumValues[i] = Math.min(leftVal, rightVal);
        }

        // Find minimum data for each subtree rooted at other level
        for(int k = (data.length / 2) - 1; k >= 1; k--){
            minimumValues[k] = Math.min(minimumValues[2 * k], minimumValues[2 * k + 1]);
        }

        // Smallest element rooted at c0 is the same as c1
        minimumValues[0] = minimumValues[1];

        return minimumValues;
    }

    public static double findPercentile(double[] data, double percentile){
        List<Double> sortedData = new ArrayList<>();
        for (double datum : data) {
            sortedData.add(datum);
        }

        Collections.sort(sortedData);

        if(percentile == 0.0){
            return sortedData.get(0);
        }

        int index = (int)Math.ceil((percentile / 100.0) * sortedData.size());

        return sortedData.get(index - 1);
    }

    public static int[] constructNzArray(double[] wavelet){
        int[] nzArray = new int[wavelet.length];

        // Populate nzArray with zero
        for(int i = 0; i < wavelet.length; i++){
            nzArray[i] = 0;
        }

        // Calculate nzArray for each subtree root, except c0
        for(int index = wavelet.length - 1; index >= 1; index--){
            if(wavelet[index] != 0.0 || nzArray[index] > 0) {
                if(wavelet[index] != 0.0)
                    nzArray[index] = nzArray[index] + 1;

                if((double)index % 2.0 == 0.0)
                    nzArray[index / 2] = nzArray[index / 2] + nzArray[index];
                else
                    nzArray[(index - 1) / 2] = nzArray[(index - 1) / 2] + nzArray[index];
            }
        }

        // Calculate non-zero coefficients rooted at c0
        if(wavelet[0] != 0.0)
            nzArray[0] = nzArray[1] + 1;
        else
            nzArray[0] = nzArray[1];

        return nzArray;
    }

    /**
     * This method will call the applyPerturbationRule method, the resulting wavelet
     * then can be used to apply the MinRelVal algorithm.
     * @param wavelet the wavelet
     * @param nzArray an array that consists of the number of nonzero coefficients rooted at each subtree
     * @param data the data
     * @param percentile percentile
     * @return An array of 'form' of each node
     */
    public static double[] perturbAndCalcNorm(double[] wavelet, int[] nzArray, double[] data, double percentile){
        double nthPercentile = findPercentile(data, percentile);
        double[] minValues = applyPerturbationRule(wavelet, nzArray, data, nthPercentile);
        double[] norm = new double[wavelet.length];

        for(int i = 0; i < wavelet.length; i++){
            norm[i] = Math.max(Math.pow(minValues[i], 2.0), Math.pow(nthPercentile, 2.0));
        }

        return norm;
    }

    public static double getOptimalNSE(double[] wavelet, double B, int q, int root,
                                       int[] nzArray, double[] norm,
                                       double[][] mValues, double[][] yValues, double[][] leftAllot,
                                       boolean[][] checked){

        int indexB = (int)Math.rint(B * q);

        if(root > wavelet.length - 1 || nzArray[root] <= B)
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
            if(wavelet[root] == 0){
                rootLeft = 0.0;
                rootRight = 0.0;
                rootSpace = 0.0;
            }
            else if(2 * root > wavelet.length - 1){
                rootLeft = (q - l) * Math.pow(wavelet[root], 2.0) / (l * norm[root]);
                rootRight = rootLeft;
                rootSpace = (double)l / (double)q;
            }
            else{
                rootLeft = (q - l) * Math.pow(wavelet[root], 2.0) / (l * norm[2 * root]);
                rootRight = (q - l) * Math.pow(wavelet[root], 2.0) / (l * norm[2 * root + 1]);
                rootSpace = (double)l / (double)q;
            }

            for(double _b = 0; _b <= B - rootSpace; _b = _b + (1.0 / (double)q)){
                double left = getOptimalNSE(wavelet, _b, q, 2 * root,
                        nzArray, norm, mValues, yValues, leftAllot, checked);
                double right = getOptimalNSE(wavelet, B - rootSpace - _b, q, 2 * root + 1,
                        nzArray, norm, mValues, yValues, leftAllot, checked);
                double max = Math.max(rootLeft + left, rootRight + right);
                if(max < mValues[root][indexB]){
                    mValues[root][indexB] = max;
                    yValues[root][indexB] = rootSpace;
                    leftAllot[root][indexB] = _b;
                }
            }
            if(wavelet[root] == 0){
                break;
            }
        }
        checked[root][indexB] = true;
        return mValues[root][indexB];
    }

    public static void callMainFunction(String pathWavelet, double b, int q, String pathData, double percentile)
            throws IOException {

        double[] wavelet = OneDHWT.fileToArrayOfDoubles(pathWavelet);
        double[] data = OneDHWT.fileToArrayOfDoubles(pathData);
        int[] nzArray = constructNzArray(wavelet);
        int length_1 = wavelet.length;
        int length_2 = (int)Math.rint(b * q) + 1;
        double[] norm = perturbAndCalcNorm(wavelet, nzArray, data, percentile);
        double[][] mValues = new double[length_1][length_2];
        double[][] yValues = new double[length_1][length_2];
        double[][] leftAllot = new double[length_1][length_2];
        boolean[][] checked = new boolean[length_1][length_2];

        double root;
        double rootSpace;

        mValues[0][(int)Math.rint(b * q)] = Double.POSITIVE_INFINITY;
        for(int l = 1; l <= q; l++){
            if(wavelet[0] == 0.0){
                root = 0;
                rootSpace = 0;
            }
            else{
                root = (q - l) * Math.pow(wavelet[0], 2.0) / (l * norm[1]);
                rootSpace = (double)l / (double)q;
            }

            for(double _b = 0; _b <= b - rootSpace; _b = _b +  (1.0 / (double)q)){
                double next = getOptimalNSE(wavelet, _b, q, 1,
                        nzArray, norm, mValues, yValues, leftAllot, checked);
                if(root + next < mValues[0][(int)Math.rint(b * q)]){
                    mValues[0][(int)Math.rint(b * q)] = root + next;
                    yValues[0][(int)Math.rint(b * q)]= rootSpace;
                    leftAllot[0][(int)Math.rint(b * q)] = _b;
                }
            }

            if(wavelet[0] == 0.0)
                break;
        }

        checked[0][(int)Math.rint(b * q)] = true;

        // Print on the console
        System.out.println();
        for(int j = 0; j < length_2; j++) {
            System.out.println(leftAllot[0][j]);
        }
    }
}
