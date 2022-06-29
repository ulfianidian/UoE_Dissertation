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
    public static void applyPerturbationRule(double[] wavelet, int[] nzArray, double[] data, double percentile){
        double nthPercentile = findPercentile(data, percentile);
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

        // delete
        for(int a = 0; a < wavelet.length; a++){
            System.out.println(nzArray[a]);
        }
        System.out.println();
        for(int a = 0; a < wavelet.length; a++){
            System.out.println(wavelet[a]);
        }
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
    public static double[] calcNorm(){

    }

    public OptimalNSE getOptimalNSE(double[] wavelet, double b, int q, int root){

    }**/
}

class OptimalNSE {
    double value;
    double yValue;
    double leftAllot;
    double b;
    boolean computed;
    int root;

    OptimalNSE(double value, double yValue, double leftAllot, double b, boolean computed, int root){
        this.value = value;
        this.yValue = yValue;
        this.leftAllot = leftAllot;
        this.b = b;
        this.computed = computed;
        this.root = root;
    }
}
