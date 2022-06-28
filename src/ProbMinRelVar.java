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
    public static void applyPerturbationRule(double[] wavelet, double[] nzArray){

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
