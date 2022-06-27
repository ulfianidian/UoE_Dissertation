import java.awt.geom.PathIterator;
import java.util.Arrays;

/**
 * @author Ulfiani Primawati
 *
 * This is an implementation of the probabilistic thresholding technique
 * that minimized the L2 (mean-squared error) error metric. This algorithm
 * deterministically keeps the most important wavelet coefficients while
 * randomly rounding the other coefficients either down to zero or up to a
 * larger value (called a "rounding value").
 *
 * This algorithm is based on a paper proposed by M. Garofalakis and
 * P. G. Gibbons in 2004, titled "Probabilistic Wavelet Synopses".
 */

public class ProbabilisticMinL2 {

    /**
     * This algorithm determines the rounding values for a given wavelet and
     * the desired b value that minimize the expected mean-squared error.
     *
     * @param wavelet the wavelet
     * @param b available space
     * @return an array of rounding values
     */
    public static double[] calcRoundingValues(double[] wavelet, int b){
        int level = 0;
        int limit = level + 1;
        double total = 0.0;
        double avail = b;
        double[] sqrtk = new double[wavelet.length];
        double[] roundingValues = new double[wavelet.length];

        for(int i = 0; i < wavelet.length; i++){
            if (i >= (int)Math.pow(2.0, limit)){
                level++;
                limit = level + 1;
            }

            sqrtk[i] = (Math.abs(wavelet[i])) / (Math.pow(2.0, (level / 2.0)));
            total = total + sqrtk[i];
        }

        // Sort the indices i in non-increasing order of sqrtk[i].
        Pair[] indexValuePair = new Pair[wavelet.length];
        for(int j = 0; j < sqrtk.length; j++){
            indexValuePair[j] = new Pair(j, sqrtk[j]);
        }
        Arrays.sort(indexValuePair);

        int  r;
        for(r = 0; r < sqrtk.length; r++){
            int idx = indexValuePair[r].getKey();
            if(sqrtk[idx] * avail < total)
                break;
            else{
                roundingValues[idx] = wavelet[idx];
                total = total - sqrtk[idx];
                avail = avail - 1;
            }
        }

        for(int k = r; k < sqrtk.length; k++){
            int idx = indexValuePair[k].getKey();
            if(wavelet[idx] != 0.0){
                roundingValues[idx] = (wavelet[idx] * total) / (sqrtk[idx] * avail);
            }
            else{
                roundingValues[idx] = 0.0; // check whether this is correct
            }
        }

        return roundingValues;
    }
}
