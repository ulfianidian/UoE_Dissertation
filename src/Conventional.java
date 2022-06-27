import java.util.Arrays;
import java.util.HashMap;

/**
 * @author Ulfiani Primawati
 *
 * This is an implementation of the deterministic wavelet thresholding technique.
 * Also known as the Conventional method. This method will greedily retain N largest
 * normalized coefficients of a wavelet. The compact wavelet then can be used to
 * reconstruct the data or to answer range aggregate queries such as range-sum and range-average.
 */
public class Conventional {

    /**
     * This method performs normalization of the wavelet using ci* = ci / sqrt(2^level(ci)) formula
     * in order to equalize the importance of each coefficient in reconstructing the data.
     * @param data The wavelet to be normalized
     * @return an array of normalized coefficients
     */
    public static double[] normalizeCoeffs(double[] data){
        double[] normalizedWavelet = new double[data.length];
        int level = 0;
        int limit = level + 1;
        for(int i = 0; i < normalizedWavelet.length; i++){
            if (i >= (int)Math.pow(2.0, limit)){
                level++;
                limit = level + 1;
            }
            normalizedWavelet[i] = data[i] / (Math.sqrt(Math.pow(2.0, level)));
        }
        return normalizedWavelet;
    }

    /**
     * This method retains b largest absolute normalized wavelet coefficients,
     * but return them in their original value.
     *
     * Example:
     * b = 8
     * Wavelet coefficients: 65, 0, 14, -15, 20, -20, 21, -21, 28, 28, 28, -28, 29, -29, -29, -29
     * Normalized coefficients: 65, 0, 14/sqrt(2), -15/sqrt (2), 10, -10, 10.5, -10.5, 14/sqrt(2),
     * 14/sqrt(2), 14/sqrt(2), -14/sqrt(2), 29/(2*sqrt(2)), -29/(2*sqrt(2)), -29/(2*sqrt(2)), -29/(2*sqrt(2))
     * Retained coefficients: 65, 0, 0, -15, 0, 0, 21, -21, 0, 0, 0, 0, 29, -29, -29, -29
     *
     * @param wavelet The wavelet before normalization
     * @param normalizedCoeffs Normalized wavelet coefficients
     * @param b the number of coefficients to be retained
     */
    public static void retainNCoeffs(double[] wavelet, double[] normalizedCoeffs, int b){
        // Populate topNValues, topNAbsNormValues, and topNIndices with the first N values, the first
        // N absolute normalized coefficients, and their corresponding indices.
        double[] topNValues = Arrays.copyOfRange(wavelet, 0, b);

        double[] topNAbsNormValues = new double[b];
        for(int k = 0; k < b; k++){
            topNAbsNormValues[k] = Math.abs(normalizedCoeffs[k]);
        }

        int[] topNIndices = new int[b];
        for(int i = 0; i < b; i++){
            topNIndices[i] = i;
        }

        // Iterate through the wavelet to get N largest coefficients.
        double nextVal;
        boolean smallestValChecked = false;

        int index = 0;
        double min = topNAbsNormValues[index];

        for(int j = b; j < wavelet.length; j++){
            nextVal = Math.abs(normalizedCoeffs[j]);

            while(!smallestValChecked) {
                for (int a = 1; a < b; a++) {
                    if (topNAbsNormValues[a] <= min) {
                        min = topNAbsNormValues[a];
                        index = a;
                    }
                }
                smallestValChecked = true;
            }

            if(nextVal > topNAbsNormValues[index]){
                topNAbsNormValues[index] = nextVal;
                topNValues[index] = wavelet[j];
                topNIndices[index] = j;
                smallestValChecked = false;
                index = 0;
                min = topNAbsNormValues[index];
            }
        }

        HashMap<Integer, Double> retainedCoeffs = storeToHashMap(topNIndices, topNValues);
        populateWavelet(retainedCoeffs, wavelet.length);
    }

    private static HashMap<Integer, Double> storeToHashMap(int[] topNIndices, double[] topNValues){
        HashMap<Integer, Double> retainedCoeffs = new HashMap<>();
        for(int i = 0; i < topNIndices.length; i++){
            retainedCoeffs.put(topNIndices[i], topNValues[i]);
        }

        return retainedCoeffs;
    }

    public static void populateWavelet(HashMap<Integer, Double> coeffsMapping, int tupleCounts){
        double[] conventionalHWT = new double[tupleCounts];
        for(int i = 0; i < tupleCounts; i++){
            if(coeffsMapping.containsKey(i))
                conventionalHWT[i] = coeffsMapping.get(i);
            else
                conventionalHWT[i] = 0.0;
        }

        //reconstruct all data
        OneDHWT.orderedFastHWTInverse(conventionalHWT);
        System.out.println();
        OneDHWT.displayData(conventionalHWT);

        //Reconstruct one element
        getOneElement(coeffsMapping, tupleCounts, 8);
    }

    /**
     * Not yet finished
     *
     * @param coeffsMapping Coefficients mapping
     * @param tupleCounts The number of tuples
     * @param nthValue Nth value
     */
    public static void getOneElement(HashMap<Integer, Double> coeffsMapping, int tupleCounts, int nthValue){
        final int NUM_OF_LOOPS = (int)(Math.log(tupleCounts) / Math.log(2.0));
        int start = 0;
        int middlePoint = tupleCounts / 2;
        int end = tupleCounts - 1;
        int i = 1;
        int level = 0;
        int increment = 0;
        double sum = coeffsMapping.get(0);

        for(int j = 0; j < NUM_OF_LOOPS; j++){
            System.out.println("start: " + start);
            System.out.println("middle point: " + middlePoint);
            System.out.println(("end: " + end));
            System.out.println("i: " + i);
            System.out.println("level: " + level);
            if(nthValue < middlePoint){
                if(coeffsMapping.containsKey(i)) {
                    sum += coeffsMapping.get(i);
                    System.out.println(coeffsMapping.get(i));
                }
                System.out.println("sum : " + sum);
                if(true){
                    i = i + increment + (int)Math.pow(2.0, (double)level);
                }
                System.out.println("next i: "+ i);
                end = middlePoint - 1;
                middlePoint = start + ( (middlePoint - start) / 2);
            }

            else{
                if(coeffsMapping.containsKey(i)) {
                    sum -= coeffsMapping.get(i);
                    System.out.println(coeffsMapping.get(i));
                }
                System.out.println("sum : " + sum);
                if(true){
                    increment++;
                    i = i + increment + (int)Math.pow(2.0, (double)level);
                }
                System.out.println("next i: "+ i);
                start = middlePoint;
                //end = start + middlePoint - 1;
                middlePoint = start + (end + 1 - start) / 2;
            }
            level++;
            System.out.println();
        }
        System.out.println(sum);
    }

    public static void calculateRangeSum(){

    }

    public static void calculateRangeAvg(){

    }

    /**
     * DO NOT USE
     *
     * @param wavelet the wavelet
     * @param b the number of retained coefficients
     */
    public static void retainNCoefficients(double[] wavelet, double[] normalizedCoeffs, int b){
        // Populate topNValues, topNAbsNormValues, and topNIndices with the first N values, the first
        // N absolute normalized coefficients, and their corresponding indices.
        double[] topNValues = Arrays.copyOfRange(wavelet, 0, b);

        double[] topNAbsNormValues = new double[b];
        for(int k = 0; k < b; k++){
            topNAbsNormValues[k] = Math.abs(normalizedCoeffs[k]);
        }

        int[] topNIndices = new int[b];
        for(int i = 0; i < b; i++){
            topNIndices[i] = i;
        }

        // Iterate through the wavelet to get N largest coefficients.
        double nextVal;
        for(int j = b; j < wavelet.length; j++){
            nextVal = Math.abs(normalizedCoeffs[j]);
            boolean isBigger = false;
            int a = 0;
            while( !isBigger && a < b){
                if(nextVal > topNAbsNormValues[a]){
                    isBigger = true;
                    topNAbsNormValues[a] = nextVal;
                    topNValues[a] = wavelet[j];
                    topNIndices[a] = j;
                }
                a++;
            }
        }

        // Print the retained coefficients and their indices on the console
        for(int x = 0; x < b; x++){
            System.out.println("index: " + topNIndices[x] + ", value: " + topNValues[x]);
        }

        HashMap<Integer, Double> retainedCoeffs = storeToHashMap(topNIndices, topNValues);
    }
}
