import java.util.*;
import cern.colt.matrix.impl.SparseObjectMatrix2D;

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
     * perturbation.
     *
     * @param wavelet the wavelet
     * @param nzArray an array that consists of the number of nonzero coefficients rooted at each subtree
     * @param data the original data
     * @param nthPercentile Nth percentile value
     * @return minimum data values rooted at each subtree
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
     * @return An array of 'norm' of each node
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

    public static double getOptimalNSE2(double[] wavelet, double B, int q, int root,
                                       int[] nzArray, double[] norm,
                                       HashMap<Coordinates, MinNSE> minNSEMap){

        int indexB = (int)Math.rint(B * q);

        if(root > wavelet.length - 1 || B - (double)nzArray[root] >= 1e-6)
            return 0;

        if(nzArray[root] > indexB)
            return Double.POSITIVE_INFINITY;

        Coordinates coordinate = new Coordinates(root, indexB);
        if(minNSEMap.containsKey(coordinate)) {
            return minNSEMap.get(coordinate).getmValue();
        }

        System.out.println("Calculate " + root + ", " + (int)Math.rint(B * q));
        MinNSE storeObj = new MinNSE(Double.POSITIVE_INFINITY, 0.0, 0.0);
        minNSEMap.put(coordinate, storeObj);

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

            for(int _b = 0; _b <= indexB - l; _b++){
                double __b = (double)_b / (double) q;
                double left = getOptimalNSE2(wavelet, __b, q, 2 * root, nzArray, norm, minNSEMap);
                double right = getOptimalNSE2(wavelet, B - rootSpace - __b, q, 2 * root + 1,
                        nzArray, norm, minNSEMap);
                double max = Math.max(rootLeft + left, rootRight + right);
                if(max < minNSEMap.get(coordinate).getmValue()){
                    minNSEMap.get(coordinate).setmValue(max);
                    minNSEMap.get(coordinate).setyValue(rootSpace);
                    minNSEMap.get(coordinate).setLeftAllot(__b);
                }
            }

            if(wavelet[root] == 0)
                break;
        }

        return minNSEMap.get(coordinate).getmValue();
    }

    public static void callMainFunction2(double[] wavelet, double b, int q, double[] data, double percentile){
        System.out.println("Construct nzArray...");
        int[] nzArray = constructNzArray(wavelet);
        System.out.println("Perform perturbation rule and calculate norm(i)...");
        double[] norm = perturbAndCalcNorm(wavelet, nzArray, data, percentile);

        System.out.println("Build mapping...");
        HashMap<Coordinates, MinNSE> minNSEMap = new HashMap<>();

        double root;
        double rootSpace;

        Coordinates coordinate = new Coordinates(0, (int)b * q);
        MinNSE storeObj = new MinNSE(Double.POSITIVE_INFINITY, 0.0, 0.0);
        minNSEMap.put(coordinate, storeObj);

        for(int l = 1; l <= q; l++){
            if(wavelet[0] == 0.0){
                root = 0.0;
                rootSpace = 0.0;
            }
            else{
                root = (q - l) * Math.pow(wavelet[0], 2.0) / (l * norm[1]);
                rootSpace = (double)l / (double)q;
            }

            for(int _b = 0; _b <= (int)Math.rint(b * q) - l; _b++){
                double __b = (double)_b / (double)q;
                double next = getOptimalNSE2(wavelet, __b, q, 1, nzArray, norm, minNSEMap);
                if(root + next < minNSEMap.get(coordinate).getmValue()){
                    minNSEMap.get(coordinate).setmValue(root + next);
                    minNSEMap.get(coordinate).setyValue(rootSpace);
                    minNSEMap.get(coordinate).setLeftAllot(__b);
                }
            }

            if(wavelet[0] == 0.0)
                break;
        }

        // Put all optimal y values for each coefficient in an array
        System.out.println("Build probability array...");
        double[] chosenY = new double[wavelet.length];
        double[] bValue = new double[wavelet.length];

        bValue[0] = b;

        chosenY[0] = minNSEMap.get(coordinate).getyValue();
        bValue[1] = minNSEMap.get(coordinate).getLeftAllot();


        for(int i = 1; i < chosenY.length; i++){
            Coordinates position = new Coordinates(i, (int)Math.rint(bValue[i] * q));
            if(minNSEMap.containsKey(position))
                chosenY[i] = minNSEMap.get(position).getyValue();
            if(i * 2 < chosenY.length){
                if(minNSEMap.containsKey(position)) {
                    bValue[i * 2] = minNSEMap.get(position).getLeftAllot();
                    bValue[i * 2 + 1] = bValue[i] - bValue[i * 2] - chosenY[i];
                }
            }
        }

        System.out.println();
        for (double v : chosenY) {
            System.out.println(v);
        }

        performCoinFlips(wavelet, chosenY);
    }

    public static double getOptimalNSE(double[] wavelet, double B, int q, int root,
                                       int[] nzArray, double[] norm,
                                       double[][] mValues, double[][] yValues, double[][] leftAllot,
                                       boolean[][] checked){

        int indexB = (int)Math.rint(B * q);

        if(root > wavelet.length - 1 || B - (double)nzArray[root] >= 1e-6) { // was  { nzArray[root] <= B }
            if(root == 526)
                System.out.println("root 526 " + B + " " + nzArray[root] + " enough space");
            return 0;
        }

        if(nzArray[root] > B * q) {
            if (root == 526)
                System.out.println("root 526 " + B + " not enough space");
            return Double.POSITIVE_INFINITY;
        }

        if(root == 526)
            System.out.println("root 526 " + B + " calculate NSE");
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

            for(int _b = 0; _b <= indexB - l; _b++){
                double __b = (double)_b / (double)q;
                double left = getOptimalNSE(wavelet, __b, q, 2 * root,
                        nzArray, norm, mValues, yValues, leftAllot, checked);
                double right = getOptimalNSE(wavelet, B - rootSpace - __b, q, 2 * root + 1,
                        nzArray, norm, mValues, yValues, leftAllot, checked);
                double max = Math.max(rootLeft + left, rootRight + right);
                if(max < mValues[root][indexB]){
                    mValues[root][indexB] = max;
                    yValues[root][indexB] = rootSpace;
                    leftAllot[root][indexB] = __b;
                }
            }

            if(wavelet[root] == 0){
                break;
            }
        }
        checked[root][indexB] = true;
        return mValues[root][indexB];
    }

    public static void callMainFunction(double[] wavelet, double b, int q, double[] data, double percentile){

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

        mValues[0][(int)(b * q)] = Double.POSITIVE_INFINITY;

        for(int l = 1; l <= q; l++){
            if(wavelet[0] == 0.0){
                root = 0;
                rootSpace = 0;
            }
            else{
                root = (q - l) * Math.pow(wavelet[0], 2.0) / (l * norm[1]);
                rootSpace = (double)l / (double)q;
            }

            for(int _b = 0; _b <= (int)Math.rint(b * q) - l; _b++){
                double __b = (double)_b / (double)q;
                double next = getOptimalNSE(wavelet, __b, q, 1,
                        nzArray, norm, mValues, yValues, leftAllot, checked);
                if(root + next < mValues[0][(int)Math.rint(b * q)]){
                    mValues[0][(int)Math.rint(b * q)] = root + next;
                    yValues[0][(int)Math.rint(b * q)]= rootSpace;
                    leftAllot[0][(int)Math.rint(b * q)] = __b;
                }
            }

            if(wavelet[0] == 0.0)
                break;
        }

        checked[0][(int)Math.rint(b * q)] = true;

        // Put all optimal y values for each coefficient in an array
        double[] chosenY = new double[length_1];
        double[] bValue = new double[length_1];

        for(int i = 0; i < length_2; i++){
            //System.out.print(i + " ");
            for(int j = 0; j < length_1; j++){
                System.out.print(mValues[j][i] + "\t");
            }
            System.out.println();
        }

        bValue[0] = b;

        chosenY[0] = yValues[0][(int)(b * q)];
        bValue[1] = leftAllot[0][(int)(b * q)];

        for(int i = 1; i < length_1; i++){
            chosenY[i] = yValues[i][(int)(Math.rint(bValue[i] * q))];
            if(i * 2 < length_1) {
                bValue[i * 2] = leftAllot[i][(int)(Math.rint(bValue[i] * q))];
                bValue[i * 2 + 1] = bValue[i] - bValue[i * 2] - chosenY[i];
            }
        }

        double sum = 0;
        for(int i = 0; i < length_1; i++){
            System.out.println(chosenY[i]);
            sum += chosenY[i];
        }

        performCoinFlips(wavelet, chosenY);
    }

    public static void performCoinFlips(double[] wavelet, double[] chosenY){
        for(int i = 0; i < wavelet.length; i++){
            if(roundUp(chosenY[i]))
                wavelet[i] = wavelet[i] / chosenY[i];

            else
                wavelet[i] = 0.0;
        }
    }

    public static boolean roundUp(double yVal){
        double randomDouble = Math.random();

        if(yVal == 0.0)
            return false;

        return randomDouble <= yVal;
    }

    public static double getOptimalNSE3(double[] wavelet, double b, int q, int root,
                                       int[] nzArray, double[] norm,
                                       SparseObjectMatrix2D objectMatrix2D){

        int indexB = (int)Math.rint(b * q);

        if(root > wavelet.length - 1 || b - (double)nzArray[root] >= 0)// was  { nzArray[root] <= B } 1e-6
            return 0;

        if(nzArray[root] > indexB)
            return Double.POSITIVE_INFINITY;

        if(!(objectMatrix2D.getQuick(indexB, root) == null))
            return ((MinNSE)objectMatrix2D.getQuick(indexB, root)).getmValue();

        MinNSE obj = new MinNSE(Double.POSITIVE_INFINITY, 0, 0);
        objectMatrix2D.setQuick(indexB, root, obj);

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

            for(int _b = 0; _b <= indexB - (int)Math.rint(rootSpace * q); _b++){
                double __b = (double)_b / (double)q;
                double left = getOptimalNSE3(wavelet, __b, q, 2 * root,
                        nzArray, norm, objectMatrix2D);
                double right = getOptimalNSE3(wavelet, b - rootSpace - __b, q, 2 * root + 1,
                        nzArray, norm, objectMatrix2D);
                double max = Math.max(rootLeft + left, rootRight + right);
                if(max < ((MinNSE)objectMatrix2D.getQuick(indexB, root)).getmValue()){
                    MinNSE newObj = new MinNSE(max, rootSpace, __b);
                    objectMatrix2D.setQuick(indexB, root, newObj);
                }
            }

            if(wavelet[root] == 0){
                break;
            }
        }
        return ((MinNSE)objectMatrix2D.getQuick(indexB, root)).getmValue();
    }

    public static void callMainFunction3(double[] wavelet, double b, int q, double[] data, double percentile){

        int[] nzArray = constructNzArray(wavelet);
        int length_1 = wavelet.length;                                                      // column (number of roots)
        int length_2 = (int)Math.rint(b * q) + 1;                                           // row
        double[] norm = perturbAndCalcNorm(wavelet, nzArray, data, percentile);
        SparseObjectMatrix2D objectMatrix2D = new SparseObjectMatrix2D(length_2, length_1);

        double root;
        double rootSpace;

        int row = (int)b * q;
        int column = 0;
        MinNSE obj = new MinNSE(Double.POSITIVE_INFINITY, 0, 0);
        objectMatrix2D.setQuick(row, column, obj);

        for(int l = 1; l <= q; l++){
            if(wavelet[0] == 0.0){
                root = 0;
                rootSpace = 0;
            }
            else{
                root = (q - l) * Math.pow(wavelet[0], 2.0) / (l * norm[1]);
                rootSpace = (double)l / (double)q;
            }

            for(int _b = 0; _b <= (int)Math.rint(b * q) - (int)Math.rint(rootSpace * q); _b++){
                double __b = (double)_b / (double)q;
                double next = getOptimalNSE3(wavelet, __b, q, 1,
                        nzArray, norm, objectMatrix2D);
                if(root + next < ((MinNSE)objectMatrix2D.getQuick(row, column)).getmValue()){
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


        for(int i = 1; i < length_1; i++){
            int ithRow = (int)(Math.rint(bValue[i] * q));
            if(objectMatrix2D.getQuick(ithRow, i) == null) {
                chosenY[i] = 0;
            }
            else
                chosenY[i] = ((MinNSE)objectMatrix2D.getQuick(ithRow, i)).getyValue();

            if(i * 2 < length_1) {
                if(objectMatrix2D.getQuick(ithRow, i) == null)
                    bValue[i * 2] = 0;
                else
                    bValue[i * 2] = ((MinNSE)objectMatrix2D.getQuick(ithRow, i)).getLeftAllot();
                bValue[i * 2 + 1] = bValue[i] - bValue[i * 2] - chosenY[i];
            }
        }

        double sum = 0;
        for(int i = 0; i < length_1; i++){
//            if (i < 1000) {
//                System.out.println(i + " " + chosenY[i]);
//            }
            sum += chosenY[i];
        }

        performCoinFlips(wavelet, chosenY);
        System.out.println("sum " + sum);
        System.out.println(nzArray[0]);
    }
}
