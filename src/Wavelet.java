import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.*;

public class Wavelet {
    public static void main(String[] args) throws IOException {

        final int Q = 5;
        final int PERCENTILE = 10;

        // Generate zipfian frequency and create the dataset
        if(args[0].equals("generate")){
           if(args[1].equals("sequential")){
                generateZipf(args[1], Integer.parseInt(args[2]), Integer.parseInt(args[3]),
                        Double.parseDouble(args[4]), Integer.parseInt(args[5]), args[6], args[7]);
           }
           else if(args[1].equals("nonconsecutive")){
               generateNonconsecutiveZipf(args[1], Integer.parseInt(args[2]),
                       Integer.parseInt(args[3]), Double.parseDouble(args[4]),
                       Integer.parseInt(args[5]), args[6], args[7], Integer.parseInt(args[8]));
           }
           else if(args[1].equals("random")){
               ZipfianGenerator uniRandom = new ZipfianGenerator(Integer.parseInt(args[2]),
                       Integer.parseInt(args[3]), Integer.parseInt(args[4]), args[5]);
               uniRandom.generateUniformData();
           }
        }

        // Build the full wavelet
        else if(args[0].equals("build") && args[1].equals("wavelet")){
            double[] wavelet = OneDHWT.orderedFastHWT(args[2]);
            OneDHWT.saveDataToFile(wavelet, args[3]);
        }

        else if(args[0].equals("wavelet")){
            double[] wavelet = OneDHWT.orderedFastHWT(args[4]);

            // Perform conventional thresholding and save the resulting wavelet coefficients in a file
            if(args[1].equals("conventional")){
                // Perform normalization
                double[] normalizedCoeffs = Conventional.normalizeCoeffs(wavelet);

                // Retain B largest coefficients
                HashMap<Integer, Double> retainedCoeffs = Conventional.retainNCoeffs(wavelet,
                        normalizedCoeffs, Integer.parseInt(args[2]));

                if(args[3].equals("full")){
                    OneDHWT.saveDataToFile(Conventional.populateWavelet(retainedCoeffs, wavelet.length), args[5]);
                }
                else if(args[3].equals("summary")){
                    int i = 0;
                    String[] mapping = new String[retainedCoeffs.size()];

                    for(Map.Entry<Integer, Double> entry : retainedCoeffs.entrySet()){
                        mapping[i] = entry.getKey() + ", " + entry.getValue();
                        i++;
                    }

                    OneDHWT.saveDataToFile(mapping, args[5]);
                }
            }

            // Perform Probabilistic MinL2 thresholding and save the resulting wavelet coefficients in a file
            else if(args[1].equals("minL2")){
                ProbabilisticMinL2.waveletMinL2(wavelet, Integer.parseInt(args[2]));
                if(args[3].equals("full")){
                    OneDHWT.saveDataToFile(wavelet, args[5]);
                }
                else if(args[3].equals("summary")){
                    // Build summary
                }
            }

            // Compute the probability distribution using minRelVar algorithm,
            // Then save the probability distribution in a file
            else if(args[1].equals("minRelVar")){
                double[] probability = ProbMinRelVar.callMainFunction2(wavelet, Integer.parseInt(args[2]), Q,
                        OneDHWT.fileToArrayOfDoubles(args[4]), PERCENTILE);
                OneDHWT.saveDataToFile(probability, args[5]);
            }

            // Compute the probability distribution using minRelBias algorithm
            // Then save the probability distribution into a file
            else if(args[1].equals("minRelBias")){
                double[] probability = ProbMinRelBias.callMainFunction1(wavelet, Integer.parseInt(args[2]), Q,
                        OneDHWT.fileToArrayOfDoubles(args[4]), PERCENTILE);
                OneDHWT.saveDataToFile(probability, args[5]);
            }
        }

        else if(args[0].equals("reconstruct")){
            double[] data = OneDHWT.orderedFastHWTInverse(args[3]);
//            double[] wavelet = OneDHWT.fileToArrayOfDoubles(args[3]);
//            double[] data = new double[wavelet.length];
//            for(int i = 0; i < wavelet.length; i++)
//                data[i] = getOneElement(wavelet, wavelet.length, i);

            getMeanRelError(data, args[5], PERCENTILE);
            getMaxRelError(data, args[5], PERCENTILE);
            getPercentileRelativeError(data, args[5], PERCENTILE);
            OneDHWT.saveDataToFile(data, args[4]);
        }

        else if(args[0].equals("flip-coin")){
            double[] probability = OneDHWT.fileToArrayOfDoubles(args[2]);
            double[] wavelet = OneDHWT.orderedFastHWT(args[3]);
            if(args[1].equals("MinRelVar")){
                ProbMinRelVar.performCoinFlips(wavelet, probability);
                OneDHWT.saveDataToFile(wavelet, args[4]);
            }
            else if (args[1].equals("MinRelBias")){
                ProbMinRelBias.performCoinFlips(wavelet, probability);
                OneDHWT.saveDataToFile(wavelet, args[4]);
            }
        }

        else if(args[0].equals("get")){

            if(args[1].equals("mean-relative-error")){
                getMeanRelError(args[2], args[3], PERCENTILE);
            }
            else if(args[1].equals("max-relative-error")){
                getMaxRelError(args[2], args[3], PERCENTILE);
            }
            else if(args[1].equals("25-percentile-relative-error")){
                getPercentileRelativeError(args[2], args[3], PERCENTILE);
            }
        }

        else if(args[0].equals("sort")){
            FileInputStream fileStream = new FileInputStream(args[1]);
            BufferedReader bufferedReader = new BufferedReader(new InputStreamReader(fileStream));
            ArrayList<Double> arrayOfDoubles = new ArrayList<>();

            //Read the file, line by line
            String stringLine;
            while((stringLine = bufferedReader.readLine()) != null){
                //Print the value on the console
                arrayOfDoubles.add(Double.valueOf(stringLine));
            }

            //Close the input stream
            bufferedReader.close();

            Collections.sort(arrayOfDoubles);

            double[] doublePrims = new double[arrayOfDoubles.size()];
            int i = 0;
            for(Double d: arrayOfDoubles){
                doublePrims[i++] = d.doubleValue();
            }

            arrayOfDoubles = null;

            OneDHWT.saveDataToFile(doublePrims, args[2]);
        }

        else if(args[0].equals("flip-and-reconstruct")){
            if(args[1].equals("MinRelVar")){
                double[] probability = OneDHWT.fileToArrayOfDoubles(args[2]);
                for(int i = 4; i < 9; i++) {
                    double[] wavelet = OneDHWT.orderedFastHWT(args[3]);
                    ProbMinRelVar.performCoinFlips(wavelet, probability);
                    OneDHWT.saveDataToFile(wavelet, args[i]);
                    OneDHWT.orderedFastHWTInverse(wavelet);
                    getMeanRelError(wavelet, args[3], PERCENTILE);
                    getMaxRelError(wavelet, args[3], PERCENTILE);
                    getPercentileRelativeError(wavelet, args[3], PERCENTILE);
                    System.out.println();
                }
            }
            else if(args[1].equals("MinRelBias")){
                double[] probability = OneDHWT.fileToArrayOfDoubles(args[2]);
                for(int i = 4; i < 9; i++) {
                    double[] wavelet = OneDHWT.orderedFastHWT(args[3]);
                    ProbMinRelBias.performCoinFlips(wavelet, probability);
                    OneDHWT.saveDataToFile(wavelet, args[i]);
                    OneDHWT.orderedFastHWTInverse(wavelet);
                    getMeanRelError(wavelet, args[3], PERCENTILE);
                    getMaxRelError(wavelet, args[3], PERCENTILE);
                    getPercentileRelativeError(wavelet, args[3], PERCENTILE);
                    System.out.println();
                }
            }
            else if(args[1].equals("MinL2")){
                for(int i = 4; i < 9; i++) {
                    double[] wavelet = OneDHWT.orderedFastHWT(args[3]);
                    ProbabilisticMinL2.waveletMinL2(wavelet, Integer.parseInt(args[2]));
                    OneDHWT.saveDataToFile(wavelet, args[i]);
                    OneDHWT.orderedFastHWTInverse(wavelet);
                    getMeanRelError(wavelet, args[3], PERCENTILE);
                    getMaxRelError(wavelet, args[3], PERCENTILE);
                    getPercentileRelativeError(wavelet, args[3], PERCENTILE);
                    System.out.println();
                }
            }
        }

        else if(args[0].equals("wavelet-recursive")){

            if(args[1].equals("conventional")){
                int i = 3;      // the first argument that stores the number of coefficients to be retained
                int j = 4;      // the first argument that stores the destination path for each b
                while(i < args.length){
                    buildConventional(args[2], args[i], args[j]);
                    i = i + 2;
                    j = j + 2;
                }
            }

            else if(args[1].equals("minRelVar")){
                int i = 3;      // the first argument that stores the number of coefficients to be retained
                int j = 4;      // the first argument that stores the destination path for each b
                while(i < args.length){
                    System.out.println("Calculate probability distribution for b: " + args[i]);
                    buildMinRelVar(args[2], args[i], Q, PERCENTILE, args[j]);
                    i = i + 2;
                    j = j + 2;
                    System.out.println();
                }
            }

            else if(args[1].equals("minRelBias")){
                int i = 3;      // the first argument that stores the number of coefficients to be retained
                int j = 4;      // the first argument that stores the destination path for each b
                while(i < args.length){
                    buildMinRelBias(args[2], args[i], Q, PERCENTILE, args[j]);
                    i = i + 2;
                    j = j + 2;
                }
            }
        }

        else if(args[0].equals("get-multiple-distribution")){
            if(args[1].equals("minRelVar")){
                int i = 2;
                int j = 3;
                int k = 4;
                while(i < args.length){
                    System.out.println("Calculate probability distribution for b: " + args[j]);
                    buildMinRelVar(args[i], args[j], Q, PERCENTILE, args[k]);
                    System.out.println("Destination path: " + k);
                    i = i + 3;
                    j = j + 3;
                    k = k + 3;
                    System.out.println();
                }
            }
            else if(args[1].equals("minRelBias")){
                int i = 2;
                int j = 3;
                int k = 4;
                while(i < args.length){
                    System.out.println("Calculate probability distribution for b: " + args[j]);
                    buildMinRelBias(args[i], args[j], Q, PERCENTILE, args[k]);
                    System.out.println("Destination path: " + k);
                    i = i + 3;
                    j = j + 3;
                    k = k + 3;
                    System.out.println();
                }
            }
        }

        else if(args[0].equals("print-full-to-summary")){
            for(int i = 1; i < args.length; i++){
                double[] wavelet = OneDHWT.fileToArrayOfDoubles(args[i]);
                HashMap<Integer, Double> summary = OneDHWT.fullToSummary(wavelet);
                summary.entrySet().forEach(entry->{
                    System.out.println(entry.getKey() + "\t" + entry.getValue());
                });
            }
            System.out.println();
        }

        else{
            System.out.println("Enter a valid command!");
        }
    }

    public static void generateZipf(String mode, int keySpace, int totalRecords, double zParam,
                                    int start, String permutation, String writeTo) throws IOException {
        ZipfianGenerator zipfianGenerator = new ZipfianGenerator(mode, keySpace, totalRecords, zParam,
                start, permutation, writeTo);
        zipfianGenerator.generateZipf();
    }

    public static void generateNonconsecutiveZipf(String mode, int keySpace, int totalRecords, double zParam,
                                    int start, String permutation, String writeTo, int gap) throws IOException {
        ZipfianGenerator zipfianGenerator = new ZipfianGenerator(mode, keySpace, totalRecords, zParam,
                start, permutation, writeTo, gap);
        zipfianGenerator.generateZipf();
    }

    public static void getMeanRelError(String dataPath, String originalDataPath, double percentile) throws IOException{
        double[] data = OneDHWT.orderedFastHWTInverse(dataPath);
        double[] originalData = OneDHWT.fileToArrayOfDoubles(originalDataPath);
        double sanityBound = findPercentile(originalData, percentile);

        double sum = 0.0;
        double curr;
        for(int i = 0; i < data.length; i++){
            curr = Math.abs(data[i] - originalData[i])/Math.max(originalData[i], sanityBound);
            sum = sum + curr;
        }

        double error = sum / data.length;
        System.out.println("Mean relative error: " + error);
    }

    public static void getMeanRelError(double[] data, String originalDataPath, double percentile) throws IOException{
        double[] originalData = OneDHWT.fileToArrayOfDoubles(originalDataPath);
        double sanityBound = findPercentile(originalData, percentile);

        double sum = 0.0;
        double curr;
        for(int i = 0; i < data.length; i++){
            curr = Math.abs(data[i] - originalData[i])/Math.max(originalData[i], sanityBound);
            sum = sum + curr;
        }

        double error = sum / data.length;
        System.out.println("Mean relative error: " + error);
    }

    public static void getMaxRelError(String dataPath, String originalDataPath, double percentile) throws IOException{
        double[] data = OneDHWT.orderedFastHWTInverse(dataPath);
        double[] originalData = OneDHWT.fileToArrayOfDoubles(originalDataPath);
        double sanityBound = findPercentile(originalData, percentile);

        double max = 0.0;
        double curr = 0.0;

        for(int i = 0; i < data.length; i++){
            curr = Math.abs(data[i] - originalData[i])/Math.max(originalData[i], sanityBound);
            if(curr > max){
                max = curr;
            }
        }

        System.out.println("Maximum relative error: " + max);
    }

    public static void getMaxRelError(double[] data, String originalDataPath, double percentile) throws IOException{
        double[] originalData = OneDHWT.fileToArrayOfDoubles(originalDataPath);
        double sanityBound = findPercentile(originalData, percentile);

        double max = 0.0;
        double curr = 0.0;

        for(int i = 0; i < data.length; i++){
            curr = Math.abs(data[i] - originalData[i])/Math.max(originalData[i], sanityBound);
            if(curr > max){
                max = curr;
            }
        }

        System.out.println("Maximum relative error: " + max);
    }

    public static void getPercentileRelativeError(String dataPath, String originalDataPath, double percentile)
            throws IOException{
        double[] data = OneDHWT.orderedFastHWTInverse(dataPath);
        double[] originalData = OneDHWT.fileToArrayOfDoubles(originalDataPath);
        double sanityBound = findPercentile(originalData, percentile);
        double curr;

        List<Double> error = new ArrayList<>();
        for(int i = 0; i < data.length; i++){
            curr = Math.abs(data[i] - originalData[i])/Math.max(originalData[i], sanityBound);
            error.add(curr);
        }

        Collections.sort(error);

        if(percentile == 0.0){
            System.out.println("25-percentile relative error: " + error.get(0));
        }

        int index = (int)Math.ceil((75 / 100.0) * error.size());

        System.out.println("25-percentile relative error: " + error.get(index - 1));
    }

    public static void getPercentileRelativeError(double[] data, String originalDataPath, double percentile)
            throws IOException{
        double[] originalData = OneDHWT.fileToArrayOfDoubles(originalDataPath);
        double sanityBound = findPercentile(originalData, percentile);
        double curr;

        List<Double> error = new ArrayList<>();
        for(int i = 0; i < data.length; i++){
            curr = Math.abs(data[i] - originalData[i])/Math.max(originalData[i], sanityBound);
            error.add(curr);
        }

        Collections.sort(error);

        if(percentile == 0.0){
            System.out.println("25-percentile relative error: " + error.get(0));
        }

        int index = (int)Math.ceil((75 / 100.0) * error.size());

        System.out.println("25-percentile relative error: " + error.get(index - 1));
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

    public static void OneDWaveletDecomp(String filePath, String writeTo) throws IOException {
        double[] wavelet = OneDHWT.orderedFastHWT2(filePath);

        //Print the wavelet on the console
        for (double w : wavelet) {
            System.out.println(w);
        }

        OneDHWT.saveDataToFile(wavelet, writeTo);
    }

    public static void ReconstructData(String filePath, String writeTo) throws IOException {
        double[] data = OneDHWT.orderedFastHWTInverse(filePath);
        for (double datum : data) {
            System.out.println(datum);
        }
        OneDHWT.saveDataToFile(data, writeTo);
    }

    public static double getOneElement(double[] wavelet, int tupleCounts, int nthValue){
        final int NUM_OF_LOOPS = (int)(Math.log(tupleCounts) / Math.log(2.0));
        int start = 0;
        int middlePoint = tupleCounts / 2;
        int end = tupleCounts - 1;
        int i = 1;
        double sum = wavelet[0];

        for(int j = 0; j < NUM_OF_LOOPS; j++){
            if(nthValue < middlePoint){
                sum += wavelet[i];
                i = 2 * i;
                end = middlePoint - 1;
                middlePoint = start + ( (middlePoint - start) / 2);
            }

            else{
                sum -= wavelet[i];
                i = 2 * i + 1;
                start = middlePoint;
                middlePoint = start + (end + 1 - start) / 2;
            }
        }
        return sum;
    }

    /**
     * Build conventional wavelet and save the resulting wavelet in a file. All the coefficients
     * are going to be saved in the file, including all the zeros.
     * @param dataPath path to the original data
     * @param b the number of coefficients to be retained
     * @param destinationPath destination path
     * @throws IOException path is not valid
     */
    public static void buildConventional(String dataPath, String b, String destinationPath) throws IOException {
        double[] wavelet = OneDHWT.orderedFastHWT(dataPath);

        // Perform normalization
        double[] normalizedCoeffs = Conventional.normalizeCoeffs(wavelet);

        // Retain B largest coefficients
        HashMap<Integer, Double> retainedCoeffs = Conventional.retainNCoeffs(wavelet,
                normalizedCoeffs, Integer.parseInt(b));

        OneDHWT.saveDataToFile(Conventional.populateWavelet(retainedCoeffs, wavelet.length), destinationPath);
    }

    /**
     * Build minRelVar probability distribution and save the resulting probability
     * distribution in a file. All the probability values are going to be saved in the file,
     * including all the zeros.
     *
     * @param dataPath path to the original data
     * @param b the number of coefficients to be retained
     * @param Q q parameter for calculating the probability distribution
     * @param PERCENTILE percentile
     * @param destinationPath destination path for storing the probability distribution
     * @throws IOException path does not exist
     */
    public static void buildMinRelVar(String dataPath, String b, int Q, int PERCENTILE, String destinationPath)
            throws IOException{
        double[] wavelet = OneDHWT.orderedFastHWT(dataPath);
        double[] probability = ProbMinRelVar.callMainFunction2(wavelet, Integer.parseInt(b), Q,
                OneDHWT.fileToArrayOfDoubles(dataPath), PERCENTILE);
        OneDHWT.saveDataToFile(probability, destinationPath);
    }

    /**
     * Build minRelBias probability distribution and save the resulting probability
     * distribution in a file. All the probability values are going to be saved in the file,
     * including all the zeros.
     *
     * @param dataPath path to the original data
     * @param b the number of coefficients to be retained
     * @param Q q parameter for calculating the probability distribution
     * @param PERCENTILE percentile
     * @param destinationPath destination path for storing the probability distribution
     * @throws IOException path does not exist
     */
    public static void buildMinRelBias(String dataPath, String b, int Q, int PERCENTILE, String destinationPath)
            throws IOException{
        double[] wavelet = OneDHWT.orderedFastHWT(dataPath);
        double[] probability = ProbMinRelBias.callMainFunction1(wavelet, Integer.parseInt(b), Q,
                OneDHWT.fileToArrayOfDoubles(dataPath), PERCENTILE);
        OneDHWT.saveDataToFile(probability, destinationPath);
    }
}
