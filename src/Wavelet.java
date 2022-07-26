import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.*;

public class Wavelet {
    public static void main(String[] args) throws IOException {

        final int Q = 10;
        final int PERCENTILE = 10;

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
        }

        else if(args[0].equals("build") && args[1].equals("wavelet")){
            double[] wavelet = OneDHWT.orderedFastHWT(args[2]);
            OneDHWT.saveDataToFile(wavelet, args[3]);
        }

        else if(args[0].equals("wavelet")){
            double[] wavelet = OneDHWT.orderedFastHWT(args[4]);

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

            else if(args[1].equals("minL2")){
                ProbabilisticMinL2.waveletMinL2(wavelet, Integer.parseInt(args[2]));
                if(args[3].equals("full")){
                    OneDHWT.saveDataToFile(wavelet, args[5]);
                }
                else if(args[3].equals("summary")){
                    // Build summary
                }
            }

            else if(args[1].equals("minRelVar")){
                ProbMinRelVar.callMainFunction(wavelet, Integer.parseInt(args[2]), Q,
                        OneDHWT.fileToArrayOfDoubles(args[4]), PERCENTILE);
                if(args[3].equals("full")){
                    OneDHWT.saveDataToFile(wavelet, args[5]);
                }
                else if(args[3].equals("summary")){
                    // Build summary
                }
            }

            else if(args[1].equals("minRelBias")){
                ProbMinRelBias.callMainFunction(wavelet, Integer.parseInt(args[2]), Q,
                        OneDHWT.fileToArrayOfDoubles(args[4]), PERCENTILE);
                if(args[3].equals("full")){
                    OneDHWT.saveDataToFile(wavelet, args[5]);
                }
                else if(args[3].equals("summary")){
                    // Build summary
                }
            }
        }

        else if(args[0].equals("reconstruct")){
            double[] data = OneDHWT.orderedFastHWTInverse(args[3]);

            OneDHWT.saveDataToFile(data, args[4]);
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
}
