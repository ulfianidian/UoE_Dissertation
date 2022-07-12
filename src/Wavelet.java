import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

public class Wavelet {
    public static void main(String[] args) throws IOException {

        final int Q = 20;
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
