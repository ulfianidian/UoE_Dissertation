import java.io.IOException;

public class Wavelet {
    public static void main(String[] args) throws IOException {

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
        else{
            System.out.println("Enter a valid command!");
        }

        /**
        ProbMinRelBias.callMainFunction("/Users/ulfianidian/Desktop/Dissertation/1DHaar/data/wavelet/wavelet1",
                8.0, 10, "/Users/ulfianidian/Desktop/Dissertation/1DHaar/data/original_data/smalldata",
                10);**/


        /**OneDWaveletDecomp("/Users/ulfianidian/Desktop/Dissertation/1DHaar/data/original_data/smalldata",
                "/Users/ulfianidian/Desktop/Dissertation/1DHaar/data/wavelet/wavelet1");**/


        /**ReconstructData("/Users/ulfianidian/Desktop/Dissertation/1DHaar/data/OneMillionWavelet",
                "/Users/ulfianidian/Desktop/Dissertation/1DHaar/data/ZipfDistribution2964476632319761509.txt");**/

        /**
         // Get lambda
         System.out.println();
         double[] wavelet = OneDHWT.fileToArrayOfDoubles("/Users/ulfianidian/Desktop/Dissertation/1DHaar/data/wavelet/wavelet1");
         ProbabilisticMinL2.waveletMinL2(wavelet, 8);**/

        /**getConventionalWavelet("/Users/ulfianidian/Desktop/Dissertation/1DHaar/data/original_data/smalldata",
                8, "/Users/ulfianidian/Desktop/Dissertation/1DHaar/data/wavelet/wavelet1");**/
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
        double[] wavelet = OneDHWT.orderedFastHWT(filePath);

        //Print the wavelet on the console
        for (double w : wavelet) {
            System.out.println(w);
        }

        OneDHWT.saveDataToFile(wavelet, writeTo);
    }

    public static void ReconstructData(String filePath, String writeTo) throws IOException {
        double[] data = OneDHWT.orderedFastHWTInverse(filePath);
        for(int i = 0; i < data.length; i++){
            System.out.println(data[i]);
        }
        OneDHWT.saveDataToFile(data, writeTo);
    }

    public static void getConventionalWavelet(String filePath, int b , String writeTo)
        throws IOException {
        double[] wavelet = OneDHWT.orderedFastHWT(filePath);

        //Perform normalization
        double[] normalizedCoeffs = Conventional.normalizeCoeffs(wavelet);

        //Retain B largest coefficients
        Conventional.retainNCoeffs(wavelet, normalizedCoeffs, b);
        //OneDHWT.saveDataToFile(normalizedCoeffs, writeTo);
    }

}
