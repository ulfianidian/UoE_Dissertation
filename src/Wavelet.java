import java.io.IOException;

public class Wavelet {
    public static void main(String[] args) throws IOException {

        if(args.length != 6 || args.length != 4){
           // Add message and then return
        }

        //generateZipfData(1000, 1048576);


        /**OneDWaveletDecomp("/Users/ulfianidian/Desktop/Dissertation/1DHaar/data/original_data/smalldata",
                "/Users/ulfianidian/Desktop/Dissertation/1DHaar/data/wavelet/wavelet1");**/


        /**ReconstructData("/Users/ulfianidian/Desktop/Dissertation/1DHaar/data/OneMillionWavelet",
                "/Users/ulfianidian/Desktop/Dissertation/1DHaar/data/ZipfDistribution2964476632319761509.txt");**/

        getConventionalWavelet("/Users/ulfianidian/Desktop/Dissertation/1DHaar/data/original_data/smalldata",
                16, "/Users/ulfianidian/Desktop/Dissertation/1DHaar/data/wavelet/wavelet1");
    }

    public static void generateZipfData(int keySpace, int tupleCounts) throws IOException {
        ZipfDataGenerator zipfDataGenerator = new ZipfDataGenerator(keySpace, tupleCounts);
        zipfDataGenerator.generateZipf();
    }

    public static void OneDWaveletDecomp(String filePath, String writeTo) throws IOException {
        double[] wavelet = OneDHWT.orderedFastHWT(filePath);

        //Print the wavelet on the console
        for(int i = 0; i < wavelet.length; i++){
            System.out.println(wavelet[i]);
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
