import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

/**
 * @author Ulfiani Primawati
 *
 * This is an implementation of 1-Dimensional Ordered Fast Haar Wavelet Transform and
 * Ordered Fast Inverse Haar Wavelet Transform.
 */

public class OneDHWT {

    /**
     * This method is used to display the data
     * @param data Array of Data
     */
    public static void displayData(double[] data){
        System.out.println("Data:");
        for (double datum : data) {
            System.out.println(datum);
        }
        System.out.println();
    }

    /**
     * Check whether n is power of 2
     * @param n The number of data in the array
     * @return True if n is power of 2, false otherwise
     */
    public static boolean isPowerOfTwo(int n){
        if(n < 1){
            return false;
        }
        else{
            double logNBase2 = (Math.log(n) / Math.log(2));   //Calculate log_2(n)
            return Math.abs(logNBase2 - (int)logNBase2) == 0;
        }
    }

    /**
     * Finds the largest power of two integer that is less than or equal to n.
     * @param n The number of data in the array
     * @return The largest power of two integer that is <= n.
     */
    public static int largestPowerOf2LessThanN(int n){
        if(isPowerOfTwo(n)){
            return n;
        }
        else{
            int curr = n-1;
            while(curr > 0){
                if(isPowerOfTwo(curr))
                    return curr;
                else
                    --curr;
            }
            return 0;
        }
    }

    /**
     * Takes i number subset of the data, where i is the largest power of two
     * integer that is less than or equal to the number of element in the array.
     * @param data The array of data
     * @return The subset of the array
     */
    public static double[] largestPowerOf2Subdata(double[] data){
        if(isPowerOfTwo(data.length)){
            return data;
        }
        else{
            int i = largestPowerOf2LessThanN(data.length);
            if(i == 0)
                return null;

            double[] subData = new double[i];
            System.arraycopy(data, 0, subData, 0, i);
            return subData;
        }
    }

    public static double[] padWithZeroes(double[] data){
        int i = largestPowerOf2LessThanN(data.length);

        int logNBase2 = (int)(Math.log(i) / Math.log(2));
        int i_2 = logNBase2 + 1;
        return Arrays.copyOf(data, (int)Math.pow(2, i_2));
    }

    /**
     * Write the data to a file
     * @param data The data
     * @param filePath File Path
     * @throws FileNotFoundException File not found exception
     * @throws IOException IO exception
     */
    public static void saveDataToFile(double[] data, String filePath)
            throws FileNotFoundException, IOException {
        FileWriter file = new FileWriter(filePath);
        BufferedWriter buffer = new BufferedWriter(file);

        for(double d: data){
            buffer.write(Double.toString(d));
            buffer.newLine();
        }

        buffer.flush();

        //Close the input stream
        buffer.close();
        file.close();
    }

    /**
     * Write the data to a file
     * @param data The data
     * @param filePath File Path
     * @throws FileNotFoundException File not found exception
     * @throws IOException IO exception
     */
    public static void saveDataToFile(String[] data, String filePath)
            throws FileNotFoundException, IOException {
        FileWriter file = new FileWriter(filePath);
        BufferedWriter buffer = new BufferedWriter(file);

        for(String d: data){
            buffer.write(d);
            buffer.newLine();
        }

        buffer.flush();

        //Close the input stream
        buffer.close();
        file.close();
    }

    /**
     * Read data from a file and return it as an array of doubles
     * @param filePath Path of the file
     * @return The data as an array of doubles
     * @throws FileNotFoundException File not found
     * @throws IOException IO exception
     */
    public static double[] fileToArrayOfDoubles(String filePath)
            throws FileNotFoundException, IOException{
        FileInputStream fileStream = new FileInputStream(filePath);
        BufferedReader bufferedReader = new BufferedReader(new InputStreamReader(fileStream));
        ArrayList<Double> arrayOfDoubles = new ArrayList<>();

        //Read the file, line by line
        String stringLine;
        while((stringLine = bufferedReader.readLine()) != null){
            //Print the value on the console
            //System.out.println(stringLine);
            arrayOfDoubles.add(Double.valueOf(stringLine));
        }

        //Close the input stream
        bufferedReader.close();

        double[] doublePrims = new double[arrayOfDoubles.size()];
        int i = 0;
        for(Double d: arrayOfDoubles){
            doublePrims[i++] = d.doubleValue();
        }

        arrayOfDoubles = null;

        return doublePrims;
    }

    /**
     * Performs the Ordered Fast Haar Wavelet Transform which takes input from a file
     * Each sample in the file is a double written on a separate line.
     * @param filePath Path of the file
     * @return HWT of the data
     * @throws IOException IO Exception
     * @throws FileNotFoundException File not found exception
     */
    public static double[] orderedFastHWT(String filePath) throws IOException, FileNotFoundException{
        double[] data = fileToArrayOfDoubles(filePath);
        final int n = data.length;

        if(!OneDHWT.isPowerOfTwo(n)){
            double[] newData;
            newData = padWithZeroes(data);
            orderedFastHWT(newData);
            return newData;
        }
        else{
            orderedFastHWT(data);
            return data;
        }
    }

    public static double[] checkLength(double[] data){
        final int n = data.length;
        double[] newData;

        if(!OneDHWT.isPowerOfTwo(n)){
            newData = padWithZeroes(data);
            orderedFastHWT(newData);
            return newData;
        }
        else{
            orderedFastHWT(data);
            return data;
        }
    }

    /**
     * Performs the Ordered Fast Haar Wavelet Transform where the input is data array
     * @param data Data array
     */
    public static void orderedFastHWT(double[] data){
        final int n = data.length;

        //If n is not an integer power of 2, then return
        if(!OneDHWT.isPowerOfTwo(n)){
            return;
        }

        //Compute the number of sweeps (e.g. if n = 8, then NUM_OF_SWEEPS is 3)
        final int NUM_OF_SWEEPS = (int)(Math.log(n) / Math.log(2.0));

        double aCoeff, cCoeff;

        //If the number of sweeps is 1
        if(NUM_OF_SWEEPS == 1){
            aCoeff = (data[0] + data[1]) / 2.0;
            cCoeff = (data[0] - data[1]) / 2.0;
            data[0] = aCoeff;
            data[1] = cCoeff;
            return;
        }

        //If the number of sweeps is greater than 1
        double[] aCoefficients;
        double[] cCoefficients;
        for(int nthSweep = 1; nthSweep < NUM_OF_SWEEPS; nthSweep++){

            //size is the number of a-coefficients and c coefficients
            //at nthSweep. E.g. if the data has 8 elements:
            //Sweep 1: 4 a-coefficients and 4 c-coefficients;
            //Sweep 2: 2 a-coefficients and 2 c-coefficients;
            //Sweep 3: 1 a-coefficient and 1 c-coefficient.
            double remainingSweeps = (double)(NUM_OF_SWEEPS - nthSweep);
            int size = (int) Math.pow(2.0, remainingSweeps);

            //endOfIndex is the index of last element in data[] at nthSweep
            double sweeps = remainingSweeps + 1.0;
            int endOfIndex = ((int)Math.pow(2.0, sweeps)) - 1;

            aCoefficients = new double[size];
            cCoefficients = new double[size];
            int index_a = 0; // index for aCoefficients
            int index_c = 0; // index for cCoefficients

            //Build an array of a-coefficients and an array of c-coefficients
            for(int i = 0; i <= endOfIndex; i += 2){
                cCoefficients[index_a++] = (data[i] - data[i + 1]) / 2.0;
                aCoefficients[index_c++] = (data[i] + data[i + 1]) / 2.0;
            }

            //Replace the first half of the data array with a-coefficients
            //and the second half with c-coefficients.
            for(int i = 0; i < size; i++){
                data[i] = aCoefficients[i];
                data[i + size] = cCoefficients[i];
            }
        }

        //Do the last sweep
        cCoeff = (data[0] - data[1]) / 2.0;
        aCoeff = (data[0] + data[1]) / 2.0;
        data[0] = aCoeff;
        data[1] = cCoeff;
    }

    /**
     * Performs the Ordered Fast Haar Wavelet Transform which takes input from a file
     * Each sample in the file is a double written on a separate line.
     * @param filePath Path of the file
     * @return HWT of the data
     * @throws IOException IO Exception
     * @throws FileNotFoundException File not found exception
     */
    public static double[] orderedFastHWT2(String filePath) throws IOException, FileNotFoundException{
        double[] data = fileToArrayOfDoubles(filePath);
        orderedFastHWT(data);
        return data;
    }

    /**
     * Performs the Ordered Fast Haar Wavelet Transform where the input is data array
     * @param data Data array
     */
    public static void orderedFastHWT2(double[] data){
        final int n = data.length;

        //If n is not an integer power of 2, then return
        if(!OneDHWT.isPowerOfTwo(n)){
            return;
        }

        //Compute the number of sweeps (e.g. if n = 8, then NUM_OF_SWEEPS is 3)
        final int NUM_OF_SWEEPS = (int)(Math.log(n) / Math.log(2.0));

        double aCoeff, cCoeff;

        //If the number of sweeps is 1
        if(NUM_OF_SWEEPS == 1){
            aCoeff = (data[0] + data[1]) / 2.0;
            cCoeff = (data[0] - data[1]) / 2.0;
            data[0] = aCoeff;
            data[1] = cCoeff;
            return;
        }

        //If the number of sweeps is greater than 1
        double[] aCoefficients;
        double[] cCoefficients;
        for(int nthSweep = 1; nthSweep < NUM_OF_SWEEPS; nthSweep++){

            //size is the number of a-coefficients and c coefficients
            //at nthSweep. E.g. if the data has 8 elements:
            //Sweep 1: 4 a-coefficients and 4 c-coefficients;
            //Sweep 2: 2 a-coefficients and 2 c-coefficients;
            //Sweep 3: 1 a-coefficient and 1 c-coefficient.
            double remainingSweeps = (double)(NUM_OF_SWEEPS - nthSweep);
            int size = (int) Math.pow(2.0, remainingSweeps);

            //endOfIndex is the index of last element in data[] at nthSweep
            double sweeps = remainingSweeps + 1.0;
            int endOfIndex = ((int)Math.pow(2.0, sweeps)) - 1;

            aCoefficients = new double[size];
            cCoefficients = new double[size];
            int index_a = 0; // index for aCoefficients
            int index_c = 0; // index for cCoefficients

            //Build an array of a-coefficients and an array of c-coefficients
            for(int i = 0; i <= endOfIndex; i += 2){
                cCoefficients[index_a++] = (data[i] - data[i + 1]) / 2.0;
                aCoefficients[index_c++] = (data[i] + data[i + 1]) / 2.0;
            }

            //Replace the first half of the data array with a-coefficients
            //and the second half with c-coefficients.
            for(int i = 0; i < size; i++){
                data[i] = aCoefficients[i];
                data[i + size] = cCoefficients[i];
            }
        }

        //Do the last sweep
        cCoeff = (data[0] - data[1]) / 2.0;
        aCoeff = (data[0] + data[1]) / 2.0;
        data[0] = aCoeff;
        data[1] = cCoeff;
    }

    /**
     * Reconstruct data based on a Haar Wavelet coefficients.
     * @param wavelet Haar Wavelet Coefficients
     */
    public static void orderedFastHWTInverse(double[] wavelet){
        int n = wavelet.length;

        //Check the validity of the wavelet1
        if(!OneDHWT.isPowerOfTwo(n) || n < 2){
            return;
        }

        final int NUM_OF_SWEEPS = (int)(Math.log(n) / Math.log(2.0));
        int gap;
        double data_0 = 0;
        double data_1 = 0;
        double[] values;

        for(int j = 1; j <= NUM_OF_SWEEPS; j++){
            values = null;
            gap = (int) (Math.pow(2.0, j - 1)); // Gap between a and its corresponding c at level j
            values = new double[2 * gap];

            for(int i = 0; i < gap; i++){
                data_0 = wavelet[i] + wavelet[i + gap];
                data_1 = wavelet[i] - wavelet[i + gap];
                values[2 * i] = data_0;
                values[2 * i + 1] = data_1;
            }

            System.arraycopy(values, 0, wavelet, 0, gap * 2);
        }
    }

    /**
     * Reconstruct data based on a Haar Wavelet coefficients.
     * This method takes input from a file where each coefficient is
     * a double written on a separate line.
     * @param filePath File path
     * @return Reconstructed data
     * @throws IOException IO Exception
     * @throws FileNotFoundException File not found exception
     */
    public static double[] orderedFastHWTInverse(String filePath) throws IOException, FileNotFoundException{
        double[] wavelet = fileToArrayOfDoubles(filePath);
        orderedFastHWTInverse(wavelet);
        return wavelet;
    }

    public static HashMap<Integer, Double> fullToSummary(double[] wavelet){
        HashMap<Integer, Double> mapping = new HashMap<>();
        for(int i = 0; i < wavelet.length; i++){
            if(wavelet[i] != 0.0){
                mapping.put(i, wavelet[i]);
            }
        }
        return mapping;
    }
}
