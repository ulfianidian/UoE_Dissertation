import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;

import org.apache.commons.math3.distribution.ZipfDistribution;
import org.apache.commons.io.FileUtils;


/**
 * @author Ulfiani Primawati
 *
 * This is a class to generate data based on a Zipfian Distribution
 */
public class ZipfianGenerator {
    protected static String MODE;
    protected static int KEY_SPACE;
    protected static int TOTAL_RECORDS;
    protected static double Z_PARAM;
    protected static int START;
    protected static String PERMUTATION;
    protected static String WRITE_TO;
    protected static int GAP;


    /**
     * Constructor method for ZipfDataGenerator class
     * @param keySpace The number of distinct values
     * @param totalRecords Desired tuple counts
     */
    public ZipfianGenerator(String mode, int keySpace, int totalRecords, double zParam,
                            int start, String permutation, String writeTo){
        MODE = mode;
        KEY_SPACE = keySpace;
        TOTAL_RECORDS = totalRecords;
        Z_PARAM = zParam;
        START = start;
        PERMUTATION = permutation;
        WRITE_TO = writeTo;
    }

    public ZipfianGenerator(String mode, int keySpace, int totalRecords, double zParam,
                            int start, String permutation, String writeTo, int gap){
        MODE = mode;
        KEY_SPACE = keySpace;
        TOTAL_RECORDS = totalRecords;
        Z_PARAM = zParam;
        START = start;
        PERMUTATION = permutation;
        WRITE_TO = writeTo;
        GAP = gap;
    }

    public void generateZipf() throws IOException{
        System.out.println("Creating test data...");
        long startingTime = System.currentTimeMillis();
        ZipfDistribution zipfDistribution = new ZipfDistribution(KEY_SPACE, Z_PARAM);
        int ind;

        File tempFile = new File(WRITE_TO);

        File file = File.createTempFile("ZipfDistribution", ".txt", tempFile);

        if(file.exists())
            FileUtils.forceDelete(file);

        FileWriter fileWriter = new FileWriter(file);

        if(MODE.equals("sequential") && START == 1 && PERMUTATION.equals("normal")) {
            try {
                for (int i = 0; i < TOTAL_RECORDS; i++) {
                    fileWriter.write(String.valueOf((double)zipfDistribution.sample()));
                    fileWriter.write('\n');
                }
            }
            finally {
                if (fileWriter != null)
                    fileWriter.close();
            }
        }

        else{
            double[] keysDouble = new double[KEY_SPACE];

            for(int j = 1; j <= KEY_SPACE; j++){
                keysDouble[j - 1] = j;
            }

            if(MODE.equals("nonconsecutive")){
                for(int k = 1; k < KEY_SPACE; k++){
                    keysDouble[k] = keysDouble[k - 1] + (double)GAP + 1.0;
                }
            }

            if(START > 1){
                for(int l = 0; l < KEY_SPACE; l++){
                    keysDouble[l] = keysDouble[l] + START - 1;
                }
            }

            if(PERMUTATION.equals("bell-shaped")){
                double[] copyKeySpace = new double[KEY_SPACE];
                System.arraycopy(keysDouble, 0, copyKeySpace, 0, KEY_SPACE);

                int i = KEY_SPACE - 1;
                for(int j = 0; j < KEY_SPACE - (KEY_SPACE / 2); j++){
                    keysDouble[i] = copyKeySpace[j];
                    i -= 2;
                }
                i = KEY_SPACE - 2;
                for(int j = KEY_SPACE - 1; i >= 0 && j >= KEY_SPACE / 2; j--){
                    keysDouble[i] = copyKeySpace[j];
                    i -= 2;
                }
            }
            else if(PERMUTATION.equals("u-shaped")){
                double[] copyKeySpace = new double[KEY_SPACE];
                System.arraycopy(keysDouble, 0, copyKeySpace, 0, KEY_SPACE);

                int i  = 0;
                for(int j = 0; i < KEY_SPACE && j <= KEY_SPACE / 2; j++){
                    keysDouble[i] = copyKeySpace[j];
                    i += 2;
                }
                i = 1;
                for(int j = KEY_SPACE - 1; i < KEY_SPACE && j >= KEY_SPACE / 2; j--){
                    keysDouble[i] = copyKeySpace[j];
                    i += 2;
                }
            }
            else if(PERMUTATION.equals("random")){
                ArrayList<Double> list = new ArrayList<>();
                for(int i = 0; i < KEY_SPACE; i++){
                    list.add(keysDouble[i]);
                }
                Collections.shuffle(list);

                for(int i = 0; i < KEY_SPACE; i++)
                    keysDouble[i] = list.get(i);
            }

            try{
                for(int i = 0; i < TOTAL_RECORDS; i++){
                    ind = zipfDistribution.sample() - 1;
                    fileWriter.write(String.valueOf(keysDouble[ind]));
                    fileWriter.write('\n');
                }
            }
            finally {
                if (fileWriter != null)
                    fileWriter.close();
            }
        }

        System.out.println("Time for creating the data: " +
                (System.currentTimeMillis() - startingTime) + "milliseconds.");
        System.out.println("Path: " + file.getAbsolutePath());
    }
}
