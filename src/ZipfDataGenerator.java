import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;

import org.apache.commons.math3.distribution.ZipfDistribution;
import org.apache.commons.io.FileUtils;

public class ZipfDataGenerator {
    protected static int KEY_SPACE;
    protected static int TOTAL_RECORDS;

    public ZipfDataGenerator(int keySpace, int totalRecords){
        KEY_SPACE = keySpace;
        TOTAL_RECORDS = totalRecords;
    }

    public String generateZipf() throws IOException{

        //double[] keysDouble;
        //keysDouble = generateKeySpace(200, 150);

        System.out.println("Creating test data...");
        long startingTime = System.currentTimeMillis();
        ZipfDistribution zipfDistribution = new ZipfDistribution(KEY_SPACE, 1.0);
        int ind;

        File tempFile = new File("/Users/ulfianidian/Desktop/Dissertation/1DHaar/data");

        File file = File.createTempFile("ZipfDistribution", ".txt", tempFile);

        if(file.exists())
            FileUtils.forceDelete(file);

        FileWriter fileWriter = new FileWriter(file);

        /**
        try{
            for(int i = 0; i < TOTAL_RECORDS; i++){
                ind = zipfDistribution.sample() - 1;
                fileWriter.write(String.valueOf(keysDouble[ind]));
                fileWriter.write('\n');
            }
        }**/

        try{
            for(int i = 0; i < TOTAL_RECORDS; i++) {
                fileWriter.write(String.valueOf((double) zipfDistribution.sample()));
                fileWriter.write('\n');
            }
        }
         finally {
            if (fileWriter != null)
                fileWriter.close();
        }

        System.out.println("Time for creating the data: " +
                (System.currentTimeMillis() - startingTime) + "milliseconds.");
        System.out.println("Path: " + file.getAbsolutePath());

        return file.getAbsolutePath();
    }

    private double[] generateKeySpace(int source, int size){
        ArrayList<Double> list = new ArrayList<Double>();
        for(int i = 1; i <= source; i++)
            list.add((double)i);
        Collections.shuffle(list);

        double[] outputArray = new double[size];
        for(int j = 0; j < size; j++){
            outputArray[j] = list.get(j);
        }
        return outputArray;
    }
}
