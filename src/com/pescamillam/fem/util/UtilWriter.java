package com.pescamillam.fem.util;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;

/**
 * Class to write to output file
 * 
 * @author Peter Escamilla (pescamilla@unab.edu.co)
 */
public class UtilWriter {
    
    private static BufferedWriter bw; 

    
    /**
     * Prints to file the received string
     * 
     * @param stringToWrite string to write to file
     */
    public static void writeToFile(String stringToWrite) {
        initializeWriterIfNecessary();
        try {
            bw.write(stringToWrite);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     * In case the writer doesn't exists this creates it
     */
    private static void initializeWriterIfNecessary() {
        if (bw == null) {
            try {
                bw = new BufferedWriter(new FileWriter("output.txt"));
            } catch (IOException e) {
                e.printStackTrace();
            } 
        }
    }
}
