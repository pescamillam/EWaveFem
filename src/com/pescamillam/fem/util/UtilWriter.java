package com.pescamillam.fem.util;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;

public class UtilWriter {
    
    private static BufferedWriter bw; 
    
    public static void writeToFile(String stringToWrite) {
        initializeWriterIfNecessary();
        try {
            System.out.println(stringToWrite);
            bw.write(stringToWrite);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

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
