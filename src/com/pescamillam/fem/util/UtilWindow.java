package com.pescamillam.fem.util;

import com.pescamillam.fem.element.Cst;
import org.apache.commons.math3.linear.FieldMatrix;
import org.apache.commons.math3.util.BigReal;

import javax.swing.*;
import java.awt.*;
import java.awt.image.BufferStrategy;
import java.math.BigDecimal;
import java.util.List;

public class UtilWindow {
    
    public static void printMovingNode(List<Cst> elements, FieldMatrix<BigReal>[] displacement,
            FieldMatrix<BigReal>[] speed, FieldMatrix<BigReal>[] acceleration, FieldMatrix<BigReal>[] force, Integer numTimes) {
        final String title = "EWaveFem node displacement";
        final int width = 1200;
        final int height = 500;

        JFrame frame = new JFrame(title);
        frame.setSize(width, height);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.setLocationRelativeTo(null);
        frame.setResizable(true);
        frame.setVisible(true);

        //Creating the canvas.
        Canvas canvas = new Canvas();

        canvas.setSize(width, height);
        canvas.setBackground(Color.BLACK);
        canvas.setVisible(true);
        canvas.setFocusable(false);
        
        frame.add(canvas);

        canvas.createBufferStrategy(3);

        BufferStrategy bufferStrategy;
        Graphics graphics = canvas.getGraphics();
        bufferStrategy = canvas.getBufferStrategy();
        graphics = bufferStrategy.getDrawGraphics();
        while (true) {
            graphics.setColor(Color.GREEN);
            for (int i = 0; i < numTimes-1; i++) {
                for (int j = 8*2+2 + ((8*2)/2+1+6); j < 8*2+2 + ((8*2)/2+1+6) + 1; j++) {
//                for (int j = 0; j < displacement[i].getData().length; j++) {
                    graphics.setColor(Color.LIGHT_GRAY);
                    int y1 = displacement[i].getData()[j][0].bigDecimalValue()
                        .multiply(new BigDecimal("10000"))
                        .intValue();
                    int y2 = displacement[i+1].getData()[j][0].bigDecimalValue()
                            .multiply(new BigDecimal("10000"))
                            .intValue();
                    graphics.drawLine(i*5, y1+250, (i+1)*5, y2+250);
                    if (acceleration[i] != null && acceleration[i+1] != null) {
                        graphics.setColor(Color.GREEN);
                        y1 = acceleration[i].getData()[j][0].bigDecimalValue()
                                .multiply(new BigDecimal("0.001"))
                                .intValue();
                        y2 = acceleration[i+1].getData()[j][0].bigDecimalValue()
                                .multiply(new BigDecimal("0.001"))
                                .intValue();
                        graphics.drawLine(i*5, y1+250, (i+1)*5, y2+250);
                    }
                    if (speed[i] != null && speed[i+1] != null) {
                        graphics.setColor(Color.YELLOW);
                        y1 = speed[i].getData()[j][0].bigDecimalValue()
                                .multiply(new BigDecimal("10"))
                                .intValue();
                        y2 = speed[i+1].getData()[j][0].bigDecimalValue()
                                .multiply(new BigDecimal("10"))
                                .intValue();
                        graphics.drawLine(i*5, y1+250, (i+1)*5, y2+250);
                    }
                    graphics.setColor(Color.CYAN);
                    y1 = force[i].getData()[j][0].bigDecimalValue()
                            .multiply(new BigDecimal("0.002"))
                            .intValue();
                    y2 = force[i+1].getData()[j][0].bigDecimalValue()
                            .multiply(new BigDecimal("0.002"))
                            .intValue();
                    graphics.drawLine(i*5, y1+250, (i+1)*5, y2+250);
                }
                
            }
            
            bufferStrategy.show();
            graphics.dispose();

            try {
                Thread.sleep(200L);
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
        }

    }

    public static void printElements(List<Cst> elements, FieldMatrix<BigReal>[] displacement,
                                      FieldMatrix<BigReal>[] speed, FieldMatrix<BigReal>[] acceleration, FieldMatrix<BigReal>[] force,
                                      Integer numX, Integer numY, Integer numTimes) {
        final String title = "EWaveFem";
        final int width = 30*(numX+2);
        final int height = 30*(numY+3);

        //Creating the frame.
        JFrame frame = new JFrame(title);

        frame.setSize(width, height);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.setLocationRelativeTo(null);
        frame.setResizable(false);
        frame.setVisible(true);

        //Creating the canvas.
        Canvas canvas = new Canvas();

        canvas.setSize(width, height);
        canvas.setBackground(Color.BLACK);
        canvas.setVisible(true);
        canvas.setFocusable(false);


        //Putting it all together.
        frame.add(canvas);

        canvas.createBufferStrategy(3);

        boolean running = true;

        BufferStrategy bufferStrategy;
        Graphics graphics;
        int i = 0;

        while (running) {
            bufferStrategy = canvas.getBufferStrategy();
            graphics = bufferStrategy.getDrawGraphics();
            graphics.clearRect(0, 0, width, height);

            graphics.setColor(Color.GREEN);
            graphics.clearRect(0, 0, 30*(numX+2), 30*(numY+3));

            for (int m = 0; m <= numY; m++) {
                for (int n = 0; n <= numX; n++) {
                    int x = n*30 + displacement[i].getData()[m*(numX+1)*2+n*2][0].bigDecimalValue()
                            //.multiply(new BigDecimal("500"))
                            .intValue();
                    int y = m*30 + displacement[i].getData()[m*(numX+1)*2+n*2+1][0].bigDecimalValue()
                            //.multiply(new BigDecimal("500"))
                            .intValue();
                    graphics.setColor(Color.GREEN);
                    graphics.drawOval(x, y, 10, 10);
//                    if (n < numX) {
//                        graphics.drawLine(x + 5, y + 5, x + 35, y + 5);
//                        if (m < numY) {
//                            graphics.drawLine(x + 5, y + 5, x + 35, y + 35);
//                            graphics.drawLine(x + 5, y + 35, x + 35, y + 35);
//                            graphics.drawLine(x + 5, y + 5, x + 5, y + 35);
//                            graphics.drawLine(x + 35, y + 5, x + 35, y + 35);
//                        }
//                    }
                    if (speed[i] != null) {
                        graphics.setColor(Color.BLUE);
                        graphics.drawLine(x + 5, y + 5,
                                x + 5 + speed[i].getData()[m*(numX+1)*2+n*2][0].bigDecimalValue()
                                        .multiply(new BigDecimal("0.01"))
                                        .intValue(),
                                y + 5 + speed[i].getData()[m*(numX+1)*2+n*2+1][0].bigDecimalValue()
                                        .multiply(new BigDecimal("0.01"))
                                        .intValue());
                    }

                    if (acceleration[i] != null) {
                        graphics.setColor(Color.RED);
                        graphics.drawLine(x + 5, y + 5,
                                x + 5 + acceleration[i].getData()[m*(numX+1)*2+n*2][0].bigDecimalValue()
                                        .multiply(new BigDecimal("0.0001"))
                                        .intValue(),
                                y + 5 + acceleration[i].getData()[m*(numX+1)*2+n*2+1][0].bigDecimalValue()
                                        .multiply(new BigDecimal("0.0001"))
                                        .intValue());
                    }

                    if (force[i] != null) {
                        graphics.setColor(Color.CYAN);
                        graphics.drawLine(x + 5, y + 5,
                                x + 5 + force[i].getData()[m*(numX+1)*2+n*2][0].bigDecimalValue()
                                        .multiply(new BigDecimal("0.001"))
                                        .intValue(),
                                y + 5 + force[i].getData()[m*(numX+1)*2+n*2+1][0].bigDecimalValue()
                                        .multiply(new BigDecimal("0.001"))
                                        .intValue());
                    }
                }
            }

            graphics.drawString("t: " + i, 100, 30*(numY+1));

            bufferStrategy.show();
            graphics.dispose();
            try {
                Thread.sleep(200L);
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
            if (i < numTimes - 1) {
                i++;
            } else {
                i = 0;
            }
        }
    }
}
