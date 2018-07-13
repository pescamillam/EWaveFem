package com.pescamillam.fem.util;

import com.pescamillam.fem.element.Cst;
import org.apache.commons.math3.linear.FieldMatrix;
import org.apache.commons.math3.util.BigReal;

import javax.swing.*;
import java.awt.*;
import java.awt.image.BufferStrategy;
import java.math.BigDecimal;
import java.text.DecimalFormat;
import java.util.List;

/**
 * Util class to show result windows
 * 
 * @author Peter Escamilla (pescamilla@unab.edu.co)
 */
public class UtilWindow {

    static final DecimalFormat DF = new DecimalFormat("0.###E0");
    
    /**
     * Shows an image with the displacement, speed, acceleration and force of a node
     * 
     * @param elements List with all the elements of the grid
     * @param displacement vector of displacement
     * @param speed vector of speed
     * @param acceleration vector of acceleration
     * @param force vector of force
     * @param numTimes number of iterations
     */
    public static void printMovingNode(List<Cst> elements, FieldMatrix<BigReal>[] displacement,
            FieldMatrix<BigReal>[] speed, FieldMatrix<BigReal>[] acceleration, FieldMatrix<BigReal>[] force, 
            Integer numTimes, int numNode) {
        //assigns the title of the window
        final String title = "EWaveFem node displacement";
        //defines the size of the window
        final int width = 1200;
        final int height = 500;

        //creates the window with the given size
        JFrame frame = new JFrame(title);
        frame.setSize(width, height);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.setLocationRelativeTo(null);
        frame.setResizable(true);
        frame.setVisible(true);

        //Creating the canvas where the vectors will be drawn
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
        //keeps drawing the vectors
        while (true) {

            BigDecimal minDisplacement = new BigDecimal(Double.MAX_VALUE);
            BigDecimal maxDisplacement = new BigDecimal(Double.MIN_VALUE);

            BigDecimal minSpeed = new BigDecimal(Double.MAX_VALUE);
            BigDecimal maxSpeed = new BigDecimal(Double.MIN_VALUE);

            BigDecimal minAcceleration = new BigDecimal(Double.MAX_VALUE);
            BigDecimal maxAcceleration = new BigDecimal(Double.MIN_VALUE);

            BigDecimal minForce = new BigDecimal(Double.MAX_VALUE);
            BigDecimal maxForce = new BigDecimal(Double.MIN_VALUE);

            //iterates for every possible time
            for (int i = 0; i < numTimes-1; i++) {
                for (int j = numNode; j < numNode + 1; j++) {
                    //Draws the displacement
                    graphics.setColor(Color.LIGHT_GRAY);
                    ((Graphics2D)graphics).setStroke(new BasicStroke(3));
                    int y1 = displacement[i].getData()[j][0].bigDecimalValue()
                        .multiply(new BigDecimal("100"))
                        .intValue();
                    int y2 = displacement[i+1].getData()[j][0].bigDecimalValue()
                            .multiply(new BigDecimal("100"))
                            .intValue();
                    graphics.drawLine(i*5, y1+250, (i+1)*5, y2+250);

                    if (displacement[i].getData()[j][0].bigDecimalValue().compareTo(maxDisplacement) > 0) {
                        maxDisplacement = displacement[i].getData()[j][0].bigDecimalValue();
                    }
                    if (displacement[i].getData()[j][0].bigDecimalValue().compareTo(minDisplacement) < 0) {
                        minDisplacement = displacement[i].getData()[j][0].bigDecimalValue();
                    }

                    //Draws the acceleration
                    if (acceleration[i] != null && acceleration[i+1] != null) {
                        graphics.setColor(Color.GREEN);
                        y1 = acceleration[i].getData()[j][0].bigDecimalValue()
                                .multiply(new BigDecimal("0.00001"))
                                .intValue();
                        y2 = acceleration[i+1].getData()[j][0].bigDecimalValue()
                                .multiply(new BigDecimal("0.00001"))
                                .intValue();
                        graphics.drawLine(i*5, y1+250, (i+1)*5, y2+250);

                        if (acceleration[i].getData()[j][0].bigDecimalValue().compareTo(maxAcceleration) > 0) {
                            maxAcceleration = acceleration[i].getData()[j][0].bigDecimalValue();
                        }
                        if (acceleration[i].getData()[j][0].bigDecimalValue().compareTo(minAcceleration) < 0) {
                            minAcceleration = acceleration[i].getData()[j][0].bigDecimalValue();
                        }
                    }
                    
                    //Draws the speed
                    if (speed[i] != null && speed[i+1] != null) {
                        graphics.setColor(Color.YELLOW);
                        y1 = speed[i].getData()[j][0].bigDecimalValue()
                                .multiply(new BigDecimal("0.06"))
                                .intValue();
                        y2 = speed[i+1].getData()[j][0].bigDecimalValue()
                                .multiply(new BigDecimal("0.06"))
                                .intValue();
                        graphics.drawLine(i*5, y1+250, (i+1)*5, y2+250);

                        if (speed[i].getData()[j][0].bigDecimalValue().compareTo(maxSpeed) > 0) {
                            maxSpeed = speed[i].getData()[j][0].bigDecimalValue();
                        }
                        if (speed[i].getData()[j][0].bigDecimalValue().compareTo(minSpeed) < 0) {
                            minSpeed = speed[i].getData()[j][0].bigDecimalValue();
                        }
                    }
                    
                    //Draws the force
                    graphics.setColor(Color.CYAN);
                    y1 = force[i].getData()[j][0].bigDecimalValue()
                            .multiply(new BigDecimal("0.002"))
                            .intValue();
                    y2 = force[i+1].getData()[j][0].bigDecimalValue()
                            .multiply(new BigDecimal("0.002"))
                            .intValue();
                    graphics.drawLine(i*5, y1+250, (i+1)*5, y2+250);

                    if (force[i].getData()[j][0].bigDecimalValue().compareTo(maxForce) > 0) {
                        maxForce = force[i].getData()[j][0].bigDecimalValue();
                    }
                    if (force[i].getData()[j][0].bigDecimalValue().compareTo(minForce) < 0) {
                        minForce = force[i].getData()[j][0].bigDecimalValue();
                    }
                    
                    //Draws a 0 line for reference
                    graphics.setColor(Color.WHITE);
                    ((Graphics2D)graphics).setStroke(new BasicStroke(1));
                    graphics.drawLine(i*5, 250, (i+1)*5, 250);
                }
                
            }

            String label = "Node " + numNode/2 + (numNode%2 == 0 ? " x " : " y ");
            graphics.setColor(Color.RED);
            graphics.drawChars(label.toCharArray(), 0, label.length(), 10, 10);

            label = "Displacement";
            graphics.setColor(Color.LIGHT_GRAY);
            graphics.drawChars(label.toCharArray(), 0, label.length(), 10, 20);
            label = "max: " + DF.format(maxDisplacement) + " inches min " + DF.format(minDisplacement) + " inches";
            graphics.drawChars(label.toCharArray(), 0, label.length(), 100, 20);

            label = "Speed";
            graphics.setColor(Color.YELLOW);
            graphics.drawChars(label.toCharArray(), 0, label.length(), 10, 30);
            label = "max: " + DF.format(maxSpeed) + " in/s min " + DF.format(minSpeed) + " in/s";
            graphics.drawChars(label.toCharArray(), 0, label.length(), 100, 30);

            label = "Acceleration";
            graphics.setColor(Color.GREEN);
            graphics.drawChars(label.toCharArray(), 0, label.length(), 10, 40);
            label = "max: " + DF.format(maxAcceleration) + " in/s2 min " + DF.format(minAcceleration) + " in/s2";
            graphics.drawChars(label.toCharArray(), 0, label.length(), 100, 40);

            label = "Force";
            graphics.setColor(Color.CYAN);
            graphics.drawChars(label.toCharArray(), 0, label.length(), 10, 50);
            label = "max: " + DF.format(maxForce) + " pound/in2 min " + DF.format(minForce) + " pound/in2";
            graphics.drawChars(label.toCharArray(), 0, label.length(), 100, 50);

            label = "0 Line";
            graphics.setColor(Color.WHITE);
            graphics.drawChars(label.toCharArray(), 0, label.length(), 10, 60);
            
            bufferStrategy.show();
            graphics.dispose();

            try {
                //keeps drawing every 0.2 seconds
                Thread.sleep(200L);
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
        }

    }

    /**
     * Shows a window with an animation of all the nodes with a line representing speed, acceleration and force for
     * each node
     * 
     * @param elements List with all the elements of the grid
     * @param displacement vector of displacement
     * @param speed vector of speed
     * @param acceleration vector of acceleration
     * @param force vector of force
     * @param numX number of elements in X
     * @param numY number of elements in Y
     * @param numTimes number of iterations
     */
    public static void printElements(List<Cst> elements, FieldMatrix<BigReal>[] displacement,
                                      FieldMatrix<BigReal>[] speed, FieldMatrix<BigReal>[] acceleration, FieldMatrix<BigReal>[] force,
                                      Integer numX, Integer numY, Integer numTimes) {
        //assigns the title of the window
        final String title = "EWaveFem";
        final int width = 60*(numX+2);
        final int height = 60*(numY+4);

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
        Graphics2D graphics;
        int i = 0;

        //keeps executing to show the animation
        while (running) {
            bufferStrategy = canvas.getBufferStrategy();
            graphics = (Graphics2D)bufferStrategy.getDrawGraphics();
            graphics.setStroke(new BasicStroke(2));
            graphics.clearRect(0, 0, width, height);

            graphics.setColor(Color.GREEN);
            graphics.clearRect(0, 0, 30*(numX+2), 30*(numY+3));

            for (int m = 0; m <= numY; m++) {
                for (int n = 0; n <= numX; n++) {

                    //places the node based on the displacement
                    int x = n*60 + displacement[i].getData()[m*(numX+1)*2+n*2][0].bigDecimalValue()
                            .multiply(new BigDecimal("30"))
                            .intValue();
                    int y = m*60 + displacement[i].getData()[m*(numX+1)*2+n*2+1][0].bigDecimalValue()
                            .multiply(new BigDecimal("30"))
                            .intValue();

                    //draws the node
                    graphics.setColor(Color.GREEN);
                    graphics.drawOval(x, y, 20, 20);

                    graphics.setStroke(new BasicStroke(4));
                    //draws the speed vector
                    if (speed[i] != null) {
                        graphics.setColor(Color.BLUE);
                        graphics.drawLine(x + 10, y + 10,
                                x + 10 + speed[i].getData()[m*(numX+1)*2+n*2][0].bigDecimalValue()
                                        .multiply(new BigDecimal("0.02"))
                                        .intValue(),
                                y + 10 + speed[i].getData()[m*(numX+1)*2+n*2+1][0].bigDecimalValue()
                                        .multiply(new BigDecimal("0.02"))
                                        .intValue());
                    }

                    //draws the acceleration vector
                    if (acceleration[i] != null) {
                        graphics.setColor(Color.RED);
                        graphics.drawLine(x + 10, y + 10,
                                x + 10 + acceleration[i].getData()[m*(numX+1)*2+n*2][0].bigDecimalValue()
                                        .multiply(new BigDecimal("0.000003"))
                                        .intValue(),
                                y + 10 + acceleration[i].getData()[m*(numX+1)*2+n*2+1][0].bigDecimalValue()
                                        .multiply(new BigDecimal("0.000003"))
                                        .intValue());
                    }

                    //draws the force vector
                    if (force[i] != null) {
                        graphics.setColor(Color.CYAN);
                        graphics.drawLine(x + 10, y + 10,
                                x + 10 + force[i].getData()[m*(numX+1)*2+n*2][0].bigDecimalValue()
                                        .multiply(new BigDecimal("0.001"))
                                        .intValue(),
                                y + 10 + force[i].getData()[m*(numX+1)*2+n*2+1][0].bigDecimalValue()
                                        .multiply(new BigDecimal("0.001"))
                                        .intValue());
                    }
                    graphics.setStroke(new BasicStroke(2));
                }
            }

            //draws the current time iteration
            graphics.drawString("t: " + i, 10, 60*(numY+1));
            graphics.setColor(Color.BLUE);
            graphics.drawString("Speed", 10, 60*(numY+1)+15);
            graphics.setStroke(new BasicStroke(4));
            graphics.drawLine(80, 60*(numY+1)+10, 90, 60*(numY+1)+10);
            graphics.setStroke(new BasicStroke(2));
            graphics.drawString(DF.format(new BigDecimal("1000")) + " in/s", 95, 60*(numY+1)+15);
            graphics.setColor(Color.RED);
            graphics.drawString("Acceleration", 10, 60*(numY+1)+30);
            graphics.setStroke(new BasicStroke(4));
            graphics.drawLine(80, 60*(numY+1)+25, 90, 60*(numY+1)+25);
            graphics.setStroke(new BasicStroke(2));
            graphics.drawString(DF.format(new BigDecimal("10000000")) + " in/s2", 95, 60*(numY+1)+30);
            graphics.setColor(Color.CYAN);
            graphics.drawString("Force", 10, 60*(numY+1)+45);
            graphics.setStroke(new BasicStroke(4));
            graphics.drawLine(80, 60*(numY+1)+40, 90, 60*(numY+1)+40);
            graphics.setStroke(new BasicStroke(2));
            graphics.drawString(DF.format(new BigDecimal("10000")) + " pound/in2", 95, 60*(numY+1)+45);

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
                //restarts the unit count when reaches the end
                i = 0;
            }
        }
    }
}
