package com.pescamillam.fem.window;

import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.*;

import com.pescamillam.fem.Application;
import com.pescamillam.fem.model.InputValues;
import com.pescamillam.fem.util.Constants;
import com.sun.org.apache.bcel.internal.classfile.ConstantString;
import com.sun.xml.internal.bind.v2.runtime.reflect.opt.Const;


/**
 * Class that shows the initial form to insert the constants to be used
 * 
 * @author Peter Escamilla (pescamilla@unab.edu.co)
 */
public class InitialFormWindow {

    /** Label for number of iterations */
    private static JLabel numTimesLabel = new JLabel("Número de iteraciones");
    /** Field text for number of iterations */
    private static JTextField numTimesText = new JTextField(String.valueOf(Constants.NUM_TIMES));

    private static JLabel widthLabel = new JLabel("Ancho del dominio (metros)");

    private static JTextField widthText = new JTextField(String.valueOf(Constants.DEFAULT_WIDTH));

    private static JLabel heightLabel = new JLabel("Alto del dominio (metros)");

    private static JTextField heightText = new JTextField(String.valueOf(Constants.DEFAULT_HEIGHT));

    /** Label for number of elements in X */
    private static JLabel numXLabel = new JLabel("Número de elementos a lo ancho");
    /** Field text for number of elements in X */
    private static JTextField numXText = new JTextField(String.valueOf(Constants.NUM_X));

    /** Label for number of elements in Y */
    private static JLabel numYLabel = new JLabel("Número de elementos a lo alto");
    /** Text field for number of elements in Y */
    private static JTextField numYText = new JTextField(String.valueOf(Constants.NUM_Y));

    /** Label for thickness field */
    private static JLabel thicknessLabel = new JLabel("Grosor (metros)");
    /** Text Field for thickness value */
    private static JTextField thicknessText = new JTextField(Constants.THICKNESS);

    /** Label for elasticity field */
    private static JLabel elasticityLabel = new JLabel("Elasticidad (N/m^2)");
    /** Text field for elasticity value */
    private static JTextField elasticityText = new JTextField(Constants.ELASTICITY);

    /** Label for density field */
    private static JLabel densityLabel = new JLabel("Densidad (Kg/m^3)");
    /** Text field for density value */
    private static JTextField densityText = new JTextField(Constants.DENSITY);

    /** Label for element area field */
    private static JLabel areaLabel = new JLabel("Area (m^2)");
    /** Text field for element area value */
    private static JTextField areaText = new JTextField(Constants.AREA);

    /** Label for poisson field */
    private static JLabel poissonLabel = new JLabel("Poisson");
    /** Text field for poisson value */
    private static JTextField poissonText = new JTextField(Constants.POISSON);

    /** Label for delta time field */
    private static JLabel deltaTimeLabel = new JLabel("Delta Time (segundos)");
    /** Text field for delta time value */
    private static JTextField deltaTimeText = new JTextField(Constants.DELTA_TIME);

    /** Label for the initial displacement field */
    private static JLabel initialDisplacementLabel = new JLabel("Desplazamiento inicial (metros)");
    /** Text field for the initial displacement value */
    private static JTextField initialDisplacementText = new JTextField(Constants.DEFAULT_INITIAL_DISPLACEMENT);

    /** Frame object where the form will be shown */
    private static JFrame frame;

    /** Creates the form with a squared layout */
    public static void createFormInitialValues() {
        //Assigns a grid layout of 12 rows and 2 columns
        JPanel p = new JPanel(new GridLayout(15, 2, 10, 10));

        addMargin();

        p.add(widthLabel);
        p.add(widthText);

        p.add(heightLabel);
        p.add(heightText);
        
        p.add(numTimesLabel);
        p.add(numTimesText);
        
        p.add(numXLabel);
        p.add(numXText);
        
        p.add(numYLabel);
        p.add(numYText);
        
        p.add(thicknessLabel);
        p.add(thicknessText);
        
        p.add(elasticityLabel);
        p.add(elasticityText);
        
        p.add(densityLabel);
        p.add(densityText);
        
        p.add(areaLabel);
        p.add(areaText);
        
        p.add(poissonLabel);
        p.add(poissonText);
        
        p.add(deltaTimeLabel);
        p.add(deltaTimeText);

        p.add(initialDisplacementLabel);
        p.add(initialDisplacementText);

        // button to start the process
        JButton button = new JButton("Ejecutar");
        p.add(button);

        // assigns the action of making the processing to the button
        button.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                button.setEnabled(false);
                //Starts the processing in background
                new Thread(new ExecutionThread()).start();
                JOptionPane.showMessageDialog(frame, "Ejecutando por favor espere");
                frame.dispose();
            }
        });
        //Create and set up the window.
        frame = new JFrame("Ingresar valores iniciales");

        //Set up the content pane.
        p.setOpaque(true);  //content panes must be opaque
        frame.setContentPane(p);

        //Display the window.
        frame.pack();
        frame.setVisible(true);
    }

    private static void addMargin() {
        widthLabel.setBorder(BorderFactory.createEmptyBorder(0,20,0,20));
        heightLabel.setBorder(BorderFactory.createEmptyBorder(0,20,0,20));
        numTimesLabel.setBorder(BorderFactory.createEmptyBorder(0,20,0,20));
        numXLabel.setBorder(BorderFactory.createEmptyBorder(0,20,0,20));
        numYLabel.setBorder(BorderFactory.createEmptyBorder(0,20,0,20));
        thicknessLabel.setBorder(BorderFactory.createEmptyBorder(0,20,0,20));
        elasticityLabel.setBorder(BorderFactory.createEmptyBorder(0,20,0,20));
        densityLabel.setBorder(BorderFactory.createEmptyBorder(0,20,0,20));
        areaLabel.setBorder(BorderFactory.createEmptyBorder(0,20,0,20));
        poissonLabel.setBorder(BorderFactory.createEmptyBorder(0,20,0,20));
        deltaTimeLabel.setBorder(BorderFactory.createEmptyBorder(0,20,0,20));
        initialDisplacementLabel.setBorder(BorderFactory.createEmptyBorder(0,20,0,20));
    }

    //Thread with the processing part to be executed in background
    private static class ExecutionThread implements Runnable {
        @Override
        public void run() {
            // Starts the processing with the entered values
            Application.executeFiniteElementProcess(new InputValues.Builder()
                    .withWidth(widthText.getText())
                    .withHeight(heightText.getText())
                    .withArea(areaText.getText())
                    .withThickness(thicknessText.getText())
                    .withDensity(densityText.getText())
                    .withElasticity(elasticityText.getText())
                    .withPoisson(poissonText.getText())
                    .withDeltaTime(deltaTimeText.getText())
                    .withNumX(numXText.getText())
                    .withNumY(numYText.getText())
                    .withNumTimes(numTimesText.getText())
                    .withInitialDisplacement(initialDisplacementText.getText())
                    .build());
        }
    }
}
