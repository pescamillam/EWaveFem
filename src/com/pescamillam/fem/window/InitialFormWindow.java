package com.pescamillam.fem.window;

import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JTextField;

import com.pescamillam.fem.Application;
import com.pescamillam.fem.model.InputValues;
import com.pescamillam.fem.util.Constants;

public class InitialFormWindow {

    private static JLabel numTimesLabel = new JLabel("Número de iteraciones");
    private static JTextField numTimesText = new JTextField(String.valueOf(Constants.NUM_TIMES));

    private static JLabel numXLabel = new JLabel("Número de elementos a lo ancho");
    private static JTextField numXText = new JTextField(String.valueOf(Constants.NUM_X));

    private static JLabel numYLabel = new JLabel("Número de elementos a lo alto");
    private static JTextField numYText = new JTextField(String.valueOf(Constants.NUM_Y));

    private static JLabel thicknessLabel = new JLabel("Grosor");
    private static JTextField thicknessText = new JTextField(Constants.THICKNESS);

    private static JLabel elasticityLabel = new JLabel("Elasticidad");
    private static JTextField elasticityText = new JTextField(Constants.ELASTICITY);

    private static JLabel densityLabel = new JLabel("Densidad");
    private static JTextField densityText = new JTextField(Constants.DENSITY);

    private static JLabel areaLabel = new JLabel("Area");
    private static JTextField areaText = new JTextField(Constants.AREA);

    private static JLabel poissonLabel = new JLabel("Poisson");
    private static JTextField poissonText = new JTextField(Constants.POISSON);

    private static JLabel deltaTimeLabel = new JLabel("Delta Time");
    private static JTextField deltaTimeText = new JTextField(Constants.DELTA_TIME);
    
    private static JFrame frame;

    public static void createFormInitialValues() {
        JPanel p = new JPanel(new GridLayout(12, 2, 10, 10));
        
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

        JButton button = new JButton("Ejecutar");
        p.add(button);

        button.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                button.setEnabled(false);
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
    
    private static class ExecutionThread implements Runnable {
        @Override
        public void run() {
            // TODO Auto-generated method stub
            Application.executeFiniteElementProcess(new InputValues.Builder()
                    .withArea(areaText.getText())
                    .withThickness(thicknessText.getText())
                    .withDensity(densityText.getText())
                    .withElasticity(elasticityText.getText())
                    .withPoisson(poissonText.getText())
                    .withDeltaTime(deltaTimeText.getText())
                    .withNumX(numXText.getText())
                    .withNumY(numYText.getText())
                    .withNumTimes(numTimesText.getText())
                    .build());
        }
    }
}
