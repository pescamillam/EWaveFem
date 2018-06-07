package com.pescamillam.fem;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import com.pescamillam.fem.util.Constants;
import com.pescamillam.fem.util.UtilWindow;
import org.apache.commons.math3.linear.FieldLUDecomposition;
import org.apache.commons.math3.linear.FieldMatrix;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.util.BigReal;

import com.pescamillam.fem.element.Cst;
import com.pescamillam.fem.element.cst.ElementOne;
import com.pescamillam.fem.element.cst.ElementTwo;
import com.pescamillam.fem.model.InputValues;
import com.pescamillam.fem.util.UtilWriter;
import com.pescamillam.fem.window.InitialFormWindow;

import static com.pescamillam.fem.util.Constants.FOUR;
import static com.pescamillam.fem.util.Constants.MINUS_ONE;
import static com.pescamillam.fem.util.Constants.TWELVE;
import static com.pescamillam.fem.util.Constants.TWO;
import static com.pescamillam.fem.util.Constants.DF;
import static org.apache.commons.math3.util.BigReal.ONE;
import static org.apache.commons.math3.util.BigReal.ZERO;

/**
 * Class that executes the application
 * 
 * @author Peter Escamilla (pescamilla@unab.edu.co)
 */
public class Application {

    /**
     * Method that starts the application
     * 
     * @param args not used
     */
    public static void main(String[] args) {
        InitialFormWindow.createFormInitialValues();
    }
    
    /**
     * Executes the finite element method with given parameters
     * 
     * @param input Initial parameters to be used by Finite element method
     */
    @SuppressWarnings("unchecked")
    public static void executeFiniteElementProcess(InputValues input) {

        //Records the initial execution time
        Long startTime = System.nanoTime();
        
        //parameters
        //thickness t
        BigReal thickness = new BigReal(input.getThickness());
        //elasticity module E
        BigReal elasticity = new BigReal(input.getElasticity());
        //density p
        BigReal density = new BigReal(input.getDensity());
        //Area A
        BigReal area = new BigReal(input.getArea());
        //poisson ratio v
        BigReal poisson = new BigReal(input.getPoisson());
        
        //delta time
        BigReal deltaTime = new BigReal(input.getDeltaTime());
        //delta time squared (used to simplify operations
        BigReal deltaTimeSquare = deltaTime.multiply(deltaTime);

        //number of horizontal elements
        Integer numX = new Integer(input.getNumX());
        //number of vertical elements
        Integer numY = new Integer(input.getNumY());
        //number of execution times
        Integer numTimes = new Integer(input.getNumTimes());

        //prints initial parameters
        writeToFile("=== Constants ===");
        writeToFile("area: " + area.bigDecimalValue().toPlainString());
        writeToFile("thickness: " + thickness.bigDecimalValue().toPlainString());
        writeToFile("elasticity: " + elasticity.bigDecimalValue().toPlainString());
        writeToFile("density: " + density.bigDecimalValue().toPlainString());
        writeToFile("poisson: " + poisson.bigDecimalValue().toPlainString());
        writeToFile("\n\n");

        //(1-2v)
        BigReal poisson1m2v = poisson.multiply(TWO).negate().add(ONE);

        //(1+v)
        BigReal poisson1pv = poisson.add(ONE);
        
        //variables
        //acceleration
        FieldMatrix<BigReal>[] acceleration = new FieldMatrix[numTimes];

        List<Cst> elements = new ArrayList<>();
        
        //Stiffness matrix
        BigReal[][] stiffnessMatrix = new BigReal[2* numX* numY+2* numX+2* numY+2][2* numX* numY+2* numX+2* numY+2];
        
        //mass matrix
        BigReal[][] massMatrix = new BigReal[2* numX* numY+2* numX+2* numY+2][2* numX* numY+2* numX+2* numY+2];

        for (BigReal[] row : stiffnessMatrix) {
            Arrays.fill(row, ZERO);
        }

        for (BigReal[] row : massMatrix) {
            Arrays.fill(row, ZERO);
        }
        

        //Stiffness matrix for local element one
        BigReal[][] localMatrixElemOne = ElementOne.getLocalStiffnessMatrix(poisson);

        //Stiffness matrix for local element two
        BigReal[][] localMatrixElemTwo = ElementTwo.getLocalStiffnessMatrix(poisson);

        //points matrix
        for (int i = 0; i < numX; i++) {
            for (int j = 0; j < numY; j++) {
                //adds local matrix to general stiffness matrix
                ElementOne.appendElementOneToStiffnessMatrix(elements, stiffnessMatrix, localMatrixElemOne, i, j, numX);
                ElementTwo.appendElementTwoToStiffnessMatrix(elements, stiffnessMatrix, localMatrixElemTwo, i, j, numX);

                //adds local mass matrix to general mass matrix
                ElementOne.appendElementOneToMassMatrix(massMatrix, i, j, numX);
                ElementTwo.appendElementTwoToMassMatrix(massMatrix, i, j, numX);
            }
        }

        //Field matrix to easily process the mass matrix
        FieldMatrix<BigReal> massFieldMatrix = MatrixUtils.createFieldMatrix(massMatrix);
        massFieldMatrix = massFieldMatrix.scalarMultiply(density.multiply(thickness).multiply(area).divide(TWELVE));

        writeToFile("=== Starting inversing mass matrix ===");
        Long startInverseTime = System.nanoTime();

        //calculates the inverse mass matrix
        FieldMatrix<BigReal> inverseMassMatrix = new FieldLUDecomposition<BigReal>(massFieldMatrix).getSolver().getInverse();
        writeToFile("=== Finished inversing mass matrix " + (System.nanoTime() - startInverseTime)/1000000000.0 + "s ===");

        FieldMatrix<BigReal> stiffnessFieldMatrix = MatrixUtils.createFieldMatrix(stiffnessMatrix);

        //multiplies the stiffness matrix with the constant values
        stiffnessFieldMatrix = stiffnessFieldMatrix.scalarMultiply(thickness.multiply(elasticity).divide(FOUR.multiply(area).multiply(poisson1pv).multiply(poisson1m2v)));
        writeToFile("=== stiffness matrix ===");
        printMatrix(stiffnessFieldMatrix);

        //creates the force vector
//        BigReal ONE_HUNDRED = new BigReal("3000");
        BigReal[] force0Vector = new BigReal[2 * numX * numY + 2 * numX + 2 * numY + 2];
        Arrays.fill(force0Vector, ZERO);
//        for (int i = 2; i < numX-2 + 1; i++) {
//            force0Vector[i*2+1] = ONE_HUNDRED;
//        }

        FieldMatrix<BigReal>[] force = new FieldMatrix[numTimes];
        for (int i = 0; i < numTimes; i++) {

            //the force vector is a sinosoidal function
//            force[i] = MatrixUtils.createColumnFieldMatrix(force0Vector)
//                    .scalarMultiply(new BigReal(Math.sin(i/9.0)));
//            if (i < 30) {
                force[i] = MatrixUtils.createColumnFieldMatrix(force0Vector);
//            } else if (i < 70) {
//                force[i] = MatrixUtils.createColumnFieldMatrix(force0Vector)
//                        .scalarMultiply(new BigReal(i-30))
//                        .scalarMultiply(Constants.MINUS_ONE);
//            } else {
//                force[i] = MatrixUtils.createColumnFieldMatrix(force0Vector)
//                        .scalarMultiply(new BigReal("20"))
//                        .scalarMultiply(Constants.MINUS_ONE);
//            }
        }

        //acceleration = mass inverse * (F0 - [K]{d0}) as {d0} = 0 
        //acceleration = mass inverse * (F0)
        writeToFile("=== acceleration 0 ===");

        //calculates acceleration for time 0
        acceleration[0] = inverseMassMatrix.multiply(force[0]);
        printMatrix(acceleration[0]);

        //calculates displacement for time -1
        FieldMatrix<BigReal> displacementM1 = acceleration[0].scalarMultiply(deltaTimeSquare.divide(TWO));

        writeToFile("=== displacement -1 ===");
        printMatrix(displacementM1);
        
        //calculates displacement for time 0
        BigReal[] displacement0 = new BigReal[2* numX* numY+2* numX+2* numY+2];
        Arrays.fill(displacement0, ZERO);
        FieldMatrix<BigReal>[] displacement = new FieldMatrix[numTimes];
        FieldMatrix<BigReal>[] speed = new FieldMatrix[numTimes];
        displacement[0] = MatrixUtils.createColumnFieldMatrix(displacement0);

        writeToFile("=== displacement 0 ===");
        printMatrix(displacement[0]);

        //calculates displacement for time 1
        displacement[1] = inverseMassMatrix.multiply(
                force[0].scalarMultiply(deltaTimeSquare)
                //.add(massFieldMatrix.scalarMultiply(TWO_BIG_REAL).add(stiffnessFieldMatrix.scalarMultiply(deltaTime.multiply(deltaTime))).multiply(displacementM1))
                .add(massFieldMatrix.multiply(displacementM1).scalarMultiply(MINUS_ONE))
                );
        
        for (int i = 1; i < numX*2+2; i+=2) {
            displacement[1].setEntry(i, 0, new BigReal(Math.sin(1/9.0)).multiply(new BigReal("0.9")));
        }

        writeToFile("=== displacement 1 ===");
        printMatrix(displacement[1]);

        for (int i = 2; i < numTimes; i++) {

            //calculates displacement for time i
            displacement[i] = inverseMassMatrix.multiply(force[i-1].scalarMultiply(deltaTimeSquare)
                                .add(massFieldMatrix.scalarMultiply(TWO).add(stiffnessFieldMatrix.scalarMultiply(deltaTimeSquare.multiply(MINUS_ONE))).multiply(displacement[i-1]))
                                .add(massFieldMatrix.multiply(displacement[i-2]).scalarMultiply(MINUS_ONE))
                                );
            
            displacement[i].setEntry(0, 0, BigReal.ZERO);
            displacement[i].setEntry(1, 0, BigReal.ZERO);
            
            displacement[i].setEntry((numX+1)*2, 0, BigReal.ZERO);
            displacement[i].setEntry((numX+1)*2+1, 0, BigReal.ZERO);
            
//            for (int j = 1; j < numX*2+2; j+=2) {
//                displacement[i].setEntry(j, 0, new BigReal(Math.sin(i/9.0)));
//            }
//            
//            for (int j = 0; j < numX*2+2; j+=2) {
//                displacement[i].setEntry(j, 0, BigReal.ZERO);
//            }

            writeToFile("==== displacement " + i + " ====");
            printMatrix(displacement[i]);

            //calculates speed for time i-1
            speed[i-1] = displacement[i].add(displacement[i-2].scalarMultiply(MINUS_ONE)).scalarMultiply(ONE.divide(TWO.multiply(deltaTime)));

            writeToFile("=== speed " + (i-1) + " ===");
            printMatrix(speed[i-1]);

            //calculates acceleration for time i
            acceleration[i] = inverseMassMatrix.multiply(force[i].add(stiffnessFieldMatrix.multiply(displacement[i]).scalarMultiply(MINUS_ONE)));
            writeToFile("=== acceleration " + i + " ===");
            printMatrix(acceleration[i]);

            writeToFile("=== force " + i + " ===");
            printMatrix(force[i]);
        }
        
        writeToFile("Total time: " + (System.nanoTime() - startTime)/1000000000.0 + "s");
        System.out.println("Total time: " + (System.nanoTime() - startTime)/1000000000.0 + "s");
        new Thread(new Runnable() {
            @Override
            public void run() {
                UtilWindow.printElements(elements, displacement, speed, acceleration, force, numX, numY, numTimes);
            }
        }).start();

        for (int i = 1; i < elements.size(); i+=10) {
            final int a = i;
            new Thread(new Runnable() {
                @Override
                public void run() {
                        UtilWindow.printMovingNode(elements, displacement, speed, acceleration, force, numTimes, a);
                }
            }).start();
        }
    }

    private static void printMatrix(FieldMatrix<BigReal> fieldMatrix) {
        
        writeToFile("====  ====");
        StringBuilder str = new StringBuilder();
        for (BigReal[] column : fieldMatrix.getData()) {
            
            for (BigReal unit : column) {
                str.append(DF.format(unit.bigDecimalValue())).append("\t");
            }
            str.append("\n");
        }
        writeToFile(str.toString());
        
    }

    private static void writeToFile(String string) {
        UtilWriter.writeToFile(string);
        UtilWriter.writeToFile("\n");
    }
}