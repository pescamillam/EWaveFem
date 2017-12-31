package com.pescamillam.fem;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

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

public class Application {

    public static void main(String[] args) {
        InitialFormWindow.createFormInitialValues();
    }
    
    @SuppressWarnings("unchecked")
    public static void executeFiniteElementProcess(InputValues input) {

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
        BigReal deltaTimeSquare = deltaTime.multiply(deltaTime);
        
        Integer numX = new Integer(input.getNumX());
        Integer numY = new Integer(input.getNumY());
        Integer numTimes = new Integer(input.getNumTimes());

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
        

        BigReal[][] localMatrixElemOne = ElementOne.getLocalStiffnessMatrix(poisson);

        BigReal[][] localMatrixElemTwo = ElementTwo.getLocalStiffnessMatrix(poisson);

        //points matrix
        for (int i = 0; i < numX; i++) {
            for (int j = 0; j < numY; j++) {
                ElementOne.appendElementOneToStiffnessMatrix(elements, stiffnessMatrix, localMatrixElemOne, i, j, numX);
                ElementTwo.appendElementTwoToStiffnessMatrix(elements, stiffnessMatrix, localMatrixElemTwo, i, j, numX);

                ElementOne.appendElementOneToMassMatrix(massMatrix, i, j, numX);
                ElementTwo.appendElementTwoToMassMatrix(massMatrix, i, j, numX);
            }
        }

        FieldMatrix<BigReal> massFieldMatrix = MatrixUtils.createFieldMatrix(massMatrix);
        massFieldMatrix = massFieldMatrix.scalarMultiply(density.multiply(thickness).multiply(area).divide(TWELVE));

        writeToFile("=== Starting inversing mass matrix ===");
        Long startInverseTime = System.nanoTime();
        FieldMatrix<BigReal> inverseMassMatrix = new FieldLUDecomposition<BigReal>(massFieldMatrix).getSolver().getInverse();
        writeToFile("=== Finished inversing mass matrix " + (System.nanoTime() - startInverseTime)/1000000000.0 + "s ===");

        FieldMatrix<BigReal> stiffnessFieldMatrix = MatrixUtils.createFieldMatrix(stiffnessMatrix);
        stiffnessFieldMatrix = stiffnessFieldMatrix.scalarMultiply(thickness.multiply(elasticity).divide(FOUR.multiply(area).multiply(poisson1pv).multiply(poisson1m2v)));
        writeToFile("=== stiffness matrix ===");
        printMatrix(stiffnessFieldMatrix);

        //force vectors
        BigReal ONE_HUNDRED = new BigReal("30000");
        BigReal[] force1Vector = new BigReal[2* numX* numY+2* numX+2* numY+2];
        Arrays.fill(force1Vector, ZERO);
        BigReal[] force0Vector = new BigReal[2* numX* numY+2* numX+2* numY+2];
        Arrays.fill(force0Vector, ZERO);
        force0Vector[numX*2+2 + ((numX*2)/2+1+6)] = ONE_HUNDRED;
        BigReal[] vector0 = new BigReal[2* numX* numY+2* numX+2* numY+2];
        Arrays.fill(vector0, ZERO);

        FieldMatrix<BigReal>[] force = new FieldMatrix[numTimes];
        for (int i = 0; i < numTimes; i++) {
//            if (i < numTimes/3) {
                force[i] = MatrixUtils.createColumnFieldMatrix(force0Vector)
                                    .scalarMultiply(new BigReal(Math.sin(i/10.0)))
                ;
//            } else {
//                force[i] = MatrixUtils.createColumnFieldMatrix(vector0);
//            }
        }

        //acceleration = mass inverse * (F0 - [K]{d0}) as {d0} = 0 
        //acceleration = mass inverse * (F0)
        writeToFile("=== acceleration 0 ===");
        acceleration[0] = inverseMassMatrix.multiply(force[0]);
        printMatrix(acceleration[0]);

        //displacement -1
        FieldMatrix<BigReal> displacementM1 = acceleration[0].scalarMultiply(deltaTimeSquare.divide(TWO));

        writeToFile("=== displacement -1 ===");
        printMatrix(displacementM1);
        
        //displacement 1
        BigReal[] displacement0 = new BigReal[2* numX* numY+2* numX+2* numY+2];
        Arrays.fill(displacement0, ZERO);
//        displacement0[5] = BigReal.ONE;
        FieldMatrix<BigReal>[] displacement = new FieldMatrix[numTimes];
        FieldMatrix<BigReal>[] speed = new FieldMatrix[numTimes];
        displacement[0] = MatrixUtils.createColumnFieldMatrix(displacement0);

        writeToFile("=== displacement 0 ===");
        printMatrix(displacement[0]);
        
        displacement[1] = inverseMassMatrix.multiply(
                force[0].scalarMultiply(deltaTimeSquare)
                //.add(massFieldMatrix.scalarMultiply(TWO_BIG_REAL).add(stiffnessFieldMatrix.scalarMultiply(deltaTime.multiply(deltaTime))).multiply(displacementM1))
                .add(massFieldMatrix.multiply(displacementM1).scalarMultiply(MINUS_ONE))
                );

        writeToFile("=== displacement 1 ===");
        printMatrix(displacement[1]);

        for (int i = 2; i < numTimes; i++) {
            
            displacement[i] = inverseMassMatrix.multiply(force[i-1].scalarMultiply(deltaTimeSquare)
                                .add(massFieldMatrix.scalarMultiply(TWO).add(stiffnessFieldMatrix.scalarMultiply(deltaTimeSquare.multiply(MINUS_ONE))).multiply(displacement[i-1]))
                                .add(massFieldMatrix.multiply(displacement[i-2]).scalarMultiply(MINUS_ONE))
                                );
            
            writeToFile("==== displacement " + i + " ====");
            printMatrix(displacement[i]);
            
            speed[i-1] = displacement[i].add(displacement[i-2].scalarMultiply(MINUS_ONE)).scalarMultiply(ONE.divide(TWO.multiply(deltaTime)));

            writeToFile("=== speed " + (i-1) + " ===");
            printMatrix(speed[i-1]);

            acceleration[i] = inverseMassMatrix.multiply(force[i].add(stiffnessFieldMatrix.multiply(displacement[i]).scalarMultiply(MINUS_ONE)));
            writeToFile("=== acceleration " + i + " ===");
            printMatrix(acceleration[i]);
            
            writeToFile("=== force " + i + " ===");
            printMatrix(force[i]);
        }
        
        writeToFile("Total time: " + (System.nanoTime() - startTime)/1000000000.0 + "s");
        UtilWindow.printElements(elements, displacement, speed, acceleration, force, numX, numY, numTimes);
        UtilWindow.printMovingNode(elements, displacement, speed, acceleration, force, numTimes);
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