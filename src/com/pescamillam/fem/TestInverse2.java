package com.pescamillam.fem;

import java.io.IOException;
import java.util.Arrays;

import org.apache.commons.math3.linear.FieldLUDecomposition;
import org.apache.commons.math3.linear.FieldMatrix;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.util.BigReal;

import com.pescamillam.fem.util.UtilFem;
import com.pescamillam.fem.util.UtilWriter;

public class TestInverse2 {
    private static final BigReal MINUS_ONE = new BigReal("-1");
    private static final BigReal TWO = new BigReal("2");

    @SuppressWarnings("unchecked")
    public static void main(String[] args) throws IOException {
        BigReal deltaTime = new BigReal("0.00025");
        BigReal deltaTimeSquare = deltaTime.multiply(deltaTime);
        
        BigReal area = new BigReal("1");
        BigReal density = new BigReal("0.00073");
        BigReal elasticity = new BigReal("30000000");
        BigReal longitude = new BigReal("100");
        writeToFile("=== Constants ===");
        writeToFile("area: " + area.bigDecimalValue().toPlainString());
        writeToFile("density: " + density.bigDecimalValue().toPlainString());
        writeToFile("elasticity: " + elasticity.bigDecimalValue().toPlainString());
        writeToFile("longitude: " + longitude.bigDecimalValue().toPlainString());
        writeToFile("\n");
        
        BigReal[][] stiffness = new BigReal[2][2];
        
        
        stiffness[0][0] = TWO;
        stiffness[0][1] = MINUS_ONE;
        
        stiffness[1][0] = MINUS_ONE;
        stiffness[1][1] = BigReal.ONE;
        FieldMatrix<BigReal> stiffnessFieldMatrix = MatrixUtils.createFieldMatrix(stiffness).scalarMultiply(area.multiply(elasticity).divide(longitude));
        writeToFile("===stiffness matrix===");
        printMatrix(stiffnessFieldMatrix);
        
        BigReal[][] massMatrix = new BigReal[2][2];
        massMatrix[1][1] = BigReal.ONE;
        massMatrix[0][1] = massMatrix[1][0] = BigReal.ZERO;
        massMatrix[0][0] = TWO;

        BigReal[] force0Vector = new BigReal[2];
        Arrays.fill(force0Vector, BigReal.ZERO);
        force0Vector[1] = new BigReal("1000");
        FieldMatrix<BigReal> force = MatrixUtils.createColumnFieldMatrix(force0Vector);
        
        FieldMatrix<BigReal> massFieldMatrix = MatrixUtils.createFieldMatrix(massMatrix).scalarMultiply(density.multiply(area).multiply(longitude).divide(TWO));
        writeToFile("===mass matrix===");
        printMatrix(massFieldMatrix);
        FieldMatrix<BigReal> inverseMassMatrix = new FieldLUDecomposition<BigReal>(massFieldMatrix).getSolver().getInverse();
        writeToFile("===inverse mass matrix===");
        printMatrix(inverseMassMatrix);
        
        
        BigReal[] displacementVector = new BigReal[2];
        Arrays.fill(displacementVector, BigReal.ZERO);
        FieldMatrix<BigReal>[] displacement = new FieldMatrix[100];
        displacement[0] = MatrixUtils.createColumnFieldMatrix(displacementVector);
        printMatrix(force);
        
        FieldMatrix<BigReal>[] acceleration = new FieldMatrix[100];
        acceleration[0] = UtilFem.getAcceleration(inverseMassMatrix, stiffnessFieldMatrix, force, displacement[0]);
        
        writeToFile("===acceleration0===");
        printMatrix(acceleration[0]);
        
        FieldMatrix<BigReal> displacementM1 = acceleration[0].scalarMultiply(deltaTimeSquare.divide(TWO));
        writeToFile("===displacement-1===");
        printMatrix(displacementM1);

        displacement[1] = inverseMassMatrix.multiply(force.scalarMultiply(deltaTimeSquare)
                            .add(massFieldMatrix.scalarMultiply(TWO).add(stiffnessFieldMatrix.scalarMultiply(deltaTimeSquare.multiply(MINUS_ONE))).multiply(displacement[0]))
                            .add(massFieldMatrix.multiply(displacementM1).scalarMultiply(MINUS_ONE))
                            );
        
        writeToFile("===displacement1===");
        printMatrix(displacement[1]);
        
        acceleration[1] = inverseMassMatrix.multiply(force.add(stiffnessFieldMatrix.multiply(displacement[1]).scalarMultiply(MINUS_ONE)));
        writeToFile("=== a1 ===");
        printMatrix(acceleration[1]);
        
        
        for (int i = 2; i < 11; i++) {
            displacement[i] = inverseMassMatrix.multiply(force.scalarMultiply(deltaTimeSquare)
                                .add(massFieldMatrix.scalarMultiply(TWO).add(stiffnessFieldMatrix.scalarMultiply(deltaTimeSquare.multiply(MINUS_ONE))).multiply(displacement[i-1]))
                                .add(massFieldMatrix.multiply(displacement[i-2]).scalarMultiply(MINUS_ONE))
                                );
            writeToFile("=== d" + i + " ===");
            printMatrix(displacement[i]);

            FieldMatrix<BigReal>[] speed = new FieldMatrix[100];
            speed[i-1] = displacement[i].add(displacement[i-2].scalarMultiply(MINUS_ONE)).scalarMultiply(BigReal.ONE.divide(TWO.multiply(deltaTime)));

            writeToFile("=== v" + (i-1) + " ===");
            printMatrix(speed[i-1]);
            
            acceleration[i] = inverseMassMatrix.multiply(force.add(stiffnessFieldMatrix.multiply(displacement[i]).scalarMultiply(MINUS_ONE)));
            writeToFile("=== a" + i + " ===");
            printMatrix(acceleration[i]);
        }
    }

    private static void writeToFile(String string) {
        UtilWriter.writeToFile(string);
        UtilWriter.writeToFile("\n");
    }

    private static void printMatrix(FieldMatrix<BigReal> fieldMatrix) throws IOException {
        
        writeToFile("====  ====");
        StringBuilder str = new StringBuilder();
        for (BigReal[] column : fieldMatrix.getData()) {
            
            for (BigReal unit : column) {
                str.append(unit.bigDecimalValue().toPlainString()).append("\t");
            }
            str.append("\n");
        }
        writeToFile(str.toString());
        
    }
}
