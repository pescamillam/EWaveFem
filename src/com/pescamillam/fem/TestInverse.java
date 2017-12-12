package com.pescamillam.fem;

import java.util.Arrays;

import org.apache.commons.math3.linear.FieldLUDecomposition;
import org.apache.commons.math3.linear.FieldMatrix;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.util.BigReal;

public class TestInverse {
    private static final BigReal MINUS_ONE = new BigReal("-1");
    private static final BigReal TWO = new BigReal("2");



    public static void main(String[] args) {
        
        BigReal deltaTime = new BigReal("0.00025");
        BigReal deltaTimeSquare = deltaTime.multiply(deltaTime);
        
        BigReal area = new BigReal("1");
        BigReal density = new BigReal("0.00073");
        BigReal elasticity = new BigReal("30000000");
        BigReal longitude = new BigReal("100");
        
        BigReal[][] stiffness = new BigReal[3][3];
        stiffness[0][0] = new BigReal("1");
        stiffness[0][1] = MINUS_ONE;
        stiffness[0][2] = new BigReal("0");
        
        stiffness[1][0] = MINUS_ONE;
        stiffness[1][1] = TWO;
        stiffness[1][2] = MINUS_ONE;
        
        stiffness[2][0] = new BigReal("0");
        stiffness[2][1] = MINUS_ONE;
        stiffness[2][2] = new BigReal("1");
        FieldMatrix<BigReal> stiffnessFieldMatrix = MatrixUtils.createFieldMatrix(stiffness).scalarMultiply(area.multiply(elasticity).divide(longitude));
        
        BigReal[][] massMatrix = new BigReal[3][3];
        massMatrix[0][0] = massMatrix[2][2] = BigReal.ONE;
        massMatrix[0][1] = massMatrix[0][2] = massMatrix[1][0] = massMatrix[1][2] = massMatrix[2][0] = massMatrix[2][1] = BigReal.ZERO;
        massMatrix[1][1] = TWO;

        BigReal[] force0Vector = new BigReal[3];
        Arrays.fill(force0Vector, BigReal.ZERO);
        force0Vector[2] = new BigReal("1000");
        FieldMatrix<BigReal> force = MatrixUtils.createColumnFieldMatrix(force0Vector);
        
        FieldMatrix<BigReal> massFieldMatrix = MatrixUtils.createFieldMatrix(massMatrix).scalarMultiply(density.multiply(area).multiply(longitude).divide(TWO));
        FieldMatrix<BigReal> inverseMassMatrix = new FieldLUDecomposition<BigReal>(massFieldMatrix).getSolver().getInverse();
        
        BigReal[] displacementVector = new BigReal[3];
        Arrays.fill(displacementVector, BigReal.ZERO);
        FieldMatrix<BigReal>[] displacement = new FieldMatrix[100];
        displacement[0] = MatrixUtils.createColumnFieldMatrix(displacementVector);
        printMatrix(force);
        
        FieldMatrix<BigReal>[] acceleration = new FieldMatrix[100];
        acceleration[0] = inverseMassMatrix.multiply(force.add(stiffnessFieldMatrix.multiply(displacement[0]).scalarMultiply(MINUS_ONE)));
        
        System.out.println("===acceleration===");
        printMatrix(acceleration[0]);
        
        FieldMatrix<BigReal> displacementM1 = acceleration[0].scalarMultiply(deltaTimeSquare.divide(TWO));
        printMatrix(displacementM1);

        displacement[1] = inverseMassMatrix.multiply(force.scalarMultiply(deltaTimeSquare)
                            .add(massFieldMatrix.scalarMultiply(TWO).add(stiffnessFieldMatrix.scalarMultiply(deltaTimeSquare.multiply(MINUS_ONE))).multiply(displacement[0]))
                            .add(massFieldMatrix.multiply(displacementM1).scalarMultiply(MINUS_ONE))
                            );
        
        printMatrix(displacement[1]);
        
        acceleration[1] = inverseMassMatrix.multiply(force.add(stiffnessFieldMatrix.multiply(displacement[1]).scalarMultiply(MINUS_ONE)));
        System.out.println("=== a1 ===");
        printMatrix(acceleration[1]);
        
        
        for (int i = 2; i < 100; i++) {
            displacement[i] = inverseMassMatrix.multiply(force.scalarMultiply(deltaTimeSquare)
                                .add(massFieldMatrix.scalarMultiply(TWO).add(stiffnessFieldMatrix.scalarMultiply(deltaTimeSquare.multiply(MINUS_ONE))).multiply(displacement[i-1]))
                                .add(massFieldMatrix.multiply(displacement[i-2]).scalarMultiply(MINUS_ONE))
                                );
            System.out.println("=== d" + i + " ===");
            printMatrix(displacement[i]);

            FieldMatrix<BigReal>[] speed = new FieldMatrix[100];
            speed[i-1] = displacement[i].add(displacement[i-2].scalarMultiply(MINUS_ONE)).scalarMultiply(BigReal.ONE.divide(TWO.multiply(deltaTime)));

            System.out.println("=== v" + (i-1) + " ===");
            printMatrix(speed[i-1]);
            
            acceleration[i] = inverseMassMatrix.multiply(force.add(stiffnessFieldMatrix.multiply(displacement[i]).scalarMultiply(MINUS_ONE)));
            System.out.println("=== a" + i + " ===");
            printMatrix(acceleration[i]);
        }
        
//        FieldMatrix<BigReal> c = new FieldLUDecomposition<BigReal>(stiffnessFieldMatrix).getSolver().getInverse();
//        
//        printMatrix(stiffnessFieldMatrix);
//        
//        printMatrix(c);
//        
//        printMatrix(stiffnessFieldMatrix.multiply(c));
    }
    


    private static void printMatrix(FieldMatrix<BigReal> fieldMatrix) {
        
        System.out.println("====  ====");
        StringBuilder str = new StringBuilder();
        for (BigReal[] column : fieldMatrix.getData()) {
            
            for (BigReal unit : column) {
                str.append(unit.bigDecimalValue().toPlainString()).append("\t");
            }
            str.append("\n");
        }
        System.out.println(str.toString());
        
    }
}
