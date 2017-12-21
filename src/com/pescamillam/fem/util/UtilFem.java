package com.pescamillam.fem.util;

import org.apache.commons.math3.linear.FieldMatrix;
import org.apache.commons.math3.util.BigReal;

public class UtilFem {
    private static final BigReal MINUS_ONE = new BigReal("-1");
    
    public static FieldMatrix<BigReal> getAcceleration(FieldMatrix<BigReal> inverseMassMatrix, FieldMatrix<BigReal> stiffnessFieldMatrix,
            FieldMatrix<BigReal> force, FieldMatrix<BigReal> displacement) {
        return inverseMassMatrix.multiply(force.add(stiffnessFieldMatrix.multiply(displacement).scalarMultiply(MINUS_ONE)));
    }
}
