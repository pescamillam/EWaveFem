package com.pescamillam.fem.element.cst;

import static com.pescamillam.fem.util.Constants.TWO;
import static org.apache.commons.math3.util.BigReal.ONE;
import static org.apache.commons.math3.util.BigReal.ZERO;

import java.math.BigDecimal;
import java.util.List;

import org.apache.commons.math3.util.BigReal;

import com.pescamillam.fem.element.Cst;
import com.pescamillam.fem.element.Point;

public class ElementOne {

    public static BigReal[][] getLocalStiffnessMatrix(BigReal poisson) {

        BigReal[][] localMatrix = new BigReal[6][6];

        //(1-v)
        BigReal poisson1mv = poisson.negate().add(new BigReal("1"));
        
        //(1-2v)/2
        BigReal poisson1m2vo2 = poisson.multiply(new BigReal("2")).negate().add(new BigReal("1")).divide(new BigReal("2"));

        //element 1
        //
        // i
        // |\
        // | \
        // |__\
        // m   j

        //beta i: y_j - y_m
        BigReal betaIel1 = new BigReal("0");
        //beta j: y_m - y_i
        BigReal betaJel1 = new BigReal("30");
        //beta m: y_i - y_j
        BigReal betaMel1 = new BigReal("-30");
        
        //gamma_i: x_m - x_j
        BigReal gammaIel1 = new BigReal("-30");
        //gamma_j: x_i - x_m
        BigReal gammaJel1 = new BigReal("0");
        //gamma_m: x_j - x_i
        BigReal gammaMel1 = new BigReal("30");
        
        //1-1
        localMatrix[0][0] = betaIel1.multiply(betaIel1).multiply(poisson1mv).add(gammaIel1.multiply(gammaIel1).multiply(poisson1m2vo2));
        //1-2
        localMatrix[0][1] = betaIel1.multiply(betaIel1).multiply(poisson).add(betaIel1.multiply(gammaIel1).multiply(poisson1m2vo2));
        //1-3
        localMatrix[0][2] = betaIel1.multiply(betaJel1).multiply(poisson1mv).add(gammaIel1.multiply(gammaJel1).multiply(poisson1m2vo2));
        //1-4
        localMatrix[0][3] = betaIel1.multiply(gammaJel1).multiply(poisson).add(betaJel1.multiply(gammaIel1).multiply(poisson1m2vo2));
        //1-5
        localMatrix[0][4] = betaIel1.multiply(betaMel1).multiply(poisson1mv).add(gammaIel1.multiply(gammaMel1).multiply(poisson1m2vo2));
        //1-6
        localMatrix[0][5] = betaIel1.multiply(gammaMel1).multiply(poisson).add(betaMel1.multiply(gammaIel1).multiply(poisson1m2vo2));
        
        //2-1
        localMatrix[1][0] = localMatrix[0][1];
        //2-2
        localMatrix[1][1] = gammaIel1.multiply(gammaIel1).multiply(poisson1mv).add(betaIel1.multiply(betaIel1).multiply(poisson1m2vo2));
        //2-3
        localMatrix[1][2] = betaJel1.multiply(gammaIel1).multiply(poisson).add(betaIel1.multiply(gammaJel1).multiply(poisson1m2vo2));
        //2-4
        localMatrix[1][3] = gammaIel1.multiply(gammaJel1).multiply(poisson1mv).add(betaIel1.multiply(betaJel1).multiply(poisson1m2vo2));
        //2-5
        localMatrix[1][4] = betaMel1.multiply(gammaIel1).multiply(poisson).add(betaIel1.multiply(gammaMel1).multiply(poisson1m2vo2));
        //2-6
        localMatrix[1][5] = gammaIel1.multiply(gammaMel1).multiply(poisson1mv).add(betaIel1.multiply(betaMel1).multiply(poisson1m2vo2));

        //3-3
        localMatrix[2][2] = betaJel1.multiply(betaJel1).multiply(poisson1mv).add(gammaJel1.multiply(gammaJel1).multiply(poisson1m2vo2));
        //3-4
        localMatrix[2][3] = betaJel1.multiply(gammaJel1).multiply(poisson).add(betaJel1.multiply(gammaJel1).multiply(poisson1m2vo2));
        //3-5
        localMatrix[2][4] = betaJel1.multiply(betaMel1).multiply(poisson1mv).add(gammaJel1.multiply(gammaMel1).multiply(poisson1m2vo2));
        //3-6
        localMatrix[2][5] = betaJel1.multiply(gammaMel1).multiply(poisson).add(betaMel1.multiply(gammaJel1).multiply(poisson1m2vo2));

        //4-3
        localMatrix[3][2] = localMatrix[2][3];
        //4-4
        localMatrix[3][3] = gammaJel1.multiply(gammaJel1).multiply(poisson1mv).add(betaJel1.multiply(betaJel1).multiply(poisson1m2vo2));
        //4-5
        localMatrix[3][4] = betaMel1.multiply(gammaJel1).multiply(poisson).add(betaJel1.multiply(gammaMel1).multiply(poisson1m2vo2));
        //4-6
        localMatrix[3][5] = gammaIel1.multiply(gammaMel1).multiply(poisson1mv).add(betaJel1.multiply(betaMel1).multiply(poisson1m2vo2));

        //5-5
        localMatrix[4][4] = betaMel1.multiply(betaMel1).multiply(poisson1mv).add(gammaMel1.multiply(gammaMel1).multiply(poisson1m2vo2));
        //5-6
        localMatrix[4][5] = betaMel1.multiply(gammaMel1).multiply(poisson).add(betaMel1.multiply(gammaMel1).multiply(poisson1m2vo2));

        //6-5
        localMatrix[5][4] = localMatrix[4][5];
        //6-6
        localMatrix[5][5] = gammaMel1.multiply(gammaMel1).multiply(poisson1mv).add(betaMel1.multiply(betaMel1).multiply(poisson1m2vo2));
        
        return localMatrix;
    }

    public static void appendElementOneToStiffnessMatrix(List<Cst> elements,
            BigReal[][] stiffnessMatrix, BigReal[][] localMatrixElementOne, int i, int j, Integer numX) {
        Cst element = new Cst(new Point(new BigDecimal(i*10), new BigDecimal(j*10)), 
                new Point(new BigDecimal(i*10+10), new BigDecimal(j*10+10)), 
                new Point(new BigDecimal(i*10), new BigDecimal(j*10+10)));
        elements.add(element);

        //top left corner
        stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*j+i)*2] = stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*j+i)*2].add(localMatrixElementOne[0][0]);
        stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*j+i)*2+1] = stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*j+i)*2+1].add(localMatrixElementOne[0][1]);
        stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*j+i)*2] = stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*j+i)*2].add(localMatrixElementOne[0][1]);
        stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*j+i)*2+1] = stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*j+i)*2+1].add(localMatrixElementOne[1][1]);

        //bottom left corner
        stiffnessMatrix[((numX+1)*(j+1)+i)*2][((numX+1)*(j+1)+i)*2] = stiffnessMatrix[((numX+1)*(j+1)+i)*2][((numX+1)*(j+1)+i)*2].add(localMatrixElementOne[4][4]);
        stiffnessMatrix[((numX+1)*(j+1)+i)*2][((numX+1)*(j+1)+i)*2+1] = stiffnessMatrix[((numX+1)*(j+1)+i)*2][((numX+1)*(j+1)+i)*2+1].add(localMatrixElementOne[4][5]);
        stiffnessMatrix[((numX+1)*(j+1)+i)*2+1][((numX+1)*(j+1)+i)*2] = stiffnessMatrix[((numX+1)*(j+1)+i)*2+1][((numX+1)*(j+1)+i)*2].add(localMatrixElementOne[5][4]);
        stiffnessMatrix[((numX+1)*(j+1)+i)*2+1][((numX+1)*(j+1)+i)*2+1] = stiffnessMatrix[((numX+1)*(j+1)+i)*2+1][((numX+1)*(j+1)+i)*2+1].add(localMatrixElementOne[5][5]);

        //bottom right corner
        stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2][((numX+1)*(j+1)+(i+1))*2] = stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2][((numX+1)*(j+1)+(i+1))*2].add(localMatrixElementOne[2][2]);
        stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2][((numX+1)*(j+1)+(i+1))*2+1] = stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2][((numX+1)*(j+1)+(i+1))*2+1].add(localMatrixElementOne[2][3]);
        stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2+1][((numX+1)*(j+1)+(i+1))*2] = stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2+1][((numX+1)*(j+1)+(i+1))*2].add(localMatrixElementOne[2][3]);
        stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2+1][((numX+1)*(j+1)+(i+1))*2+1] = stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2+1][((numX+1)*(j+1)+(i+1))*2+1].add(localMatrixElementOne[3][3]);

        //ij
        stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*(j+1)+(i+1))*2] = stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*(j+1)+(i+1))*2].add(localMatrixElementOne[0][2]);
        stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*(j+1)+(i+1))*2+1] = stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*(j+1)+(i+1))*2+1].add(localMatrixElementOne[0][3]);
        stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*(j+1)+(i+1))*2] = stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*(j+1)+(i+1))*2].add(localMatrixElementOne[1][2]);
        stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*(j+1)+(i+1))*2+1] = stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*(j+1)+(i+1))*2+1].add(localMatrixElementOne[1][3]);

        stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2][((numX+1)*j+i)*2] = stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*(j+1)+(i+1))*2];
        stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2+1][((numX+1)*j+i)*2] = stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*(j+1)+(i+1))*2+1];
        stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2][((numX+1)*j+i)*2+1] = stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*(j+1)+(i+1))*2];
        stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2+1][((numX+1)*j+i)*2+1] = stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*(j+1)+(i+1))*2+1];

        //im
        stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*(j+1)+i)*2] = stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*(j+1)+i)*2].add(localMatrixElementOne[0][4]);
        stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*(j+1)+i)*2+1] = stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*(j+1)+i)*2+1].add(localMatrixElementOne[0][5]);
        stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*(j+1)+i)*2] = stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*(j+1)+i)*2].add(localMatrixElementOne[1][4]);
        stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*(j+1)+i)*2+1] = stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*(j+1)+i)*2+1].add(localMatrixElementOne[1][5]);

        stiffnessMatrix[((numX+1)*(j+1)+i)*2][((numX+1)*j+i)*2] = stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*(j+1)+i)*2];
        stiffnessMatrix[((numX+1)*(j+1)+i)*2+1][((numX+1)*j+i)*2] = stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*(j+1)+i)*2+1];
        stiffnessMatrix[((numX+1)*(j+1)+i)*2][((numX+1)*j+i)*2+1] = stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*(j+1)+i)*2];
        stiffnessMatrix[((numX+1)*(j+1)+i)*2+1][((numX+1)*j+i)*2+1] = stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*(j+1)+i)*2+1];
        
        //jm
        stiffnessMatrix[((numX+1)*(j+1)+i)*2][((numX+1)*(j+1)+(i+1))*2] = stiffnessMatrix[((numX+1)*(j+1)+i)*2][((numX+1)*(j+1)+(i+1))*2].add(localMatrixElementOne[2][4]);
        stiffnessMatrix[((numX+1)*(j+1)+i)*2+1][((numX+1)*(j+1)+(i+1))*2] = stiffnessMatrix[((numX+1)*(j+1)+i)*2+1][((numX+1)*(j+1)+(i+1))*2].add(localMatrixElementOne[3][4]);
        stiffnessMatrix[((numX+1)*(j+1)+i)*2][((numX+1)*(j+1)+(i+1))*2+1] = stiffnessMatrix[((numX+1)*(j+1)+i)*2][((numX+1)*(j+1)+(i+1))*2+1].add(localMatrixElementOne[2][5]);
        stiffnessMatrix[((numX+1)*(j+1)+i)*2+1][((numX+1)*(j+1)+(i+1))*2+1] = stiffnessMatrix[((numX+1)*(j+1)+i)*2+1][((numX+1)*(j+1)+(i+1))*2+1].add(localMatrixElementOne[3][5]);

        stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2][((numX+1)*(j+1)+i)*2] = stiffnessMatrix[((numX+1)*(j+1)+i)*2][((numX+1)*(j+1)+(i+1))*2];
        stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2][((numX+1)*(j+1)+i)*2+1] = stiffnessMatrix[((numX+1)*(j+1)+i)*2+1][((numX+1)*(j+1)+(i+1))*2];
        stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2+1][((numX+1)*(j+1)+i)*2] = stiffnessMatrix[((numX+1)*(j+1)+i)*2][((numX+1)*(j+1)+(i+1))*2+1];
        stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2+1][((numX+1)*(j+1)+i)*2+1] = stiffnessMatrix[((numX+1)*(j+1)+i)*2+1][((numX+1)*(j+1)+(i+1))*2+1];
    }

    public static void appendElementOneToMassMatrix(BigReal[][] massMatrix, int i, int j, Integer numX) {
        //element 1
        //
        // i
        // |\
        // | \
        // |__\
        // m   j

        //top left corner
        massMatrix[((numX+1)*j+i)*2][((numX+1)*j+i)*2] = massMatrix[((numX+1)*j+i)*2][((numX+1)*j+i)*2].add(TWO);
        massMatrix[((numX+1)*j+i)*2][((numX+1)*j+i)*2+1] = massMatrix[((numX+1)*j+i)*2][((numX+1)*j+i)*2+1].add(ZERO);
        massMatrix[((numX+1)*j+i)*2+1][((numX+1)*j+i)*2] = massMatrix[((numX+1)*j+i)*2+1][((numX+1)*j+i)*2].add(ZERO);
        massMatrix[((numX+1)*j+i)*2+1][((numX+1)*j+i)*2+1] = massMatrix[((numX+1)*j+i)*2+1][((numX+1)*j+i)*2+1].add(TWO);

        //bottom left corner
        massMatrix[((numX+1)*(j+1)+i)*2][((numX+1)*(j+1)+i)*2] = massMatrix[((numX+1)*(j+1)+i)*2][((numX+1)*(j+1)+i)*2].add(TWO);
        massMatrix[((numX+1)*(j+1)+i)*2][((numX+1)*(j+1)+i)*2+1] = massMatrix[((numX+1)*(j+1)+i)*2][((numX+1)*(j+1)+i)*2+1].add(ZERO);
        massMatrix[((numX+1)*(j+1)+i)*2+1][((numX+1)*(j+1)+i)*2] = massMatrix[((numX+1)*(j+1)+i)*2+1][((numX+1)*(j+1)+i)*2].add(ZERO);
        massMatrix[((numX+1)*(j+1)+i)*2+1][((numX+1)*(j+1)+i)*2+1] = massMatrix[((numX+1)*(j+1)+i)*2+1][((numX+1)*(j+1)+i)*2+1].add(TWO);

        //bottom right corner
        massMatrix[((numX+1)*(j+1)+(i+1))*2][((numX+1)*(j+1)+(i+1))*2] = massMatrix[((numX+1)*(j+1)+(i+1))*2][((numX+1)*(j+1)+(i+1))*2].add(TWO);
        massMatrix[((numX+1)*(j+1)+(i+1))*2][((numX+1)*(j+1)+(i+1))*2+1] = massMatrix[((numX+1)*(j+1)+(i+1))*2][((numX+1)*(j+1)+(i+1))*2+1].add(ZERO);
        massMatrix[((numX+1)*(j+1)+(i+1))*2+1][((numX+1)*(j+1)+(i+1))*2] = massMatrix[((numX+1)*(j+1)+(i+1))*2+1][((numX+1)*(j+1)+(i+1))*2].add(ZERO);
        massMatrix[((numX+1)*(j+1)+(i+1))*2+1][((numX+1)*(j+1)+(i+1))*2+1] = massMatrix[((numX+1)*(j+1)+(i+1))*2+1][((numX+1)*(j+1)+(i+1))*2+1].add(TWO);

        //ij
        massMatrix[((numX+1)*j+i)*2][((numX+1)*(j+1)+(i+1))*2] = massMatrix[((numX+1)*j+i)*2][((numX+1)*(j+1)+(i+1))*2].add(ONE);
        massMatrix[((numX+1)*j+i)*2][((numX+1)*(j+1)+(i+1))*2+1] = massMatrix[((numX+1)*j+i)*2][((numX+1)*(j+1)+(i+1))*2+1].add(BigReal.ZERO);
        massMatrix[((numX+1)*j+i)*2+1][((numX+1)*(j+1)+(i+1))*2] = massMatrix[((numX+1)*j+i)*2+1][((numX+1)*(j+1)+(i+1))*2].add(BigReal.ZERO);
        massMatrix[((numX+1)*j+i)*2+1][((numX+1)*(j+1)+(i+1))*2+1] = massMatrix[((numX+1)*j+i)*2+1][((numX+1)*(j+1)+(i+1))*2+1].add(ONE);

        massMatrix[((numX+1)*(j+1)+(i+1))*2][((numX+1)*j+i)*2] = massMatrix[((numX+1)*j+i)*2][((numX+1)*(j+1)+(i+1))*2];
        massMatrix[((numX+1)*(j+1)+(i+1))*2+1][((numX+1)*j+i)*2] = massMatrix[((numX+1)*j+i)*2][((numX+1)*(j+1)+(i+1))*2+1];
        massMatrix[((numX+1)*(j+1)+(i+1))*2][((numX+1)*j+i)*2+1] = massMatrix[((numX+1)*j+i)*2+1][((numX+1)*(j+1)+(i+1))*2];
        massMatrix[((numX+1)*(j+1)+(i+1))*2+1][((numX+1)*j+i)*2+1] = massMatrix[((numX+1)*j+i)*2+1][((numX+1)*(j+1)+(i+1))*2+1];

        //im
        massMatrix[((numX+1)*j+i)*2][((numX+1)*(j+1)+i)*2] = massMatrix[((numX+1)*j+i)*2][((numX+1)*(j+1)+i)*2].add(ONE);
        massMatrix[((numX+1)*j+i)*2][((numX+1)*(j+1)+i)*2+1] = massMatrix[((numX+1)*j+i)*2][((numX+1)*(j+1)+i)*2+1].add(ZERO);
        massMatrix[((numX+1)*j+i)*2+1][((numX+1)*(j+1)+i)*2] = massMatrix[((numX+1)*j+i)*2+1][((numX+1)*(j+1)+i)*2].add(ZERO);
        massMatrix[((numX+1)*j+i)*2+1][((numX+1)*(j+1)+i)*2+1] = massMatrix[((numX+1)*j+i)*2+1][((numX+1)*(j+1)+i)*2+1].add(ONE);

        massMatrix[((numX+1)*(j+1)+i)*2][((numX+1)*j+i)*2] = massMatrix[((numX+1)*j+i)*2][((numX+1)*(j+1)+i)*2];
        massMatrix[((numX+1)*(j+1)+i)*2+1][((numX+1)*j+i)*2] = massMatrix[((numX+1)*j+i)*2][((numX+1)*(j+1)+i)*2+1];
        massMatrix[((numX+1)*(j+1)+i)*2][((numX+1)*j+i)*2+1] = massMatrix[((numX+1)*j+i)*2+1][((numX+1)*(j+1)+i)*2];
        massMatrix[((numX+1)*(j+1)+i)*2+1][((numX+1)*j+i)*2+1] = massMatrix[((numX+1)*j+i)*2+1][((numX+1)*(j+1)+i)*2+1];
        
        //jm
        massMatrix[((numX+1)*(j+1)+i)*2][((numX+1)*(j+1)+(i+1))*2] = massMatrix[((numX+1)*(j+1)+i)*2][((numX+1)*(j+1)+(i+1))*2].add(ONE);
        massMatrix[((numX+1)*(j+1)+i)*2+1][((numX+1)*(j+1)+(i+1))*2] = massMatrix[((numX+1)*(j+1)+i)*2+1][((numX+1)*(j+1)+(i+1))*2].add(ZERO);
        massMatrix[((numX+1)*(j+1)+i)*2][((numX+1)*(j+1)+(i+1))*2+1] = massMatrix[((numX+1)*(j+1)+i)*2][((numX+1)*(j+1)+(i+1))*2+1].add(ZERO);
        massMatrix[((numX+1)*(j+1)+i)*2+1][((numX+1)*(j+1)+(i+1))*2+1] = massMatrix[((numX+1)*(j+1)+i)*2+1][((numX+1)*(j+1)+(i+1))*2+1].add(ONE);

        massMatrix[((numX+1)*(j+1)+(i+1))*2][((numX+1)*(j+1)+i)*2] = massMatrix[((numX+1)*(j+1)+i)*2][((numX+1)*(j+1)+(i+1))*2];
        massMatrix[((numX+1)*(j+1)+(i+1))*2][((numX+1)*(j+1)+i)*2+1] = massMatrix[((numX+1)*(j+1)+i)*2+1][((numX+1)*(j+1)+(i+1))*2];
        massMatrix[((numX+1)*(j+1)+(i+1))*2+1][((numX+1)*(j+1)+i)*2] = massMatrix[((numX+1)*(j+1)+i)*2][((numX+1)*(j+1)+(i+1))*2+1];
        massMatrix[((numX+1)*(j+1)+(i+1))*2+1][((numX+1)*(j+1)+i)*2+1] = massMatrix[((numX+1)*(j+1)+i)*2+1][((numX+1)*(j+1)+(i+1))*2+1];
    }
}
