package com.pescamillam.fem.element.cst;

import static com.pescamillam.fem.util.Constants.TWO;
import static org.apache.commons.math3.util.BigReal.ONE;
import static org.apache.commons.math3.util.BigReal.ZERO;

import java.math.BigDecimal;
import java.util.List;

import org.apache.commons.math3.util.BigReal;

import com.pescamillam.fem.element.Cst;
import com.pescamillam.fem.element.Point;

/**
 * Util class to manage the element two type
 * 
 * @author Peter Escamilla (pescamilla@unab.edu.co)
 */
public class ElementTwo {

    /**
     * Gets the local stiffness matrix for an element two
     * @param poisson poisson module as given by input values
     * @return local stiffness matrix of the element two
     */
    public static BigReal[][] getLocalStiffnessMatrix(BigReal poisson) {

        BigReal[][] localMatrix = new BigReal[6][6];

        //(1-v)
        BigReal poisson1mv = poisson.negate().add(new BigReal("1"));
        
        //(1-2v)/2
        BigReal poisson1m2vo2 = poisson.multiply(new BigReal("2")).negate().add(new BigReal("1")).divide(new BigReal("2"));

        //element 2
        //
        // i___j
        //  \  |
        //   \ |
        //    \|
        //     m
        
        //beta i: y_j - y_m
        BigReal betaIel2 = new BigReal("-30");
        //beta j: y_m - y_i
        BigReal betaJel2 = new BigReal("30");
        //beta m: y_i - y_j
        BigReal betaMel2 = new BigReal("0");
        
        //gamma_i: x_m - x_j
        BigReal gammaIel2 = new BigReal("0");
        //gamma_j: x_i - x_m
        BigReal gammaJel2 = new BigReal("-30");
        //gamma_m: x_j - x_i
        BigReal gammaMel2 = new BigReal("30");
        
        //1-1
        localMatrix[0][0] = betaIel2.multiply(betaIel2).multiply(poisson1mv).add(gammaIel2.multiply(gammaIel2).multiply(poisson1m2vo2));
        //1-2
        localMatrix[0][1] = betaIel2.multiply(betaIel2).multiply(poisson).add(betaIel2.multiply(gammaIel2).multiply(poisson1m2vo2));
        //1-3
        localMatrix[0][2] = betaIel2.multiply(betaJel2).multiply(poisson1mv).add(gammaIel2.multiply(gammaJel2).multiply(poisson1m2vo2));
        //1-4
        localMatrix[0][3] = betaIel2.multiply(gammaJel2).multiply(poisson).add(betaJel2.multiply(gammaIel2).multiply(poisson1m2vo2));
        //1-5
        localMatrix[0][4] = betaIel2.multiply(betaMel2).multiply(poisson1mv).add(gammaIel2.multiply(gammaMel2).multiply(poisson1m2vo2));
        //1-6
        localMatrix[0][5] = betaIel2.multiply(gammaMel2).multiply(poisson).add(betaMel2.multiply(gammaIel2).multiply(poisson1m2vo2));

        //2-2
        localMatrix[1][1] = gammaIel2.multiply(gammaIel2).multiply(poisson1mv).add(betaIel2.multiply(betaIel2).multiply(poisson1m2vo2));
        //2-3
        localMatrix[1][2] = betaJel2.multiply(gammaIel2).multiply(poisson).add(betaIel2.multiply(gammaJel2).multiply(poisson1m2vo2));
        //2-4
        localMatrix[1][3] = gammaIel2.multiply(gammaJel2).multiply(poisson1mv).add(betaIel2.multiply(betaJel2).multiply(poisson1m2vo2));
        //2-5
        localMatrix[1][4] = betaMel2.multiply(gammaIel2).multiply(poisson).add(betaIel2.multiply(gammaMel2).multiply(poisson1m2vo2));
        //2-6
        localMatrix[1][5] = gammaIel2.multiply(gammaMel2).multiply(poisson1mv).add(betaIel2.multiply(betaMel2).multiply(poisson1m2vo2));

        //3-3
        localMatrix[2][2] = betaJel2.multiply(betaJel2).multiply(poisson1mv).add(gammaJel2.multiply(gammaJel2).multiply(poisson1m2vo2));
        //3-4
        localMatrix[2][3] = betaJel2.multiply(gammaJel2).multiply(poisson).add(betaJel2.multiply(gammaJel2).multiply(poisson1m2vo2));
        //3-5
        localMatrix[2][4] = betaJel2.multiply(betaMel2).multiply(poisson1mv).add(gammaJel2.multiply(gammaMel2).multiply(poisson1m2vo2));
        //3-6
        localMatrix[2][5] = betaJel2.multiply(gammaMel2).multiply(poisson).add(betaMel2.multiply(gammaJel2).multiply(poisson1m2vo2));

        //4-4
        localMatrix[3][3] = gammaJel2.multiply(gammaJel2).multiply(poisson1mv).add(betaJel2.multiply(betaJel2).multiply(poisson1m2vo2));
        //4-5
        localMatrix[3][4] = betaMel2.multiply(gammaJel2).multiply(poisson).add(betaJel2.multiply(gammaMel2).multiply(poisson1m2vo2));
        //4-6
        localMatrix[3][5] = gammaIel2.multiply(gammaMel2).multiply(poisson1mv).add(betaJel2.multiply(betaMel2).multiply(poisson1m2vo2));

        //5-5
        localMatrix[4][4] = betaMel2.multiply(betaMel2).multiply(poisson1mv).add(gammaMel2.multiply(gammaMel2).multiply(poisson1m2vo2));
        //5-6
        localMatrix[4][5] = betaMel2.multiply(gammaMel2).multiply(poisson).add(betaMel2.multiply(gammaMel2).multiply(poisson1m2vo2));

        //6-6
        localMatrix[5][5] = gammaMel2.multiply(gammaMel2).multiply(poisson1mv).add(betaMel2.multiply(betaMel2).multiply(poisson1m2vo2));
        return localMatrix;
    }

    /**
     * Adds a local matrix element two to the global stiffness matrix
     * 
     * @param elements list of elements
     * @param stiffnessMatrix global stiffness matrix
     * @param localMatrixElementTwo local stiffness matrix to add
     * @param i element placement in X
     * @param j element placement in Y
     * @param numX number of elements in X
     */
    public static void appendElementTwoToStiffnessMatrix(List<Cst> elements,
            BigReal[][] stiffnessMatrix, BigReal[][] localMatrixElementTwo, int i, int j, Integer numX) {
        Cst element2 = new Cst(new Point(new BigDecimal(i*10), new BigDecimal(j*10)), 
                new Point(new BigDecimal(i*10+10), new BigDecimal(j*10)), 
                new Point(new BigDecimal(i*10+10), new BigDecimal(j*10+10)));
        elements.add(element2);
        //top left corner
        stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*j+i)*2] = stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*j+i)*2].add(localMatrixElementTwo[0][0]);
        stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*j+i)*2+1] = stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*j+i)*2+1].add(localMatrixElementTwo[0][1]);
        stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*j+i)*2] = stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*j+i)*2].add(localMatrixElementTwo[0][1]);
        stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*j+i)*2+1] = stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*j+i)*2+1].add(localMatrixElementTwo[1][1]);

        //top right corner
        stiffnessMatrix[((numX+1)*j+(i+1))*2][((numX+1)*j+(i+1))*2] = stiffnessMatrix[((numX+1)*j+(i+1))*2][((numX+1)*j+(i+1))*2].add(localMatrixElementTwo[2][2]);
        stiffnessMatrix[((numX+1)*j+(i+1))*2][((numX+1)*j+(i+1))*2+1] = stiffnessMatrix[((numX+1)*j+(i+1))*2][((numX+1)*j+(i+1))*2+1].add(localMatrixElementTwo[2][3]);
        stiffnessMatrix[((numX+1)*j+(i+1))*2+1][((numX+1)*j+(i+1))*2] = stiffnessMatrix[((numX+1)*j+(i+1))*2+1][((numX+1)*j+(i+1))*2].add(localMatrixElementTwo[2][3]);
        stiffnessMatrix[((numX+1)*j+(i+1))*2+1][((numX+1)*j+(i+1))*2+1] = stiffnessMatrix[((numX+1)*j+(i+1))*2+1][((numX+1)*j+(i+1))*2+1].add(localMatrixElementTwo[3][3]);

        //bottom right corner
        stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2][((numX+1)*(j+1)+(i+1))*2] = stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2][((numX+1)*(j+1)+(i+1))*2].add(localMatrixElementTwo[4][4]);
        stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2][((numX+1)*(j+1)+(i+1))*2+1] = stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2][((numX+1)*(j+1)+(i+1))*2+1].add(localMatrixElementTwo[4][5]);
        stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2+1][((numX+1)*(j+1)+(i+1))*2] = stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2+1][((numX+1)*(j+1)+(i+1))*2].add(localMatrixElementTwo[4][5]);
        stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2+1][((numX+1)*(j+1)+(i+1))*2+1] = stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2+1][((numX+1)*(j+1)+(i+1))*2+1].add(localMatrixElementTwo[5][5]);
        
        //ij
        stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*j+(i+1))*2] = stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*j+(i+1))*2].add(localMatrixElementTwo[0][2]);
        stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*j+(i+1))*2+1] = stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*j+(i+1))*2+1].add(localMatrixElementTwo[0][3]);
        stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*j+(i+1))*2] = stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*j+(i+1))*2].add(localMatrixElementTwo[1][2]);
        stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*j+(i+1))*2+1] = stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*j+(i+1))*2+1].add(localMatrixElementTwo[1][3]);
        

        stiffnessMatrix[((numX+1)*j+(i+1))*2][((numX+1)*j+i)*2] = stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*j+(i+1))*2];
        stiffnessMatrix[((numX+1)*j+(i+1))*2+1][((numX+1)*j+i)*2] = stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*j+(i+1))*2+1];
        stiffnessMatrix[((numX+1)*j+(i+1))*2][((numX+1)*j+i)*2+1] = stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*j+(i+1))*2];
        stiffnessMatrix[((numX+1)*j+(i+1))*2+1][((numX+1)*j+i)*2+1] = stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*j+(i+1))*2+1];
        
        //im
        stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*(j+1)+(i+1))*2] = stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*(j+1)+(i+1))*2].add(localMatrixElementTwo[0][4]);
        stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*(j+1)+(i+1))*2+1] = stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*(j+1)+(i+1))*2+1].add(localMatrixElementTwo[0][5]);
        stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*(j+1)+(i+1))*2] = stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*(j+1)+(i+1))*2].add(localMatrixElementTwo[1][4]);
        stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*(j+1)+(i+1))*2+1] = stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*(j+1)+(i+1))*2+1].add(localMatrixElementTwo[1][5]);

        stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2][((numX+1)*j+i)*2] = stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*(j+1)+(i+1))*2];
        stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2+1][((numX+1)*j+i)*2] = stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*(j+1)+(i+1))*2+1];
        stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2][((numX+1)*j+i)*2+1] = stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*(j+1)+(i+1))*2];
        stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2+1][((numX+1)*j+i)*2+1] = stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*(j+1)+(i+1))*2+1];
        
        //jm
        stiffnessMatrix[((numX+1)*j+(i+1))*2][((numX+1)*(j+1)+(i+1))*2] = stiffnessMatrix[((numX+1)*j+(i+1))*2][((numX+1)*(j+1)+(i+1))*2].add(localMatrixElementTwo[2][4]);
        stiffnessMatrix[((numX+1)*j+(i+1))*2][((numX+1)*(j+1)+(i+1))*2+1] = stiffnessMatrix[((numX+1)*j+(i+1))*2][((numX+1)*(j+1)+(i+1))*2+1].add(localMatrixElementTwo[2][5]);
        stiffnessMatrix[((numX+1)*j+(i+1))*2+1][((numX+1)*(j+1)+(i+1))*2] = stiffnessMatrix[((numX+1)*j+(i+1))*2+1][((numX+1)*(j+1)+(i+1))*2].add(localMatrixElementTwo[3][4]);
        stiffnessMatrix[((numX+1)*j+(i+1))*2+1][((numX+1)*(j+1)+(i+1))*2+1] = stiffnessMatrix[((numX+1)*j+(i+1))*2+1][((numX+1)*(j+1)+(i+1))*2+1].add(localMatrixElementTwo[3][5]);

        stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2][((numX+1)*j+(i+1))*2] = stiffnessMatrix[((numX+1)*j+(i+1))*2][((numX+1)*(j+1)+(i+1))*2];
        stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2+1][((numX+1)*j+(i+1))*2] = stiffnessMatrix[((numX+1)*j+(i+1))*2][((numX+1)*(j+1)+(i+1))*2+1];
        stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2][((numX+1)*j+(i+1))*2+1] = stiffnessMatrix[((numX+1)*j+(i+1))*2+1][((numX+1)*(j+1)+(i+1))*2];
        stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2+1][((numX+1)*j+(i+1))*2+1] = stiffnessMatrix[((numX+1)*j+(i+1))*2+1][((numX+1)*(j+1)+(i+1))*2+1];
    }

    /**
     * Adds a local mass matrix element two to the global mass matrix
     * 
     * @param massMatrix global mass matrix
     * @param i element placement in x
     * @param j element placement in y
     * @param numX total number of elements in X
     */
    public static void appendElementTwoToMassMatrix(BigReal[][] massMatrix, int i, int j, Integer numX) {
        //element 2
        //
        // i___j
        //  \  |
        //   \ |
        //    \|
        //     m

        //top left corner
        massMatrix[((numX+1)*j+i)*2][((numX+1)*j+i)*2] = massMatrix[((numX+1)*j+i)*2][((numX+1)*j+i)*2].add(TWO);
        massMatrix[((numX+1)*j+i)*2][((numX+1)*j+i)*2+1] = massMatrix[((numX+1)*j+i)*2][((numX+1)*j+i)*2+1].add(ZERO);
        massMatrix[((numX+1)*j+i)*2+1][((numX+1)*j+i)*2] = massMatrix[((numX+1)*j+i)*2+1][((numX+1)*j+i)*2].add(ZERO);
        massMatrix[((numX+1)*j+i)*2+1][((numX+1)*j+i)*2+1] = massMatrix[((numX+1)*j+i)*2+1][((numX+1)*j+i)*2+1].add(TWO);

        //top right corner
        massMatrix[((numX+1)*j+(i+1))*2][((numX+1)*j+(i+1))*2] = massMatrix[((numX+1)*j+(i+1))*2][((numX+1)*j+(i+1))*2].add(TWO);
        massMatrix[((numX+1)*j+(i+1))*2][((numX+1)*j+(i+1))*2+1] = massMatrix[((numX+1)*j+(i+1))*2][((numX+1)*j+(i+1))*2+1].add(ZERO);
        massMatrix[((numX+1)*j+(i+1))*2+1][((numX+1)*j+(i+1))*2] = massMatrix[((numX+1)*j+(i+1))*2+1][((numX+1)*j+(i+1))*2].add(ZERO);
        massMatrix[((numX+1)*j+(i+1))*2+1][((numX+1)*j+(i+1))*2+1] = massMatrix[((numX+1)*j+(i+1))*2+1][((numX+1)*j+(i+1))*2+1].add(TWO);

        //bottom right corner
        massMatrix[((numX+1)*(j+1)+(i+1))*2][((numX+1)*(j+1)+(i+1))*2] = massMatrix[((numX+1)*(j+1)+(i+1))*2][((numX+1)*(j+1)+(i+1))*2].add(TWO);
        massMatrix[((numX+1)*(j+1)+(i+1))*2][((numX+1)*(j+1)+(i+1))*2+1] = massMatrix[((numX+1)*(j+1)+(i+1))*2][((numX+1)*(j+1)+(i+1))*2+1].add(ZERO);
        massMatrix[((numX+1)*(j+1)+(i+1))*2+1][((numX+1)*(j+1)+(i+1))*2] = massMatrix[((numX+1)*(j+1)+(i+1))*2+1][((numX+1)*(j+1)+(i+1))*2].add(ZERO);
        massMatrix[((numX+1)*(j+1)+(i+1))*2+1][((numX+1)*(j+1)+(i+1))*2+1] = massMatrix[((numX+1)*(j+1)+(i+1))*2+1][((numX+1)*(j+1)+(i+1))*2+1].add(TWO);
        
        //ij
        massMatrix[((numX+1)*j+i)*2][((numX+1)*j+(i+1))*2] = massMatrix[((numX+1)*j+i)*2][((numX+1)*j+(i+1))*2].add(ONE);
        massMatrix[((numX+1)*j+i)*2][((numX+1)*j+(i+1))*2+1] = massMatrix[((numX+1)*j+i)*2][((numX+1)*j+(i+1))*2+1].add(ZERO);
        massMatrix[((numX+1)*j+i)*2+1][((numX+1)*j+(i+1))*2] = massMatrix[((numX+1)*j+i)*2+1][((numX+1)*j+(i+1))*2].add(ZERO);
        massMatrix[((numX+1)*j+i)*2+1][((numX+1)*j+(i+1))*2+1] = massMatrix[((numX+1)*j+i)*2+1][((numX+1)*j+(i+1))*2+1].add(ONE);
        

        massMatrix[((numX+1)*j+(i+1))*2][((numX+1)*j+i)*2] = massMatrix[((numX+1)*j+i)*2][((numX+1)*j+(i+1))*2];
        massMatrix[((numX+1)*j+(i+1))*2+1][((numX+1)*j+i)*2] = massMatrix[((numX+1)*j+i)*2][((numX+1)*j+(i+1))*2+1];
        massMatrix[((numX+1)*j+(i+1))*2][((numX+1)*j+i)*2+1] = massMatrix[((numX+1)*j+i)*2+1][((numX+1)*j+(i+1))*2];
        massMatrix[((numX+1)*j+(i+1))*2+1][((numX+1)*j+i)*2+1] = massMatrix[((numX+1)*j+i)*2+1][((numX+1)*j+(i+1))*2+1];
        
        //im
        massMatrix[((numX+1)*j+i)*2][((numX+1)*(j+1)+(i+1))*2] = massMatrix[((numX+1)*j+i)*2][((numX+1)*(j+1)+(i+1))*2].add(ONE);
        massMatrix[((numX+1)*j+i)*2][((numX+1)*(j+1)+(i+1))*2+1] = massMatrix[((numX+1)*j+i)*2][((numX+1)*(j+1)+(i+1))*2+1].add(ZERO);
        massMatrix[((numX+1)*j+i)*2+1][((numX+1)*(j+1)+(i+1))*2] = massMatrix[((numX+1)*j+i)*2+1][((numX+1)*(j+1)+(i+1))*2].add(ZERO);
        massMatrix[((numX+1)*j+i)*2+1][((numX+1)*(j+1)+(i+1))*2+1] = massMatrix[((numX+1)*j+i)*2+1][((numX+1)*(j+1)+(i+1))*2+1].add(ONE);

        massMatrix[((numX+1)*(j+1)+(i+1))*2][((numX+1)*j+i)*2] = massMatrix[((numX+1)*j+i)*2][((numX+1)*(j+1)+(i+1))*2];
        massMatrix[((numX+1)*(j+1)+(i+1))*2+1][((numX+1)*j+i)*2] = massMatrix[((numX+1)*j+i)*2][((numX+1)*(j+1)+(i+1))*2+1];
        massMatrix[((numX+1)*(j+1)+(i+1))*2][((numX+1)*j+i)*2+1] = massMatrix[((numX+1)*j+i)*2+1][((numX+1)*(j+1)+(i+1))*2];
        massMatrix[((numX+1)*(j+1)+(i+1))*2+1][((numX+1)*j+i)*2+1] = massMatrix[((numX+1)*j+i)*2+1][((numX+1)*(j+1)+(i+1))*2+1];
        
        //jm
        massMatrix[((numX+1)*j+(i+1))*2][((numX+1)*(j+1)+(i+1))*2] = massMatrix[((numX+1)*j+(i+1))*2][((numX+1)*(j+1)+(i+1))*2].add(ONE);
        massMatrix[((numX+1)*j+(i+1))*2][((numX+1)*(j+1)+(i+1))*2+1] = massMatrix[((numX+1)*j+(i+1))*2][((numX+1)*(j+1)+(i+1))*2+1].add(ZERO);
        massMatrix[((numX+1)*j+(i+1))*2+1][((numX+1)*(j+1)+(i+1))*2] = massMatrix[((numX+1)*j+(i+1))*2+1][((numX+1)*(j+1)+(i+1))*2].add(ZERO);
        massMatrix[((numX+1)*j+(i+1))*2+1][((numX+1)*(j+1)+(i+1))*2+1] = massMatrix[((numX+1)*j+(i+1))*2+1][((numX+1)*(j+1)+(i+1))*2+1].add(ONE);

        massMatrix[((numX+1)*(j+1)+(i+1))*2][((numX+1)*j+(i+1))*2] = massMatrix[((numX+1)*j+(i+1))*2][((numX+1)*(j+1)+(i+1))*2];
        massMatrix[((numX+1)*(j+1)+(i+1))*2+1][((numX+1)*j+(i+1))*2] = massMatrix[((numX+1)*j+(i+1))*2][((numX+1)*(j+1)+(i+1))*2+1];
        massMatrix[((numX+1)*(j+1)+(i+1))*2][((numX+1)*j+(i+1))*2+1] = massMatrix[((numX+1)*j+(i+1))*2+1][((numX+1)*(j+1)+(i+1))*2];
        massMatrix[((numX+1)*(j+1)+(i+1))*2+1][((numX+1)*j+(i+1))*2+1] = massMatrix[((numX+1)*j+(i+1))*2+1][((numX+1)*(j+1)+(i+1))*2+1];
    }

}
