package com.pescamillam.fem;

import java.math.BigDecimal;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import com.pescamillam.fem.util.UtilWindow;
import org.apache.commons.math3.linear.FieldLUDecomposition;
import org.apache.commons.math3.linear.FieldMatrix;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.util.BigReal;

import com.pescamillam.fem.element.Cst;
import com.pescamillam.fem.element.Point;
import com.pescamillam.fem.util.UtilWriter;

import static com.pescamillam.fem.util.Constants.FOUR;
import static com.pescamillam.fem.util.Constants.MINUS_ONE;
import static com.pescamillam.fem.util.Constants.NUM_TIMES;
import static com.pescamillam.fem.util.Constants.NUM_X;
import static com.pescamillam.fem.util.Constants.NUM_Y;
import static com.pescamillam.fem.util.Constants.TWELVE;
import static com.pescamillam.fem.util.Constants.TWO;
import static org.apache.commons.math3.util.BigReal.ONE;
import static org.apache.commons.math3.util.BigReal.ZERO;

public class Application {
    
    private static DecimalFormat df;
    
    @SuppressWarnings("unchecked")
    public static void main(String[] args) {
        Long startTime = System.nanoTime();
        setDecimalFormat();
        
        //parameters
        //thickness t
        BigReal thickness = new BigReal("1");
        //elasticity module E
        BigReal elasticity = new BigReal("3000000");
        //density p
        BigReal density = new BigReal("0.00073");
        //Area A
        BigReal area = new BigReal("450");
        //poisson ratio v
        BigReal poisson = new BigReal("0.3");
        
        //delta time
        BigReal deltaTime = new BigReal("0.00005");
        BigReal deltaTimeSquare = deltaTime.multiply(deltaTime); 

        writeToFile("=== Constants ===");
        writeToFile("area: " + area.bigDecimalValue().toPlainString());
        writeToFile("thickness: " + thickness.bigDecimalValue().toPlainString());
        writeToFile("elasticity: " + elasticity.bigDecimalValue().toPlainString());
        writeToFile("density: " + density.bigDecimalValue().toPlainString());
        writeToFile("poisson: " + poisson.bigDecimalValue().toPlainString());
        writeToFile("\n\n");
        
        //(1-v)
        BigReal poisson1mv = poisson.negate().add(new BigReal("1"));
        
        //(1-2v)/2
        BigReal poisson1m2vo2 = poisson.multiply(new BigReal("2")).negate().add(new BigReal("1")).divide(new BigReal("2"));
        
        //(1-2v)
        BigReal poisson1m2v = poisson.multiply(new BigReal("2")).negate().add(new BigReal("1"));

        //(1+v)
        BigReal poisson1pv = poisson.add(new BigReal("1"));
        
        //variables
        //acceleration
        FieldMatrix<BigReal>[] acceleration = new FieldMatrix[NUM_TIMES];

        List<Cst> elements = new ArrayList<>();
        
        //Stiffness matrix
        BigReal[][] stiffnessMatrix = new BigReal[2* NUM_X* NUM_Y+2* NUM_X+2* NUM_Y+2][2* NUM_X* NUM_Y+2* NUM_X+2* NUM_Y+2];
        
        //mass matrix
        BigReal[][] massMatrix = new BigReal[2* NUM_X* NUM_Y+2* NUM_X+2* NUM_Y+2][2* NUM_X* NUM_Y+2* NUM_X+2* NUM_Y+2];

        for (BigReal[] row : stiffnessMatrix) {
            Arrays.fill(row, ZERO);
        }

        for (BigReal[] row : massMatrix) {
            Arrays.fill(row, ZERO);
        }
        

        
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
        BigReal m11el1 = betaIel1.multiply(betaIel1).multiply(poisson1mv).add(gammaIel1.multiply(gammaIel1).multiply(poisson1m2vo2));
        //1-2
        BigReal m12el1 = betaIel1.multiply(betaIel1).multiply(poisson).add(betaIel1.multiply(gammaIel1).multiply(poisson1m2vo2));
        //1-3
        BigReal m13el1 = betaIel1.multiply(betaJel1).multiply(poisson1mv).add(gammaIel1.multiply(gammaJel1).multiply(poisson1m2vo2));
        //1-4
        BigReal m14el1 = betaIel1.multiply(gammaJel1).multiply(poisson).add(betaJel1.multiply(gammaIel1).multiply(poisson1m2vo2));
        //1-5
        BigReal m15el1 = betaIel1.multiply(betaMel1).multiply(poisson1mv).add(gammaIel1.multiply(gammaMel1).multiply(poisson1m2vo2));
        //1-6
        BigReal m16el1 = betaIel1.multiply(gammaMel1).multiply(poisson).add(betaMel1.multiply(gammaIel1).multiply(poisson1m2vo2));
        
        //2-1
        BigReal m21el1 = m12el1;
        //2-2
        BigReal m22el1 = gammaIel1.multiply(gammaIel1).multiply(poisson1mv).add(betaIel1.multiply(betaIel1).multiply(poisson1m2vo2));
        //2-3
        BigReal m23el1 = betaJel1.multiply(gammaIel1).multiply(poisson).add(betaIel1.multiply(gammaJel1).multiply(poisson1m2vo2));
        //2-4
        BigReal m24el1 = gammaIel1.multiply(gammaJel1).multiply(poisson1mv).add(betaIel1.multiply(betaJel1).multiply(poisson1m2vo2));
        //2-5
        BigReal m25el1 = betaMel1.multiply(gammaIel1).multiply(poisson).add(betaIel1.multiply(gammaMel1).multiply(poisson1m2vo2));
        //2-6
        BigReal m26el1 = gammaIel1.multiply(gammaMel1).multiply(poisson1mv).add(betaIel1.multiply(betaMel1).multiply(poisson1m2vo2));
        
        //3-1
//        BigReal m31el1 = m13el1;
        //3-2
//        BigReal m32 = m23el1;
        //3-3
        BigReal m33el1 = betaJel1.multiply(betaJel1).multiply(poisson1mv).add(gammaJel1.multiply(gammaJel1).multiply(poisson1m2vo2));
        //3-4
        BigReal m34el1 = betaJel1.multiply(gammaJel1).multiply(poisson).add(betaJel1.multiply(gammaJel1).multiply(poisson1m2vo2));
        //3-5
        BigReal m35el1 = betaJel1.multiply(betaMel1).multiply(poisson1mv).add(gammaJel1.multiply(gammaMel1).multiply(poisson1m2vo2));
        //3-6
        BigReal m36el1 = betaJel1.multiply(gammaMel1).multiply(poisson).add(betaMel1.multiply(gammaJel1).multiply(poisson1m2vo2));
        
        //4-1
//        BigReal m41el1 = m14el1;
        //4-2
//        BigReal m42el1 = m24el1;
        //4-3
        BigReal m43el1 = m34el1;
        //4-4
        BigReal m44el1 = gammaJel1.multiply(gammaJel1).multiply(poisson1mv).add(betaJel1.multiply(betaJel1).multiply(poisson1m2vo2));
        //4-5
        BigReal m45el1 = betaMel1.multiply(gammaJel1).multiply(poisson).add(betaJel1.multiply(gammaMel1).multiply(poisson1m2vo2));
        //4-6
        BigReal m46el1 = gammaIel1.multiply(gammaMel1).multiply(poisson1mv).add(betaJel1.multiply(betaMel1).multiply(poisson1m2vo2));
        
        //5-1
//        BigReal m51el1 = m15el1;
        //5-2
//        BigReal m52el1 = m25el1;
        //5-3
//        BigReal m53el1 = m35el1;
        //5-4
//        BigReal m54el1 = m45el1;
        //5-5
        BigReal m55el1 = betaMel1.multiply(betaMel1).multiply(poisson1mv).add(gammaMel1.multiply(gammaMel1).multiply(poisson1m2vo2));
        //5-6
        BigReal m56el1 = betaMel1.multiply(gammaMel1).multiply(poisson).add(betaMel1.multiply(gammaMel1).multiply(poisson1m2vo2));
        
        //6-1
//        BigReal m61el1 = m16el1;
        //6-2
//        BigReal m62el1 = m26el1;
        //6-3
//        BigReal m63el1 = m36el1;
        //6-4
//        BigReal m64el1 = m46el1;
        //6-5
        BigReal m65el1 = m56el1;
        //6-6
        BigReal m66el1 = gammaMel1.multiply(gammaMel1).multiply(poisson1mv).add(betaMel1.multiply(betaMel1).multiply(poisson1m2vo2));
        

        
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
        BigReal m11el2 = betaIel2.multiply(betaIel2).multiply(poisson1mv).add(gammaIel2.multiply(gammaIel2).multiply(poisson1m2vo2));
        //1-2
        BigReal m12el2 = betaIel2.multiply(betaIel2).multiply(poisson).add(betaIel2.multiply(gammaIel2).multiply(poisson1m2vo2));
        //1-3
        BigReal m13el2 = betaIel2.multiply(betaJel2).multiply(poisson1mv).add(gammaIel2.multiply(gammaJel2).multiply(poisson1m2vo2));
        //1-4
        BigReal m14el2 = betaIel2.multiply(gammaJel2).multiply(poisson).add(betaJel2.multiply(gammaIel2).multiply(poisson1m2vo2));
        //1-5
        BigReal m15el2 = betaIel2.multiply(betaMel2).multiply(poisson1mv).add(gammaIel2.multiply(gammaMel2).multiply(poisson1m2vo2));
        //1-6
        BigReal m16el2 = betaIel2.multiply(gammaMel2).multiply(poisson).add(betaMel2.multiply(gammaIel2).multiply(poisson1m2vo2));
        
        //2-1
        BigReal m21el2 = m12el2;
        //2-2
        BigReal m22el2 = gammaIel2.multiply(gammaIel2).multiply(poisson1mv).add(betaIel2.multiply(betaIel2).multiply(poisson1m2vo2));
        //2-3
        BigReal m23el2 = betaJel2.multiply(gammaIel2).multiply(poisson).add(betaIel2.multiply(gammaJel2).multiply(poisson1m2vo2));
        //2-4
        BigReal m24el2 = gammaIel2.multiply(gammaJel2).multiply(poisson1mv).add(betaIel2.multiply(betaJel2).multiply(poisson1m2vo2));
        //2-5
        BigReal m25el2 = betaMel2.multiply(gammaIel2).multiply(poisson).add(betaIel2.multiply(gammaMel2).multiply(poisson1m2vo2));
        //2-6
        BigReal m26el2 = gammaIel2.multiply(gammaMel2).multiply(poisson1mv).add(betaIel2.multiply(betaMel2).multiply(poisson1m2vo2));
        
        //3-1
//        BigReal m31el2 = m13el2;
        //3-2
//        BigReal m32el2 = m23el2;
        //3-3
        BigReal m33el2 = betaJel2.multiply(betaJel2).multiply(poisson1mv).add(gammaJel2.multiply(gammaJel2).multiply(poisson1m2vo2));
        //3-4
        BigReal m34el2 = betaJel2.multiply(gammaJel2).multiply(poisson).add(betaJel2.multiply(gammaJel2).multiply(poisson1m2vo2));
        //3-5
        BigReal m35el2 = betaJel2.multiply(betaMel2).multiply(poisson1mv).add(gammaJel2.multiply(gammaMel2).multiply(poisson1m2vo2));
        //3-6
        BigReal m36el2 = betaJel2.multiply(gammaMel2).multiply(poisson).add(betaMel2.multiply(gammaJel2).multiply(poisson1m2vo2));
        
        //4-1
//        BigReal m41el2 = m14el2;
        //4-2
//        BigReal m42el2 = m24el2;
        //4-3
        BigReal m43el2 = m34el2;
        //4-4
        BigReal m44el2 = gammaJel2.multiply(gammaJel2).multiply(poisson1mv).add(betaJel2.multiply(betaJel2).multiply(poisson1m2vo2));
        //4-5
        BigReal m45el2 = betaMel2.multiply(gammaJel2).multiply(poisson).add(betaJel2.multiply(gammaMel2).multiply(poisson1m2vo2));
        //4-6
        BigReal m46el2 = gammaIel2.multiply(gammaMel2).multiply(poisson1mv).add(betaJel2.multiply(betaMel2).multiply(poisson1m2vo2));
        
        //5-1
//        BigReal m51el2 = m15el2;
        //5-2
//        BigReal m52el2 = m25el2;
//        5-3
//        BigReal m53el2 = m35el2;
        //5-4
//        BigReal m54el2 = m45el2;
        //5-5
        BigReal m55el2 = betaMel2.multiply(betaMel2).multiply(poisson1mv).add(gammaMel2.multiply(gammaMel2).multiply(poisson1m2vo2));
        //5-6
        BigReal m56el2 = betaMel2.multiply(gammaMel2).multiply(poisson).add(betaMel2.multiply(gammaMel2).multiply(poisson1m2vo2));
        
        //6-1
//        BigReal m61el2 = m16el2;
        //6-2
//        BigReal m62el2 = m26el2;
        //6-3
//        BigReal m63el2 = m36el2;
        //6-4
//        BigReal m64el2 = m46el2;
        //6-5
        BigReal m65el2 = m56el2;
        //6-6
        BigReal m66el2 = gammaMel2.multiply(gammaMel2).multiply(poisson1mv).add(betaMel2.multiply(betaMel2).multiply(poisson1m2vo2));


        //points matrix
        for (int i = 0; i < NUM_X; i++) {
            for (int j = 0; j < NUM_Y; j++) {
                setElementOneStiffnessMatrix(elements, stiffnessMatrix, m11el1, m12el1, m13el1,
                        m14el1, m15el1, m16el1, m21el1, m22el1, m23el1, m24el1, m25el1, m26el1,
                        m33el1, m34el1, m35el1, m36el1, m43el1, m44el1, m45el1, m46el1, m55el1,
                        m56el1, m65el1, m66el1, i, j);
                
                setElementTwoStiffnessMatrix(elements, stiffnessMatrix, m11el2, m12el2, m13el2,
                        m14el2, m15el2, m16el2, m21el2, m22el2, m23el2, m24el2, m25el2, m26el2,
                        m33el2, m34el2, m35el2, m36el2, m43el2, m44el2, m45el2, m46el2, m55el2,
                        m56el2, m65el2, m66el2, i, j);

                setElementOneMassMatrix(massMatrix, i, j, TWO, ZERO, ONE);

                setElementTwoMassMatrix(massMatrix, i, j, TWO, ZERO, ONE);
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
        BigReal[] force1Vector = new BigReal[2* NUM_X* NUM_Y+2* NUM_X+2* NUM_Y+2];
        Arrays.fill(force1Vector, ZERO);
//        force1Vector[4] = ONE_HUNDRED;
//        force1Vector[5] = ONE_HUNDRED;
        BigReal[] force0Vector = new BigReal[2* NUM_X* NUM_Y+2* NUM_X+2* NUM_Y+2];
        Arrays.fill(force0Vector, ZERO);
        force0Vector[NUM_X*2+2 + ((NUM_X*2)/2+1+6)] = ONE_HUNDRED;
        BigReal[] vector0 = new BigReal[2* NUM_X* NUM_Y+2* NUM_X+2* NUM_Y+2];
        Arrays.fill(vector0, ZERO);


//        for (int i = 0; i < Constants.NUM_X + 1; i++) {
//            f
        FieldMatrix<BigReal>[] force = new FieldMatrix[NUM_TIMES];
        for (int i = 0; i < NUM_TIMES; i++) {
            if (i < NUM_TIMES/3) {
                force[i] = MatrixUtils.createColumnFieldMatrix(force0Vector)
                //                    .scalarMultiply(new BigReal(Math.sin(i/10.0)))
                ;
            } else {
                force[i] = MatrixUtils.createColumnFieldMatrix(vector0);
            }
        }

        //acceleration = mass inverse * (F0 - [K]{d0}) as {d0} = 0 
        //acceleration = mass inverse * (F0)
        writeToFile("=== acceleration 0 ===");
        acceleration[0] = inverseMassMatrix.multiply(force[0]);
        printMatrix(acceleration[0]);

        //displacement -1
        FieldMatrix<BigReal> displacementM1 = acceleration[0].scalarMultiply(deltaTimeSquare.divide(new BigReal("2")));

        writeToFile("=== displacement -1 ===");
        printMatrix(displacementM1);
        
        //displacement 1
        BigReal[] displacement0 = new BigReal[2* NUM_X* NUM_Y+2* NUM_X+2* NUM_Y+2];
        Arrays.fill(displacement0, ZERO);
//        displacement0[5] = BigReal.ONE;
        FieldMatrix<BigReal>[] displacement = new FieldMatrix[NUM_TIMES];
        FieldMatrix<BigReal>[] speed = new FieldMatrix[NUM_TIMES];
        displacement[0] = MatrixUtils.createColumnFieldMatrix(displacement0);

        writeToFile("=== displacement 0 ===");
        printMatrix(displacement[0]);
        
        displacement[1] = inverseMassMatrix.multiply(
                force[0].scalarMultiply(deltaTimeSquare)
                //.add(massFieldMatrix.scalarMultiply(TWO_BIG_REAL).add(stiffnessFieldMatrix.scalarMultiply(deltaTime.multiply(deltaTime))).multiply(displacementM1))
                .add(massFieldMatrix.multiply(displacementM1).scalarMultiply(MINUS_ONE))
                );

//        writeToFile("=== displacement 1 ===");
//        printMatrix(displacement[1]);

        for (int i = 2; i < NUM_TIMES; i++) {
            
            displacement[i] = inverseMassMatrix.multiply(force[i-1].scalarMultiply(deltaTimeSquare)
                                .add(massFieldMatrix.scalarMultiply(TWO).add(stiffnessFieldMatrix.scalarMultiply(deltaTimeSquare.multiply(MINUS_ONE))).multiply(displacement[i-1]))
                                .add(massFieldMatrix.multiply(displacement[i-2]).scalarMultiply(MINUS_ONE))
                                );
            
//            writeToFile("==== displacement " + i + " ====");
//            printMatrix(displacement[i]);
            
            speed[i-1] = displacement[i].add(displacement[i-2].scalarMultiply(MINUS_ONE)).scalarMultiply(ONE.divide(TWO.multiply(deltaTime)));

//            writeToFile("=== speed " + (i-1) + " ===");
//            printMatrix(speed[i-1]);
//
            acceleration[i] = inverseMassMatrix.multiply(force[i].add(stiffnessFieldMatrix.multiply(displacement[i]).scalarMultiply(MINUS_ONE)));
//            writeToFile("=== acceleration " + i + " ===");
//            printMatrix(acceleration[i]);
            
            writeToFile("=== force " + i + " ===");
            printMatrix(force[i]);
        }
        
        writeToFile("Total time: " + (System.nanoTime() - startTime)/1000000000.0 + "s");
        UtilWindow.printElements(elements, displacement, speed, acceleration, force);
        
    }

    private static void setElementTwoStiffnessMatrix(List<Cst> elements,
            BigReal[][] stiffnessMatrix, BigReal m11el2, BigReal m12el2, BigReal m13el2,
            BigReal m14el2, BigReal m15el2, BigReal m16el2, BigReal m21el2, BigReal m22el2,
            BigReal m23el2, BigReal m24el2, BigReal m25el2, BigReal m26el2, BigReal m33el2,
            BigReal m34el2, BigReal m35el2, BigReal m36el2, BigReal m43el2, BigReal m44el2,
            BigReal m45el2, BigReal m46el2, BigReal m55el2, BigReal m56el2, BigReal m65el2,
            BigReal m66el2, int i, int j) {
        Cst element2 = new Cst(new Point(new BigDecimal(i*10), new BigDecimal(j*10)), 
                new Point(new BigDecimal(i*10+10), new BigDecimal(j*10)), 
                new Point(new BigDecimal(i*10+10), new BigDecimal(j*10+10)));
        elements.add(element2);
        //top left corner
        stiffnessMatrix[((NUM_X+1)*j+i)*2][((NUM_X+1)*j+i)*2] = stiffnessMatrix[((NUM_X+1)*j+i)*2][((NUM_X+1)*j+i)*2].add(m11el2);
        stiffnessMatrix[((NUM_X+1)*j+i)*2][((NUM_X+1)*j+i)*2+1] = stiffnessMatrix[((NUM_X+1)*j+i)*2][((NUM_X+1)*j+i)*2+1].add(m12el2);
        stiffnessMatrix[((NUM_X+1)*j+i)*2+1][((NUM_X+1)*j+i)*2] = stiffnessMatrix[((NUM_X+1)*j+i)*2+1][((NUM_X+1)*j+i)*2].add(m21el2);
        stiffnessMatrix[((NUM_X+1)*j+i)*2+1][((NUM_X+1)*j+i)*2+1] = stiffnessMatrix[((NUM_X+1)*j+i)*2+1][((NUM_X+1)*j+i)*2+1].add(m22el2);

        //top right corner
        stiffnessMatrix[((NUM_X+1)*j+(i+1))*2][((NUM_X+1)*j+(i+1))*2] = stiffnessMatrix[((NUM_X+1)*j+(i+1))*2][((NUM_X+1)*j+(i+1))*2].add(m33el2);
        stiffnessMatrix[((NUM_X+1)*j+(i+1))*2][((NUM_X+1)*j+(i+1))*2+1] = stiffnessMatrix[((NUM_X+1)*j+(i+1))*2][((NUM_X+1)*j+(i+1))*2+1].add(m34el2);
        stiffnessMatrix[((NUM_X+1)*j+(i+1))*2+1][((NUM_X+1)*j+(i+1))*2] = stiffnessMatrix[((NUM_X+1)*j+(i+1))*2+1][((NUM_X+1)*j+(i+1))*2].add(m43el2);
        stiffnessMatrix[((NUM_X+1)*j+(i+1))*2+1][((NUM_X+1)*j+(i+1))*2+1] = stiffnessMatrix[((NUM_X+1)*j+(i+1))*2+1][((NUM_X+1)*j+(i+1))*2+1].add(m44el2);

        //bottom right corner
        stiffnessMatrix[((NUM_X+1)*(j+1)+(i+1))*2][((NUM_X+1)*(j+1)+(i+1))*2] = stiffnessMatrix[((NUM_X+1)*(j+1)+(i+1))*2][((NUM_X+1)*(j+1)+(i+1))*2].add(m55el2);
        stiffnessMatrix[((NUM_X+1)*(j+1)+(i+1))*2][((NUM_X+1)*(j+1)+(i+1))*2+1] = stiffnessMatrix[((NUM_X+1)*(j+1)+(i+1))*2][((NUM_X+1)*(j+1)+(i+1))*2+1].add(m56el2);
        stiffnessMatrix[((NUM_X+1)*(j+1)+(i+1))*2+1][((NUM_X+1)*(j+1)+(i+1))*2] = stiffnessMatrix[((NUM_X+1)*(j+1)+(i+1))*2+1][((NUM_X+1)*(j+1)+(i+1))*2].add(m65el2);
        stiffnessMatrix[((NUM_X+1)*(j+1)+(i+1))*2+1][((NUM_X+1)*(j+1)+(i+1))*2+1] = stiffnessMatrix[((NUM_X+1)*(j+1)+(i+1))*2+1][((NUM_X+1)*(j+1)+(i+1))*2+1].add(m66el2);
        
        //ij
        stiffnessMatrix[((NUM_X+1)*j+i)*2][((NUM_X+1)*j+(i+1))*2] = stiffnessMatrix[((NUM_X+1)*j+i)*2][((NUM_X+1)*j+(i+1))*2].add(m13el2);
        stiffnessMatrix[((NUM_X+1)*j+i)*2][((NUM_X+1)*j+(i+1))*2+1] = stiffnessMatrix[((NUM_X+1)*j+i)*2][((NUM_X+1)*j+(i+1))*2+1].add(m14el2);
        stiffnessMatrix[((NUM_X+1)*j+i)*2+1][((NUM_X+1)*j+(i+1))*2] = stiffnessMatrix[((NUM_X+1)*j+i)*2+1][((NUM_X+1)*j+(i+1))*2].add(m23el2);
        stiffnessMatrix[((NUM_X+1)*j+i)*2+1][((NUM_X+1)*j+(i+1))*2+1] = stiffnessMatrix[((NUM_X+1)*j+i)*2+1][((NUM_X+1)*j+(i+1))*2+1].add(m24el2);
        

        stiffnessMatrix[((NUM_X+1)*j+(i+1))*2][((NUM_X+1)*j+i)*2] = stiffnessMatrix[((NUM_X+1)*j+i)*2][((NUM_X+1)*j+(i+1))*2];
        stiffnessMatrix[((NUM_X+1)*j+(i+1))*2+1][((NUM_X+1)*j+i)*2] = stiffnessMatrix[((NUM_X+1)*j+i)*2][((NUM_X+1)*j+(i+1))*2+1];
        stiffnessMatrix[((NUM_X+1)*j+(i+1))*2][((NUM_X+1)*j+i)*2+1] = stiffnessMatrix[((NUM_X+1)*j+i)*2+1][((NUM_X+1)*j+(i+1))*2];
        stiffnessMatrix[((NUM_X+1)*j+(i+1))*2+1][((NUM_X+1)*j+i)*2+1] = stiffnessMatrix[((NUM_X+1)*j+i)*2+1][((NUM_X+1)*j+(i+1))*2+1];
        
        //im
        stiffnessMatrix[((NUM_X+1)*j+i)*2][((NUM_X+1)*(j+1)+(i+1))*2] = stiffnessMatrix[((NUM_X+1)*j+i)*2][((NUM_X+1)*(j+1)+(i+1))*2].add(m15el2);
        stiffnessMatrix[((NUM_X+1)*j+i)*2][((NUM_X+1)*(j+1)+(i+1))*2+1] = stiffnessMatrix[((NUM_X+1)*j+i)*2][((NUM_X+1)*(j+1)+(i+1))*2+1].add(m16el2);
        stiffnessMatrix[((NUM_X+1)*j+i)*2+1][((NUM_X+1)*(j+1)+(i+1))*2] = stiffnessMatrix[((NUM_X+1)*j+i)*2+1][((NUM_X+1)*(j+1)+(i+1))*2].add(m25el2);
        stiffnessMatrix[((NUM_X+1)*j+i)*2+1][((NUM_X+1)*(j+1)+(i+1))*2+1] = stiffnessMatrix[((NUM_X+1)*j+i)*2+1][((NUM_X+1)*(j+1)+(i+1))*2+1].add(m26el2);

        stiffnessMatrix[((NUM_X+1)*(j+1)+(i+1))*2][((NUM_X+1)*j+i)*2] = stiffnessMatrix[((NUM_X+1)*j+i)*2][((NUM_X+1)*(j+1)+(i+1))*2];
        stiffnessMatrix[((NUM_X+1)*(j+1)+(i+1))*2+1][((NUM_X+1)*j+i)*2] = stiffnessMatrix[((NUM_X+1)*j+i)*2][((NUM_X+1)*(j+1)+(i+1))*2+1];
        stiffnessMatrix[((NUM_X+1)*(j+1)+(i+1))*2][((NUM_X+1)*j+i)*2+1] = stiffnessMatrix[((NUM_X+1)*j+i)*2+1][((NUM_X+1)*(j+1)+(i+1))*2];
        stiffnessMatrix[((NUM_X+1)*(j+1)+(i+1))*2+1][((NUM_X+1)*j+i)*2+1] = stiffnessMatrix[((NUM_X+1)*j+i)*2+1][((NUM_X+1)*(j+1)+(i+1))*2+1];
        
        //jm
        stiffnessMatrix[((NUM_X+1)*j+(i+1))*2][((NUM_X+1)*(j+1)+(i+1))*2] = stiffnessMatrix[((NUM_X+1)*j+(i+1))*2][((NUM_X+1)*(j+1)+(i+1))*2].add(m35el2);
        stiffnessMatrix[((NUM_X+1)*j+(i+1))*2][((NUM_X+1)*(j+1)+(i+1))*2+1] = stiffnessMatrix[((NUM_X+1)*j+(i+1))*2][((NUM_X+1)*(j+1)+(i+1))*2+1].add(m36el2);
        stiffnessMatrix[((NUM_X+1)*j+(i+1))*2+1][((NUM_X+1)*(j+1)+(i+1))*2] = stiffnessMatrix[((NUM_X+1)*j+(i+1))*2+1][((NUM_X+1)*(j+1)+(i+1))*2].add(m45el2);
        stiffnessMatrix[((NUM_X+1)*j+(i+1))*2+1][((NUM_X+1)*(j+1)+(i+1))*2+1] = stiffnessMatrix[((NUM_X+1)*j+(i+1))*2+1][((NUM_X+1)*(j+1)+(i+1))*2+1].add(m46el2);

        stiffnessMatrix[((NUM_X+1)*(j+1)+(i+1))*2][((NUM_X+1)*j+(i+1))*2] = stiffnessMatrix[((NUM_X+1)*j+(i+1))*2][((NUM_X+1)*(j+1)+(i+1))*2];
        stiffnessMatrix[((NUM_X+1)*(j+1)+(i+1))*2+1][((NUM_X+1)*j+(i+1))*2] = stiffnessMatrix[((NUM_X+1)*j+(i+1))*2][((NUM_X+1)*(j+1)+(i+1))*2+1];
        stiffnessMatrix[((NUM_X+1)*(j+1)+(i+1))*2][((NUM_X+1)*j+(i+1))*2+1] = stiffnessMatrix[((NUM_X+1)*j+(i+1))*2+1][((NUM_X+1)*(j+1)+(i+1))*2];
        stiffnessMatrix[((NUM_X+1)*(j+1)+(i+1))*2+1][((NUM_X+1)*j+(i+1))*2+1] = stiffnessMatrix[((NUM_X+1)*j+(i+1))*2+1][((NUM_X+1)*(j+1)+(i+1))*2+1];
    }

    private static void setElementOneStiffnessMatrix(List<Cst> elements,
            BigReal[][] stiffnessMatrix, BigReal m11el1, BigReal m12el1, BigReal m13el1,
            BigReal m14el1, BigReal m15el1, BigReal m16el1, BigReal m21el1, BigReal m22el1,
            BigReal m23el1, BigReal m24el1, BigReal m25el1, BigReal m26el1, BigReal m33el1,
            BigReal m34el1, BigReal m35el1, BigReal m36el1, BigReal m43el1, BigReal m44el1,
            BigReal m45el1, BigReal m46el1, BigReal m55el1, BigReal m56el1, BigReal m65el1,
            BigReal m66el1, int i, int j) {
        Cst element = new Cst(new Point(new BigDecimal(i*10), new BigDecimal(j*10)), 
                new Point(new BigDecimal(i*10+10), new BigDecimal(j*10+10)), 
                new Point(new BigDecimal(i*10), new BigDecimal(j*10+10)));
        elements.add(element);

        //top left corner
        stiffnessMatrix[((NUM_X+1)*j+i)*2][((NUM_X+1)*j+i)*2] = stiffnessMatrix[((NUM_X+1)*j+i)*2][((NUM_X+1)*j+i)*2].add(m11el1);
        stiffnessMatrix[((NUM_X+1)*j+i)*2][((NUM_X+1)*j+i)*2+1] = stiffnessMatrix[((NUM_X+1)*j+i)*2][((NUM_X+1)*j+i)*2+1].add(m12el1);
        stiffnessMatrix[((NUM_X+1)*j+i)*2+1][((NUM_X+1)*j+i)*2] = stiffnessMatrix[((NUM_X+1)*j+i)*2+1][((NUM_X+1)*j+i)*2].add(m21el1);
        stiffnessMatrix[((NUM_X+1)*j+i)*2+1][((NUM_X+1)*j+i)*2+1] = stiffnessMatrix[((NUM_X+1)*j+i)*2+1][((NUM_X+1)*j+i)*2+1].add(m22el1);

        //bottom left corner
        stiffnessMatrix[((NUM_X+1)*(j+1)+i)*2][((NUM_X+1)*(j+1)+i)*2] = stiffnessMatrix[((NUM_X+1)*(j+1)+i)*2][((NUM_X+1)*(j+1)+i)*2].add(m55el1);
        stiffnessMatrix[((NUM_X+1)*(j+1)+i)*2][((NUM_X+1)*(j+1)+i)*2+1] = stiffnessMatrix[((NUM_X+1)*(j+1)+i)*2][((NUM_X+1)*(j+1)+i)*2+1].add(m56el1);
        stiffnessMatrix[((NUM_X+1)*(j+1)+i)*2+1][((NUM_X+1)*(j+1)+i)*2] = stiffnessMatrix[((NUM_X+1)*(j+1)+i)*2+1][((NUM_X+1)*(j+1)+i)*2].add(m65el1);
        stiffnessMatrix[((NUM_X+1)*(j+1)+i)*2+1][((NUM_X+1)*(j+1)+i)*2+1] = stiffnessMatrix[((NUM_X+1)*(j+1)+i)*2+1][((NUM_X+1)*(j+1)+i)*2+1].add(m66el1);

        //bottom right corner
        stiffnessMatrix[((NUM_X+1)*(j+1)+(i+1))*2][((NUM_X+1)*(j+1)+(i+1))*2] = stiffnessMatrix[((NUM_X+1)*(j+1)+(i+1))*2][((NUM_X+1)*(j+1)+(i+1))*2].add(m33el1);
        stiffnessMatrix[((NUM_X+1)*(j+1)+(i+1))*2][((NUM_X+1)*(j+1)+(i+1))*2+1] = stiffnessMatrix[((NUM_X+1)*(j+1)+(i+1))*2][((NUM_X+1)*(j+1)+(i+1))*2+1].add(m34el1);
        stiffnessMatrix[((NUM_X+1)*(j+1)+(i+1))*2+1][((NUM_X+1)*(j+1)+(i+1))*2] = stiffnessMatrix[((NUM_X+1)*(j+1)+(i+1))*2+1][((NUM_X+1)*(j+1)+(i+1))*2].add(m43el1);
        stiffnessMatrix[((NUM_X+1)*(j+1)+(i+1))*2+1][((NUM_X+1)*(j+1)+(i+1))*2+1] = stiffnessMatrix[((NUM_X+1)*(j+1)+(i+1))*2+1][((NUM_X+1)*(j+1)+(i+1))*2+1].add(m44el1);

        //ij
        stiffnessMatrix[((NUM_X+1)*j+i)*2][((NUM_X+1)*(j+1)+(i+1))*2] = stiffnessMatrix[((NUM_X+1)*j+i)*2][((NUM_X+1)*(j+1)+(i+1))*2].add(m13el1);
        stiffnessMatrix[((NUM_X+1)*j+i)*2][((NUM_X+1)*(j+1)+(i+1))*2+1] = stiffnessMatrix[((NUM_X+1)*j+i)*2][((NUM_X+1)*(j+1)+(i+1))*2+1].add(m14el1);
        stiffnessMatrix[((NUM_X+1)*j+i)*2+1][((NUM_X+1)*(j+1)+(i+1))*2] = stiffnessMatrix[((NUM_X+1)*j+i)*2+1][((NUM_X+1)*(j+1)+(i+1))*2].add(m23el1);
        stiffnessMatrix[((NUM_X+1)*j+i)*2+1][((NUM_X+1)*(j+1)+(i+1))*2+1] = stiffnessMatrix[((NUM_X+1)*j+i)*2+1][((NUM_X+1)*(j+1)+(i+1))*2+1].add(m24el1);

        stiffnessMatrix[((NUM_X+1)*(j+1)+(i+1))*2][((NUM_X+1)*j+i)*2] = stiffnessMatrix[((NUM_X+1)*j+i)*2][((NUM_X+1)*(j+1)+(i+1))*2];
        stiffnessMatrix[((NUM_X+1)*(j+1)+(i+1))*2+1][((NUM_X+1)*j+i)*2] = stiffnessMatrix[((NUM_X+1)*j+i)*2][((NUM_X+1)*(j+1)+(i+1))*2+1];
        stiffnessMatrix[((NUM_X+1)*(j+1)+(i+1))*2][((NUM_X+1)*j+i)*2+1] = stiffnessMatrix[((NUM_X+1)*j+i)*2+1][((NUM_X+1)*(j+1)+(i+1))*2];
        stiffnessMatrix[((NUM_X+1)*(j+1)+(i+1))*2+1][((NUM_X+1)*j+i)*2+1] = stiffnessMatrix[((NUM_X+1)*j+i)*2+1][((NUM_X+1)*(j+1)+(i+1))*2+1];

        //im
        stiffnessMatrix[((NUM_X+1)*j+i)*2][((NUM_X+1)*(j+1)+i)*2] = stiffnessMatrix[((NUM_X+1)*j+i)*2][((NUM_X+1)*(j+1)+i)*2].add(m15el1);
        stiffnessMatrix[((NUM_X+1)*j+i)*2][((NUM_X+1)*(j+1)+i)*2+1] = stiffnessMatrix[((NUM_X+1)*j+i)*2][((NUM_X+1)*(j+1)+i)*2+1].add(m16el1);
        stiffnessMatrix[((NUM_X+1)*j+i)*2+1][((NUM_X+1)*(j+1)+i)*2] = stiffnessMatrix[((NUM_X+1)*j+i)*2+1][((NUM_X+1)*(j+1)+i)*2].add(m25el1);
        stiffnessMatrix[((NUM_X+1)*j+i)*2+1][((NUM_X+1)*(j+1)+i)*2+1] = stiffnessMatrix[((NUM_X+1)*j+i)*2+1][((NUM_X+1)*(j+1)+i)*2+1].add(m26el1);

        stiffnessMatrix[((NUM_X+1)*(j+1)+i)*2][((NUM_X+1)*j+i)*2] = stiffnessMatrix[((NUM_X+1)*j+i)*2][((NUM_X+1)*(j+1)+i)*2];
        stiffnessMatrix[((NUM_X+1)*(j+1)+i)*2+1][((NUM_X+1)*j+i)*2] = stiffnessMatrix[((NUM_X+1)*j+i)*2][((NUM_X+1)*(j+1)+i)*2+1];
        stiffnessMatrix[((NUM_X+1)*(j+1)+i)*2][((NUM_X+1)*j+i)*2+1] = stiffnessMatrix[((NUM_X+1)*j+i)*2+1][((NUM_X+1)*(j+1)+i)*2];
        stiffnessMatrix[((NUM_X+1)*(j+1)+i)*2+1][((NUM_X+1)*j+i)*2+1] = stiffnessMatrix[((NUM_X+1)*j+i)*2+1][((NUM_X+1)*(j+1)+i)*2+1];
        
        //jm
        stiffnessMatrix[((NUM_X+1)*(j+1)+i)*2][((NUM_X+1)*(j+1)+(i+1))*2] = stiffnessMatrix[((NUM_X+1)*(j+1)+i)*2][((NUM_X+1)*(j+1)+(i+1))*2].add(m35el1);
        stiffnessMatrix[((NUM_X+1)*(j+1)+i)*2+1][((NUM_X+1)*(j+1)+(i+1))*2] = stiffnessMatrix[((NUM_X+1)*(j+1)+i)*2+1][((NUM_X+1)*(j+1)+(i+1))*2].add(m45el1);
        stiffnessMatrix[((NUM_X+1)*(j+1)+i)*2][((NUM_X+1)*(j+1)+(i+1))*2+1] = stiffnessMatrix[((NUM_X+1)*(j+1)+i)*2][((NUM_X+1)*(j+1)+(i+1))*2+1].add(m36el1);
        stiffnessMatrix[((NUM_X+1)*(j+1)+i)*2+1][((NUM_X+1)*(j+1)+(i+1))*2+1] = stiffnessMatrix[((NUM_X+1)*(j+1)+i)*2+1][((NUM_X+1)*(j+1)+(i+1))*2+1].add(m46el1);

        stiffnessMatrix[((NUM_X+1)*(j+1)+(i+1))*2][((NUM_X+1)*(j+1)+i)*2] = stiffnessMatrix[((NUM_X+1)*(j+1)+i)*2][((NUM_X+1)*(j+1)+(i+1))*2];
        stiffnessMatrix[((NUM_X+1)*(j+1)+(i+1))*2][((NUM_X+1)*(j+1)+i)*2+1] = stiffnessMatrix[((NUM_X+1)*(j+1)+i)*2+1][((NUM_X+1)*(j+1)+(i+1))*2];
        stiffnessMatrix[((NUM_X+1)*(j+1)+(i+1))*2+1][((NUM_X+1)*(j+1)+i)*2] = stiffnessMatrix[((NUM_X+1)*(j+1)+i)*2][((NUM_X+1)*(j+1)+(i+1))*2+1];
        stiffnessMatrix[((NUM_X+1)*(j+1)+(i+1))*2+1][((NUM_X+1)*(j+1)+i)*2+1] = stiffnessMatrix[((NUM_X+1)*(j+1)+i)*2+1][((NUM_X+1)*(j+1)+(i+1))*2+1];
    }

    private static void setElementTwoMassMatrix(BigReal[][] massMatrix, int i, int j, BigReal two,
            BigReal zero, BigReal one) {
        //element 2
        //
        // i___j
        //  \  |
        //   \ |
        //    \|
        //     m

        //top left corner
        massMatrix[((NUM_X+1)*j+i)*2][((NUM_X+1)*j+i)*2] = massMatrix[((NUM_X+1)*j+i)*2][((NUM_X+1)*j+i)*2].add(two);
        massMatrix[((NUM_X+1)*j+i)*2][((NUM_X+1)*j+i)*2+1] = massMatrix[((NUM_X+1)*j+i)*2][((NUM_X+1)*j+i)*2+1].add(zero);
        massMatrix[((NUM_X+1)*j+i)*2+1][((NUM_X+1)*j+i)*2] = massMatrix[((NUM_X+1)*j+i)*2+1][((NUM_X+1)*j+i)*2].add(zero);
        massMatrix[((NUM_X+1)*j+i)*2+1][((NUM_X+1)*j+i)*2+1] = massMatrix[((NUM_X+1)*j+i)*2+1][((NUM_X+1)*j+i)*2+1].add(two);

        //top right corner
        massMatrix[((NUM_X+1)*j+(i+1))*2][((NUM_X+1)*j+(i+1))*2] = massMatrix[((NUM_X+1)*j+(i+1))*2][((NUM_X+1)*j+(i+1))*2].add(two);
        massMatrix[((NUM_X+1)*j+(i+1))*2][((NUM_X+1)*j+(i+1))*2+1] = massMatrix[((NUM_X+1)*j+(i+1))*2][((NUM_X+1)*j+(i+1))*2+1].add(zero);
        massMatrix[((NUM_X+1)*j+(i+1))*2+1][((NUM_X+1)*j+(i+1))*2] = massMatrix[((NUM_X+1)*j+(i+1))*2+1][((NUM_X+1)*j+(i+1))*2].add(zero);
        massMatrix[((NUM_X+1)*j+(i+1))*2+1][((NUM_X+1)*j+(i+1))*2+1] = massMatrix[((NUM_X+1)*j+(i+1))*2+1][((NUM_X+1)*j+(i+1))*2+1].add(two);

        //bottom right corner
        massMatrix[((NUM_X+1)*(j+1)+(i+1))*2][((NUM_X+1)*(j+1)+(i+1))*2] = massMatrix[((NUM_X+1)*(j+1)+(i+1))*2][((NUM_X+1)*(j+1)+(i+1))*2].add(two);
        massMatrix[((NUM_X+1)*(j+1)+(i+1))*2][((NUM_X+1)*(j+1)+(i+1))*2+1] = massMatrix[((NUM_X+1)*(j+1)+(i+1))*2][((NUM_X+1)*(j+1)+(i+1))*2+1].add(zero);
        massMatrix[((NUM_X+1)*(j+1)+(i+1))*2+1][((NUM_X+1)*(j+1)+(i+1))*2] = massMatrix[((NUM_X+1)*(j+1)+(i+1))*2+1][((NUM_X+1)*(j+1)+(i+1))*2].add(zero);
        massMatrix[((NUM_X+1)*(j+1)+(i+1))*2+1][((NUM_X+1)*(j+1)+(i+1))*2+1] = massMatrix[((NUM_X+1)*(j+1)+(i+1))*2+1][((NUM_X+1)*(j+1)+(i+1))*2+1].add(two);
        
        //ij
        massMatrix[((NUM_X+1)*j+i)*2][((NUM_X+1)*j+(i+1))*2] = massMatrix[((NUM_X+1)*j+i)*2][((NUM_X+1)*j+(i+1))*2].add(one);
        massMatrix[((NUM_X+1)*j+i)*2][((NUM_X+1)*j+(i+1))*2+1] = massMatrix[((NUM_X+1)*j+i)*2][((NUM_X+1)*j+(i+1))*2+1].add(zero);
        massMatrix[((NUM_X+1)*j+i)*2+1][((NUM_X+1)*j+(i+1))*2] = massMatrix[((NUM_X+1)*j+i)*2+1][((NUM_X+1)*j+(i+1))*2].add(zero);
        massMatrix[((NUM_X+1)*j+i)*2+1][((NUM_X+1)*j+(i+1))*2+1] = massMatrix[((NUM_X+1)*j+i)*2+1][((NUM_X+1)*j+(i+1))*2+1].add(one);
        

        massMatrix[((NUM_X+1)*j+(i+1))*2][((NUM_X+1)*j+i)*2] = massMatrix[((NUM_X+1)*j+i)*2][((NUM_X+1)*j+(i+1))*2];
        massMatrix[((NUM_X+1)*j+(i+1))*2+1][((NUM_X+1)*j+i)*2] = massMatrix[((NUM_X+1)*j+i)*2][((NUM_X+1)*j+(i+1))*2+1];
        massMatrix[((NUM_X+1)*j+(i+1))*2][((NUM_X+1)*j+i)*2+1] = massMatrix[((NUM_X+1)*j+i)*2+1][((NUM_X+1)*j+(i+1))*2];
        massMatrix[((NUM_X+1)*j+(i+1))*2+1][((NUM_X+1)*j+i)*2+1] = massMatrix[((NUM_X+1)*j+i)*2+1][((NUM_X+1)*j+(i+1))*2+1];
        
        //im
        massMatrix[((NUM_X+1)*j+i)*2][((NUM_X+1)*(j+1)+(i+1))*2] = massMatrix[((NUM_X+1)*j+i)*2][((NUM_X+1)*(j+1)+(i+1))*2].add(one);
        massMatrix[((NUM_X+1)*j+i)*2][((NUM_X+1)*(j+1)+(i+1))*2+1] = massMatrix[((NUM_X+1)*j+i)*2][((NUM_X+1)*(j+1)+(i+1))*2+1].add(zero);
        massMatrix[((NUM_X+1)*j+i)*2+1][((NUM_X+1)*(j+1)+(i+1))*2] = massMatrix[((NUM_X+1)*j+i)*2+1][((NUM_X+1)*(j+1)+(i+1))*2].add(zero);
        massMatrix[((NUM_X+1)*j+i)*2+1][((NUM_X+1)*(j+1)+(i+1))*2+1] = massMatrix[((NUM_X+1)*j+i)*2+1][((NUM_X+1)*(j+1)+(i+1))*2+1].add(one);

        massMatrix[((NUM_X+1)*(j+1)+(i+1))*2][((NUM_X+1)*j+i)*2] = massMatrix[((NUM_X+1)*j+i)*2][((NUM_X+1)*(j+1)+(i+1))*2];
        massMatrix[((NUM_X+1)*(j+1)+(i+1))*2+1][((NUM_X+1)*j+i)*2] = massMatrix[((NUM_X+1)*j+i)*2][((NUM_X+1)*(j+1)+(i+1))*2+1];
        massMatrix[((NUM_X+1)*(j+1)+(i+1))*2][((NUM_X+1)*j+i)*2+1] = massMatrix[((NUM_X+1)*j+i)*2+1][((NUM_X+1)*(j+1)+(i+1))*2];
        massMatrix[((NUM_X+1)*(j+1)+(i+1))*2+1][((NUM_X+1)*j+i)*2+1] = massMatrix[((NUM_X+1)*j+i)*2+1][((NUM_X+1)*(j+1)+(i+1))*2+1];
        
        //jm
        massMatrix[((NUM_X+1)*j+(i+1))*2][((NUM_X+1)*(j+1)+(i+1))*2] = massMatrix[((NUM_X+1)*j+(i+1))*2][((NUM_X+1)*(j+1)+(i+1))*2].add(one);
        massMatrix[((NUM_X+1)*j+(i+1))*2][((NUM_X+1)*(j+1)+(i+1))*2+1] = massMatrix[((NUM_X+1)*j+(i+1))*2][((NUM_X+1)*(j+1)+(i+1))*2+1].add(zero);
        massMatrix[((NUM_X+1)*j+(i+1))*2+1][((NUM_X+1)*(j+1)+(i+1))*2] = massMatrix[((NUM_X+1)*j+(i+1))*2+1][((NUM_X+1)*(j+1)+(i+1))*2].add(zero);
        massMatrix[((NUM_X+1)*j+(i+1))*2+1][((NUM_X+1)*(j+1)+(i+1))*2+1] = massMatrix[((NUM_X+1)*j+(i+1))*2+1][((NUM_X+1)*(j+1)+(i+1))*2+1].add(one);

        massMatrix[((NUM_X+1)*(j+1)+(i+1))*2][((NUM_X+1)*j+(i+1))*2] = massMatrix[((NUM_X+1)*j+(i+1))*2][((NUM_X+1)*(j+1)+(i+1))*2];
        massMatrix[((NUM_X+1)*(j+1)+(i+1))*2+1][((NUM_X+1)*j+(i+1))*2] = massMatrix[((NUM_X+1)*j+(i+1))*2][((NUM_X+1)*(j+1)+(i+1))*2+1];
        massMatrix[((NUM_X+1)*(j+1)+(i+1))*2][((NUM_X+1)*j+(i+1))*2+1] = massMatrix[((NUM_X+1)*j+(i+1))*2+1][((NUM_X+1)*(j+1)+(i+1))*2];
        massMatrix[((NUM_X+1)*(j+1)+(i+1))*2+1][((NUM_X+1)*j+(i+1))*2+1] = massMatrix[((NUM_X+1)*j+(i+1))*2+1][((NUM_X+1)*(j+1)+(i+1))*2+1];
    }

    private static void setElementOneMassMatrix(BigReal[][] massMatrix, int i, int j, BigReal two,
            BigReal zero, BigReal one) {
        //element 1
        //
        // i
        // |\
        // | \
        // |__\
        // m   j

        //top left corner
        massMatrix[((NUM_X+1)*j+i)*2][((NUM_X+1)*j+i)*2] = massMatrix[((NUM_X+1)*j+i)*2][((NUM_X+1)*j+i)*2].add(two);
        massMatrix[((NUM_X+1)*j+i)*2][((NUM_X+1)*j+i)*2+1] = massMatrix[((NUM_X+1)*j+i)*2][((NUM_X+1)*j+i)*2+1].add(zero);
        massMatrix[((NUM_X+1)*j+i)*2+1][((NUM_X+1)*j+i)*2] = massMatrix[((NUM_X+1)*j+i)*2+1][((NUM_X+1)*j+i)*2].add(zero);
        massMatrix[((NUM_X+1)*j+i)*2+1][((NUM_X+1)*j+i)*2+1] = massMatrix[((NUM_X+1)*j+i)*2+1][((NUM_X+1)*j+i)*2+1].add(two);

        //bottom left corner
        massMatrix[((NUM_X+1)*(j+1)+i)*2][((NUM_X+1)*(j+1)+i)*2] = massMatrix[((NUM_X+1)*(j+1)+i)*2][((NUM_X+1)*(j+1)+i)*2].add(two);
        massMatrix[((NUM_X+1)*(j+1)+i)*2][((NUM_X+1)*(j+1)+i)*2+1] = massMatrix[((NUM_X+1)*(j+1)+i)*2][((NUM_X+1)*(j+1)+i)*2+1].add(zero);
        massMatrix[((NUM_X+1)*(j+1)+i)*2+1][((NUM_X+1)*(j+1)+i)*2] = massMatrix[((NUM_X+1)*(j+1)+i)*2+1][((NUM_X+1)*(j+1)+i)*2].add(zero);
        massMatrix[((NUM_X+1)*(j+1)+i)*2+1][((NUM_X+1)*(j+1)+i)*2+1] = massMatrix[((NUM_X+1)*(j+1)+i)*2+1][((NUM_X+1)*(j+1)+i)*2+1].add(two);

        //bottom right corner
        massMatrix[((NUM_X+1)*(j+1)+(i+1))*2][((NUM_X+1)*(j+1)+(i+1))*2] = massMatrix[((NUM_X+1)*(j+1)+(i+1))*2][((NUM_X+1)*(j+1)+(i+1))*2].add(two);
        massMatrix[((NUM_X+1)*(j+1)+(i+1))*2][((NUM_X+1)*(j+1)+(i+1))*2+1] = massMatrix[((NUM_X+1)*(j+1)+(i+1))*2][((NUM_X+1)*(j+1)+(i+1))*2+1].add(zero);
        massMatrix[((NUM_X+1)*(j+1)+(i+1))*2+1][((NUM_X+1)*(j+1)+(i+1))*2] = massMatrix[((NUM_X+1)*(j+1)+(i+1))*2+1][((NUM_X+1)*(j+1)+(i+1))*2].add(zero);
        massMatrix[((NUM_X+1)*(j+1)+(i+1))*2+1][((NUM_X+1)*(j+1)+(i+1))*2+1] = massMatrix[((NUM_X+1)*(j+1)+(i+1))*2+1][((NUM_X+1)*(j+1)+(i+1))*2+1].add(two);

        //ij
        massMatrix[((NUM_X+1)*j+i)*2][((NUM_X+1)*(j+1)+(i+1))*2] = massMatrix[((NUM_X+1)*j+i)*2][((NUM_X+1)*(j+1)+(i+1))*2].add(one);
        massMatrix[((NUM_X+1)*j+i)*2][((NUM_X+1)*(j+1)+(i+1))*2+1] = massMatrix[((NUM_X+1)*j+i)*2][((NUM_X+1)*(j+1)+(i+1))*2+1].add(zero);
        massMatrix[((NUM_X+1)*j+i)*2+1][((NUM_X+1)*(j+1)+(i+1))*2] = massMatrix[((NUM_X+1)*j+i)*2+1][((NUM_X+1)*(j+1)+(i+1))*2].add(zero);
        massMatrix[((NUM_X+1)*j+i)*2+1][((NUM_X+1)*(j+1)+(i+1))*2+1] = massMatrix[((NUM_X+1)*j+i)*2+1][((NUM_X+1)*(j+1)+(i+1))*2+1].add(one);

        massMatrix[((NUM_X+1)*(j+1)+(i+1))*2][((NUM_X+1)*j+i)*2] = massMatrix[((NUM_X+1)*j+i)*2][((NUM_X+1)*(j+1)+(i+1))*2];
        massMatrix[((NUM_X+1)*(j+1)+(i+1))*2+1][((NUM_X+1)*j+i)*2] = massMatrix[((NUM_X+1)*j+i)*2][((NUM_X+1)*(j+1)+(i+1))*2+1];
        massMatrix[((NUM_X+1)*(j+1)+(i+1))*2][((NUM_X+1)*j+i)*2+1] = massMatrix[((NUM_X+1)*j+i)*2+1][((NUM_X+1)*(j+1)+(i+1))*2];
        massMatrix[((NUM_X+1)*(j+1)+(i+1))*2+1][((NUM_X+1)*j+i)*2+1] = massMatrix[((NUM_X+1)*j+i)*2+1][((NUM_X+1)*(j+1)+(i+1))*2+1];

        //im
        massMatrix[((NUM_X+1)*j+i)*2][((NUM_X+1)*(j+1)+i)*2] = massMatrix[((NUM_X+1)*j+i)*2][((NUM_X+1)*(j+1)+i)*2].add(one);
        massMatrix[((NUM_X+1)*j+i)*2][((NUM_X+1)*(j+1)+i)*2+1] = massMatrix[((NUM_X+1)*j+i)*2][((NUM_X+1)*(j+1)+i)*2+1].add(zero);
        massMatrix[((NUM_X+1)*j+i)*2+1][((NUM_X+1)*(j+1)+i)*2] = massMatrix[((NUM_X+1)*j+i)*2+1][((NUM_X+1)*(j+1)+i)*2].add(zero);
        massMatrix[((NUM_X+1)*j+i)*2+1][((NUM_X+1)*(j+1)+i)*2+1] = massMatrix[((NUM_X+1)*j+i)*2+1][((NUM_X+1)*(j+1)+i)*2+1].add(one);

        massMatrix[((NUM_X+1)*(j+1)+i)*2][((NUM_X+1)*j+i)*2] = massMatrix[((NUM_X+1)*j+i)*2][((NUM_X+1)*(j+1)+i)*2];
        massMatrix[((NUM_X+1)*(j+1)+i)*2+1][((NUM_X+1)*j+i)*2] = massMatrix[((NUM_X+1)*j+i)*2][((NUM_X+1)*(j+1)+i)*2+1];
        massMatrix[((NUM_X+1)*(j+1)+i)*2][((NUM_X+1)*j+i)*2+1] = massMatrix[((NUM_X+1)*j+i)*2+1][((NUM_X+1)*(j+1)+i)*2];
        massMatrix[((NUM_X+1)*(j+1)+i)*2+1][((NUM_X+1)*j+i)*2+1] = massMatrix[((NUM_X+1)*j+i)*2+1][((NUM_X+1)*(j+1)+i)*2+1];
        
        //jm
        massMatrix[((NUM_X+1)*(j+1)+i)*2][((NUM_X+1)*(j+1)+(i+1))*2] = massMatrix[((NUM_X+1)*(j+1)+i)*2][((NUM_X+1)*(j+1)+(i+1))*2].add(one);
        massMatrix[((NUM_X+1)*(j+1)+i)*2+1][((NUM_X+1)*(j+1)+(i+1))*2] = massMatrix[((NUM_X+1)*(j+1)+i)*2+1][((NUM_X+1)*(j+1)+(i+1))*2].add(zero);
        massMatrix[((NUM_X+1)*(j+1)+i)*2][((NUM_X+1)*(j+1)+(i+1))*2+1] = massMatrix[((NUM_X+1)*(j+1)+i)*2][((NUM_X+1)*(j+1)+(i+1))*2+1].add(zero);
        massMatrix[((NUM_X+1)*(j+1)+i)*2+1][((NUM_X+1)*(j+1)+(i+1))*2+1] = massMatrix[((NUM_X+1)*(j+1)+i)*2+1][((NUM_X+1)*(j+1)+(i+1))*2+1].add(one);

        massMatrix[((NUM_X+1)*(j+1)+(i+1))*2][((NUM_X+1)*(j+1)+i)*2] = massMatrix[((NUM_X+1)*(j+1)+i)*2][((NUM_X+1)*(j+1)+(i+1))*2];
        massMatrix[((NUM_X+1)*(j+1)+(i+1))*2][((NUM_X+1)*(j+1)+i)*2+1] = massMatrix[((NUM_X+1)*(j+1)+i)*2+1][((NUM_X+1)*(j+1)+(i+1))*2];
        massMatrix[((NUM_X+1)*(j+1)+(i+1))*2+1][((NUM_X+1)*(j+1)+i)*2] = massMatrix[((NUM_X+1)*(j+1)+i)*2][((NUM_X+1)*(j+1)+(i+1))*2+1];
        massMatrix[((NUM_X+1)*(j+1)+(i+1))*2+1][((NUM_X+1)*(j+1)+i)*2+1] = massMatrix[((NUM_X+1)*(j+1)+i)*2+1][((NUM_X+1)*(j+1)+(i+1))*2+1];
    }

    private static void printMatrix(FieldMatrix<BigReal> fieldMatrix) {
        
        writeToFile("====  ====");
        StringBuilder str = new StringBuilder();
        for (BigReal[] column : fieldMatrix.getData()) {
            
            for (BigReal unit : column) {
                str.append(df.format(unit.bigDecimalValue())).append("\t");
            }
            str.append("\n");
        }
        writeToFile(str.toString());
        
    }

    private static void writeToFile(String string) {
        UtilWriter.writeToFile(string);
        UtilWriter.writeToFile("\n");
    }

    private static void setDecimalFormat() {
        df = new DecimalFormat();
        df.setMaximumFractionDigits(10);
        df.setMinimumFractionDigits(0);
        df.setGroupingUsed(false);
    }
}