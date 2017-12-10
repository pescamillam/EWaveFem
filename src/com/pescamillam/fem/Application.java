package com.pescamillam.fem;

import java.awt.Canvas;
import java.awt.Color;
import java.awt.Graphics;
import java.awt.image.BufferStrategy;
import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import javax.swing.JFrame;

import org.apache.commons.math3.linear.FieldLUDecomposition;
import org.apache.commons.math3.linear.FieldMatrix;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.util.BigReal;

import com.pescamillam.fem.element.Cst;
import com.pescamillam.fem.element.Point;

public class Application {
    
    private static final BigReal TWO_BIG_REAL = new BigReal("2");
    private static final BigReal MINUS_ONE = new BigReal("-1");

    private static final int numX = 15;
    private static final int numY = 15;
    
    public static void main(String[] args) {
        Long startTime = System.nanoTime();
        
        //parameters
        //thickness t
        BigReal thickness = new BigReal("1");
        //elasticity module E
        BigReal elasticity = new BigReal("30000000");
        //density p
        BigReal density = new BigReal(0.00073);
        //Area A
        BigReal area = new BigReal("1");
        //poisson ratio v
        BigReal poisson = new BigReal("0.3");
        
        //delta time
        BigReal deltaTime = new BigReal("0.025");
        BigReal deltaTimeSquare = deltaTime.multiply(deltaTime); 
        
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
        FieldMatrix<BigReal> acceleration;

        List<Cst> elements = new ArrayList<>();
        
        //Stiffness matrix
        BigReal[][] stiffnessMatrix = new BigReal[2*numX*numY+2*numX+2*numY+2][2*numX*numY+2*numX+2*numY+2];
        
        //mass matrix
        BigReal[][] massMatrix = new BigReal[2*numX*numY+2*numX+2*numY+2][2*numX*numY+2*numX+2*numY+2];

        for (BigReal[] row : stiffnessMatrix) {
            Arrays.fill(row, BigReal.ZERO);
        }

        for (BigReal[] row : massMatrix) {
            Arrays.fill(row, BigReal.ZERO);
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
        BigReal betaJel1 = new BigReal("10");
        //beta m: y_i - y_j
        BigReal betaMel1 = new BigReal("-10");
        
        //gamma_i: x_m - x_j
        BigReal gammaIel1 = new BigReal("-10");
        //gamma_j: x_i - x_m
        BigReal gammaJel1 = new BigReal("0");
        //gamma_m: x_j - x_i
        BigReal gammaMel1 = new BigReal("10");
        
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
        BigReal m31el1 = m13el1;
        //3-2
        BigReal m32 = m23el1;
        //3-3
        BigReal m33el1 = betaJel1.multiply(betaJel1).multiply(poisson1mv).add(gammaJel1.multiply(gammaJel1).multiply(poisson1m2vo2));
        //3-4
        BigReal m34el1 = betaJel1.multiply(gammaJel1).multiply(poisson).add(betaJel1.multiply(gammaJel1).multiply(poisson1m2vo2));
        //3-5
        BigReal m35el1 = betaJel1.multiply(betaMel1).multiply(poisson1mv).add(gammaJel1.multiply(gammaMel1).multiply(poisson1m2vo2));
        //3-6
        BigReal m36el1 = betaJel1.multiply(gammaMel1).multiply(poisson).add(betaMel1.multiply(gammaJel1).multiply(poisson1m2vo2));
        
        //4-1
        BigReal m41el1 = m14el1;
        //4-2
        BigReal m42el1 = m24el1;
        //4-3
        BigReal m43el1 = m34el1;
        //4-4
        BigReal m44el1 = gammaJel1.multiply(gammaJel1).multiply(poisson1mv).add(betaJel1.multiply(betaJel1).multiply(poisson1m2vo2));
        //4-5
        BigReal m45el1 = betaMel1.multiply(gammaJel1).multiply(poisson).add(betaJel1.multiply(gammaMel1).multiply(poisson1m2vo2));
        //4-6
        BigReal m46el1 = gammaIel1.multiply(gammaMel1).multiply(poisson1mv).add(betaJel1.multiply(betaMel1).multiply(poisson1m2vo2));
        
        //5-1
        BigReal m51el1 = m15el1;
        //5-2
        BigReal m52el1 = m25el1;
        //5-3
        BigReal m53el1 = m35el1;
        //5-4
        BigReal m54el1 = m45el1;
        //5-5
        BigReal m55el1 = betaMel1.multiply(betaMel1).multiply(poisson1mv).add(gammaMel1.multiply(gammaMel1).multiply(poisson1m2vo2));
        //5-6
        BigReal m56el1 = betaMel1.multiply(gammaMel1).multiply(poisson).add(betaMel1.multiply(gammaMel1).multiply(poisson1m2vo2));
        
        //6-1
        BigReal m61el1 = m16el1;
        //6-2
        BigReal m62el1 = m26el1;
        //6-3
        BigReal m63el1 = m36el1;
        //6-4
        BigReal m64el1 = m46el1;
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
        BigReal betaIel2 = new BigReal("-10");
        //beta j: y_m - y_i
        BigReal betaJel2 = new BigReal("10");
        //beta m: y_i - y_j
        BigReal betaMel2 = new BigReal("0");
        
        //gamma_i: x_m - x_j
        BigReal gammaIel2 = new BigReal("0");
        //gamma_j: x_i - x_m
        BigReal gammaJel2 = new BigReal("-10");
        //gamma_m: x_j - x_i
        BigReal gammaMel2 = new BigReal("10");
        
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
        BigReal m31el2 = m13el2;
        //3-2
        BigReal m32el2 = m23el2;
        //3-3
        BigReal m33el2 = betaJel2.multiply(betaJel2).multiply(poisson1mv).add(gammaJel2.multiply(gammaJel2).multiply(poisson1m2vo2));
        //3-4
        BigReal m34el2 = betaJel2.multiply(gammaJel2).multiply(poisson).add(betaJel2.multiply(gammaJel2).multiply(poisson1m2vo2));
        //3-5
        BigReal m35el2 = betaJel2.multiply(betaMel2).multiply(poisson1mv).add(gammaJel2.multiply(gammaMel2).multiply(poisson1m2vo2));
        //3-6
        BigReal m36el2 = betaJel2.multiply(gammaMel2).multiply(poisson).add(betaMel2.multiply(gammaJel2).multiply(poisson1m2vo2));
        
        //4-1
        BigReal m41el2 = m14el2;
        //4-2
        BigReal m42el2 = m24el2;
        //4-3
        BigReal m43el2 = m34el2;
        //4-4
        BigReal m44el2 = gammaJel2.multiply(gammaJel2).multiply(poisson1mv).add(betaJel2.multiply(betaJel2).multiply(poisson1m2vo2));
        //4-5
        BigReal m45el2 = betaMel2.multiply(gammaJel2).multiply(poisson).add(betaJel2.multiply(gammaMel2).multiply(poisson1m2vo2));
        //4-6
        BigReal m46el2 = gammaIel2.multiply(gammaMel2).multiply(poisson1mv).add(betaJel2.multiply(betaMel2).multiply(poisson1m2vo2));
        
        //5-1
        BigReal m51el2 = m15el2;
        //5-2
        BigReal m52el2 = m25el2;
        //5-3
        BigReal m53el2 = m35el2;
        //5-4
        BigReal m54el2 = m45el2;
        //5-5
        BigReal m55el2 = betaMel2.multiply(betaMel2).multiply(poisson1mv).add(gammaMel2.multiply(gammaMel2).multiply(poisson1m2vo2));
        //5-6
        BigReal m56el2 = betaMel2.multiply(gammaMel2).multiply(poisson).add(betaMel2.multiply(gammaMel2).multiply(poisson1m2vo2));
        
        //6-1
        BigReal m61el2 = m16el2;
        //6-2
        BigReal m62el2 = m26el2;
        //6-3
        BigReal m63el2 = m36el2;
        //6-4
        BigReal m64el2 = m46el2;
        //6-5
        BigReal m65el2 = m56el2;
        //6-6
        BigReal m66el2 = gammaMel2.multiply(gammaMel2).multiply(poisson1mv).add(betaMel2.multiply(betaMel2).multiply(poisson1m2vo2));


        //points matrix
        for (int i = 0; i < numX; i++) {
            for (int j = 0; j < numY; j++) {
                Cst element = new Cst(new Point(new BigDecimal(i*10), new BigDecimal(j*10)), 
                        new Point(new BigDecimal(i*10+10), new BigDecimal(j*10+10)), 
                        new Point(new BigDecimal(i*10), new BigDecimal(j*10+10)));
                elements.add(element);

//                System.out.println("==== 1 elem ====");
//                System.out.println(m11 + "\t" + m12 + "\t" + m13 + "\t" + m14 + "\t" + m15 + "\t" + m16);
//                System.out.println(m21 + "\t" + m22 + "\t" + m23 + "\t" + m24 + "\t" + m25 + "\t" + m26);
//                System.out.println(m31 + "\t" + m32 + "\t" + m33 + "\t" + m34 + "\t" + m35 + "\t" + m36);
//                System.out.println(m41 + "\t" + m42 + "\t" + m43 + "\t" + m44 + "\t" + m45 + "\t" + m46);
//                System.out.println(m51 + "\t" + m52 + "\t" + m53 + "\t" + m54 + "\t" + m55 + "\t" + m56);
//                System.out.println(m61 + "\t" + m62 + "\t" + m63 + "\t" + m64 + "\t" + m65 + "\t" + m66);
                

                //top left corner
                stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*j+i)*2] = stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*j+i)*2].add(m11el1);
                stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*j+i)*2+1] = stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*j+i)*2+1].add(m12el1);
                stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*j+i)*2] = stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*j+i)*2].add(m21el1);
                stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*j+i)*2+1] = stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*j+i)*2+1].add(m22el1);

                //bottom left corner
                stiffnessMatrix[((numX+1)*(j+1)+i)*2][((numX+1)*(j+1)+i)*2] = stiffnessMatrix[((numX+1)*(j+1)+i)*2][((numX+1)*(j+1)+i)*2].add(m55el1);
                stiffnessMatrix[((numX+1)*(j+1)+i)*2][((numX+1)*(j+1)+i)*2+1] = stiffnessMatrix[((numX+1)*(j+1)+i)*2][((numX+1)*(j+1)+i)*2+1].add(m56el1);
                stiffnessMatrix[((numX+1)*(j+1)+i)*2+1][((numX+1)*(j+1)+i)*2] = stiffnessMatrix[((numX+1)*(j+1)+i)*2+1][((numX+1)*(j+1)+i)*2].add(m65el1);
                stiffnessMatrix[((numX+1)*(j+1)+i)*2+1][((numX+1)*(j+1)+i)*2+1] = stiffnessMatrix[((numX+1)*(j+1)+i)*2+1][((numX+1)*(j+1)+i)*2+1].add(m66el1);

                //bottom right corner
                stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2][((numX+1)*(j+1)+(i+1))*2] = stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2][((numX+1)*(j+1)+(i+1))*2].add(m33el1);
                stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2][((numX+1)*(j+1)+(i+1))*2+1] = stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2][((numX+1)*(j+1)+(i+1))*2+1].add(m34el1);
                stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2+1][((numX+1)*(j+1)+(i+1))*2] = stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2+1][((numX+1)*(j+1)+(i+1))*2].add(m43el1);
                stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2+1][((numX+1)*(j+1)+(i+1))*2+1] = stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2+1][((numX+1)*(j+1)+(i+1))*2+1].add(m44el1);

                //ij
                stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*(j+1)+(i+1))*2] = stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*(j+1)+(i+1))*2].add(m13el1);
                stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*(j+1)+(i+1))*2+1] = stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*(j+1)+(i+1))*2+1].add(m14el1);
                stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*(j+1)+(i+1))*2] = stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*(j+1)+(i+1))*2].add(m23el1);
                stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*(j+1)+(i+1))*2+1] = stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*(j+1)+(i+1))*2+1].add(m24el1);

                stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2][((numX+1)*j+i)*2] = stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*(j+1)+(i+1))*2];
                stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2+1][((numX+1)*j+i)*2] = stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*(j+1)+(i+1))*2+1];
                stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2][((numX+1)*j+i)*2+1] = stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*(j+1)+(i+1))*2];
                stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2+1][((numX+1)*j+i)*2+1] = stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*(j+1)+(i+1))*2+1];

                //im
                stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*(j+1)+i)*2] = stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*(j+1)+i)*2].add(m15el1);
                stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*(j+1)+i)*2+1] = stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*(j+1)+i)*2+1].add(m16el1);
                stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*(j+1)+i)*2] = stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*(j+1)+i)*2].add(m25el1);
                stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*(j+1)+i)*2+1] = stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*(j+1)+i)*2+1].add(m26el1);

                stiffnessMatrix[((numX+1)*(j+1)+i)*2][((numX+1)*j+i)*2] = stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*(j+1)+i)*2];
                stiffnessMatrix[((numX+1)*(j+1)+i)*2+1][((numX+1)*j+i)*2] = stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*(j+1)+i)*2+1];
                stiffnessMatrix[((numX+1)*(j+1)+i)*2][((numX+1)*j+i)*2+1] = stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*(j+1)+i)*2];
                stiffnessMatrix[((numX+1)*(j+1)+i)*2+1][((numX+1)*j+i)*2+1] = stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*(j+1)+i)*2+1];
                
                //jm
                stiffnessMatrix[((numX+1)*(j+1)+i)*2][((numX+1)*(j+1)+(i+1))*2] = stiffnessMatrix[((numX+1)*(j+1)+i)*2][((numX+1)*(j+1)+(i+1))*2].add(m35el1);
                stiffnessMatrix[((numX+1)*(j+1)+i)*2+1][((numX+1)*(j+1)+(i+1))*2] = stiffnessMatrix[((numX+1)*(j+1)+i)*2+1][((numX+1)*(j+1)+(i+1))*2].add(m45el1);
                stiffnessMatrix[((numX+1)*(j+1)+i)*2][((numX+1)*(j+1)+(i+1))*2+1] = stiffnessMatrix[((numX+1)*(j+1)+i)*2][((numX+1)*(j+1)+(i+1))*2+1].add(m36el1);
                stiffnessMatrix[((numX+1)*(j+1)+i)*2+1][((numX+1)*(j+1)+(i+1))*2+1] = stiffnessMatrix[((numX+1)*(j+1)+i)*2+1][((numX+1)*(j+1)+(i+1))*2+1].add(m46el1);

                stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2][((numX+1)*(j+1)+i)*2] = stiffnessMatrix[((numX+1)*(j+1)+i)*2][((numX+1)*(j+1)+(i+1))*2];
                stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2][((numX+1)*(j+1)+i)*2+1] = stiffnessMatrix[((numX+1)*(j+1)+i)*2+1][((numX+1)*(j+1)+(i+1))*2];
                stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2+1][((numX+1)*(j+1)+i)*2] = stiffnessMatrix[((numX+1)*(j+1)+i)*2][((numX+1)*(j+1)+(i+1))*2+1];
                stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2+1][((numX+1)*(j+1)+i)*2+1] = stiffnessMatrix[((numX+1)*(j+1)+i)*2+1][((numX+1)*(j+1)+(i+1))*2+1];
                
                Cst element2 = new Cst(new Point(new BigDecimal(i*10), new BigDecimal(j*10)), 
                        new Point(new BigDecimal(i*10+10), new BigDecimal(j*10)), 
                        new Point(new BigDecimal(i*10+10), new BigDecimal(j*10+10)));
                elements.add(element2);
                //top left corner
                stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*j+i)*2] = stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*j+i)*2].add(m11el2);
                stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*j+i)*2+1] = stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*j+i)*2+1].add(m12el2);
                stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*j+i)*2] = stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*j+i)*2].add(m21el2);
                stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*j+i)*2+1] = stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*j+i)*2+1].add(m22el2);

                //top right corner
                stiffnessMatrix[((numX+1)*j+(i+1))*2][((numX+1)*j+(i+1))*2] = stiffnessMatrix[((numX+1)*j+(i+1))*2][((numX+1)*j+(i+1))*2].add(m33el2);
                stiffnessMatrix[((numX+1)*j+(i+1))*2][((numX+1)*j+(i+1))*2+1] = stiffnessMatrix[((numX+1)*j+(i+1))*2][((numX+1)*j+(i+1))*2+1].add(m34el2);
                stiffnessMatrix[((numX+1)*j+(i+1))*2+1][((numX+1)*j+(i+1))*2] = stiffnessMatrix[((numX+1)*j+(i+1))*2+1][((numX+1)*j+(i+1))*2].add(m43el2);
                stiffnessMatrix[((numX+1)*j+(i+1))*2+1][((numX+1)*j+(i+1))*2+1] = stiffnessMatrix[((numX+1)*j+(i+1))*2+1][((numX+1)*j+(i+1))*2+1].add(m44el2);

                //bottom right corner
                stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2][((numX+1)*(j+1)+(i+1))*2] = stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2][((numX+1)*(j+1)+(i+1))*2].add(m55el2);
                stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2][((numX+1)*(j+1)+(i+1))*2+1] = stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2][((numX+1)*(j+1)+(i+1))*2+1].add(m56el2);
                stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2+1][((numX+1)*(j+1)+(i+1))*2] = stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2+1][((numX+1)*(j+1)+(i+1))*2].add(m65el2);
                stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2+1][((numX+1)*(j+1)+(i+1))*2+1] = stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2+1][((numX+1)*(j+1)+(i+1))*2+1].add(m66el2);
                
                //ij
                stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*j+(i+1))*2] = stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*j+(i+1))*2].add(m13el2);
                stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*j+(i+1))*2+1] = stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*j+(i+1))*2+1].add(m14el2);
                stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*j+(i+1))*2] = stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*j+(i+1))*2].add(m23el2);
                stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*j+(i+1))*2+1] = stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*j+(i+1))*2+1].add(m24el2);
                

                stiffnessMatrix[((numX+1)*j+(i+1))*2][((numX+1)*j+i)*2] = stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*j+(i+1))*2];
                stiffnessMatrix[((numX+1)*j+(i+1))*2+1][((numX+1)*j+i)*2] = stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*j+(i+1))*2+1];
                stiffnessMatrix[((numX+1)*j+(i+1))*2][((numX+1)*j+i)*2+1] = stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*j+(i+1))*2];
                stiffnessMatrix[((numX+1)*j+(i+1))*2+1][((numX+1)*j+i)*2+1] = stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*j+(i+1))*2+1];
                
                //im
                stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*(j+1)+(i+1))*2] = stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*(j+1)+(i+1))*2].add(m15el2);
                stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*(j+1)+(i+1))*2+1] = stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*(j+1)+(i+1))*2+1].add(m16el2);
                stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*(j+1)+(i+1))*2] = stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*(j+1)+(i+1))*2].add(m25el2);
                stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*(j+1)+(i+1))*2+1] = stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*(j+1)+(i+1))*2+1].add(m26el2);

                stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2][((numX+1)*j+i)*2] = stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*(j+1)+(i+1))*2];
                stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2+1][((numX+1)*j+i)*2] = stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*(j+1)+(i+1))*2+1];
                stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2][((numX+1)*j+i)*2+1] = stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*(j+1)+(i+1))*2];
                stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2+1][((numX+1)*j+i)*2+1] = stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*(j+1)+(i+1))*2+1];
                
                //jm
                stiffnessMatrix[((numX+1)*j+(i+1))*2][((numX+1)*(j+1)+(i+1))*2] = stiffnessMatrix[((numX+1)*j+(i+1))*2][((numX+1)*(j+1)+(i+1))*2].add(m35el2);
                stiffnessMatrix[((numX+1)*j+(i+1))*2][((numX+1)*(j+1)+(i+1))*2+1] = stiffnessMatrix[((numX+1)*j+(i+1))*2][((numX+1)*(j+1)+(i+1))*2+1].add(m36el2);
                stiffnessMatrix[((numX+1)*j+(i+1))*2+1][((numX+1)*(j+1)+(i+1))*2] = stiffnessMatrix[((numX+1)*j+(i+1))*2+1][((numX+1)*(j+1)+(i+1))*2].add(m45el2);
                stiffnessMatrix[((numX+1)*j+(i+1))*2+1][((numX+1)*(j+1)+(i+1))*2+1] = stiffnessMatrix[((numX+1)*j+(i+1))*2+1][((numX+1)*(j+1)+(i+1))*2+1].add(m46el2);

                stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2][((numX+1)*j+(i+1))*2] = stiffnessMatrix[((numX+1)*j+(i+1))*2][((numX+1)*(j+1)+(i+1))*2];
                stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2+1][((numX+1)*j+(i+1))*2] = stiffnessMatrix[((numX+1)*j+(i+1))*2][((numX+1)*(j+1)+(i+1))*2+1];
                stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2][((numX+1)*j+(i+1))*2+1] = stiffnessMatrix[((numX+1)*j+(i+1))*2+1][((numX+1)*(j+1)+(i+1))*2];
                stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2+1][((numX+1)*j+(i+1))*2+1] = stiffnessMatrix[((numX+1)*j+(i+1))*2+1][((numX+1)*(j+1)+(i+1))*2+1];

//                System.out.println("==== 2 elem ====");
//                System.out.println(m11 + " " + m12 + " " + m13 + " " + m14 + " " + m15 + " " + m16);
//                System.out.println(m21 + " " + m22 + " " + m23 + " " + m24 + " " + m25 + " " + m26);
//                System.out.println(m31 + " " + m32 + " " + m33 + " " + m34 + " " + m35 + " " + m36);
//                System.out.println(m41 + " " + m42 + " " + m43 + " " + m44 + " " + m45 + " " + m46);
//                System.out.println(m51 + " " + m52 + " " + m53 + " " + m54 + " " + m55 + " " + m56);
//                System.out.println(m61 + " " + m62 + " " + m63 + " " + m64 + " " + m65 + " " + m66);

//                BigDecimal constant = density.multiply(area).multiply(thickness);
                BigReal constant = new BigReal("1");
                BigReal two = new BigReal(2).multiply(constant);
                BigReal zero = new BigReal(0);
                BigReal one = constant;

                //element 1
                //
                // i
                // |\
                // | \
                // |__\
                // m   j

                //top left corner
                massMatrix[((numX+1)*j+i)*2][((numX+1)*j+i)*2] = massMatrix[((numX+1)*j+i)*2][((numX+1)*j+i)*2].add(two);
                massMatrix[((numX+1)*j+i)*2][((numX+1)*j+i)*2+1] = massMatrix[((numX+1)*j+i)*2][((numX+1)*j+i)*2+1].add(zero);
                massMatrix[((numX+1)*j+i)*2+1][((numX+1)*j+i)*2] = massMatrix[((numX+1)*j+i)*2+1][((numX+1)*j+i)*2].add(zero);
                massMatrix[((numX+1)*j+i)*2+1][((numX+1)*j+i)*2+1] = massMatrix[((numX+1)*j+i)*2+1][((numX+1)*j+i)*2+1].add(two);

                //bottom left corner
                massMatrix[((numX+1)*(j+1)+i)*2][((numX+1)*(j+1)+i)*2] = massMatrix[((numX+1)*(j+1)+i)*2][((numX+1)*(j+1)+i)*2].add(two);
                massMatrix[((numX+1)*(j+1)+i)*2][((numX+1)*(j+1)+i)*2+1] = massMatrix[((numX+1)*(j+1)+i)*2][((numX+1)*(j+1)+i)*2+1].add(zero);
                massMatrix[((numX+1)*(j+1)+i)*2+1][((numX+1)*(j+1)+i)*2] = massMatrix[((numX+1)*(j+1)+i)*2+1][((numX+1)*(j+1)+i)*2].add(zero);
                massMatrix[((numX+1)*(j+1)+i)*2+1][((numX+1)*(j+1)+i)*2+1] = massMatrix[((numX+1)*(j+1)+i)*2+1][((numX+1)*(j+1)+i)*2+1].add(two);

                //bottom right corner
                massMatrix[((numX+1)*(j+1)+(i+1))*2][((numX+1)*(j+1)+(i+1))*2] = massMatrix[((numX+1)*(j+1)+(i+1))*2][((numX+1)*(j+1)+(i+1))*2].add(two);
                massMatrix[((numX+1)*(j+1)+(i+1))*2][((numX+1)*(j+1)+(i+1))*2+1] = massMatrix[((numX+1)*(j+1)+(i+1))*2][((numX+1)*(j+1)+(i+1))*2+1].add(zero);
                massMatrix[((numX+1)*(j+1)+(i+1))*2+1][((numX+1)*(j+1)+(i+1))*2] = massMatrix[((numX+1)*(j+1)+(i+1))*2+1][((numX+1)*(j+1)+(i+1))*2].add(zero);
                massMatrix[((numX+1)*(j+1)+(i+1))*2+1][((numX+1)*(j+1)+(i+1))*2+1] = massMatrix[((numX+1)*(j+1)+(i+1))*2+1][((numX+1)*(j+1)+(i+1))*2+1].add(two);

                //ij
                massMatrix[((numX+1)*j+i)*2][((numX+1)*(j+1)+(i+1))*2] = massMatrix[((numX+1)*j+i)*2][((numX+1)*(j+1)+(i+1))*2].add(one);
                massMatrix[((numX+1)*j+i)*2][((numX+1)*(j+1)+(i+1))*2+1] = massMatrix[((numX+1)*j+i)*2][((numX+1)*(j+1)+(i+1))*2+1].add(zero);
                massMatrix[((numX+1)*j+i)*2+1][((numX+1)*(j+1)+(i+1))*2] = massMatrix[((numX+1)*j+i)*2+1][((numX+1)*(j+1)+(i+1))*2].add(zero);
                massMatrix[((numX+1)*j+i)*2+1][((numX+1)*(j+1)+(i+1))*2+1] = massMatrix[((numX+1)*j+i)*2+1][((numX+1)*(j+1)+(i+1))*2+1].add(one);

                massMatrix[((numX+1)*(j+1)+(i+1))*2][((numX+1)*j+i)*2] = massMatrix[((numX+1)*j+i)*2][((numX+1)*(j+1)+(i+1))*2];
                massMatrix[((numX+1)*(j+1)+(i+1))*2+1][((numX+1)*j+i)*2] = massMatrix[((numX+1)*j+i)*2][((numX+1)*(j+1)+(i+1))*2+1];
                massMatrix[((numX+1)*(j+1)+(i+1))*2][((numX+1)*j+i)*2+1] = massMatrix[((numX+1)*j+i)*2+1][((numX+1)*(j+1)+(i+1))*2];
                massMatrix[((numX+1)*(j+1)+(i+1))*2+1][((numX+1)*j+i)*2+1] = massMatrix[((numX+1)*j+i)*2+1][((numX+1)*(j+1)+(i+1))*2+1];

                //im
                massMatrix[((numX+1)*j+i)*2][((numX+1)*(j+1)+i)*2] = massMatrix[((numX+1)*j+i)*2][((numX+1)*(j+1)+i)*2].add(one);
                massMatrix[((numX+1)*j+i)*2][((numX+1)*(j+1)+i)*2+1] = massMatrix[((numX+1)*j+i)*2][((numX+1)*(j+1)+i)*2+1].add(zero);
                massMatrix[((numX+1)*j+i)*2+1][((numX+1)*(j+1)+i)*2] = massMatrix[((numX+1)*j+i)*2+1][((numX+1)*(j+1)+i)*2].add(zero);
                massMatrix[((numX+1)*j+i)*2+1][((numX+1)*(j+1)+i)*2+1] = massMatrix[((numX+1)*j+i)*2+1][((numX+1)*(j+1)+i)*2+1].add(one);

                massMatrix[((numX+1)*(j+1)+i)*2][((numX+1)*j+i)*2] = massMatrix[((numX+1)*j+i)*2][((numX+1)*(j+1)+i)*2];
                massMatrix[((numX+1)*(j+1)+i)*2+1][((numX+1)*j+i)*2] = massMatrix[((numX+1)*j+i)*2][((numX+1)*(j+1)+i)*2+1];
                massMatrix[((numX+1)*(j+1)+i)*2][((numX+1)*j+i)*2+1] = massMatrix[((numX+1)*j+i)*2+1][((numX+1)*(j+1)+i)*2];
                massMatrix[((numX+1)*(j+1)+i)*2+1][((numX+1)*j+i)*2+1] = massMatrix[((numX+1)*j+i)*2+1][((numX+1)*(j+1)+i)*2+1];
                
                //jm
                massMatrix[((numX+1)*(j+1)+i)*2][((numX+1)*(j+1)+(i+1))*2] = massMatrix[((numX+1)*(j+1)+i)*2][((numX+1)*(j+1)+(i+1))*2].add(one);
                massMatrix[((numX+1)*(j+1)+i)*2+1][((numX+1)*(j+1)+(i+1))*2] = massMatrix[((numX+1)*(j+1)+i)*2+1][((numX+1)*(j+1)+(i+1))*2].add(zero);
                massMatrix[((numX+1)*(j+1)+i)*2][((numX+1)*(j+1)+(i+1))*2+1] = massMatrix[((numX+1)*(j+1)+i)*2][((numX+1)*(j+1)+(i+1))*2+1].add(zero);
                massMatrix[((numX+1)*(j+1)+i)*2+1][((numX+1)*(j+1)+(i+1))*2+1] = massMatrix[((numX+1)*(j+1)+i)*2+1][((numX+1)*(j+1)+(i+1))*2+1].add(one);

                massMatrix[((numX+1)*(j+1)+(i+1))*2][((numX+1)*(j+1)+i)*2] = massMatrix[((numX+1)*(j+1)+i)*2][((numX+1)*(j+1)+(i+1))*2];
                massMatrix[((numX+1)*(j+1)+(i+1))*2][((numX+1)*(j+1)+i)*2+1] = massMatrix[((numX+1)*(j+1)+i)*2+1][((numX+1)*(j+1)+(i+1))*2];
                massMatrix[((numX+1)*(j+1)+(i+1))*2+1][((numX+1)*(j+1)+i)*2] = massMatrix[((numX+1)*(j+1)+i)*2][((numX+1)*(j+1)+(i+1))*2+1];
                massMatrix[((numX+1)*(j+1)+(i+1))*2+1][((numX+1)*(j+1)+i)*2+1] = massMatrix[((numX+1)*(j+1)+i)*2+1][((numX+1)*(j+1)+(i+1))*2+1];

                //element 2
                //
                // i___j
                //  \  |
                //   \ |
                //    \|
                //     m

                //top left corner
                massMatrix[((numX+1)*j+i)*2][((numX+1)*j+i)*2] = massMatrix[((numX+1)*j+i)*2][((numX+1)*j+i)*2].add(two);
                massMatrix[((numX+1)*j+i)*2][((numX+1)*j+i)*2+1] = massMatrix[((numX+1)*j+i)*2][((numX+1)*j+i)*2+1].add(zero);
                massMatrix[((numX+1)*j+i)*2+1][((numX+1)*j+i)*2] = massMatrix[((numX+1)*j+i)*2+1][((numX+1)*j+i)*2].add(zero);
                massMatrix[((numX+1)*j+i)*2+1][((numX+1)*j+i)*2+1] = massMatrix[((numX+1)*j+i)*2+1][((numX+1)*j+i)*2+1].add(two);

                //top right corner
                massMatrix[((numX+1)*j+(i+1))*2][((numX+1)*j+(i+1))*2] = massMatrix[((numX+1)*j+(i+1))*2][((numX+1)*j+(i+1))*2].add(two);
                massMatrix[((numX+1)*j+(i+1))*2][((numX+1)*j+(i+1))*2+1] = massMatrix[((numX+1)*j+(i+1))*2][((numX+1)*j+(i+1))*2+1].add(zero);
                massMatrix[((numX+1)*j+(i+1))*2+1][((numX+1)*j+(i+1))*2] = massMatrix[((numX+1)*j+(i+1))*2+1][((numX+1)*j+(i+1))*2].add(zero);
                massMatrix[((numX+1)*j+(i+1))*2+1][((numX+1)*j+(i+1))*2+1] = massMatrix[((numX+1)*j+(i+1))*2+1][((numX+1)*j+(i+1))*2+1].add(two);

                //bottom right corner
                massMatrix[((numX+1)*(j+1)+(i+1))*2][((numX+1)*(j+1)+(i+1))*2] = massMatrix[((numX+1)*(j+1)+(i+1))*2][((numX+1)*(j+1)+(i+1))*2].add(two);
                massMatrix[((numX+1)*(j+1)+(i+1))*2][((numX+1)*(j+1)+(i+1))*2+1] = massMatrix[((numX+1)*(j+1)+(i+1))*2][((numX+1)*(j+1)+(i+1))*2+1].add(zero);
                massMatrix[((numX+1)*(j+1)+(i+1))*2+1][((numX+1)*(j+1)+(i+1))*2] = massMatrix[((numX+1)*(j+1)+(i+1))*2+1][((numX+1)*(j+1)+(i+1))*2].add(zero);
                massMatrix[((numX+1)*(j+1)+(i+1))*2+1][((numX+1)*(j+1)+(i+1))*2+1] = massMatrix[((numX+1)*(j+1)+(i+1))*2+1][((numX+1)*(j+1)+(i+1))*2+1].add(two);
                
                //ij
                massMatrix[((numX+1)*j+i)*2][((numX+1)*j+(i+1))*2] = massMatrix[((numX+1)*j+i)*2][((numX+1)*j+(i+1))*2].add(one);
                massMatrix[((numX+1)*j+i)*2][((numX+1)*j+(i+1))*2+1] = massMatrix[((numX+1)*j+i)*2][((numX+1)*j+(i+1))*2+1].add(zero);
                massMatrix[((numX+1)*j+i)*2+1][((numX+1)*j+(i+1))*2] = massMatrix[((numX+1)*j+i)*2+1][((numX+1)*j+(i+1))*2].add(zero);
                massMatrix[((numX+1)*j+i)*2+1][((numX+1)*j+(i+1))*2+1] = massMatrix[((numX+1)*j+i)*2+1][((numX+1)*j+(i+1))*2+1].add(one);
                

                massMatrix[((numX+1)*j+(i+1))*2][((numX+1)*j+i)*2] = massMatrix[((numX+1)*j+i)*2][((numX+1)*j+(i+1))*2];
                massMatrix[((numX+1)*j+(i+1))*2+1][((numX+1)*j+i)*2] = massMatrix[((numX+1)*j+i)*2][((numX+1)*j+(i+1))*2+1];
                massMatrix[((numX+1)*j+(i+1))*2][((numX+1)*j+i)*2+1] = massMatrix[((numX+1)*j+i)*2+1][((numX+1)*j+(i+1))*2];
                massMatrix[((numX+1)*j+(i+1))*2+1][((numX+1)*j+i)*2+1] = massMatrix[((numX+1)*j+i)*2+1][((numX+1)*j+(i+1))*2+1];
                
                //im
                massMatrix[((numX+1)*j+i)*2][((numX+1)*(j+1)+(i+1))*2] = massMatrix[((numX+1)*j+i)*2][((numX+1)*(j+1)+(i+1))*2].add(one);
                massMatrix[((numX+1)*j+i)*2][((numX+1)*(j+1)+(i+1))*2+1] = massMatrix[((numX+1)*j+i)*2][((numX+1)*(j+1)+(i+1))*2+1].add(zero);
                massMatrix[((numX+1)*j+i)*2+1][((numX+1)*(j+1)+(i+1))*2] = massMatrix[((numX+1)*j+i)*2+1][((numX+1)*(j+1)+(i+1))*2].add(zero);
                massMatrix[((numX+1)*j+i)*2+1][((numX+1)*(j+1)+(i+1))*2+1] = massMatrix[((numX+1)*j+i)*2+1][((numX+1)*(j+1)+(i+1))*2+1].add(one);

                massMatrix[((numX+1)*(j+1)+(i+1))*2][((numX+1)*j+i)*2] = massMatrix[((numX+1)*j+i)*2][((numX+1)*(j+1)+(i+1))*2];
                massMatrix[((numX+1)*(j+1)+(i+1))*2+1][((numX+1)*j+i)*2] = massMatrix[((numX+1)*j+i)*2][((numX+1)*(j+1)+(i+1))*2+1];
                massMatrix[((numX+1)*(j+1)+(i+1))*2][((numX+1)*j+i)*2+1] = massMatrix[((numX+1)*j+i)*2+1][((numX+1)*(j+1)+(i+1))*2];
                massMatrix[((numX+1)*(j+1)+(i+1))*2+1][((numX+1)*j+i)*2+1] = massMatrix[((numX+1)*j+i)*2+1][((numX+1)*(j+1)+(i+1))*2+1];
                
                //jm
                massMatrix[((numX+1)*j+(i+1))*2][((numX+1)*(j+1)+(i+1))*2] = massMatrix[((numX+1)*j+(i+1))*2][((numX+1)*(j+1)+(i+1))*2].add(one);
                massMatrix[((numX+1)*j+(i+1))*2][((numX+1)*(j+1)+(i+1))*2+1] = massMatrix[((numX+1)*j+(i+1))*2][((numX+1)*(j+1)+(i+1))*2+1].add(zero);
                massMatrix[((numX+1)*j+(i+1))*2+1][((numX+1)*(j+1)+(i+1))*2] = massMatrix[((numX+1)*j+(i+1))*2+1][((numX+1)*(j+1)+(i+1))*2].add(zero);
                massMatrix[((numX+1)*j+(i+1))*2+1][((numX+1)*(j+1)+(i+1))*2+1] = massMatrix[((numX+1)*j+(i+1))*2+1][((numX+1)*(j+1)+(i+1))*2+1].add(one);

                massMatrix[((numX+1)*(j+1)+(i+1))*2][((numX+1)*j+(i+1))*2] = massMatrix[((numX+1)*j+(i+1))*2][((numX+1)*(j+1)+(i+1))*2];
                massMatrix[((numX+1)*(j+1)+(i+1))*2+1][((numX+1)*j+(i+1))*2] = massMatrix[((numX+1)*j+(i+1))*2][((numX+1)*(j+1)+(i+1))*2+1];
                massMatrix[((numX+1)*(j+1)+(i+1))*2][((numX+1)*j+(i+1))*2+1] = massMatrix[((numX+1)*j+(i+1))*2+1][((numX+1)*(j+1)+(i+1))*2];
                massMatrix[((numX+1)*(j+1)+(i+1))*2+1][((numX+1)*j+(i+1))*2+1] = massMatrix[((numX+1)*j+(i+1))*2+1][((numX+1)*(j+1)+(i+1))*2+1];
            }
            
        }
        System.out.println("Phase 1 time: " + (System.nanoTime() - startTime)/1000000000.0);
//        System.out.println("==== Stiffness ====");
//        for (BigReal[] column : stiffnessMatrix) {
//            for (BigReal unit : column) {
//                System.out.print(unit + "\t");
//            }
//            System.out.println();
//        }
//        
//        System.out.println("==== mass ====");
//        for (BigReal[] column : massMatrix) {
//            for (BigReal unit : column) {
//                System.out.print(unit.bigDecimalValue() + "\t");
//            }
//            System.out.println();
//        }

        FieldMatrix<BigReal> massFieldMatrix = MatrixUtils.createFieldMatrix(massMatrix);
        System.out.println("Phase 2 time: " + (System.nanoTime() - startTime)/1000000000.0);
        FieldMatrix<BigReal> inverseMassMatrix = new FieldLUDecomposition<BigReal>(MatrixUtils.createFieldMatrix(massMatrix)).getSolver().getInverse();
        System.out.println("Phase 3 time: " + (System.nanoTime() - startTime)/1000000000.0);
        FieldMatrix<BigReal> stiffnessFieldMatrix = MatrixUtils.createFieldMatrix(stiffnessMatrix);
        System.out.println("Phase 4 time: " + (System.nanoTime() - startTime)/1000000000.0);

        
//        System.out.println("==== mass inverse ====");
//        for (BigReal[] column : inverseMassMatrix.getData()) {
//            for (BigReal unit : column) {
//                System.out.print(unit.bigDecimalValue() + "\t");
//            }
//            System.out.println();
//        }
        
        //acceleration = mass inverse * (F0 - [K]{d0}) as {d0} = 0 
        //acceleration = mass inverse * (F0)
        
        //force vectors
        BigReal ONE_HUNDRED = new BigReal("100");
        BigReal[] force1Vector = new BigReal[2*numX*numY+2*numX+2*numY+2];
        Arrays.fill(force1Vector, BigReal.ZERO);
        force1Vector[5] = ONE_HUNDRED;
        force1Vector[4] = ONE_HUNDRED;
        force1Vector[6] = ONE_HUNDRED;
        force1Vector[7] = ONE_HUNDRED;
        BigReal[] force0Vector = new BigReal[2*numX*numY+2*numX+2*numY+2];
        Arrays.fill(force0Vector, BigReal.ZERO);
        FieldMatrix<BigReal>[] force = new FieldMatrix[100];
        force[0] = MatrixUtils.createColumnFieldMatrix(force1Vector);
        for (int i = 1; i < 100; i++) {
            force[i] = MatrixUtils.createColumnFieldMatrix(force0Vector);
        }
        
        acceleration = inverseMassMatrix.multiply(force[0]);
        

//        System.out.println("==== acceleration 0 ====");
//        for (BigReal[] column : acceleration.getData()) {
//            for (BigReal unit : column) {
//                System.out.print(unit.bigDecimalValue() + "\t");
//            }
//            System.out.println();
//        }

        //displacement -1
        FieldMatrix<BigReal> displacementM1 = acceleration.scalarMultiply(deltaTime.multiply(deltaTime).divide(new BigReal("2")));

//        System.out.println("==== displacement -1 ====");
//        for (BigReal[] column : displacementM1.getData()) {
//            for (BigReal unit : column) {
//                System.out.print(unit.bigDecimalValue() + "\t");
//            }
//            System.out.println();
//        }
        
        //displacement 1
        BigReal[] displacement0 = new BigReal[2*numX*numY+2*numX+2*numY+2];
        Arrays.fill(displacement0, BigReal.ZERO);
        FieldMatrix<BigReal>[] displacement = new FieldMatrix[100];
        displacement[0] = MatrixUtils.createColumnFieldMatrix(displacement0);

        displacement[1] = inverseMassMatrix.multiply(
                force[0].scalarMultiply(deltaTime.multiply(deltaTime))
                //.add(massFieldMatrix.scalarMultiply(TWO_BIG_REAL).add(stiffnessFieldMatrix.scalarMultiply(deltaTime.multiply(deltaTime))).multiply(displacementM1))
                .add(massFieldMatrix.multiply(displacementM1).scalarMultiply(MINUS_ONE))
                );
        for (int i = 2; i < 100; i++) {
            displacement[i] = inverseMassMatrix.multiply(
                        force[i-1].scalarMultiply(deltaTime.multiply(deltaTime)).add(
                                massFieldMatrix.scalarMultiply(TWO_BIG_REAL).add(stiffnessFieldMatrix.scalarMultiply(deltaTime.multiply(deltaTime).multiply(MINUS_ONE))).multiply(displacement[i-1])
                                ).add(massFieldMatrix.multiply(displacement[i-2]).scalarMultiply(MINUS_ONE))
                    );
            System.out.println(i);
        }
        
        System.out.println("Total time: " + (System.nanoTime() - startTime)/1000000000.0 + "s");
        printElements(elements, displacement);
        
        System.out.println(elasticity);
    }

//    private static void printMatrix(FieldMatrix<BigReal> fieldMatrix) {
//        System.out.println("====  ====");
//        for (BigReal[] column : fieldMatrix.getData()) {
//            for (BigReal unit : column) {
//                System.out.print(unit.bigDecimalValue() + "\t");
//            }
//            System.out.println();
//        }
//        
//    }

    private static void printElements(List<Cst> elements, FieldMatrix<BigReal>[] displacement) {
        final String title = "Test Window";
        final int width = 1200;
        final int height = width / 16 * 9;

        //Creating the frame.
        JFrame frame = new JFrame(title);

        frame.setSize(width, height);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.setLocationRelativeTo(null);
        frame.setResizable(false);
        frame.setVisible(true);

        //Creating the canvas.
        Canvas canvas = new Canvas();

        canvas.setSize(width, height);
        canvas.setBackground(Color.BLACK);
        canvas.setVisible(true);
        canvas.setFocusable(false);


        //Putting it all together.
        frame.add(canvas);

        canvas.createBufferStrategy(3);

        boolean running = true;

        BufferStrategy bufferStrategy;
        Graphics graphics;
        int i = 0;

        while (running) {
            bufferStrategy = canvas.getBufferStrategy();
            graphics = bufferStrategy.getDrawGraphics();
            graphics.clearRect(0, 0, width, height);

            graphics.setColor(Color.GREEN);
            graphics.clearRect(0, 0, 1000, 1000);

            for (int m = 0; m <= numY; m++) {
                for (int n = 0; n <= numX; n++) {
                    graphics.drawOval(n*30 + displacement[i].getData()[m*numX+n][0].bigDecimalValue().multiply(new BigDecimal("300")).intValue(), m*30 + displacement[i].getData()[m*numX+n+1][0].bigDecimalValue().multiply(new BigDecimal("300")).intValue(), 10, 10);
                }
            }
            
            graphics.drawString("t: " + i, 100, 500);
            
//            for (Cst element : elements) {
//                graphics.drawLine(element.getI().x.intValue(), element.getI().y.intValue()
//                        , element.getJ().x.intValue(), element.getJ().y.intValue());
//                graphics.drawLine(element.getI().x.intValue(), element.getI().y.intValue()
//                        , element.getM().x.intValue(), element.getM().y.intValue());
//                graphics.drawLine(element.getM().x.intValue(), element.getM().y.intValue()
//                        , element.getJ().x.intValue(), element.getJ().y.intValue());
//            }

            bufferStrategy.show();
            graphics.dispose();
            try {
                Thread.sleep(200L);
            } catch (InterruptedException e) {
                // TODO Auto-generated catch block
                e.printStackTrace();
            }
            if (i < 99) {
                i++;
            } else {
                i = 0;
            }
        }
    }
}