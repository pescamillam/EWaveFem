package com.pescamillam.fem;

import java.awt.Canvas;
import java.awt.Color;
import java.awt.Graphics;
import java.awt.Window;
import java.awt.image.BufferStrategy;
import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import javax.swing.JFrame;

import com.pescamillam.fem.element.Cst;
import com.pescamillam.fem.element.Point;

public class Application {
    
    
    
    public static void main(String[] args) {
        System.out.println("hola");
        
        //parameters
        //thickness t
        BigDecimal thickness = new BigDecimal("1");
        //elasticity module E
        BigDecimal elasticity = new BigDecimal("30000000");
        //density p
        BigDecimal density = new BigDecimal(0.00073);
        //Area A
        BigDecimal area = new BigDecimal("1");
        //poisson ratio v
        BigDecimal poisson = new BigDecimal("0.3");
        
        //(1-v)
        BigDecimal poisson1mv = poisson.negate().add(new BigDecimal("1"));
        
        //(1-2v)/2
        BigDecimal poisson1m2vo2 = poisson.multiply(new BigDecimal("2")).negate().add(new BigDecimal("1")).divide(new BigDecimal("2"));
        
        //(1-2v)
        BigDecimal poisson1m2v = poisson.multiply(new BigDecimal("2")).negate().add(new BigDecimal("1"));

        //(1+v)
        BigDecimal poisson1pv = poisson.add(new BigDecimal("1"));
        
        //variables
        //acceleration
        BigDecimal[] acceleration = new BigDecimal[100];

        int numX = 5;
        int numY = 3;

        List<Cst> elements = new ArrayList<>();
        
        //Stiffness matrix
        BigDecimal[][] stiffnessMatrix = new BigDecimal[2*numX*numY+2*numX+2*numY+2][2*numX*numY+2*numX+2*numY+2];
        
        //mass matrix
        BigDecimal[][] massMatrix = new BigDecimal[2*numX*numY+2*numX+2*numY+2][2*numX*numY+2*numX+2*numY+2];

        for (BigDecimal[] row : stiffnessMatrix) {
            Arrays.fill(row, BigDecimal.ZERO);
        }

        for (BigDecimal[] row : massMatrix) {
            Arrays.fill(row, BigDecimal.ZERO);
        }

        //points matrix
        for (int i = 0; i < numX; i++) {
            for (int j = 0; j < numY; j++) {
                Cst element = new Cst(new Point(new BigDecimal(i*10), new BigDecimal(j*10)), 
                        new Point(new BigDecimal(i*10+10), new BigDecimal(j*10+10)), 
                        new Point(new BigDecimal(i*10), new BigDecimal(j*10+10)));
                elements.add(element);
                
                //element 1
                //
                // i
                // |\
                // | \
                // |__\
                // m   j
                
                //beta i: y_j - y_m
                BigDecimal betaI = new BigDecimal("0");
                //beta j: y_m - y_i
                BigDecimal betaJ = new BigDecimal("10");
                //beta m: y_i - y_j
                BigDecimal betaM = new BigDecimal("-10");
                
                //gamma_i: x_m - x_j
                BigDecimal gammaI = new BigDecimal("-10");
                //gamma_j: x_i - x_m
                BigDecimal gammaJ = new BigDecimal("0");
                //gamma_m: x_j - x_i
                BigDecimal gammaM = new BigDecimal("10");
                
                //1-1
                BigDecimal m11 = betaI.multiply(betaI).multiply(poisson1mv).add(gammaI.multiply(gammaI).multiply(poisson1m2vo2));
                //1-2
                BigDecimal m12 = betaI.multiply(betaI).multiply(poisson).add(betaI.multiply(gammaI).multiply(poisson1m2vo2));
                //1-3
                BigDecimal m13 = betaI.multiply(betaJ).multiply(poisson1mv).add(gammaI.multiply(gammaJ).multiply(poisson1m2vo2));
                //1-4
                BigDecimal m14 = betaI.multiply(gammaJ).multiply(poisson).add(betaJ.multiply(gammaI).multiply(poisson1m2vo2));
                //1-5
                BigDecimal m15 = betaI.multiply(betaM).multiply(poisson1mv).add(gammaI.multiply(gammaM).multiply(poisson1m2vo2));
                //1-6
                BigDecimal m16 = betaI.multiply(gammaM).multiply(poisson).add(betaM.multiply(gammaI).multiply(poisson1m2vo2));
                
                //2-1
                BigDecimal m21 = m12;
                //2-2
                BigDecimal m22 = gammaI.multiply(gammaI).multiply(poisson1mv).add(betaI.multiply(betaI).multiply(poisson1m2vo2));
                //2-3
                BigDecimal m23 = betaJ.multiply(gammaI).multiply(poisson).add(betaI.multiply(gammaJ).multiply(poisson1m2vo2));
                //2-4
                BigDecimal m24 = gammaI.multiply(gammaJ).multiply(poisson1mv).add(betaI.multiply(betaJ).multiply(poisson1m2vo2));
                //2-5
                BigDecimal m25 = betaM.multiply(gammaI).multiply(poisson).add(betaI.multiply(gammaM).multiply(poisson1m2vo2));
                //2-6
                BigDecimal m26 = gammaI.multiply(gammaM).multiply(poisson1mv).add(betaI.multiply(betaM).multiply(poisson1m2vo2));
                
                //3-1
                BigDecimal m31 = m13;
                //3-2
                BigDecimal m32 = m23;
                //3-3
                BigDecimal m33 = betaJ.multiply(betaJ).multiply(poisson1mv).add(gammaJ.multiply(gammaJ).multiply(poisson1m2vo2));
                //3-4
                BigDecimal m34 = betaJ.multiply(gammaJ).multiply(poisson).add(betaJ.multiply(gammaJ).multiply(poisson1m2vo2));
                //3-5
                BigDecimal m35 = betaJ.multiply(betaM).multiply(poisson1mv).add(gammaJ.multiply(gammaM).multiply(poisson1m2vo2));
                //3-6
                BigDecimal m36 = betaJ.multiply(gammaM).multiply(poisson).add(betaM.multiply(gammaJ).multiply(poisson1m2vo2));
                
                //4-1
                BigDecimal m41 = m14;
                //4-2
                BigDecimal m42 = m24;
                //4-3
                BigDecimal m43 = m34;
                //4-4
                BigDecimal m44 = gammaJ.multiply(gammaJ).multiply(poisson1mv).add(betaJ.multiply(betaJ).multiply(poisson1m2vo2));
                //4-5
                BigDecimal m45 = betaM.multiply(gammaJ).multiply(poisson).add(betaJ.multiply(gammaM).multiply(poisson1m2vo2));
                //4-6
                BigDecimal m46 = gammaI.multiply(gammaM).multiply(poisson1mv).add(betaJ.multiply(betaM).multiply(poisson1m2vo2));
                
                //5-1
                BigDecimal m51 = m15;
                //5-2
                BigDecimal m52 = m25;
                //5-3
                BigDecimal m53 = m35;
                //5-4
                BigDecimal m54 = m45;
                //5-5
                BigDecimal m55 = betaM.multiply(betaM).multiply(poisson1mv).add(gammaM.multiply(gammaM).multiply(poisson1m2vo2));
                //5-6
                BigDecimal m56 = betaM.multiply(gammaM).multiply(poisson).add(betaM.multiply(gammaM).multiply(poisson1m2vo2));
                
                //6-1
                BigDecimal m61 = m16;
                //6-2
                BigDecimal m62 = m26;
                //6-3
                BigDecimal m63 = m36;
                //6-4
                BigDecimal m64 = m46;
                //6-5
                BigDecimal m65 = m56;
                //6-6
                BigDecimal m66 = gammaM.multiply(gammaM).multiply(poisson1mv).add(betaM.multiply(betaM).multiply(poisson1m2vo2));

                System.out.println("==== 1 elem ====");
                System.out.println(m11 + "\t" + m12 + "\t" + m13 + "\t" + m14 + "\t" + m15 + "\t" + m16);
                System.out.println(m21 + "\t" + m22 + "\t" + m23 + "\t" + m24 + "\t" + m25 + "\t" + m26);
                System.out.println(m31 + "\t" + m32 + "\t" + m33 + "\t" + m34 + "\t" + m35 + "\t" + m36);
                System.out.println(m41 + "\t" + m42 + "\t" + m43 + "\t" + m44 + "\t" + m45 + "\t" + m46);
                System.out.println(m51 + "\t" + m52 + "\t" + m53 + "\t" + m54 + "\t" + m55 + "\t" + m56);
                System.out.println(m61 + "\t" + m62 + "\t" + m63 + "\t" + m64 + "\t" + m65 + "\t" + m66);
                

                //top left corner
                stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*j+i)*2] = stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*j+i)*2].add(m11);
                stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*j+i)*2+1] = stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*j+i)*2+1].add(m12);
                stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*j+i)*2] = stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*j+i)*2].add(m21);
                stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*j+i)*2+1] = stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*j+i)*2+1].add(m22);

                //bottom left corner
                stiffnessMatrix[((numX+1)*(j+1)+i)*2][((numX+1)*(j+1)+i)*2] = stiffnessMatrix[((numX+1)*(j+1)+i)*2][((numX+1)*(j+1)+i)*2].add(m55);
                stiffnessMatrix[((numX+1)*(j+1)+i)*2][((numX+1)*(j+1)+i)*2+1] = stiffnessMatrix[((numX+1)*(j+1)+i)*2][((numX+1)*(j+1)+i)*2+1].add(m56);
                stiffnessMatrix[((numX+1)*(j+1)+i)*2+1][((numX+1)*(j+1)+i)*2] = stiffnessMatrix[((numX+1)*(j+1)+i)*2+1][((numX+1)*(j+1)+i)*2].add(m65);
                stiffnessMatrix[((numX+1)*(j+1)+i)*2+1][((numX+1)*(j+1)+i)*2+1] = stiffnessMatrix[((numX+1)*(j+1)+i)*2+1][((numX+1)*(j+1)+i)*2+1].add(m66);

                //bottom right corner
                stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2][((numX+1)*(j+1)+(i+1))*2] = stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2][((numX+1)*(j+1)+(i+1))*2].add(m33);
                stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2][((numX+1)*(j+1)+(i+1))*2+1] = stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2][((numX+1)*(j+1)+(i+1))*2+1].add(m34);
                stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2+1][((numX+1)*(j+1)+(i+1))*2] = stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2+1][((numX+1)*(j+1)+(i+1))*2].add(m43);
                stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2+1][((numX+1)*(j+1)+(i+1))*2+1] = stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2+1][((numX+1)*(j+1)+(i+1))*2+1].add(m44);

                //ij
                stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*(j+1)+(i+1))*2] = stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*(j+1)+(i+1))*2].add(m13);
                stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*(j+1)+(i+1))*2+1] = stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*(j+1)+(i+1))*2+1].add(m14);
                stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*(j+1)+(i+1))*2] = stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*(j+1)+(i+1))*2].add(m23);
                stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*(j+1)+(i+1))*2+1] = stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*(j+1)+(i+1))*2+1].add(m24);

                stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2][((numX+1)*j+i)*2] = stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*(j+1)+(i+1))*2];
                stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2+1][((numX+1)*j+i)*2] = stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*(j+1)+(i+1))*2+1];
                stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2][((numX+1)*j+i)*2+1] = stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*(j+1)+(i+1))*2];
                stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2+1][((numX+1)*j+i)*2+1] = stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*(j+1)+(i+1))*2+1];

                //im
                stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*(j+1)+i)*2] = stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*(j+1)+i)*2].add(m15);
                stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*(j+1)+i)*2+1] = stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*(j+1)+i)*2+1].add(m16);
                stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*(j+1)+i)*2] = stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*(j+1)+i)*2].add(m25);
                stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*(j+1)+i)*2+1] = stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*(j+1)+i)*2+1].add(m26);

                stiffnessMatrix[((numX+1)*(j+1)+i)*2][((numX+1)*j+i)*2] = stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*(j+1)+i)*2];
                stiffnessMatrix[((numX+1)*(j+1)+i)*2+1][((numX+1)*j+i)*2] = stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*(j+1)+i)*2+1];
                stiffnessMatrix[((numX+1)*(j+1)+i)*2][((numX+1)*j+i)*2+1] = stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*(j+1)+i)*2];
                stiffnessMatrix[((numX+1)*(j+1)+i)*2+1][((numX+1)*j+i)*2+1] = stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*(j+1)+i)*2+1];
                
                //jm
                stiffnessMatrix[((numX+1)*(j+1)+i)*2][((numX+1)*(j+1)+(i+1))*2] = stiffnessMatrix[((numX+1)*(j+1)+i)*2][((numX+1)*(j+1)+(i+1))*2].add(m35);
                stiffnessMatrix[((numX+1)*(j+1)+i)*2+1][((numX+1)*(j+1)+(i+1))*2] = stiffnessMatrix[((numX+1)*(j+1)+i)*2+1][((numX+1)*(j+1)+(i+1))*2].add(m45);
                stiffnessMatrix[((numX+1)*(j+1)+i)*2][((numX+1)*(j+1)+(i+1))*2+1] = stiffnessMatrix[((numX+1)*(j+1)+i)*2][((numX+1)*(j+1)+(i+1))*2+1].add(m36);
                stiffnessMatrix[((numX+1)*(j+1)+i)*2+1][((numX+1)*(j+1)+(i+1))*2+1] = stiffnessMatrix[((numX+1)*(j+1)+i)*2+1][((numX+1)*(j+1)+(i+1))*2+1].add(m46);

                stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2][((numX+1)*(j+1)+i)*2] = stiffnessMatrix[((numX+1)*(j+1)+i)*2][((numX+1)*(j+1)+(i+1))*2];
                stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2][((numX+1)*(j+1)+i)*2+1] = stiffnessMatrix[((numX+1)*(j+1)+i)*2+1][((numX+1)*(j+1)+(i+1))*2];
                stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2+1][((numX+1)*(j+1)+i)*2] = stiffnessMatrix[((numX+1)*(j+1)+i)*2][((numX+1)*(j+1)+(i+1))*2+1];
                stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2+1][((numX+1)*(j+1)+i)*2+1] = stiffnessMatrix[((numX+1)*(j+1)+i)*2+1][((numX+1)*(j+1)+(i+1))*2+1];

//                stiffnessMatrix[i*j*2+i*2+j*2][i*j*2+i*2+j*2] = stiffnessMatrix[i*j*2+i*2+j*2][i*j*2+i*2+j*2].add(m55);
//                stiffnessMatrix[i*j*2+i*2+j*2][i*j*2+i*2+j*2+1] = stiffnessMatrix[i*j*2+i*2+j*2][i*j*2+i*2+j*2+1].add(m56);
//                stiffnessMatrix[i*j*2+i*2+j*2+1][i*j*2+i*2+j*2] = stiffnessMatrix[i*j*2+i*2+j*2+1][i*j*2+i*2+j*2].add(m65);
//                stiffnessMatrix[i*j*2+i*2+j*2+1][i*j*2+i*2+j*2+1] = stiffnessMatrix[i*j*2+i*2+j*2+1][i*j*2+i*2+j*2+1].add(m66);
//                
//                stiffnessMatrix[i*2][j*2+2] = stiffnessMatrix[i*2][j*2+2].add(m53);
//                stiffnessMatrix[i*2][j*2+3] = stiffnessMatrix[i*2][j*2+3].add(m54);
//                stiffnessMatrix[i*2+1][j*2+2] = stiffnessMatrix[i*2+1][j*2+2].add(m63);
//                stiffnessMatrix[i*2+1][j*2+3] = stiffnessMatrix[i*2+1][j*2+3].add(m64);
//                
//                stiffnessMatrix[i*2][j*2+4] = stiffnessMatrix[i*2][j*2+4].add(m51);
//                stiffnessMatrix[i*2][j*2+5] = stiffnessMatrix[i*2][j*2+5].add(m52);
//                stiffnessMatrix[i*2+1][j*2+4] = stiffnessMatrix[i*2+1][j*2+4].add(m61);
//                stiffnessMatrix[i*2+1][j*2+5] = stiffnessMatrix[i*2+1][j*2+5].add(m62);
//
//                stiffnessMatrix[i*2+2][j*2] = stiffnessMatrix[i*2+2][j*2].add(m35);
//                stiffnessMatrix[i*2+2][j*2+1] = stiffnessMatrix[i*2+2][j*2+1].add(m36);
//                stiffnessMatrix[i*2+3][j*2] = stiffnessMatrix[i*2+3][j*2].add(m45);
//                stiffnessMatrix[i*2+3][j*2+1] = stiffnessMatrix[i*2+3][j*2+1].add(m46);
//                
//                stiffnessMatrix[i*2+2][j*2+2] = stiffnessMatrix[i*2+2][j*2+2].add(m33);
//                stiffnessMatrix[i*2+2][j*2+3] = stiffnessMatrix[i*2+2][j*2+3].add(m34);
//                stiffnessMatrix[i*2+3][j*2+2] = stiffnessMatrix[i*2+3][j*2+2].add(m43);
//                stiffnessMatrix[i*2+3][j*2+3] = stiffnessMatrix[i*2+3][j*2+3].add(m44);
//                
//                stiffnessMatrix[i*2+2][j*2+4] = stiffnessMatrix[i*2+2][j*2+4].add(m31);
//                stiffnessMatrix[i*2+2][j*2+5] = stiffnessMatrix[i*2+2][j*2+5].add(m32);
//                stiffnessMatrix[i*2+3][j*2+4] = stiffnessMatrix[i*2+3][j*2+4].add(m41);
//                stiffnessMatrix[i*2+3][j*2+5] = stiffnessMatrix[i*2+3][j*2+5].add(m42);
//
//                stiffnessMatrix[i*2+4][j*2]   = stiffnessMatrix[i*2+4][j*2].add(m15);
//                stiffnessMatrix[i*2+4][j*2+1] = stiffnessMatrix[i*2+4][j*2+1].add(m16);
//                stiffnessMatrix[i*2+5][j*2]   = stiffnessMatrix[i*2+5][j*2].add(m25);
//                stiffnessMatrix[i*2+5][j*2+1] = stiffnessMatrix[i*2+5][j*2+1].add(m26);
//                
//                stiffnessMatrix[i*2+4][j*2+2] = stiffnessMatrix[i*2+4][j*2+2].add(m13);
//                stiffnessMatrix[i*2+4][j*2+3] = stiffnessMatrix[i*2+4][j*2+3].add(m14);
//                stiffnessMatrix[i*2+5][j*2+2] = stiffnessMatrix[i*2+5][j*2+2].add(m23);
//                stiffnessMatrix[i*2+5][j*2+3] = stiffnessMatrix[i*2+5][j*2+3].add(m24);
//                
//                stiffnessMatrix[i*2+4][j*2+4] = stiffnessMatrix[i*2+4][j*2+4].add(m11);
//                stiffnessMatrix[i*2+4][j*2+5] = stiffnessMatrix[i*2+4][j*2+5].add(m12);
//                stiffnessMatrix[i*2+5][j*2+4] = stiffnessMatrix[i*2+5][j*2+4].add(m21);
//                stiffnessMatrix[i*2+5][j*2+5] = stiffnessMatrix[i*2+5][j*2+5].add(m22);
                
                Cst element2 = new Cst(new Point(new BigDecimal(i*10), new BigDecimal(j*10)), 
                        new Point(new BigDecimal(i*10+10), new BigDecimal(j*10)), 
                        new Point(new BigDecimal(i*10+10), new BigDecimal(j*10+10)));
                elements.add(element2);
                
                //element 2
                //
                // i___j
                //  \  |
                //   \ |
                //    \|
                //     m
                
                //beta i: y_j - y_m
                betaI = new BigDecimal("-10");
                //beta j: y_m - y_i
                betaJ = new BigDecimal("10");
                //beta m: y_i - y_j
                betaM = new BigDecimal("0");
                
                //gamma_i: x_m - x_j
                gammaI = new BigDecimal("0");
                //gamma_j: x_i - x_m
                gammaJ = new BigDecimal("-10");
                //gamma_m: x_j - x_i
                gammaM = new BigDecimal("10");
                
                //1-1
                m11 = betaI.multiply(betaI).multiply(poisson1mv).add(gammaI.multiply(gammaI).multiply(poisson1m2vo2));
                //1-2
                m12 = betaI.multiply(betaI).multiply(poisson).add(betaI.multiply(gammaI).multiply(poisson1m2vo2));
                //1-3
                m13 = betaI.multiply(betaJ).multiply(poisson1mv).add(gammaI.multiply(gammaJ).multiply(poisson1m2vo2));
                //1-4
                m14 = betaI.multiply(gammaJ).multiply(poisson).add(betaJ.multiply(gammaI).multiply(poisson1m2vo2));
                //1-5
                m15 = betaI.multiply(betaM).multiply(poisson1mv).add(gammaI.multiply(gammaM).multiply(poisson1m2vo2));
                //1-6
                m16 = betaI.multiply(gammaM).multiply(poisson).add(betaM.multiply(gammaI).multiply(poisson1m2vo2));
                
                //2-1
                m21 = m12;
                //2-2
                m22 = gammaI.multiply(gammaI).multiply(poisson1mv).add(betaI.multiply(betaI).multiply(poisson1m2vo2));
                //2-3
                m23 = betaJ.multiply(gammaI).multiply(poisson).add(betaI.multiply(gammaJ).multiply(poisson1m2vo2));
                //2-4
                m24 = gammaI.multiply(gammaJ).multiply(poisson1mv).add(betaI.multiply(betaJ).multiply(poisson1m2vo2));
                //2-5
                m25 = betaM.multiply(gammaI).multiply(poisson).add(betaI.multiply(gammaM).multiply(poisson1m2vo2));
                //2-6
                m26 = gammaI.multiply(gammaM).multiply(poisson1mv).add(betaI.multiply(betaM).multiply(poisson1m2vo2));
                
                //3-1
                m31 = m13;
                //3-2
                m32 = m23;
                //3-3
                m33 = betaJ.multiply(betaJ).multiply(poisson1mv).add(gammaJ.multiply(gammaJ).multiply(poisson1m2vo2));
                //3-4
                m34 = betaJ.multiply(gammaJ).multiply(poisson).add(betaJ.multiply(gammaJ).multiply(poisson1m2vo2));
                //3-5
                m35 = betaJ.multiply(betaM).multiply(poisson1mv).add(gammaJ.multiply(gammaM).multiply(poisson1m2vo2));
                //3-6
                m36 = betaJ.multiply(gammaM).multiply(poisson).add(betaM.multiply(gammaJ).multiply(poisson1m2vo2));
                
                //4-1
                m41 = m14;
                //4-2
                m42 = m24;
                //4-3
                m43 = m34;
                //4-4
                m44 = gammaJ.multiply(gammaJ).multiply(poisson1mv).add(betaJ.multiply(betaJ).multiply(poisson1m2vo2));
                //4-5
                m45 = betaM.multiply(gammaJ).multiply(poisson).add(betaJ.multiply(gammaM).multiply(poisson1m2vo2));
                //4-6
                m46 = gammaI.multiply(gammaM).multiply(poisson1mv).add(betaJ.multiply(betaM).multiply(poisson1m2vo2));
                
                //5-1
                m51 = m15;
                //5-2
                m52 = m25;
                //5-3
                m53 = m35;
                //5-4
                m54 = m45;
                //5-5
                m55 = betaM.multiply(betaM).multiply(poisson1mv).add(gammaM.multiply(gammaM).multiply(poisson1m2vo2));
                //5-6
                m56 = betaM.multiply(gammaM).multiply(poisson).add(betaM.multiply(gammaM).multiply(poisson1m2vo2));
                
                //6-1
                m61 = m16;
                //6-2
                m62 = m26;
                //6-3
                m63 = m36;
                //6-4
                m64 = m46;
                //6-5
                m65 = m56;
                //6-6
                m66 = gammaM.multiply(gammaM).multiply(poisson1mv).add(betaM.multiply(betaM).multiply(poisson1m2vo2));

                //top left corner
                stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*j+i)*2] = stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*j+i)*2].add(m11);
                stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*j+i)*2+1] = stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*j+i)*2+1].add(m12);
                stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*j+i)*2] = stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*j+i)*2].add(m21);
                stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*j+i)*2+1] = stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*j+i)*2+1].add(m22);

                //top right corner
                stiffnessMatrix[((numX+1)*j+(i+1))*2][((numX+1)*j+(i+1))*2] = stiffnessMatrix[((numX+1)*j+(i+1))*2][((numX+1)*j+(i+1))*2].add(m33);
                stiffnessMatrix[((numX+1)*j+(i+1))*2][((numX+1)*j+(i+1))*2+1] = stiffnessMatrix[((numX+1)*j+(i+1))*2][((numX+1)*j+(i+1))*2+1].add(m34);
                stiffnessMatrix[((numX+1)*j+(i+1))*2+1][((numX+1)*j+(i+1))*2] = stiffnessMatrix[((numX+1)*j+(i+1))*2+1][((numX+1)*j+(i+1))*2].add(m43);
                stiffnessMatrix[((numX+1)*j+(i+1))*2+1][((numX+1)*j+(i+1))*2+1] = stiffnessMatrix[((numX+1)*j+(i+1))*2+1][((numX+1)*j+(i+1))*2+1].add(m44);

                //bottom right corner
                stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2][((numX+1)*(j+1)+(i+1))*2] = stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2][((numX+1)*(j+1)+(i+1))*2].add(m55);
                stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2][((numX+1)*(j+1)+(i+1))*2+1] = stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2][((numX+1)*(j+1)+(i+1))*2+1].add(m56);
                stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2+1][((numX+1)*(j+1)+(i+1))*2] = stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2+1][((numX+1)*(j+1)+(i+1))*2].add(m65);
                stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2+1][((numX+1)*(j+1)+(i+1))*2+1] = stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2+1][((numX+1)*(j+1)+(i+1))*2+1].add(m66);
                
                //ij
                stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*j+(i+1))*2] = stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*j+(i+1))*2].add(m13);
                stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*j+(i+1))*2+1] = stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*j+(i+1))*2+1].add(m14);
                stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*j+(i+1))*2] = stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*j+(i+1))*2].add(m23);
                stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*j+(i+1))*2+1] = stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*j+(i+1))*2+1].add(m24);
                

                stiffnessMatrix[((numX+1)*j+(i+1))*2][((numX+1)*j+i)*2] = stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*j+(i+1))*2];
                stiffnessMatrix[((numX+1)*j+(i+1))*2+1][((numX+1)*j+i)*2] = stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*j+(i+1))*2+1];
                stiffnessMatrix[((numX+1)*j+(i+1))*2][((numX+1)*j+i)*2+1] = stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*j+(i+1))*2];
                stiffnessMatrix[((numX+1)*j+(i+1))*2+1][((numX+1)*j+i)*2+1] = stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*j+(i+1))*2+1];
                
                //im
                stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*(j+1)+(i+1))*2] = stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*(j+1)+(i+1))*2].add(m15);
                stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*(j+1)+(i+1))*2+1] = stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*(j+1)+(i+1))*2+1].add(m16);
                stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*(j+1)+(i+1))*2] = stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*(j+1)+(i+1))*2].add(m25);
                stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*(j+1)+(i+1))*2+1] = stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*(j+1)+(i+1))*2+1].add(m26);

                stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2][((numX+1)*j+i)*2] = stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*(j+1)+(i+1))*2];
                stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2+1][((numX+1)*j+i)*2] = stiffnessMatrix[((numX+1)*j+i)*2][((numX+1)*(j+1)+(i+1))*2+1];
                stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2][((numX+1)*j+i)*2+1] = stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*(j+1)+(i+1))*2];
                stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2+1][((numX+1)*j+i)*2+1] = stiffnessMatrix[((numX+1)*j+i)*2+1][((numX+1)*(j+1)+(i+1))*2+1];
                
                //jm
                stiffnessMatrix[((numX+1)*j+(i+1))*2][((numX+1)*(j+1)+(i+1))*2] = stiffnessMatrix[((numX+1)*j+(i+1))*2][((numX+1)*(j+1)+(i+1))*2].add(m35);
                stiffnessMatrix[((numX+1)*j+(i+1))*2][((numX+1)*(j+1)+(i+1))*2+1] = stiffnessMatrix[((numX+1)*j+(i+1))*2][((numX+1)*(j+1)+(i+1))*2+1].add(m36);
                stiffnessMatrix[((numX+1)*j+(i+1))*2+1][((numX+1)*(j+1)+(i+1))*2] = stiffnessMatrix[((numX+1)*j+(i+1))*2+1][((numX+1)*(j+1)+(i+1))*2].add(m45);
                stiffnessMatrix[((numX+1)*j+(i+1))*2+1][((numX+1)*(j+1)+(i+1))*2+1] = stiffnessMatrix[((numX+1)*j+(i+1))*2+1][((numX+1)*(j+1)+(i+1))*2+1].add(m46);

                stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2][((numX+1)*j+(i+1))*2] = stiffnessMatrix[((numX+1)*j+(i+1))*2][((numX+1)*(j+1)+(i+1))*2];
                stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2+1][((numX+1)*j+(i+1))*2] = stiffnessMatrix[((numX+1)*j+(i+1))*2][((numX+1)*(j+1)+(i+1))*2+1];
                stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2][((numX+1)*j+(i+1))*2+1] = stiffnessMatrix[((numX+1)*j+(i+1))*2+1][((numX+1)*(j+1)+(i+1))*2];
                stiffnessMatrix[((numX+1)*(j+1)+(i+1))*2+1][((numX+1)*j+(i+1))*2+1] = stiffnessMatrix[((numX+1)*j+(i+1))*2+1][((numX+1)*(j+1)+(i+1))*2+1];

                System.out.println("==== 2 elem ====");
                System.out.println(m11 + " " + m12 + " " + m13 + " " + m14 + " " + m15 + " " + m16);
                System.out.println(m21 + " " + m22 + " " + m23 + " " + m24 + " " + m25 + " " + m26);
                System.out.println(m31 + " " + m32 + " " + m33 + " " + m34 + " " + m35 + " " + m36);
                System.out.println(m41 + " " + m42 + " " + m43 + " " + m44 + " " + m45 + " " + m46);
                System.out.println(m51 + " " + m52 + " " + m53 + " " + m54 + " " + m55 + " " + m56);
                System.out.println(m61 + " " + m62 + " " + m63 + " " + m64 + " " + m65 + " " + m66);

//                BigDecimal constant = density.multiply(area).multiply(thickness);
                BigDecimal constant = new BigDecimal("1");
                BigDecimal two = new BigDecimal(2).multiply(constant);
                BigDecimal zero = new BigDecimal(0);
                BigDecimal one = constant;

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
        System.out.println("==== Stiffness ====");
        for (BigDecimal[] column : stiffnessMatrix) {
            for (BigDecimal unit : column) {
                System.out.print(unit + "\t");
            }
            System.out.println();
        }
        
        System.out.println("==== mass ====");
        for (BigDecimal[] column : massMatrix) {
            for (BigDecimal unit : column) {
                System.out.print(unit + "\t");
            }
            System.out.println();
        }

        printElements(elements);
        
        System.out.println(elasticity);
    }

    private static void printElements(List<Cst> elements) {
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

        while (running) {
            bufferStrategy = canvas.getBufferStrategy();
            graphics = bufferStrategy.getDrawGraphics();
            graphics.clearRect(0, 0, width, height);

            graphics.setColor(Color.GREEN);

            for (Cst element : elements) {
                graphics.drawLine(element.getI().x.intValue(), element.getI().y.intValue()
                        , element.getJ().x.intValue(), element.getJ().y.intValue());
                graphics.drawLine(element.getI().x.intValue(), element.getI().y.intValue()
                        , element.getM().x.intValue(), element.getM().y.intValue());
                graphics.drawLine(element.getM().x.intValue(), element.getM().y.intValue()
                        , element.getJ().x.intValue(), element.getJ().y.intValue());
            }

            bufferStrategy.show();
            graphics.dispose();
            try {
                Thread.sleep(1000L);
            } catch (InterruptedException e) {
                // TODO Auto-generated catch block
                e.printStackTrace();
            }
        }
    }
}