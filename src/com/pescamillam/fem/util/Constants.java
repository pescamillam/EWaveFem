package com.pescamillam.fem.util;

import java.text.DecimalFormat;

import org.apache.commons.math3.util.BigReal;

/**
 * Constants with initial values and default decimal format
 * 
 * @author Peter Escamilla (pescamilla@unab.edu.co)
 */
public class Constants {
    public static final BigReal FOUR = new BigReal("4");
    public static final BigReal TWELVE = new BigReal("12");
    public static final BigReal TWO = new BigReal("2");
    public static final BigReal MINUS_ONE = new BigReal("-1");

    public static final String THICKNESS = "1";
    public static final String ELASTICITY = "3000000";
    public static final String DENSITY = "0.00073";
    public static final String AREA = "450";
    public static final String POISSON = "0.3";
    public static final String DELTA_TIME = "0.00005";
    public static final String NUM_X = "8";
    public static final String NUM_Y = "8";
    public static final String NUM_TIMES = "200";

    public static DecimalFormat DF = new DecimalFormat();
    static {
        DF.setMaximumFractionDigits(10);
        DF.setMinimumFractionDigits(0);
        DF.setGroupingUsed(false);
    }
}
