package com.spbpu;

import static com.spbpu.VectorOperations.VectorVectorMult;

/**
 * Set of functions and their gradients
 */
public class FunctionSet {

    public static double f(double x[])
    {
        return x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
    }
    public static double[] gradf(double[] x)
    {
        return new double[]{2*x[0], 2*x[1], 2*x[2]};
    }
    /*public static double f(double x[])
    {
        return Ec[0] * x[0] + Ec[1] * x[1];
    }
    public static double[] gradf(double[] x)
    {
        return new double[]{Ec[0], Ec[1]};
    }

    // PARAMETERS OF THE TASK
    private static double[] Ec = new double[]{6, 8};
    private static double[][] d = new double[][]{{3, 5}, {1, 6}};  *//*matrix of borders of allowed values of x
                                                                    (example: d[1][1] <= x[1] <= d[1][2]) *//*
    private static double[] p = new double[]{0.98, 0.98}; // vector of probability levels (alpha_i)
    private static double[][] EA = new double[][]{{12, 14}, {16, 12}}; // expectancy of coefficient matrix A
    private static double[][] DA = new double[][]{{3, 4}, {5, 4}}; // dispersion of coefficient matrix A
    private static double[] Eb = new double[]{140, 160}; // expectancy of vector b
    private static double[] Db = new double[]{8, 10}; // dispersion of vector b
    private static double r = 2; // coefficient of penalty function*/

    /*private static double[] Ec = new double[]{5, 8};
    private static double[][] d = new double[][]{{2, 6}, {3, 9}}; *//* matrix of borders of allowed values of x
                                                                    (example: d[1][1] <= x[1] <= d[1][2]) *//*
    private static double[] p = new double[]{0.994, 0.994}; // vector of probability levels (alpha_i)
    private static double[][] EA = new double[][]{{-10, -15}, {-20, -14}}; // expectancy of coefficient matrix A
    private static double[][] DA = new double[][]{{2, 3}, {6, 4}}; // dispersion of coefficient matrix A
    private static double[] Eb = new double[]{-50, -90}; // expectancy of vector b
    private static double[] Db = new double[]{9, 12}; // dispersion of vector b*/
    private static double r = 2; // coefficient near penalty function

    /**
     * Multiply r by 2
     */
    public static void IncreaseR() {
        r *= 2;
    }

    /**
     * Calculate values of limit functions in provided point
     * @param x - point
     * @return vector with values of limit functions
     */
    public static double[] LimitFunctions(double[] x) {
        double g1 = x[0]*x[0] + x[1]*x[1] - 20;
        double g2 = x[0] * x[1] - 10;
        double g3 = x[0] * x[2] - 8;
        double g4 = x[0] + x[1] + 3;
        double g5 = -g4;
        double g6 = x[1] + x[2] + 4;
        double g7 = -g6;
        return new double[]{g1, g2, g3, g4, g5, g6, g7};
    }
    /*public static double[] LimitFunctions(double[] x) {
        double t1 = ReverseNormalFunction(p[0]);
        double t2 = ReverseNormalFunction(p[1]);
        double g1 = VectorVectorMult(EA[0], x) - Eb[0] +
                t1 * Math.sqrt(Math.pow(DA[0][0] * x[0],2) + Math.pow(DA[0][1] * x[1],2) + Math.pow(Db[0],2));
        double g2 = VectorVectorMult(EA[1], x) - Eb[1] +
                t2 * Math.sqrt(Math.pow(DA[1][0] * x[0],2) + Math.pow(DA[1][1] * x[1],2) + Math.pow(Db[1],2));
        double g3 = d[0][0] - x[0], g4 = x[0] - d[0][1];
        double g5 = d[1][0] - x[1], g6 = x[1] - d[1][1];
        return new double[]{g1, g2, g3, g4, g5, g6};
    }*/

    /**
     * Calculate value of new objective function in provided point
     * @param x - point
     * @return value of new objective function
     */
    public static double NewFunction(double[] x) {
        double[] G = LimitFunctions(x);
        double limitFunctionsSum  = 0;
        for (double GItem: G) {
            if (GItem >= 0)
                limitFunctionsSum += GItem * GItem;
        }
        return f(x) + r * limitFunctionsSum;
    }

    public static double[] NewGradient(double[] x) {
        double[] gradFunction = gradf(x);
        double grad0 = gradFunction[0];
        double grad1 = gradFunction[1];
        double grad2 = gradFunction[2];
        double[] G = LimitFunctions(x);
        if(G[0] > 0) {
            grad0 += 2*x[0] * r * 2*G[0];
            grad1 += 2*x[1] * r * 2*G[0];
        }
        if (G[1] > 0) {
            grad0 += x[1] * r * 2*G[1];
            grad1 += x[0] * r * 2*G[1];
        }
        if (G[2] > 0) {
            grad0 += x[2] * r * 2*G[2];
            grad2 += x[0] * r * 2*G[2];
        }
        if (G[3] > 0) {
            grad0 += 1 * r * 2*G[3];
            grad1 += 1 * r * 2*G[3];
        }
        if (G[4] > 0) {
            grad0 += -1 * r * 2*G[4];
            grad1 += -1 * r * 2*G[4];
        }
        if (G[5] > 0) {
            grad1 += 1 * r * 2*G[5];
            grad2 += 1 * r * 2*G[5];
        }
        /*if (G[6] > 0) {
            grad0 += -r * 2*G[6];
            grad2 += -r * 2*G[6];
        }*/
        return new double[]{grad0, grad1, grad2};
    }
    /*public static double[] NewGradient(double[] x) {
        double[] gradfunction = gradf(x);
        double grad0 = 0;
        double grad1 = 0;
        double[] G = LimitFunctions(x);
        if(G[0] >= 0) {
            grad0 += r * 2 * G[0] * EA[0][0] + ReverseNormalFunction(p[0]) * Math.pow(DA[0][0], 2) * x[0] /
                    Math.sqrt(Math.pow(DA[0][0] * x[0], 2) + Math.pow(DA[0][1] * x[1], 2) + Math.pow(Db[0], 2));
            grad1 += r * 2 * G[0] * EA[0][1] + ReverseNormalFunction(p[0]) * Math.pow(DA[0][1], 2) * x[1] /
                    Math.sqrt(Math.pow(DA[0][0] * x[0], 2) + Math.pow(DA[0][1] * x[1], 2) + Math.pow(Db[0], 2));
        }
        if (G[1] >= 0) {
            grad0 += r * 2 * G[1] * EA[1][0] + ReverseNormalFunction(p[1]) * Math.pow(DA[1][0], 2) * x[0] /
                    Math.sqrt(Math.pow(DA[1][0] * x[0], 2) + Math.pow(DA[1][1] * x[1], 2) + Math.pow(Db[1], 2));
            grad1 += r * 2 * G[1] * EA[1][1] + ReverseNormalFunction(p[1]) * Math.pow(DA[1][1], 2) * x[1] /
                    Math.sqrt(Math.pow(DA[1][0] * x[0], 2) + Math.pow(DA[1][1] * x[1], 2) + Math.pow(Db[1], 2));
        }
        if (G[2] >= 0)
            grad0 += -r * 2 * G[2];
        if (G[3] >= 0)
            grad0 += r * 2 * G[3];
        if (G[4] >= 0)
            grad1 += -r * 2 * G[4];
        if (G[5] >= 0)
            grad1 += r * 2 * G[5];
        grad0 *= -1;
        grad1 *= -1;
        grad0 += gradfunction[0];
        grad1 += gradfunction[1];
        return new double[]{grad0, grad1};
    }*/

    /**
     * Calculate reverse normal function (normal distribution) from tables in point p
     * @param p - point
     * @return value of reverse normal function
     */
    private static double ReverseNormalFunction(double p) {
        double eps = 0.0000001;
        if (Math.abs(p - 0.5) < eps) return 0;
        if (Math.abs(p - 0.6) < eps) return 0.25;
        if (Math.abs(p - 0.7) < eps) return 0.5;
        if (Math.abs(p - 0.77) < eps) return 0.75;
        if (Math.abs(p - 0.84) < eps) return 1.0;
        if (Math.abs(p - 0.89) < eps) return 1.25;
        if (Math.abs(p - 0.93) < eps) return 1.5;
        if (Math.abs(p - 0.96) < eps) return 1.75;
        if (Math.abs(p - 0.98) < eps) return 2.0;
        if (Math.abs(p - 0.987) < eps) return 2.25;
        if (Math.abs(p - 0.994) < eps) return 2.5;
        return -1;
    }
}
