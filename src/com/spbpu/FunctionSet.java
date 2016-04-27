package com.spbpu;

import static com.spbpu.VectorOperations.VectorVectorMult;

/**
 * Set of functions and their gradients
 */
public class FunctionSet {

    /*private static double a1 = 1, a2 = 2, a3 = 3;

    public static double g(double[] x)
    {
        return a1*x[0] + x[1] + 4*Math.sqrt(1 + a2*x[0]*x[0] + a3*x[1]*x[1]);
    }
    public static double[] gradg(double[] x)
    {
        double g1 = a1 + 4*a2*x[0] / (Math.sqrt(1 + a2*x[0]*x[0] + a3*x[1]*x[1]));
        double g2 = 1 + 4*a3*x[1] / (Math.sqrt(1 + a2*x[0]*x[0] + a3*x[1]*x[1]));
        return new double[]{g1, g2};
    }*/

    private static double[] Ec = new double[]{6, 8};

    public static double f(double x[])
    {
        return Ec[0] * x[0] + Ec[1] * x[1];
    }
    public static double[] gradf(double[] x)
    {
        return new double[]{Ec[0], Ec[1]};
    }

    // PARAMETERS OF THE TASK
    private static double[][] d = new double[][]{{3, 5}, {1, 6}}; /* matrix of borders of allowed values of x
                                                                    (example: d[1][1] <= x[1] <= d[1][2]) */
    private static double[] p = new double[]{0.5, 0.5}; // vector of probability levels (alpha_i)
    private static double[][] EA = new double[][]{{12, 14}, {16, 12}}; // expectancy of coefficient matrix A
    private static double[][] DA = new double[][]{{3, 4}, {5, 4}}; // dispersion of coefficient matrix A
    private static double[] Eb = new double[]{140, 160}; // expectancy of vector b
    private static double[] Db = new double[]{8, 10}; // dispersion of vector b
    private static double r = 1; // coefficient of penalty function

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
        double t1 = ReverseNormalFunction(p[0]);
        double t2 = ReverseNormalFunction(p[1]);
        double g1 = VectorVectorMult(EA[0], x) - Eb[0] +
                t1 * Math.sqrt(Math.pow(DA[0][0] * x[0],2) + Math.pow(DA[0][1] * x[1],2) + Math.pow(Db[0],2));
        double g2 = VectorVectorMult(EA[1], x) - Eb[1] +
                t2 * Math.sqrt(Math.pow(DA[1][0] * x[0],2) + Math.pow(DA[1][1] * x[1],2) + Math.pow(Db[1],2));
        double g3 = d[0][0] - x[0], g4 = x[0] - d[0][1];
        double g5 = d[1][0] - x[1], g6 = x[1] - d[1][1];
        return new double[]{g1, g2, g3, g4, g5, g6};
    }

    /**
     * Calculate value of new objective function in provided point
     * @param x - point
     * @return value of new objective function
     */
    public static double NewFunction(double[] x) {
        double[] G = LimitFunctions(x);
        double limitFunctionsSum  = 0;
        for (double GItem: G) {
            if (GItem > 0)
                limitFunctionsSum += GItem;
        }
        return f(x) + r * limitFunctionsSum;
    }

    public static double[] NewGradient(double[] x) {
        double[] gradfunction = gradf(x);
        double grad0 = gradfunction[0];
        double grad1 = gradfunction[1];
        double[] G = LimitFunctions(x);
        if(G[0] >= 0) {
            grad0 += r * EA[0][0] + ReverseNormalFunction(p[0]) * Math.pow(DA[0][0], 2) * x[0] /
                    Math.sqrt(Math.pow(DA[0][0] * x[0], 2) + Math.pow(DA[0][1] * x[1], 2) + Math.pow(Db[0], 2));
            grad1 += r * EA[0][1] + ReverseNormalFunction(p[0]) * Math.pow(DA[0][1], 2) * x[1] /
                    Math.sqrt(Math.pow(DA[0][0] * x[0], 2) + Math.pow(DA[0][1] * x[1], 2) + Math.pow(Db[0], 2));
        }
        if (G[1] >= 0) {
            grad0 += r * EA[1][0] + ReverseNormalFunction(p[1]) * Math.pow(DA[1][0], 2) * x[0] /
                    Math.sqrt(Math.pow(DA[1][0] * x[0], 2) + Math.pow(DA[1][1] * x[1], 2) + Math.pow(Db[1], 2));
            grad1 += r * EA[1][1] + ReverseNormalFunction(p[1]) * Math.pow(DA[1][1], 2) * x[1] /
                    Math.sqrt(Math.pow(DA[1][0] * x[0], 2) + Math.pow(DA[1][1] * x[1], 2) + Math.pow(Db[1], 2));
        }
        if (G[2] >= 0)
            grad0 += -r;
        if (G[3] >= 0)
            grad0 += r;
        if (G[4] >= 0)
            grad1 += -r;
        if (G[5] >= 0)
            grad1 += r;
        return new double[]{grad0, grad1};
    }

    /**
     * Calculate reverse normal function (normal distribution) from tables in point p
     * @param p - point
     * @return value of reverse normal function
     */
    private static double ReverseNormalFunction(double p) {
        double eps = 0.00001;
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
