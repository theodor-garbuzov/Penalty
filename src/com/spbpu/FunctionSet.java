package com.spbpu;

import static com.spbpu.VectorOperations.VectorVectorMult;

/**
 * Set of functions and their gradients
 */
public class FunctionSet {

    private static double a1 = 1, a2 = 2, a3 = 3;

    public static double f(double x[])
    {
        return x[0]*x[0] + a1*x[1]*x[1] + Math.exp(a2*x[0] + c*x[1]) - x[0] + 2*x[1];
    }
    public static double[] gradf(double[] x)
    {
        double g1 = 2*x[0] + a2*Math.exp(a2*x[0] + a3*x[1]) - 1;
        double g2 = 2*a1*x[1] + c*Math.exp(a2*x[0] + a3*x[1]) + 2;
        return new double[]{g1, g2};
    }

    public double g(double[] x)
    {
        return a1*x[0] + x[1] + 4*Math.sqrt(1 + a2*x[0]*x[0] + a3*x[1]*x[1]);
    }
    public double[] gradg(double[] x)
    {
        double g1 = a1 + 4*a2*x[0] / (Math.sqrt(1 + a2*x[0]*x[0] + a3*x[1]*x[1]));
        double g2 = 1 + 4*c*x[1] / (Math.sqrt(1 + a2*x[0]*x[0] + a3*x[1]*x[1]));
        return new double[]{g1, g2};
    }

    // PARAMETERS OF THE TASK
    private static double[][] d = new double[][]{{0, 1}, {0, 1}}; /* matrix of borders of allowed values of x
                                                                    (example: d[1][1] <= x[1] <= d[1][2]) */
    private static double[] p = new double[]{0.5, 0.5}; // vector of probability levels (alpha_i)
    private static double[][] EA = new double[][]{{0, 1}, {0, 1}}; // expectancy of coefficient matrix A
    private static double[][] DA = new double[][]{{0, 1}, {0, 1}}; // dispersion of coefficient matrix A
    private static double[] Eb = new double[]{0, 1}; // expectancy of vector b
    private static double[] Db = new double[]{0, 1}; // dispersion of vector b
    private static double c = 1; // coefficient of penalty function

    /**
     * Multiply c by 2
     */
    public static void IncreaseC() {
        c *= 2;
    }

    /**
     * Calculate values of limit functions in provided point
     * @param x - point
     * @return vector with values of limit functions
     */
    public static double[] LimitFunctions(double[] x) {
        double t1 = ReverseNormalFunction(p[1]);
        double t2 = ReverseNormalFunction(p[2]);
        double g1 = VectorVectorMult(EA[1], x) - Eb[1] +
                t1 * Math.sqrt(Math.pow(DA[1][1] * x[1],2) + Math.pow(DA[1][2] * x[2],2) + Math.pow(Db[1],2));
        double g2 = VectorVectorMult(EA[2], x) - Eb[2] +
                t2 * Math.sqrt(Math.pow(DA[2][1] * x[1],2) + Math.pow(DA[2][2] * x[2],2) + Math.pow(Db[2],2));
        double g3 = d[1][1] - x[1], g4 = x[1] - d[1][2];
        double g5 = d[2][1] - x[2], g6 = x[2] - d[2][2];
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
            limitFunctionsSum += GItem;
        }
        return f(x) + c * limitFunctionsSum;
    }

    public static double[] NewGradient(double[] x) {
        double[] gradfunction = gradf(x);
        double grad1 = gradfunction[1];
        double grad2 = gradfunction[2];
        double[] G = LimitFunctions(x);
        if(G[1] >= 0) {
            grad1 += c * EA[1][1] + ReverseNormalFunction(p[1]) * Math.pow(DA[1][1], 2) * x[1] /
                    Math.sqrt(Math.pow(DA[1][1] * x[1], 2) + Math.pow(DA[1][2] * x[2], 2) + Math.pow(Db[1], 2));
            grad2 += c * EA[1][2] + ReverseNormalFunction(p[1]) * Math.pow(DA[1][2], 2) * x[2] /
                    Math.sqrt(Math.pow(DA[1][1] * x[1], 2) + Math.pow(DA[1][2] * x[2], 2) + Math.pow(Db[1], 2));
        }
        if (G[2] >= 0) {
            grad1 += c * EA[2][1] + ReverseNormalFunction(p[2]) * Math.pow(DA[2][1], 2) * x[1] /
                    Math.sqrt(Math.pow(DA[2][1] * x[1], 2) + Math.pow(DA[2][2] * x[2], 2) + Math.pow(Db[2], 2));
            grad2 += c * EA[2][2] + ReverseNormalFunction(p[2]) * Math.pow(DA[2][2], 2) * x[2] /
                    Math.sqrt(Math.pow(DA[2][1] * x[1], 2) + Math.pow(DA[2][2] * x[2], 2) + Math.pow(Db[2], 2));
        }
        if (G[3] >= 0)
            grad1 += -c;
        if (G[4] >= 0)
            grad1 += c;
        if (G[5] >= 0)
            grad2 += -c;
        if (G[6] >= 0)
            grad2 += c;
        return new double[]{grad1, grad2};
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
