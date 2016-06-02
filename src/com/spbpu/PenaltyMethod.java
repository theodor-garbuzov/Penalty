package com.spbpu;

public class PenaltyMethod {

    public static int Optimize(double[] x, Function function, Function newFunction, Gradient newGradient,
                               double eps, double eps1, double eps2) {
        int iCount = 0;
        double penalty = newFunction.f(x) - function.f(x);
        while (penalty > eps && iCount < 100) {
            System.out.println("Штраф: " + penalty + "\n");
            GradientDescent.Optimize(x, newFunction, newGradient, eps1, eps2);
            ++iCount;
            FunctionSet.IncreaseR();
            penalty = newFunction.f(x) - function.f(x);
        }
        return iCount;
    }



}
