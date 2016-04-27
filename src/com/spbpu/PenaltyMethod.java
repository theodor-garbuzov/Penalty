package com.spbpu;

import static com.spbpu.VectorOperations.*;

public class PenaltyMethod {

    public static int Optimize(double[] x, Function function, Function newFunction, Gradient newGradient,
                               double eps, double eps1, double eps2) {
        int iCount = 0;
        double tmp = newFunction.f(x) - function.f(x);
        while (tmp > eps) {
            GradientDescent.Optimize(x, newFunction, newGradient, eps1, eps2);
            ++iCount;
            FunctionSet.IncreaseR();
            tmp = newFunction.f(x) - function.f(x);
        }
        return iCount;
    }



}
