package com.spbpu;

import static com.spbpu.VectorOperations.*;

public class PenaltyMethod {

    public static int Optimize(double[] x, Function function, Function NewFunction, Gradient NewGradient,
                               double eps, double eps1, double eps2) {
        int iCount = 0;
        while (NewFunction.f(x) - function.f(x) > eps) {
            GradientDescent.Optimize(x, NewFunction, NewGradient, eps1, eps2);
            ++iCount;
            FunctionSet.IncreaseC();
        }
        return iCount;
    }



}
