package com.spbpu;

import static com.spbpu.VectorOperations.*;

public class GradientDescent {

    /**
     * Optimize function starting from a point x with provided exit parameter
     *
     * @param x        - start point, at the end it will be optimum point
     * @param function - function to optimize
     * @param gradient - function's gradient
     * @param eps      - exit parameter
     * @param eps2     - exit parameter of one-dimensional optimization
     * @return iCount -  number of method iterations
     */
    public static int Optimize(double[] x, Function function, Gradient gradient, double eps, double eps2) {
        double[] x_new = new double[x.length];
        double[] grad = gradient.gradf(x);
        double[] x1 = new double[x.length], x2 = new double[x.length]; // step vectors
        double step;
        int iCount = 0;

        while (norm2(grad) > eps) {
            // get new point x_new = x - step * grad
            step = GetNewPointSplitting(x, grad, x_new, function, eps2);
            System.arraycopy(VectorSum(x, NumberVectorMult(-step, grad)), 0, x_new, 0, x_new.length);

            // x = x_new
            System.arraycopy(x_new, 0, x, 0, x_new.length);
            grad = gradient.gradf(x);
            /*System.out.print(iCount); System.out.print("-й шаг. Градиент: "); System.out.println(norm2(grad));
            System.out.print(x[0]); System.out.print(" "); System.out.println(x[1]);*/
        }
        return iCount;
    }

    private static double GetNewStepDichotomy(double[] x, double[] grad, double[] x_new, Function function, double eps, int[] fCount) {
        double left = 0, right = 5; // начальный интервал неопределённости шага
        double step1, step2;
        double delta = eps / 20;
        while (right - left > eps) {
            step1 = (left + right) / 2 - delta;
            step2 = (left + right) / 2 + delta;
            fCount[0] += 2;
            if (function.f(VectorSum(x, NumberVectorMult(-step1, grad))) > function.f(VectorSum(x, NumberVectorMult(-step2, grad))))
                left = step1;
            else
                right = step2;
            assert (left < right);
        }
        return (left + right) / 2;
    }

    private static double GetNewPointSplitting(double[] x, double[] grad, double[] x_new, Function function, double eps) {
        boolean exitCondition1, exitCondition2;
        double step = 2;
        do {
            step /= 2;
            System.arraycopy(VectorSum(x, NumberVectorMult(-step, grad)), 0, x_new, 0, x_new.length);
            exitCondition1 = function.f(x_new) - function.f(x) < -eps * norm2(grad) * step;
            exitCondition2 = Math.abs(-eps * norm2(grad) * step) < Math.pow(10, -15);
        } while (!exitCondition1 && !exitCondition2);
        return step;
    }
}
