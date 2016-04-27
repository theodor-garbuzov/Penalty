package com.spbpu;

public class Main {

    public static void main(String[] args) {
        double[] x = new double[]{0, 0};
        int iCount;
        double eps = 0.0001; // exit parameter of penalty method
        double eps1 = 0.01; // exit parameter of gradient method
        double eps2 = 0.01; // parameter of one-dimensional optimization

        System.out.println("ЭПСИЛОН " + eps);
        iCount = PenaltyMethod.Optimize(x, FunctionSet::f, FunctionSet::NewFunction, FunctionSet::NewGradient, eps, eps1, eps2);
        //iCount = GradientDescent.Optimize(x, FunctionSet::g, FunctionSet::gradg, eps1, eps2);
        System.out.println("Количество шагов: " + iCount);
        System.out.println("Точка оптимума: " + x[0] + " " + x[1]);
        //System.arraycopy(x0, 0, x, 0, x.length);

    }
}
