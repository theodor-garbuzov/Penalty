package com.spbpu;

public class Main {

    public static void main(String[] args) {
        double[] x = new double[]{10, 10, 10};
        int iCount;
        double eps = 0.1; // exit parameter of penalty method
        double eps1 = 0.001; // exit parameter of gradient method
        double eps2 = 0.00001; // parameter of one-dimensional optimization

        System.out.println("ЭПСИЛОН " + eps);
        iCount = PenaltyMethod.Optimize(x, FunctionSet::f, FunctionSet::NewFunction, FunctionSet::NewGradient, eps, eps1, eps2);
        //iCount = GradientDescent.Optimize(x, FunctionSet::g, FunctionSet::gradg, eps1, eps2);
        System.out.println("Количество шагов: " + iCount);
        System.out.println("Точка оптимума: " + x[0] + " " + x[1] + " " + x[2]);
        //System.arraycopy(x0, 0, x, 0, x.length);

    }
}
