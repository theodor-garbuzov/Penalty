package com.spbpu;

/**
 * Interface to send a function gradient to the optimization function
 */
public interface Gradient {
    double[] gradf(double[] x);
}
