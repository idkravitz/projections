#!/usr/bin/env octave
addpath("../src");

X = eye(3)
epsilon = 1.e-8;
maxit = 1e6;
[z, reps, iter, lmb, kvec, R1, info] = SimPro(X, maxit, epsilon, [1 1], [], [])
