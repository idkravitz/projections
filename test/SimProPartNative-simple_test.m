#!/usr/bin/env octave
addpath("../lib");

function [B, C] = splitEvenOdd(A)
    B = A(:, 1:2:end);
    C = A(:, 2:2:end);
endfunction

function [B, C] = splitMatrix(A, p)
    criteria = p'*A - dot(p, p);
    B = A(:, criteria <= 0 );
    C = A(:, criteria > 0);
endfunction

A = testgen(4);

sol2 = SimProPartNative(A, @splitMatrix)
norm(sol2)
sol = SimPro(A, 1e5, 1.e-8, [-1 -1], [], [])
norm(sol)
# vim:ft=octave
