# teh benchmark of parallel and serial versions
addpath("../lib");

function [B, C] = splitMatrix(A, p)
    epsilon = 1e-8;
    criteria = p'*A - dot(p, p);
    B = A(:, criteria <= epsilon);
    C = A(:, criteria > -epsilon);
endfunction

benchmarkSizes = 10:10:1000;

function w = genInitialSplit(X, splitFunc)
    cols = columns(X);
    w = ones(1, cols) / cols;
endfunction

printf("size,serial_time,parallel_time,native_time\n");
for sz = benchmarkSizes
    printf('%d,', sz);
    A = testgen(sz);
    w = genInitialSplit(A, @splitMatrix);
    tic;
    # do calc
    s1 = SimProPartSerial(A, w, @splitMatrix);
    printf('%f,', toc);
    tic;
    # do calc
    s2 = SimProPart(A, w, @splitMatrix);
    printf('%f,', toc);
    tic;
    # do calc
    s3 = SimProPartNative(A, w, @splitMatrix);
    printf('%f\n', toc);
endfor;

# vim:ft=octave
