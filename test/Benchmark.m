# teh benchmark of parallel and serial versions
addpath("../lib");

function [B, C] = splitMatrix(A, p)
    criteria = p'*A - dot(p, p);
    B = A(:, criteria <= 0 );
    C = A(:, criteria > 0);
endfunction

benchmarkSizes = [3 5 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 100];

function w = genInitialSplit(X, splitFunc)
    do % TODO: find a better a way to force non empty decomposition
        w = rand(1, columns(X));
        w = w / norm(w);
        [X1, X2] = splitFunc(X, X*w');
    until (columns(X1) > 0 && columns(X2) > 0)
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
