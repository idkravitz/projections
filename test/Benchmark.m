# teh benchmark of parallel and serial versions
addpath("../lib");

function [B, C] = splitMatrix(A, p)
    criteria = p'*A - dot(p, p);
    B = A(:, criteria <= 0 );
    C = A(:, criteria > 0);
endfunction

benchmarkSizes = [3 10 50 100 150];

printf("[\n");
for sz = benchmarkSizes
    printf('{\n  "size": %d,\n', sz);
    A = testgen(sz);
    tic;
    # do calc
    s1 = SimProPartSerial(A, @splitMatrix);
    printf('  "serial_time": %f,\n', toc);
    tic;
    # do calc
    s2 = SimProPart(A, @splitMatrix);
    printf('  "parallel_time": %f\n},\n', toc);
endfor;
printf("]");

# vim:ft=octave
