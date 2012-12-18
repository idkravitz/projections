% X         - source matrix 
% splitFunc - function handler for some matrix decomposition method, accepts X and initial vector,
%   which formed as convex combinations of vectors from X

function z = SimProPart(X, splitFunc)
    epsilon = 1.e-8;
    do % TODO: find a better a way to force non empty decomposition
        w = rand(1, columns(X));
        w = w / norm(w);
        [X1, X2] = splitFunc(X, X*w');
    until (columns(X1) > 0 && columns(X2) > 0)

    x1 = SimPro(X1, 1e5, epsilon, [-1 -1], [], []); % TODO: don't repeat parameters lists
    x2 = SimPro(X2, 1e5, epsilon, [-1 -1], [], []);

    z = SimPro([x1 x2], 1e5, epsilon, [-1 -1], [], [])
    it = 1;
    while any(z'*X - sumsq(z) < -epsilon) && ( it < 1000)
        [X1, X2] = splitFunc(X, z);
        x1 = SimPro([X1 z], 1e5, epsilon, [-1 -1], [], []);
        x2 = SimPro([X2 z], 1e5, epsilon, [-1 -1], [], []);
        [z, reps, iter, lmb, kvec, R1, info] = SimPro([x1 x2], 1e5, epsilon, [-1 -1], [], []);
        printf(" it %4d norm(z) %20.12e\n", it, norm(z));
        it++;
    endwhile
endfunction
