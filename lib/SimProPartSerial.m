% X         - source matrix 
% splitFunc - function handler for some matrix decomposition method, accepts X and initial vector,
%   which formed as convex combinations of vectors from X

function z = SimProPartSerial(X, w, splitFunc)
    epsilon = 1.e-8;
    [X1, X2] = splitFunc(X, X*w');

    SimProAdjusted = @(x) (SimPro(x, 1e5, epsilon, [-1 -1], [], []));

    x1 = SimProAdjusted(X1);
    x2 = SimProAdjusted(X2);

    z = SimPro([x1 x2], 1e5, epsilon, [-1 -1], [], []);
    it = 1;
    while any(z'*X - sumsq(z) < -epsilon) && ( it < 5000)
        [X1, X2] = splitFunc(X, z);
        Y = cellfun(@(x) (SimPro(x, 1e5, epsilon, [-1 -1], [], [])), {[X1 z], [X2 z]}, 'UniformOutput', false);
        [z, reps, iter, lmb, kvec, R1, info] = SimProAdjusted([Y{1} Y{2}]);
#        printf(" it %4d norm(z) %20.12e\n", it, norm(z));
        it++;
    endwhile
endfunction
