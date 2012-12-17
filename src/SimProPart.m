% X         - source matrix 
% splitFunc - function handler for some matrix decomposition method, accepts X and initial vector,
%   which formed as convex combinations of vectors from X

function z = SimProPart(X, splitFunc)
    epsi = 1.e-3;
    do % TODO: find a better a way to force non empty decomposition
        V = rand(1, columns(X));
        V = V / norm(V);
        [X1, X2] = splitFunc(X, X * V', 2);
    until (columns(X1) > 0 && columns(X2) > 0)

    x1 = SimPro(X1, 1e5, 1.e-8, [-1 -1], [], []); % TODO: don't repeat parameters lists
    x2 = SimPro(X2, 1e5, 1.e-8, [-1 -1], [], []);

    z = SimPro([x1 x2], 1e5, 1.e-8, [-1 -1], [], []);

    while any((dot(X, repmat(z, 1, columns(X))) - dot(z, z)) < -epsi)
        x1 = SimPro([X1 z], 1e5, 1.e-8, [-1 -1], [], []);
        x2 = SimPro([X2 z], 1e5, 1.e-8, [-1 -1], [], []);

        z = SimPro([x1 x2], 1e5, 1.e-8, [-1 -1], [], []);
    endwhile
endfunction


