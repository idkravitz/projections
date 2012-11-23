function X = testgen(dim)

eps = .01;
X = diag(100 * rand(1, (dim - 1)));
X = [X zeros(dim-1, 1) ];
c = sum(X, 2) / dim;
X = X - c * ones(1, dim);
X = [ X; eps * ones(1, dim) ];
X = [ X eps * [ (rand(1, dim - 1) - 0.5) 2]'];

endfunction
