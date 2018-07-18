% This is an example script demonstrating how PARDISO works on a small,
% sparse, real symmetric matrix. It computes the m solutions X to the
% collection of linear systems
%
%    A * X = B
%
% using the PARDISO solver, where A is a symmetric n x n matrix, B is an
% n x m matrix, and X is another n x m matrix.
clear
verbose = false;

n = 4;  % The number of equations.
m = 3;  % The number of right-hand sides.

A = sparse([ 1  0 -2  3
             0  5  1  2
            -2  1  4 -7
             3  2 -7  5 ]);

% Generate a random collection of right-hand sides.
B = randn(n,m);

% Initialize the PARDISO internal data structures. We've told PARDISO to
% handle real symmetric matrices using the sparse direct solver.
% info = pardisoinit(-2,0);
info = pardisoinit(-2,0);

% Analyze the matrix and compute a symbolic factorization.
p    = randperm(n);
info = pardisoreorder(tril(A),info,verbose,p);
fprintf('The factors have %d nonzero entries.\n',info.iparm(18));

% Compute the numeric factorization.
info = pardisofactor(tril(A),info,verbose);
fprintf('The matrix has %d positive and %d negative eigenvalues.\n',...
	info.iparm(22),info.iparm(23));

% Compute the solutions X using the symbolic factorization.
[X info] = pardisosolve(tril(A),B,info,verbose);
fprintf('PARDISO performed %d iterative refinement steps.\n',info.iparm(7));

% Compute the residuals.
R = max(abs(A*X - B));
fprintf('The maximum residual for the solution X is %0.3g.\n',max(R(:)));

% Free the PARDISO data structures.
pardisofree(info);
clear info


