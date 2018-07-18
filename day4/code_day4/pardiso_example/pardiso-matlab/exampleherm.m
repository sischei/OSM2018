% This is an example script demonstrating how PARDISO works on a
% medium-sized Hermitian positive definite matrix.
clear

% Script parameters.
% -----------------
verbose = false;
n       = 100;
lambda  = 3;

% Create the Hermitian positive definite matrix A and the vector b in the
% linear system Ax = b.
e = ones(n,1);
A = spdiags([ i*e lambda*e -i*e ],-1:1,n,n);
b = randn(n,1);

% Compute solution to linear system
% ---------------------------------
% Initialize the PARDISO internal data structures. We've told PARDISO to
% handle Hermitian positive definite matrices using the sparse direct
% solver. 
info = pardisoinit(4,0);

% Analyze the matrix and compute a symbolic factorization.
info = pardisoreorder(tril(A),info,verbose);
fprintf('The factors have %d nonzero entries.\n',info.iparm(18));

% Compute the numeric factorization.
info = pardisofactor(tril(A),info,verbose);

% Compute the solution v using the symbolic factorization.
[x info] = pardisosolve(tril(A),b,info,verbose);
fprintf('PARDISO performed %d iterative refinement steps.\n',info.iparm(7));

% Compute the residuals.
r = abs(A*x - b);
fprintf('The maximum residual for the solution is %0.3g.\n',max(r));

% Free the PARDISO data structures.
pardisofree(info);
clear info
