% This is an example script demonstrating how PARDISO works on a fairly
% large sparse, symmetric and complex matrix.
clear

% Script parameters.
% -----------------
verbose = false;
l       = 2;
k       = l * pi;
r0      = 0.5;
R       = 2;
M       = 15 * l;
N       = 120 * l; 

fprintf('Grid size: %d x %d\n',M,N);

% Construct the complex matrix C.
% ------------------------------
n  = M*N; 
r  = linspace(r0,R,M);
dr = (R - r0)/(M - 1);

dtheta = 2 * pi/N;
theta  = 0 : dtheta : (2*pi - dtheta);
A      = sparse([],[],[],n,n,5*n);
b      = zeros(n,1);

ratio = ones(M-1,1);
for j = 2:M-1
  alpha0 = 1/(dr^2) - 1/(2*dr*r(j));
  alpha1 = 1/(r(j)*dtheta)^2;
  alpha2 = k^2 - 2/(dr^2) - 2/(r(j)^2 * dtheta^2);

  if j > 2
    ratio(j-1) = ratio(j-2) * alpha0 / alpha3;
  end

  alpha3 = 1/(dr^2) + 1/(2*dr*r(j));

  for j2 = (j-2)*N+1:(j-1)*N
     A(j2+N,j2) = alpha0;
  end
  
  j2 = (j-1)*N + 1;
  A(j2:j2+1,j2) = [ alpha2
                    alpha1 ];

  for j2 = (j-1)*N+2:j*N-1		
    A(j2-1:j2+1,j2) = [ alpha1
                        alpha2
                        alpha1 ];
  end					
  
  j2=j*N;				
  A(j2-1:j2,j2) = [ alpha1
                    alpha2 ];	
  
  for j2 = j*N+1:(j+1)*N
    A(j2-N,j2) = alpha3;
  end
  
  A((j-1)*N+1,j*N) = alpha1; 
  A(j*N,(j-1)*N+1) = alpha1; 
end

for j = 1:N
    A(j,j)=1;
end

for j = 1:N
  b(j) = -exp(i*k*r(1)*cos(theta(j)));
end

beta0      = 2/(dr^2);
beta1      = 1/(r(end)*dtheta)^2;
beta2      = k^2 - 2/(dr^2) - 2/(r(end)*dtheta)^2 + ...
             2*i*k*dr*(1/dr^2 + 1/(2*dr*r(end)));
ratio(M-1) = ratio(M-2) * beta0 / alpha3;
	
for j = (M-2)*N+1:(M-1)*N
  A(j+N,j) = beta0;
end

j2            = (M-1)*N+1;      
A(j2:j2+1,j2) = [ beta2
                  beta1 ];	

for j2 = (M-1)*N+2:M*N-1		
  A(j2-1:j2+1,j2) = [ beta1
                      beta2
                      beta1 ];
end		

j2                = M*N;
A(j2-1:j2,j2)     = [ beta1
                      beta2 ];
A((M-1)*N+1, M*N) = beta1;
A(M*N,(M-1)*N+1)  = beta1;

B        = A;
B(1:N,:) = [];
B(:,1:N) = [];
c        = b;
c(1:N)   = [];

for j = 1:N
  c(1:N) = c(1:N) - A(N+1:2*N,j) * b(j);
end

dL = sqrt(abs(1./ratio));
dR = sqrt(abs(ratio)) .* sign(ratio);
DL = kron(spdiags(dL,0,M-1,M-1),speye(N));
DR = kron(spdiags(dR,0,M-1,M-1),speye(N));
C  = DL * B * DR;
C  = (C + C.')/2;
c  = DL * c;

% Solve for v in linear system C * v = c.
% --------------------------------------
% Initialize the PARDISO internal data structures. We've told PARDISO to
% handle complex matrices using the sparse direct solver.
info = pardisoinit(6,0);

% Analyze the matrix and compute a symbolic factorization.
info = pardisoreorder(tril(C),info,verbose);
fprintf('The factors have %d nonzero entries.\n',info.iparm(18));

% Compute the numeric factorization.
info = pardisofactor(tril(C),info,verbose);

% Compute the solution v using the symbolic factorization.
[v info] = pardisosolve(tril(C),c,info,verbose);
fprintf('PARDISO performed %d iterative refinement steps.\n',info.iparm(7));

% Compute the residuals.
r = abs(C*v - c);
fprintf('The maximum residual for the solution is %0.3g.\n',max(r));

% Free the PARDISO data structures.
pardisofree(info);
clear info
