#include "mex.h"
#include "matlabmatrix.h"
#include "sparsematrix.h"
#include "pardisoinfo.h"

void mexFunction (int nlhs, mxArray* plhs[], 
		  int nrhs, const mxArray* prhs[]) {
  const mxArray* ptr;

  // Check to see if we have the correct number of input and output
  // arguments.
  if (nrhs != 4)
    mexErrMsgTxt("Incorrect number of input arguments");
  if (nlhs != 2)
    mexErrMsgTxt("Incorrect number of output arguments");

  // Get the third input, which contains PARDISO's internal data
  // structures.
  ptr = prhs[2];
  if (!PardisoInfo::isValid(ptr))
    mexErrMsgTxt("The third input must a STRUCT initialized with \
PARDISOINIT");
  PardisoInfo info(ptr);

  // Convert the first input to PARDISO's sparse matrix format.
  ptr = prhs[0];
  if (!SparseMatrix::isValid(ptr))
    mexErrMsgTxt("The first input must be a sparse, square matrix");
  SparseMatrix A(ptr,useComplexNumbers(info),useConjTranspose(info));

  // Process the second input, the matrix of right-hand sides B. Each
  // column of B corresponds to a single right-hand side vector.
  ptr = prhs[1];
  if (!MatlabMatrix::isValid(ptr))
    mexErrMsgTxt("The second input must be a matrix in DOUBLE precision");
  MatlabMatrix B(ptr,useComplexNumbers(info));

  // Report an error if we are using the iterative solver, and there
  // is more than one right-hand side (i.e. the width of the matrix B
  // is greater than 1).
  if (useIterativeSolver(info) && width(B) > 1)
    mexErrMsgTxt("The iterative solver can only compute the solution \
for a single right-hand side");

  // Get the fourth input, the level of verbosity.
  ptr = prhs[3];
  if (!mxIsLogicalScalar(ptr))
    mexErrMsgTxt("The fourth input must be either TRUE or FALSE");
  bool verbose = mxIsLogicalScalarTrue(ptr);

  // Ask PARDISO to solve the system(s) of equations.
  MatlabMatrix X(height(B),width(B),useComplexNumbers(info));
  info.solve(A,B,X,verbose);

  // The first ouput is the matrix of solutions X, and the second
  // output is the updated PARDISO internal data structures.
  plhs[0] = X;
  plhs[1] = info;
}

