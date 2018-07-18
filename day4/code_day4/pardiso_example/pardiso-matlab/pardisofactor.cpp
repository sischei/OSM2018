#include "mex.h"
#include "sparsematrix.h"
#include "pardisoinfo.h"

void mexFunction (int nlhs, mxArray* plhs[], 
		  int nrhs, const mxArray* prhs[]) {
  const mxArray* ptr;

  // Check to see if we have the correct number of input and output
  // arguments.
  if (nrhs != 3)
    mexErrMsgTxt("Incorrect number of input arguments");
  if (nlhs != 1)
    mexErrMsgTxt("Incorrect number of output arguments");

  // Get the second input, which contains PARDISO's internal data
  // structures.
  ptr = prhs[1];
  if (!PardisoInfo::isValid(ptr))
    mexErrMsgTxt("The second input must a STRUCT initialized with \
PARDISOINIT");
  PardisoInfo info(ptr);

  // Convert the first input, a matrix from MATLAB's sparse matrix
  // format, to PARDISO's sparse matrix format.
  ptr = prhs[0];
  if (!SparseMatrix::isValid(ptr))
    mexErrMsgTxt("The first input must be a sparse, square matrix");
  SparseMatrix A(ptr,useComplexNumbers(info),useConjTranspose(info));

  // Get the third input, the level of verbosity.
  ptr = prhs[2];
  if (!mxIsLogicalScalar(ptr))
    mexErrMsgTxt("The third input must be either TRUE or FALSE");
  bool verbose = mxIsLogicalScalarTrue(ptr);

  // Ask PARDISO to factorize the matrix.
  info.factor(A,verbose);

  // Return the modified PARDISO internal data structures to the user.
  plhs[0] = info;
}
