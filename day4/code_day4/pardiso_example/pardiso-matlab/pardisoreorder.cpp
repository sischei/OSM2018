#include "mex.h"
#include "sparsematrix.h"
#include "pardisoinfo.h"

void mexFunction (int nlhs, mxArray* plhs[],
		  int nrhs, const mxArray* prhs[]) {
  const mxArray* ptr;

  // Check to see if we have the correct number of input and output
  // arguments.
  if (nrhs < 3 || nrhs > 4)
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
  
  // Convert the first input to PARDISO's sparse matrix format.
  ptr = prhs[0];
  if (!SparseMatrix::isValid(ptr))
    mexErrMsgTxt("The first input must be a sparse, square matrix");
  SparseMatrix A(ptr,useComplexNumbers(info),useConjTranspose(info));

  // Get the size of the matrix.
  int n = size(A);

  // Get the third input, the level of verbosity.
  ptr = prhs[2];
  if (!mxIsLogicalScalar(ptr))
    mexErrMsgTxt("The third input must be either TRUE or FALSE");
  bool verbose = mxIsLogicalScalarTrue(ptr);

  // Get the fourth (optional) input, a user-supplied reordering of the
  // rows and columns of the system.
  int* perm = 0;
  if (nrhs == 4) {
    ptr = prhs[3];
    if (!mxIsEmpty(ptr)) {
      if (!mxIsDouble(ptr) || mxGetNumberOfElements(ptr) != n)
	mexErrMsgTxt("The permutation must be a double-precision array \
of length equal to the size of the sparse matrix");
      
      // Copy the elements of the array.
      double* p = mxGetPr(ptr);
      perm = new int[n];
      for (int i = 0; i < n; i++)
	perm[i] = (int) p[i];
    }
  }
  
  // Ask PARDISO to analyze the matrix A.
  info.reorder(A,perm,verbose);

  // Return the modified PARDISO internal data structures to the user.
  plhs[0] = info;

  // Free the dynamically allocated memory.
  delete[] perm;
}
