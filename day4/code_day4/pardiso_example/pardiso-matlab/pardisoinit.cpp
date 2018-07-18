#include "mex.h"
#include "common.h"
#include "pardisoinfo.h"

void mexFunction (int nlhs, mxArray* plhs[], 
		  int nrhs, const mxArray* prhs[]) { 
  int            mtype;   // The matrix type.
  int            solver;  // Which solver to use.
  const mxArray* ptr;     // Pointer to an input argument.

  // Check to see if we have the correct number of input and output
  // arguments.
  if (nrhs != 2)
    mexErrMsgTxt("Incorrect number of input arguments");
  if (nlhs != 1)
    mexErrMsgTxt("Incorrect number of output arguments");
  
  // The first input specifies the matrix type.
  ptr = prhs[0];
  if (mxIsDoubleScalar(ptr))
    mtype = (int) mxGetScalar(ptr);
  else
    mexErrMsgTxt("The first input must be a number specifying the matrix \
type (see the PARDISO manual for more details).");

  // The second input specifies which solver to use.
  ptr = prhs[1];
  if (mxIsDoubleScalar(ptr))
    solver = (int) mxGetScalar(ptr);
  else
    mexErrMsgTxt("The second input must be a number specifying which solver \
to use (see the PARDISO manual for more details).");

  // Create the single output containing PARDISO's internal data
  // structures.
  PardisoInfo info(mtype,solver);
  plhs[0] = info;
}
