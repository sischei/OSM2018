#include "mex.h"
#include "pardisoinfo.h"

void mexFunction (int nlhs, mxArray* plhs[],
		  int nrhs, const mxArray* prhs[]) {
  const mxArray* ptr;

  // Check to see if we have the correct number of input and output
  // arguments.
  if (nrhs != 1)
    mexErrMsgTxt("Incorrect number of input arguments");
  if (nlhs != 0)
    mexErrMsgTxt("Incorrect number of output arguments");

  // Get the input argument containing PARDISO's internal data
  // structures.
  ptr = prhs[0];
  if (!PardisoInfo::isValid(ptr))
    mexErrMsgTxt("The input must a STRUCT initialized with PARDISOINIT");
  PardisoInfo info(ptr);

  // Ask PARDISO release memory associated with all internal
  // structures.
  info.free();
}
