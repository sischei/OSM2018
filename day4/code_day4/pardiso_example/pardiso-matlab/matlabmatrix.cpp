#include "matlabmatrix.h"
#include "common.h"

// Function definitions for class MatlabMatrix.
// -----------------------------------------------------------------
bool MatlabMatrix::isValid (const mxArray* ptr) {
  return mxIsDouble(ptr) && mxGetNumberOfDimensions(ptr) <= 2;
}

MatlabMatrix::MatlabMatrix (int h, int w, bool iscomplex)
  : iscomplex(iscomplex), h(h), w(w), ar(0), ac(0) {
  initialize();
}

MatlabMatrix::MatlabMatrix (const mxArray* ptr, bool iscomplex) 
  : iscomplex(iscomplex), h(mxGetM(ptr)), w(mxGetN(ptr)), ar(0), ac(0) {
  int     n  = h * w;  // The total number of matrix entries.
  double* pr = mxGetPr(ptr);
  double* pi = mxGetPi(ptr);
  initialize();

  // Copy the entries of the MATLAB array.
  if (!mxIsComplex(ptr)) {

    // Copy the real matrix entries.
    if (!iscomplex)
      copymemory<double>(pr,ar,h*w); 
    else
      for (int i = 0; i < n; i++)
	ac[i].re = pr[i];
  }
  else {
    
    // Copy the complex matrix entries.
    if (iscomplex) 
      for (int i = 0; i < n; i++) {
	ac[i].re = pr[i];
	ac[i].im = pi[i];
      }
    else
      mexErrMsgTxt("PARDISO was initialized to handle real matrices, but \
you provided a matrix with complex entries");
  }
}

MatlabMatrix::~MatlabMatrix() {
  if (ar) delete[] ar;
  if (ac) delete[] ac;
}

MatlabMatrix::operator mxArray*() const {
  mxArray* ptr;       // The return value.
  int      n = h * w; // Number of entries in the matrix.

  if (!iscomplex) {

    // Allocate storage for the MATLAB array containing real-valued
    // entries, then copy the entries into the MATLAB array.
    ptr = mxCreateDoubleMatrix(h,w,mxREAL);
    copymemory<double>(ar,mxGetPr(ptr),n);
  }
  else {

    // Allocate storage for the MATLAB array containing complex
    // numbers as entries.
    ptr = mxCreateDoubleMatrix(h,w,mxCOMPLEX);

    // Copy the complex entries into the MATLAB array.
    double* pr = mxGetPr(ptr);
    double* pi = mxGetPi(ptr);
    for (int i = 0; i < n; i++) {
      pr[i] = ac[i].re;
      pi[i] = ac[i].im;
    }
  }

  return ptr;
}

void MatlabMatrix::initialize() {
  int n = h * w;  // The number of matrix entries.

  if (!iscomplex) {

    // Allocate memory for the matrix entries.
    ar = new double[h*w];

    // Initialize the entries of the matrix to zero.
    for (int i = 0; i < n; i++)
      ar[i] = 0;
  } else {

    // Allocate memory for the matrix entries.
    ac = new complex[h*w];

    // Initialize the complex entries of the matrix to zero.
    for (int i = 0; i < n; i++) {
      ac[i].re = 0;
      ac[i].im = 0;
    }
  }
}
