#include "sparsematrix.h"
#include "common.h"

// Function definitions for class SparseMatrix.
// -----------------------------------------------------------------
bool SparseMatrix::isValid (const mxArray* ptr) {
  return mxIsSparse(ptr) && (mxGetM(ptr) == mxGetN(ptr));
}

SparseMatrix::SparseMatrix (const mxArray* ptr, bool iscomplex, 
			    bool useconjtrans) 
  : iscomplex(iscomplex), ia(0), ja(0), ar(0), ac(0) {
  double* pr = mxGetPr(ptr);
  double* pi = mxGetPi(ptr);

  // Get the height and width of the sparse matrix, the number of
  // nonzeros, and the values of the nonzero entries.
  n   = mxGetM(ptr);
  nnz = mxGetNzmax(ptr);

  // Allocate memory for the nonzero entries, then initialize the
  // complex entries to zero.
  if (!iscomplex)
    ar = new double[nnz];
  else {
    ac = new complex[nnz];
    for (int i = 0; i < nnz; i++) {
      ac[i].re = 0;
      ac[i].im = 0;
    }
  }

  // Copy the values of the nonzero entries.
  if (!mxIsComplex(ptr)) {

    // Copy the real matrix entries.
    if (!iscomplex)
      copymemory<double>(pr,ar,nnz);
    else
      for (int i = 0; i < nnz; i++)
	ac[i].re = pr[i];
  } else {
    if (iscomplex) {

      // Copy the real matrix entries.
      for (int i = 0; i < nnz; i++)
	ac[i].re = pr[i];

      // Copy the complex matrix entries. Note that sometimes we take
      // the complex conjugate of each entry.
      if (useconjtrans)
	for (int i = 0; i < nnz; i++)
	  ac[i].im = -pi[i];
      else
	for (int i = 0; i < nnz; i++)
	  ac[i].im = pi[i];
    } else
      mexErrMsgTxt("PARDISO was initialized to handle real matrices, but \
you provided a sparse matrix with complex entries");
  }

  // Copy the row and column indices of the nonzero entries,
  // converting the matrix from 0-based C++ notation to Fortran
  // 1-based notation. Note that in the process of doing this, we
  // transpose the matrix.
  mwIndex* jc = mxGetJc(ptr);
  mwIndex* ir = mxGetIr(ptr);
  ia          = new int[n+1];
  ja          = new int[nnz];
  for (int i = 0; i < n + 1; i++)
    ia[i] = (int) jc[i] + 1;
  for (int i = 0; i < nnz; i++)
    ja[i] = (int) ir[i] + 1;
}

SparseMatrix::~SparseMatrix() {
  if (ar) delete[] ar;
  if (ac) delete[] ac;
  if (ia) delete[] ia;
  if (ja) delete[] ja;
}

// Function definitions for friends of class SparseMatrix.
// -----------------------------------------------------------------
bool isUpperTriangular (const SparseMatrix& A) {
  bool good = true;
  int  c;

  // Repeat for each row. Note that we follow the Fortran convention
  // whereby the row and column indices start at 1.
  for (int r = 1, i = 1; r <= A.n; r++)
    
    // Repeat for each nonzero entry in the row.
    for ( ; i < A.ia[r]; i++) {
      c    = A.ja[i-1];
      good = good && (c >= r);
    }

  return good;
}

bool diagonalIsPresent (const SparseMatrix& A) {
  bool valid = true;
  bool present;
  int  c;

  // Repeat for each row. Note that we follow the Fortran convention
  // whereby the row and column indices start at 1.
  for (int r = 1, i = 1; r <= A.n; r++) {

    // Check whether the diagonal entry for the ith row is present.
    present = false;
    for ( ; i < A.ia[r]; i++) {
      c       = A.ja[i-1];
      present = present || (r == c);
    }

    valid = valid && present;
  }

  return valid;
}

