#ifndef INCLUDE_SPARSEMATRIX
#define INCLUDE_SPARSEMATRIX

#include "mex.h"
#include "matlabmatrix.h"

// Class declarations.
// -----------------------------------------------------------------
class PardisoInfo;

// Class SparseMatrix.
// -----------------------------------------------------------------
// This structure stores all the information we need to know about a
// sparse symmetric matrix with double-precision real or complex
// entries. It stores the sparse matrix in PARDISO's sparse matrix
// format.
class SparseMatrix {
  friend class PardisoInfo;

public:

  // Convert the sparse matrix from MATLAB's internal sparse matrix
  // format to the PARDISO's sparse matrix format. The input should
  // point to a matrix stored in MATLAB's sparse matrix format. Note
  // that the (CONJUGATE) TRANSPOSE of the original matrix is stored.
  SparseMatrix (const mxArray* ptr, bool iscomplex, bool useconjtrans);

  // The destructor.
  ~SparseMatrix();

  // The first function returns true if and only if the matrix is
  // upper triangular. The second function returns true if and only if
  // all the diagonal entries are present (even if they are equal to
  // zero).
  friend bool isUpperTriangular (const SparseMatrix& A);
  friend bool diagonalIsPresent (const SparseMatrix& A);

  // Return true if and only if the MATLAB array is a valid square
  // sparse matrix.
  static bool isValid (const mxArray* ptr);

  // Return the size (height and width) of the matrix.
  friend int size (const SparseMatrix& A) { return A.n; };

protected:
  bool     iscomplex;
  int      n;    // Height and width of the matrix.
  int      nnz;  // The number of nonzero entries.
  int*     ia;   // An array of length n + 1 in which the ith array entry
		 // points to the index of the first nonzero element of
		 // the ith row. Note that ia[0] = 1 following Fortran's
		 // convention.
  int*     ja;   // An array of length equal to the number of nonzero
		 // entries, in which ja(i) is the column of the ith
		 // nonzero element in the sparse matrix. Note that the
		 // column indices start at 1 following Fortran's
		 // convention.
  double*  ar;   // The values of the nonzero entries when we have a
		 // real matrix.
  complex* ac;   // The values of the nonzero entries when we have a
		 // complex matrix.
};

#endif
