#ifndef INCLUDE_MATLABMATRIX
#define INCLUDE_MATLABMATRIX

#include "mex.h"

// Type definitions.
// -----------------------------------------------------------------
// Storage for a single complex number.
typedef struct {
  double re; 
  double im;
} complex;

// Class declarations.
// -----------------------------------------------------------------
class PardisoInfo;

// Class MatlabMatrix.
// -----------------------------------------------------------------
// A matrix object stores its elements in column-major format, as in
// Fortran and Matlab. This means that columns are stored one after
// another. For example, the matrix
//
//   1  2  3
//   4  5  6
//
// is stored in memory as
//
//   1  4  2  5  3  6
//
// The elements may either be real numbers or complex numbers.
class MatlabMatrix {
  friend class PardisoInfo;

public:

  // Return true if and only if the MATLAB array is a valid matrix.
  static bool isValid (const mxArray* ptr);

  // This constructor allocates memory for matrix (with either real or
  // complex entries) of the specified height and width. The entries
  // of the matrix are initialized to zero.
  MatlabMatrix (int h, int w, bool iscomplex);

  // This constructor retrieves a matrix from a Matlab array.
  MatlabMatrix (const mxArray* ptr, bool iscomplex);

  // The destructor.
  ~MatlabMatrix();

  // Return the height and width of the matrix.
  friend int height (const MatlabMatrix& A) { return A.h; };
  friend int width  (const MatlabMatrix& A) { return A.w; };

  // Convert the matrix object to a MATLAB array.
  operator mxArray*() const;

protected:
  bool     iscomplex;
  int      h;   // Height of the matrix.
  int      w;   // Width of the matrix.
  double*  ar;  // The matrix entries if it is a real matrix.
  complex* ac;  // The matrix entries if it is a complex matrix.

  // Helper function used by the class constructors.
  void initialize();
};

#endif
