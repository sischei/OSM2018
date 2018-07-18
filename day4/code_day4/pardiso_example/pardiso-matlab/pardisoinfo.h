#ifndef INCLUDE_PARDISOINFO
#define INCLUDE_PARDISOINFO

#include "mex.h"
#include "matlabmatrix.h"
#include "sparsematrix.h"

// Macros and definitions.
// -----------------------------------------------------------------
#define real_structsym         1
#define real_symmetric_pdef    2
#define real_symmetric_indef  -2
#define cmpx_structsym         3
#define cmpx_hermitian_pdef    4
#define cmpx_hermitian_indef  -4
#define cmpx_symmetric         6
#define real_nonsym           11
#define cmpx_nonsym           13

// External Fortran routines defined in the PARDISO library.
extern "C" int pardisoinit_ (void *, int *, int *, int *, double *, int *);
extern "C" int pardiso_     (void *, int *, int *, int *, int *, int *, 
			     void *, int *, int *, int *, int *, int *, 
			     int *, void *, void *, int *, double *);

// Function declarations.
// -----------------------------------------------------------------
// Reports a PARDISO error to MATLAB.
void reportPardisoError (int error);

// Class PardisoInfo.
// -----------------------------------------------------------------
// An object of this class stores all the parameters and internal data
// structures that are common to all PARDISO subroutines.
class PardisoInfo { 
public:

  // Return true if and only if the MATLAB array is a structure with
  // the right field names.
  static bool isValid (const mxArray* ptr);

  // This constructor initializes the internal structures for PARDISO
  // and allocates storage for these structures within a MATLAB
  // array. Note that the MATLAB array should not be deallocated if
  // you want to pass the structure back to MATLAB. See the PARDISO
  // manual for more information on the two inputs.
  PardisoInfo (int mtype, int solver);
  
  // This constructor ,akes a copy of the PARDISO internal data from
  // the MATLAB structure array passed as input to this function.
  explicit PardisoInfo (const mxArray* ptr);

  // Note that the destructor does not deallocate the memory
  // associated with PARDISO's internal structures.
  ~PardisoInfo() { };

  // Get the pointer to the MATLAB structure storing the object info.
  operator mxArray*() const { return matlabptr; };

  // The first function returns 1 if we are using the iterative
  // solver; otherwise returns 0. The second function returns true if
  // and only if the user specified the use of complex matrices. The
  // third function returns true if we should be using the complex
  // conjugate transpose (this is CTRANSPOSE in MATLAB). This will
  // depend on the matrix type. The fourth function returns true if
  // and only if we are performing computations with symmetric or
  // Hermitian matrices.
  friend bool useIterativeSolver   (const PardisoInfo& info);
  friend bool useComplexNumbers    (const PardisoInfo& info);
  friend bool useConjTranspose     (const PardisoInfo& info);
  friend bool useSymmetricMatrices (const PardisoInfo& info);

  // Checks whether a sparse symmetric matrix is upper triangular (in
  // other words, lower triangular in MATLAB).
  void checkMatrix (const SparseMatrix& A) const;

  // Ask PARDISO to create a symbolic factorization for the sparsity
  // structure given by the sparse matrix A. Optionally, the user may
  // provide a pre-defined permutation (with indices starting at 1) in
  // the input array "perm".
  void reorder (SparseMatrix& A, int* perm, bool verbose);

  // Ask PARDISO to create a numeric factorizization of the sparse matrix.
  void factor (SparseMatrix& A, bool verbose);

  // Ask PARDISO to solve the system(s) of equations AX = B, returning
  // the solution(s) in the matrix X.
  void solve (SparseMatrix& A, MatlabMatrix& B, MatlabMatrix& X, 
	      bool verbose);

  // Ask PARDISO release memory associated with all internal structures. 
  void free();

protected:
  int      mtype;     // The matrix type.
  void*    pt;        // PARDISO's internal data address pointer.
  int*     iparm;     // Miscellaneous information and parameters.
  double*  dparm;     // More miscellaneous information.
  mxArray* matlabptr; // Pointer to the structure in MATLAB.

  // These are routines used by some of the class methods above. The
  // first routine is used by the "free" method. The second routine is
  // used by the "reorder" and "factor" methods. And the third routine
  // is used by the "solve" method.
  void callpardiso (int phase, int msglvl);
  void callpardiso (int phase, int msglvl, SparseMatrix& A, int* perm = 0);
  void callpardiso (int phase, int msglvl, SparseMatrix& A, 
		    MatlabMatrix& B, MatlabMatrix& X);
};

#endif
