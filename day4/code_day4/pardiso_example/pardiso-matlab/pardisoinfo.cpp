#include "pardisoinfo.h"

// Global variables.
// -----------------------------------------------------------------
int maxfct = 1;
int mnum   = 1;

// Function definitions.
// -----------------------------------------------------------------
void reportPardisoError (int error) {
  char str[256];
  sprintf(str,"PARDISO reports an error of type %d",error);
  mexErrMsgTxt(str);
}

// Function definitions for class PardisoInfo.
// -----------------------------------------------------------------
bool PardisoInfo::isValid (const mxArray* ptr) {
  if (!mxIsStruct(ptr))
    return false;
  else
    return mxGetField(ptr,0,"mtype") && 
	   mxGetField(ptr,0,"pt")    &&
	   mxGetField(ptr,0,"iparm") &&
	   mxGetField(ptr,0,"dparm");
}

PardisoInfo::PardisoInfo (int mtype, int solver) 
  : mtype(mtype), pt(0), iparm(0), dparm(0), matlabptr(0) {
  mxArray*  ptr;    // Pointer to a Matlab array.
  int       error;  // Error reported by PARDISO.

  // Create the Matlab structure which will store all the global 
  // information used by PARDISO.
  const char* fieldnames[] = { "mtype", "pt", "iparm", "dparm" };
  matlabptr = mxCreateStructMatrix(1,1,4,fieldnames);

  // Initialize storage for the matrix type.
  ptr = mxCreateDoubleScalar(mtype);
  mxSetField(matlabptr,0,"mtype",ptr);

  // Initialize storage for PARDISO's internal data addresses.
  ptr = mxCreateNumericMatrix(1,64,mxINT64_CLASS,mxREAL);
  pt  = (void*) mxGetPr(ptr);
  mxSetField(matlabptr,0,"pt",ptr);
  
  // Initialize storage for the "iparm" miscellaneous information.
  // Notice that we need to make sure that the individual elements of
  // the Matlab array are the same size as the integer data type in
  // C++.
  ptr = mxCreateNumericMatrix(1,64,mxINT32_CLASS,mxREAL);
  if (mxGetElementSize(ptr) != sizeof(int))
    mexErrMsgTxt("The Matlab integer array does not have elements of the \
appropriate size");
  iparm = (int*) mxGetPr(ptr);
  mxSetField(matlabptr,0,"iparm",ptr);

  // Initialize storage for the "dparm" miscellaneous information.
  ptr   = mxCreateDoubleMatrix(1,64,mxREAL);
  dparm = mxGetPr(ptr);
  mxSetField(matlabptr,0,"dparm",ptr);

  // Initialize the internal solver addresses and set parameters to
  // their default values.
  pardisoinit_(pt,&mtype,&solver,iparm,dparm,&error);
  if (error)
    reportPardisoError(error);
  iparm[2] = 1;  // Set the number of processors.
}

PardisoInfo::PardisoInfo (const mxArray* ptr) {

  // Make a copy of the MATLAB structure array.
  matlabptr = mxDuplicateArray(ptr);

  ptr   = mxGetField(matlabptr,0,"mtype");
  mtype = (int) mxGetScalar(ptr);
 
  // Get PARDISO's internal data address data.
  ptr = mxGetField(matlabptr,0,"pt");
  pt  = mxGetPr(ptr);

  // Get the "iparm" miscellaneous information structures.
  ptr   = mxGetField(matlabptr,0,"iparm");
  iparm = (int*) mxGetPr(ptr);

  // Get the "dparm" miscellaneous information structures.
  ptr   = mxGetField(matlabptr,0,"dparm");
  dparm = mxGetPr(ptr);
}

void PardisoInfo::checkMatrix (const SparseMatrix& A) const {

  // Check to make sure a sparse symmetric matrix is upper triangular.
  if (useSymmetricMatrices(*this) && !isUpperTriangular(A)) 
    mexErrMsgTxt("You must only provide the lower triangular portion of \
a symmetric matrix");

  // Check to make sure a sparse symmetric matrix has all its diagonal
  // elements included.
  if (useSymmetricMatrices(*this) && !diagonalIsPresent(A))
    mexErrMsgTxt("All the diagonal entries must be included in a \
sparse symmetric matrix");
}


void PardisoInfo::reorder (SparseMatrix& A, int* perm, bool verbose) {

  // Check whether the user has provided a permutation vector.
  if (perm)
    iparm[4] = 1;
  else
    iparm[4] = 0;

  // Call PARDISO.
  callpardiso(11,(int) verbose,A,perm);
}

void PardisoInfo::factor (SparseMatrix& A, bool verbose) {
  
  // Call PARDISO.
  callpardiso(22,(int) verbose,A);
}

void PardisoInfo::solve (SparseMatrix& A, MatlabMatrix& B, MatlabMatrix& X, 
			 bool verbose) {

  // Make sure we have a valid sparse matrix.
  checkMatrix(A);

  // Call PARDISO. Note that we need to solve the system using the
  // transpose of the sparse matrix A.
  int& t   = iparm[11];
  t        = 1 - t;
  iparm[5] = 0;  // Don't write the solution to X.
  callpardiso(33,(int) verbose,A,B,X);
  t = 1 - t;
}

void PardisoInfo::free() {
  callpardiso(-1,0);
}

void PardisoInfo::callpardiso (int phase, int msglvl) {
  int     error = 0;
  int     d     = 0; // Integer placeholder.
  complex x;         // Double placeholder.

  pardiso_(pt,&maxfct,&mnum,&mtype,&phase,&d,&x,&d,&d,&d,&d,
	   iparm,&msglvl,&x,&x,&error,dparm);
  if (error)
    reportPardisoError(error);
}

void PardisoInfo::callpardiso (int phase, int msglvl, SparseMatrix& A, 
			       int* perm) {
  int d     = 0;  // Integer placeholder.
  int error = 0;

  // Make sure we have a valid sparse matrix.
  checkMatrix(A);

  if (!A.iscomplex) {
    double x;  // Double placeholder.

    pardiso_(pt,&maxfct,&mnum,&mtype,&phase,
	     &A.n,A.ar,A.ia,A.ja,perm,&d,
	     iparm,&msglvl,&x,&x,&error,dparm);    
  }
  else {
    complex x;  // Complex placeholder.

    pardiso_(pt,&maxfct,&mnum,&mtype,&phase,
	     &A.n,A.ac,A.ia,A.ja,perm,&d,
	     iparm,&msglvl,&x,&x,&error,dparm);        
  }

  if (error)
    reportPardisoError(error);
}

void PardisoInfo::callpardiso (int phase, int msglvl, SparseMatrix& A, 
			       MatlabMatrix& B, MatlabMatrix& X) {
  int d     = 0;  // Integer placeholder.
  int error = 0;

  // Make sure we have a valid sparse matrix.
  checkMatrix(A);

  // Call PARDISO.
  if (!A.iscomplex)
    pardiso_(pt,&maxfct,&mnum,&mtype,&phase,
	     &A.n,A.ar,A.ia,A.ja,&d,&B.w,
	     iparm,&msglvl,B.ar,X.ar,&error,dparm);
  else
    pardiso_(pt,&maxfct,&mnum,&mtype,&phase,
	     &A.n,A.ac,A.ia,A.ja,&d,&B.w,
	     iparm,&msglvl,B.ac,X.ac,&error,dparm);

  if (error)
    reportPardisoError(error);
}

// Function definitions for friends of class PardisoInfo.
// -----------------------------------------------------------------
bool useIterativeSolver (const PardisoInfo& info) { 
  return info.iparm[31]; 
}

bool useConjTranspose (const PardisoInfo& info) {
  return info.mtype == cmpx_hermitian_pdef || 
         info.mtype == cmpx_hermitian_indef;
}

bool useComplexNumbers (const PardisoInfo& info) {
  return info.mtype == cmpx_structsym       || 
         info.mtype == cmpx_hermitian_pdef  || 
         info.mtype == cmpx_hermitian_indef || 
         info.mtype == cmpx_symmetric       || 
         info.mtype == cmpx_nonsym;
}

bool useSymmetricMatrices (const PardisoInfo& info) {
  return (info.mtype % 2) == 0;
}
