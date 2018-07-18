#ifndef INCLUDE_COMMON
#define INCLUDE_COMMON

#include "mex.h"
#include <string.h>

// Function declarations.
// -----------------------------------------------------------------
// Return true if the MATLAB array is a scalar in double precision.
bool mxIsDoubleScalar (const mxArray* ptr);

// Function definitions.
// -----------------------------------------------------------------
template <class Type> void copymemory (const Type* source, Type* dest, 
				       int n) {
  memcpy(dest,source,sizeof(Type) * n);
}

#endif
