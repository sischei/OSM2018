#include "common.h"

// Function definitions.
// -----------------------------------------------------------------
bool mxIsDoubleScalar (const mxArray* ptr) {
  return mxIsDouble(ptr) && mxGetNumberOfElements(ptr);
}
