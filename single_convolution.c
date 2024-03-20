#include <stdio.h>
#include <math.h>

#include "pan_tompkins.h"

dataType single_convolution(dataType *x, int n, dataType *h, int nc)
{
  int k = 0;

  dataType result = 0;
  
  for (k = 0; k < MIN(n, nc); k++)  // It updates the buffersize with #define MIN
  {
    result += x[n-k] * h[k];
  } 
  return result;
}
