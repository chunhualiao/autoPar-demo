//#include <A++.h>
#include "simpleA++.h"
#include "omp.h" 

void interpolate1D(class floatArray &fineGrid,class floatArray &coarseGrid)
{
  int _var_0;
  int i;
  int _var_1;
  int fineGridSize = fineGrid .  length (0);
  int coarseGridSize = coarseGrid .  length (0);
// Interior fine points
  class Range If(2,_var_1 = fineGridSize - 2,2);
  class Range Ic(1,coarseGridSize - 1,1);
#if 0
#else  
  
#pragma omp parallel for private (i) firstprivate (_var_0)
  for (i = 1; i <= _var_0 - 1; i += 1) {
    fineGrid .  elem (i) = fineGrid .  elem (i) + 1;
  }
#endif  
}
