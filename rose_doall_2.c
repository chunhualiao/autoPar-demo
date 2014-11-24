#include "omp.h" 
int i;
int j;
int a[100][100];

void foo()
{
  
#pragma omp parallel for private (i,j)
  for (i = 0; i <= 99; i += 1) {
    
#pragma omp parallel for private (j)
    for (j = 0; j <= 99; j += 1) {
      a[i][j] = a[i][j] + 1;
    }
  }
}
