#include <stdio.h>
#include <math.h>

#define IMAX 2
#define JMAX 2

#define MAX(a, b) ((a) > (b) ? (a) : (b))

double L2Norm(double array[IMAX+2][JMAX+2]) 
{
  double sum = 0;
  int i, j;
  for (i = 1; i <= IMAX; i++) {
    for (j = 1; j <= JMAX; j++) {
      sum += array[i][j]*array[i][j];
    }
  }
  return sqrt(sum/(IMAX*JMAX));
}

int approxEqual(const double a, const double b)
{
  double ave = (a + b) * 0.5;
  double tol = MAX(ave*1.e-10, 1.e-10);
  return (fabs(a-b) < tol ? 1 : 0);
}

int main()
{
  double soln[IMAX+2][JMAX+2], result;

  soln[1][1] = 1;
  soln[2][1] = 3;
  soln[1][2] = 3;
  soln[2][2] = 9;

  result = L2Norm(soln);

  printf("L2 norm is %f\n", result);
  if (approxEqual(result, 5)) {
    printf("Success!\n");
  }
  else {
    printf("Failure!\n");
  }
  return(0);
}
