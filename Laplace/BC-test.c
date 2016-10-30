#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define MAX(a, b) ((a) > (b) ? (a) : (b))

double BCfunc(const double coord) 
{
  return cos(M_PI*coord);
}

int approxEqual(const double a, const double b)
{
  double ave = (a + b) * 0.5;
  double tol = MAX(ave*1.e-10, 1.e-10);
  return (fabs(a-b) < tol ? 1 : 0);
}

int BCtest(const double coord, const double interiorValue,
	    const double ghostValue)
{
  double wallValue = BCfunc(coord);
  double aveValue = (interiorValue + ghostValue) * 0.5;
  return approxEqual(wallValue, aveValue);
}

#define IMAX 20
#define JMAX 20

int testAllBCs(double soln[IMAX+2][JMAX+2]) 
{
  double dx = 1./IMAX;
  int i;
  int OK = 1;
  for (i = 1; i <= IMAX && OK; i++) {
    double x = (i-0.5) * dx;
    OK = BCtest(x, soln[i][JMAX], soln[i][JMAX+1]);
  }

  if (i == IMAX+1 && OK == 1) {
    printf("BC test passed\n");
    return 1;
  }
  else {
    printf("BC test failed\n");
    return 0;
  }
}

int main() {
  double soln[IMAX+2][JMAX+2];
  double dx = 1./IMAX;
  double dy = 1./JMAX;
  int i, j;
  
  for (i = 1; i <= IMAX; i++) {
    for (j = 1; j <= JMAX; j++) {
      soln[i][j] = 0;
    }
  }

  for (i = 1; i <= IMAX; i++) {
    double x = (i-0.5) * dx;
    soln[i][JMAX+1] = 2*BCfunc(x) - soln[i][JMAX];
  }

  testAllBCs(soln);

  exit(0);
}
