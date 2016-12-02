/* This definition must be changed if you want to use a mesh with more than 
   100 control volumes. */ 
#define NMAX 12
#define maxx 12
#define maxy 12


struct carray{
//arrays

double f1 [maxx][maxy];//second stage solution mesh
double T1 [maxx][maxy];//first stage and solution mesh
double v1 [maxx][maxy];//first stage and solution mesh
double u1 [maxx][maxy];
//array attributes
int sizex = maxx;
int sizey = maxy;
double DIMx = 0.0;
double DIMy = 0.0;
//scheme
int scheme = 0;
double ctime = 0;
};

struct crow{
double LHS [maxx][3];
double RHS [maxx];
};
/* Solve a tri-diagonal system Ax = b.
 * 
 * Uses the Thomas algorithm, which is Gauss elimination and back
 * substitution specialized for a tri-diagonal matrix.
 *
 * Input:
 *  LHS: The matrix A.  The three columns are the three non-zero diagonals.
 *     (0: below main diag; 1: on main diag; 2: above main diag)
 *  RHS: A vector containing b.
 *  iSize: Number of equations.  For situations with BC's, those are included
 *     in the count.
 *
 * Output:
 *  LHS: Garbled.
 *  RHS: The solution x.
 */
static void SolveThomas(crow & r,
			const int iSize)
{
  int i;
  /* This next line actually has no effect, but it -does- make clear that
     the values in those locations have no impact. */
  r.LHS[0][0] = r.LHS[iSize-1][2] = 0;
  /* Forward elimination */
  for (i = 0; i < iSize-1; i++) {
    r.LHS[i][2] /= r.LHS[i][1];
    r.RHS[i] /= r.LHS[i][1];
    r.LHS[i+1][1] -= r.LHS[i][2]*r.LHS[i+1][0];
    r.RHS[i+1] -= r.LHS[i+1][0]*r.RHS[i];
  }
  /* Last line of elimination */
  r.RHS[iSize-1] /= r.LHS[iSize-1][1];

  /* Back-substitution */
  for (i = iSize-2; i >= 0; i--) {
    r.RHS[i] -= r.RHS[i+1]*r.LHS[i][2];
  }
}

/* A test program to confirm proper behavior. */
#include <stdio.h>
int main()
{

  crow r;
  //double LHS[12][3], RHS[12];
  int i;
  
  for (i = 0; i <= 11; i++) {
    r.LHS[i][0] = 1;
    r.LHS[i][1] = 2+i;
    r.LHS[i][2] = 3;
    r.RHS[i]   = i;
  }

  SolveThomas(r, 12);
      
  for (i = 0; i <= 11; i++) {
    double result;
    if (i == 0) result = RHS[0]*2 + RHS[1]*3;
    else if (i == 11) result = RHS[10]*1 + RHS[11]*(2+11);
    else result = RHS[i-1]*1 + RHS[i]*(2+i) + RHS[i+1]*3;
    printf("%2d Soln: %10.6f RHS recomputed: %10.6f  Error: %10.6G\n", i, RHS[i], result, result - i);
  }

  return(0);
}
