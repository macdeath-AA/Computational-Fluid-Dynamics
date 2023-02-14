#include "functions.h"
#include "linsolvers.h"
#include "utilities.h"
#include<stdlib.h>

int main()
{

  int nx, ny;
  double *xf, *xc, *dxf, *dxc, *yf, *yc, *dyc, *dyf;
  double ** aP, **aE, **aW, **aN, **aS, **b, **Sp;
  double **T, **Tex, **kdiff, **wrk1, **wrk2;
  double xst, xen, yst, yen;
  double *Tleft, *Tright, *Ttop, *qbot;
  int i, j, max_iter;
  double relax_T, tol, l2err;
  FILE* fp;  

  // read inputs
  fp = fopen("input.in", "r");
  fscanf(fp, "%d %d\n", &nx, &ny);
  fscanf(fp, "%lf %lf\n", &xst, &xen);
  fscanf(fp, "%lf %lf\n", &yst, &yen);
  fclose(fp);

  printf("Inputs are: %d %d %lf %lf %lf %lf\n", nx, ny, xst, xen, yst, yen);

  // allocate memory
  printf("\n > Allocating Memory -- \n");
  xc  = (double *)malloc( nx * sizeof(double));      // CV centers (interiors)
  xf  = (double *)malloc( (nx+1) * sizeof(double));      // CV faces (interiors)
  dxc = (double *)malloc( (nx+1) * sizeof(double));      
  dxf = (double *)malloc( nx * sizeof(double));      
  

  yc  = (double *)malloc( ny    * sizeof(double));   // CV centers (interiors)
  yf = (double *)malloc((ny+1) * sizeof(double));        
  dyc = (double *)malloc((ny+1) * sizeof(double));        
  dyf = (double *)malloc(ny * sizeof(double));       
  
  Tleft  = (double *)malloc( ny * sizeof(double));   // left   Dirichlet condition
  Tright = (double *)malloc( ny * sizeof(double));   // right  Dirichlet condition
  Ttop   = (double *)malloc( nx * sizeof(double));   // top    Dirichlet condition
  qbot   = (double *)malloc( nx * sizeof(double));   // bottom Neumann   condition
  printf("   >> Done allocating 1D arrays -- \n");


  // allocate 2D arrays dynamically
  // -- for T --
  T = (double **)malloc(nx*sizeof(double *));
  Tex = (double **)malloc(nx * sizeof(double *));
  aP = (double **)malloc(nx * sizeof(double *));
  aE = (double **)malloc(nx * sizeof(double *));
  aW = (double **)malloc(nx * sizeof(double *));
  aN = (double **)malloc(nx * sizeof(double *));
  aS = (double **)malloc(nx * sizeof(double *));
  b = (double **)malloc(nx * sizeof(double *));
  Sp = (double **)malloc(nx * sizeof(double *));
  kdiff = (double **)malloc(nx * sizeof(double *));

  for(i=0; i<nx; i++)
    T[i] = (double *)malloc(ny*sizeof(double));
    Tex[i] = (double *)malloc(ny * sizeof(double));
    aP[i] = (double *)malloc(ny * sizeof(double));
    aE[i] = (double *)malloc(ny * sizeof(double));
    aW[i] = (double *)malloc(ny * sizeof(double));
    aN[i] = (double *)malloc(ny * sizeof(double));
    aS[i] = (double *)malloc(ny * sizeof(double));
    b[i] = (double *)malloc(ny * sizeof(double));
    Sp[i] = (double *)malloc(ny * sizeof(double));
    kdiff[i] = (double *)malloc(ny * sizeof(double));

  wrk1 = (double **)malloc((nx+2) * sizeof(double *));
  wrk2 = (double **)malloc((nx + 2) * sizeof(double *));

  for (i = 0; i < (nx + 2); i++)
    wrk1[i] = (double *)malloc((ny+2) * sizeof(double));
    wrk2[i] = (double *)malloc((ny + 2) * sizeof(double));

  // -- for work arrays wrk1, wrk2 :: note : these should be padded arrays--

  printf("   >> Done allocating 2D arrays -- \n");
  printf(" > Done allocating memory -------- \n");

  // initialize the grid
  grid(nx, xst, xen, xc, xf, dxc, dxf); // -- along x --
  grid(ny, yst, yen, yc, yf, dyc, dyf); // -- along y --
  printf("\n > Done setting up grid ---------- \n");

  set_initial_guess(nx, ny, xc, yc, T); // initial condition
  printf("\n > Done setting up initial guess -- \n");

  get_fv_coeffs(nx, ny, xc, yc, xf, yf, dxc, dxf, dyc, dyf, // grid vars
                aP, aE, aW, aN, aS, b, Sp, kdiff, T,        // coefficients
                Tleft, Tright, qbot, Ttop);                 // BC vars
  printf("\n > Done calculating fv coeffs ----- \n");

  printf("\n > Solving for T ------------- \n");
  max_iter = 100000;
  tol = 1.0e-10;
  relax_T = 1.0;
  solve_gssor(nx, ny, aP, aE, aW, aN, aS, b, T, wrk1, wrk2, max_iter, tol, relax_T);
  printf(" > Done solving for T ------------- \n");

  get_exact_soln(nx, ny, xc, yc, Tex);
  char str1[20] = "numer";
  output_soln(nx, ny, 1, xc, yc, T, str1);
  char str2[20] = "exact";
  output_soln(nx, ny, 1, xc, yc, Tex, str2);

  l2err = get_l2err_norm(nx, ny, T, Tex);
  printf("%d %d %9.5e", nx, ny, l2err);

  // free memory
  // ----1D arrays ---
  free(xc);
  free(xf);
  free(dxc);
  free(dxf);
  free(yc);
  free(yc);
  free(dyc);
  free(dyf);
  free(Tleft);
  free(Tright);
  free(Ttop);
  free(qbot);

  // --- Done 1D arrays ---

  // ----2D arrays ---
  for (i = 0; i < nx; i++)
  {
    free(T[i]);
    free(Tex[i]);  
    free(aP[i]);
    free(aE[i]);
    free(aW[i]);
    free(aN[i]);
    free(aS[i]);
    free(b[i]);
    free(Sp[i]);
    free(kdiff[i]);
  }

  free(T);
  free(Tex);
  free(aP);
  free(aE);
  free(aW);
  free(aN);
  free(aS);
  free(b);
  free(Sp);
  free(kdiff);

  for (i = 0; i < nx + 2; i++)
  {  
    free(wrk1[i]); 
    free(wrk2[i]);

  }
  free(wrk1);
  free(wrk2);
  // --- Done 2D arrays ---
  printf("\n > Done freeing up memory --------- \n");
  return 0;
}

