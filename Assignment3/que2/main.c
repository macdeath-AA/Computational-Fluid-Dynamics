#include "utilities.h"
#include "functions.h"
#include "linsolvers.h"
#include<stdlib.h>

int main()
{

  int nx, ny;
  double *xf, *xc, *dxf, *dxc, *yf, *yc, *dyc, *dyf;
  double ** aP, **aE, **aW, **aN, **aS, **b, **Sp;
  double **T, **Tex, **kdiff, **wrk1, **wrk2, **rho, **u, **v;
  double xst, xen, yst, yen, Pe;
  double *Tleft, *Tright, *Ttop, *Tbot;
  double *kleft, *kright, *ktop, *kbot;
  double *Fleft, *Fright, *Ftop, *Fbot;
  int i, j, max_iter;
  double relax_T, tol, l2err;
  FILE* fp;  

  // read inputs
  fp = fopen("input.in", "r");
  fscanf(fp, "%d %d\n", &nx, &ny);
  fscanf(fp, "%lf %lf\n", &xst, &xen);
  fscanf(fp, "%lf %lf\n", &yst, &yen);
  fscanf(fp, "%lf %lf %lf %lf\n", &rho_const, &u_const, &v_const, &Pe);
  fscanf(fp, "%lf %lf\n", &T_high, &T_low);
  fclose(fp);

  printf("Inputs are: %d %d\n %lf %lf\n %lf %lf\n %lf %lf %lf %lf\n %lf %lf\n", nx, ny, xst, xen, yst, yen, rho_const, u_const, v_const, Pe, T_high, T_low);
  // calculate k_const
  k_const = rho_const*u_const*(xen-xst)/(double)nx/Pe;
  printf("Calculated k_const = %lf", k_const);

  // allocate memory
  printf("\n > Allocating Memory -- \n");
  xc  = (double *)malloc( nx    * sizeof(double));   // CV centers (interiors)
  xf =  (double *)malloc((nx + 1) * sizeof(double));  // CV faces
  dxc = (double *)malloc((nx + 1) * sizeof(double)); // spacing betw centers (interiors)
  dxf = (double *)malloc(nx * sizeof(double));       // spacing betw faces

  yc =  (double *)malloc(ny * sizeof(double));        // CV centers (interiors)
  yf =  (double *)malloc((ny + 1) * sizeof(double));  // CV faces
  dyc = (double *)malloc((ny + 1) * sizeof(double)); // spacing betw centers (interiors)
  dyf = (double *)malloc(ny * sizeof(double));

  Tleft = (double *)malloc(ny * sizeof(double));  // left   Dirichlet condition
  Tright = (double *)malloc(ny * sizeof(double)); // right  Dirichlet condition
  Ttop = (double *)malloc(nx * sizeof(double));   // top    Dirichlet condition
  Tbot = (double *)malloc(nx * sizeof(double));   // bottom Dirchlet   condition

  Fleft = (double *)malloc(ny * sizeof(double));     // left   diffusivity array
  Fright = (double *)malloc(ny * sizeof(double));    // right  diffusivity array
  Ftop = (double *)malloc(nx * sizeof(double));      // top    diffusivity array
  Fbot = (double *)malloc(nx * sizeof(double));      // bottom diffusivity array

  kleft = (double *)malloc(ny * sizeof(double));  // left   diffusivity array
  kright = (double *)malloc(ny * sizeof(double)); // right  diffusivity array
  ktop = (double *)malloc(nx * sizeof(double));   // top    diffusivity array
  kbot = (double *)malloc(nx * sizeof(double));   // bottom diffusivity array
  printf("   >> Done allocating 1D arrays -- \n");


  // allocate 2D arrays dynamically
 
  T = (double **)malloc(nx*sizeof(double *));
  Tex = (double **)malloc(nx * sizeof(double *));
  rho = (double **)malloc(nx * sizeof(double *));
  u = (double **)malloc(nx * sizeof(double *));
  v = (double **)malloc(nx * sizeof(double *));
  aP = (double **)malloc(nx * sizeof(double *));
  aE = (double **)malloc(nx * sizeof(double *));
  aW = (double **)malloc(nx * sizeof(double *));
  aN = (double **)malloc(nx * sizeof(double *));
  aS = (double **)malloc(nx * sizeof(double *));
  Sp = (double **)malloc(nx * sizeof(double *));
  kdiff = (double **)malloc(nx * sizeof(double *));
  b = (double **)malloc(nx * sizeof(double *));

  for(i=0; i<nx; i++){
    T[i] = (double *)malloc(ny * sizeof(double));
    Tex[i] = (double *)malloc(ny * sizeof(double));
    rho[i] = (double *)malloc(ny * sizeof(double));
    u[i] = (double *)malloc(ny * sizeof(double));
    v[i] = (double *)malloc(ny * sizeof(double));
    aP[i] = (double *)malloc(ny * sizeof(double));
    aE[i] = (double *)malloc(ny * sizeof(double));
    aW[i] = (double *)malloc(ny * sizeof(double));
    aN[i] = (double *)malloc(ny * sizeof(double));
    aS[i] = (double *)malloc(ny * sizeof(double));
    Sp[i] = (double *)malloc(ny * sizeof(double));
    kdiff[i] = (double *)malloc(ny * sizeof(double));
    b[i] = (double *)malloc(ny * sizeof(double));
  }

  // -- for work arrays wrk1, wrk2 --
  wrk1 = (double **)malloc((nx+2)*sizeof(double *));
  for(i=0; i<nx+2; i++)
    wrk1[i] = (double *)malloc((ny+2)*sizeof(double));

  wrk2 = (double **)malloc((nx+2)*sizeof(double *));
  for(i=0; i<nx+2; i++)
    wrk2[i] = (double *)malloc((ny+2)*sizeof(double));


  printf("   >> Done allocating 2D arrays -- \n");
  printf(" > Done allocating memory -------- \n");
  
  // initialize the grid
  grid(nx, xst, xen, xc, xf, dxc, dxf);  // -- along x --
  grid(ny, yst, yen, yc, yf, dyc, dyf);  // -- along y --
  printf("\n > Done setting up grid ---------- \n");

  set_initial_guess(nx,ny,xc,yc,T,rho,u,v);  // initial condition
  printf("\n > Done setting up initial guess -- \n");

  get_fv_coeffs(nx, ny, xc, yc, xf, yf, dxc, dxf, dyc, dyf,            // grid vars
                aP, aE, aW, aN, aS, b, Sp, kdiff, T, rho, u, v,        // coefficients
                Tleft, Tright, Tbot, Ttop, 
                Fleft, Fright, Fbot, Ftop, 
                kleft, kright, kbot, ktop); // BC vars
  printf("\n > Done calculating fv coeffs ----- \n");

  printf("\n > Solving for T ------------- \n");
  max_iter = 100000; tol = 1.0e-15; relax_T = 1.0;
  solve_gssor(nx, ny, aP, aE, aW, aN, aS, b, T, wrk1, wrk2, max_iter, tol, relax_T);
  printf(" > Done solving for T ------------- \n");

  get_exact_soln(nx, ny, xc, yc, Tex);
  char str1[20] = "numer";  output_soln(nx, ny, 0, xc, yc, T, str1);
  char str2[20] = "exact";  output_soln(nx, ny, 0, xc, yc, Tex, str2);

  l2err = get_l2err_norm(nx, ny, T, Tex);
  printf("%d %d %9.5e", nx, ny, l2err);

  // free memory
   // ----1D arrays ---
  free(yf);
  free(yc);
  free(dyf);
  free(dyc);
  free(xf);
  free(xc);
  free(dxf);
  free(dxc);
  free(Tleft);
  free(Tright);
  free(Ttop);
  free(Tbot);
  free(Fleft);
  free(Fright);
  free(Ftop);
  free(Fbot);
  free(kleft);
  free(kright);
  free(ktop);
  free(kbot);

  // --- Done 1D arrays ---

  // ----2D arrays ---
  for (i = 0; i < nx; i++){
    free(rho[i]);
    free(u[i]);
    free(v[i]);
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

  free(rho);
  free(u);
  free(v);
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

  // --- Done 2D arrays ---
  printf("\n > Done freeing up memory --------- \n");
  return 0;
}

