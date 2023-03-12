#include "functions.h"

void grid(int nx, double xst, double xen, double *xc, double *xf, double *dxc, double *dxf)
{
  int i;
  double dxunif;

  // uniform mesh for now; 
  // can use stretching factors to place xc nodes later
  dxunif = (xen-xst)/(double)nx;

  // xc[i] s are inside the CVs
  for(i=0; i<nx; i++)
    xc[i] = ((double)i + 0.5) * dxunif;

  // xf[i] s are at the faces; mid-way betwee xc[i]s
  for (i = 0; i <= nx; i++)
  {
    xf[i] = (double)i * dxunif;
  }

  // dxc[i] s are spacing between adjacent xcs; ends are different 
  dxc[0] = xc[0]-xf[0];
  for (i = 1; i < nx; i++)
  {
    dxc[i] = xc[i] - xc[i - 1];
  }
  dxc[nx] = xf[nx] - xc[nx - 1];

  // dxf[i] s are spacing between adjacent xes
  for (i = 0; i < nx; i++)
  {
    dxf[i] = xf[i + 1] - xf[i];
  }

}

void set_initial_guess(int nx, int ny, double *xc, double *yc, double **T, double **rho, double **u, double **v)
{
  int i, j;

  // set rho to be rho_const everywhere
  for(i=0; i<nx; i++)
    for(j=0; j<ny; j++)
      rho[i][j] = rho_const;

  // set u to be u_const everywhere
  for(i=0; i<nx; i++)
    for(j=0; j<ny; j++)
      u[i][j] = u_const;

  // set v to be v_const everywhere
  for(i=0; i<nx; i++)
    for(j=0; j<ny; j++)
      v[i][j] = v_const;

  // set T to be 0 everywhere
  for(i=0; i<nx; i++)
    for(j=0; j<ny; j++)
      T[i][j] = 0.0;
}

void calc_diffusivity(int nx, int ny, double *xc, double *yc, double **T, double **kdiff, double *kleft, double *kright, double *ktop, double *kbot)
{
  int i, j;

  // calculate diffusivity (may be dependent on T)
  for(i=0; i<nx; i++)
   for(j=0; j<ny; j++)
   {
     kdiff[i][j] = k_const; // set to uniform for now
   }

  // left boundary -- Non-homogeneous Dirichlet
  for(j=0; j<ny; j++)
    kleft[j] = k_const;

  // right boundary -- Homogeneous Dirichlet
  for(j=0; j<ny; j++)
    kright[j] = k_const;

  // top boundary -- Non-homogeneous Dirichlet
  for(i=0; i<nx; i++)
    ktop[i] = k_const;

  // bottom boundary -- Homogeneous Dirichlet
  for(i=0; i<nx; i++)
    kbot[i] = k_const;
}

void calc_sources(int nx, int ny, double *xc, double *yc, double *dxf, double *dyf, double **T, double **b, double **Sp)
{
  int i, j;
  double AA = 1.0, BB = 0.8; 

  // calculate Su source (may be dependent on T)
  for(i=0; i<nx; i++)
   for(j=0; j<ny; j++)
   {
     // set source function here
     b[i][j]  = 0.0;

     // multiply by volume of CV
     b[i][j] *= dxf[i]*dyf[j];
   }

  // calculate Sp source (may be dependent on T)
  for(i=0; i<nx; i++)
   for(j=0; j<ny; j++)
   {
     Sp[i][j]  = 0.0;           // set source function here
     Sp[i][j] *= dxf[i]*dyf[j]; // multiply by volume of CV
   }
}

void set_boundary_conditions(int nx, int ny, double *xc, double *yc, double *xf, double *yf, double *dxc, double *dyc, double *dxf, double *dyf, double **T, double *Tleft, double *Tright, double *Ttop, double *Tbot, double *Fleft, double *Fright, double *Ftop, double *Fbot)
{

  int i, j;
  double bcE, bcW, bcN, bcS, kw, ke, kn, ks;

  // left boundary -- Non-homogeneous Dirichlet
  for(j=0; j<ny; j++)
    Tleft[j] = T_high;

  // right boundary -- Homogeneous Dirichlet
  for(j=0; j<ny; j++)
    Tright[j] = T_low;

  // top boundary -- Non-homogeneous Dirichlet
  for(i=0; i<nx; i++)
    Ttop[i] = T_high;

  // bottom boundary -- Homogeneous Dirichlet
  for(i=0; i<nx; i++)
    Tbot[i] = T_low;

  // now set convective fluxes at the boundaries
  // rho = rho_const everywhere; u = u_const; v = v_const everywhere
  // left boundary
  for(j=0; j<ny; j++)
    Fleft[j] = rho_const * u_const; 

  // right boundary
  for(j=0; j<ny; j++)
    Fright[j] = rho_const * u_const;

  // top boundary
  for(i=0; i<nx; i++)
    Ftop[i] = rho_const * v_const;

  // bottom boundary
  for(i=0; i<nx; i++)
    Fbot[i] = rho_const * v_const;
}

void get_exact_soln(int nx, int ny, double *xc, double *yc, double **Tex)
{
  int i, j;

  // note :: this exact solution is valid for Pe = infinity
  for(i=0; i<nx; i++)
   for(j=0; j<ny; j++)
   {
     if(fabs(yc[j] - xc[i]) < 1.0e-12)
         Tex[i][j] =  0.5*(T_high + T_low);   // on the diagonal
     else if(yc[j] > xc[i])
         Tex[i][j] = T_high;   // above diagonal
     else if(yc[j] < xc[i])
         Tex[i][j] =   T_low;   // below diagonal
   }
}

void get_fv_coeffs(int nx, int ny, double *xc, double *yc, double *xf, double *yf, double *dxc, double *dxf, double *dyc, double *dyf, double **aP, double **aE, double **aW, double **aN, double **aS, double **b, double **Sp, double **kdiff, double **T, double **rho, double **u, double **v, double *Tleft, double *Tright, double *Tbot, double *Ttop, double *Fleft, double *Fright, double *Fbot, double *Ftop, double *kleft, double *kright, double *kbot, double *ktop)
{

  int i, j;
  double bcE, bcW, bcN, bcS, kw, ke, kn, ks, ifac, Fe, Fw, Fn, Fs;
  double De, Dw, Dn, Ds, rhoe, rhow, rhon, rhos, ue, uw, vn, vs;
  int convec_scheme;

  // calculate diffusivity at [xc, yc]; may be dependent on T
  calc_diffusivity(nx, ny, xc, yc, T, kdiff, kleft, kright, ktop, kbot);

  // calculate sources Su, Sp at [xc, yc]; may be dependent on T
  calc_sources(nx, ny, xc, yc, dxf, dyf, T, b, Sp);

  // populate values in BC arrays
  set_boundary_conditions(nx, ny, xc, yc, xf, yf, dxc, dyc, dxf, dyf, 
             T, Tleft, Tright, Ttop, Tbot, Fleft, Fright, Ftop, Fbot);

  // start populating the coefficients

  // ------ Step 1 :: interior points ------
  for(i=1; i<nx-1; i++)
   for(j=1; j<ny-1; j++)
   {
      // east-west -- diffusion
      ke = 0.5*(kdiff[i][j] + kdiff[i+1][j]);  De = ke * dyf[j] / dxc[i+1];
      kw = 0.5*(kdiff[i][j] + kdiff[i-1][j]);  Dw = kw * dyf[j] / dxc[i];

      // north-south -- diffusion
      kn = 0.5*(kdiff[i][j] + kdiff[i][j+1]);  Dn = kn * dxf[i] / dyc[j+1];
      ks = 0.5*(kdiff[i][j] + kdiff[i][j-1]);  Ds = ks * dxf[i] / dyc[j];

      // east-west -- convection
      rhoe = 0.5*(rho[i][j] + rho[i+1][j]);  ue = 0.5*(u[i][j] + u[i+1][j]);  Fe = rhoe * ue * dyf[j];
      rhow = 0.5*(rho[i][j] + rho[i-1][j]);  uw = 0.5*(u[i][j] + u[i-1][j]);  Fw = rhow * ue * dyf[j];

      // north-south -- convection
      rhon = 0.5*(rho[i][j] + rho[i][j+1]);  vn = 0.5*(v[i][j] + v[i][j+1]);  Fn = rhon * vn * dxf[i];
      rhos = 0.5*(rho[i][j] + rho[i][j-1]);  vs = 0.5*(v[i][j] + v[i][j-1]);  Fs = rhos * vs * dxf[i];

      #include "convective_scheme_coeffs.h"

      // source --  has already been populated
   }

  // ------ Step 2  :: 1-boundary points ------
  // ------ Step 2a :: left boundary ----------
  i = 0;
   for(j=1; j<ny-1; j++)
   {
      // east-west
      ke = 0.5*(kdiff[i][j] + kdiff[i+1][j]);  De = ke * dyf[j] / dxc[i+1];
      kw = 0.0;                                Dw = 0.0;

      // north-south
      kn = 0.5*(kdiff[i][j] + kdiff[i][j+1]);  Dn = kn * dxf[i] / dyc[j+1];
      ks = 0.5*(kdiff[i][j] + kdiff[i][j-1]);  Ds = ks * dxf[i] / dyc[j];

      // east-west -- convection
      rhoe = 0.5*(rho[i][j] + rho[i+1][j]);  ue = 0.5*(u[i][j] + u[i+1][j]);  Fe = rhoe * ue * dyf[j];
      rhow = 0.0;                            uw = 0.0;                        Fw = 0.0;

      // north-south -- convection
      rhon = 0.5*(rho[i][j] + rho[i][j+1]);  vn = 0.5*(v[i][j] + v[i][j+1]);  Fn = rhon * vn * dxf[i];
      rhos = 0.5*(rho[i][j] + rho[i][j-1]);  vs = 0.5*(v[i][j] + v[i][j-1]);  Fs = rhos * vs * dxf[i];

      #include "convective_scheme_coeffs.h"

      // Dirichlet BC on the left; 
       // note Fleft is added to source sometimes, but is always subtracted from aP
       // set ifac, boundary coefficient (bcW) and modify b, aP as needed -- 
       ifac = 1.0; //default
       if(convec_scheme==1 && Fleft[j] < 0.0){ // UW
          ifac = 0.0;
       }
       bcW = (kleft[j]/dxc[i] + ifac * Fleft[j]) * dyf[j];
       b[i][j]  += bcW * Tleft[j];        aP[i][j] += bcW - Fleft[j]*dyf[j];
   }
  // ------ Step 2a :: left boundary done ---

  // ------ Step 2b :: right boundary ----------
  i = nx-1;
   for(j=1; j<ny-1; j++)
   {
      // east-west
       ke = 0.0;
       De = 0.0;
       kw = 0.5 * (kdiff[i][j] + kdiff[i - 1][j]);
       Dw = kw * dyf[j] / dxc[i];

       // north-south
       kn = 0.5 * (kdiff[i][j] + kdiff[i][j + 1]);
       Dn = kn * dxf[i] / dyc[j + 1];
       ks = 0.5 * (kdiff[i][j] + kdiff[i][j - 1]);
       Ds = ks * dxf[i] / dyc[j];

       // east-west -- convection
       rhoe = 0.0;
       ue = 0.0;
       Fe = 0.0;
       rhow = 0.5 * (rho[i][j] + rho[i - 1][j]);
       uw = 0.5 * (u[i][j] + u[i - 1][j]);
       Fw = rhow * uw * dyf[j];

       // north-south -- convection
       rhon = 0.5 * (rho[i][j] + rho[i][j + 1]);
       vn = 0.5 * (v[i][j] + v[i][j + 1]);
       Fn = rhon * vn * dxf[i];
       rhos = 0.5 * (rho[i][j] + rho[i][j - 1]);
       vs = 0.5 * (v[i][j] + v[i][j - 1]);
       Fs = rhos * vs * dxf[i];

#include "convective_scheme_coeffs.h"

      // Dirichlet BC on the right; 
       // note Fright is subtracted from source sometimes, but is always added to aP
       // set ifac, boundary coefficient (bcE) and modify b, aP as needed -- 
       ifac = 1.0; //default
       if (convec_scheme == 1 && Fright[j] > 0.0){ // UW
          ifac = 0.0;
       }
       bcE = (kright[j] / dxc[i] - ifac * Fright[j]) * dyf[j];
       b[i][j] += bcE * Tright[j];       aP[i][j] += bcE + Fright[j] * dyf[j];
   }
  // ------ Step 2b :: right boundary done ---

  // ------ Step 2c :: bottom boundary ----------
  j = 0;
  for(i=1; i<nx-1; i++)
   {
      // east-west
       ke = 0.5 * (kdiff[i][j] + kdiff[i + 1][j]);
       De = ke * dyf[j] / dxc[i + 1];
       kw = 0.5 * (kdiff[i][j] + kdiff[i - 1][j]);
       Dw = kw * dyf[j] / dxc[i];

       // north-south
       kn = 0.5 * (kdiff[i][j] + kdiff[i][j + 1]);
       Dn = kn * dxf[i] / dyc[j + 1];
       ks = 0.0;
       Ds = 0.0;

       // east-west -- convection
       rhoe = 0.5 * (rho[i][j] + rho[i + 1][j]);
       ue = 0.5 * (u[i][j] + u[i + 1][j]);
       Fe = rhoe * ue * dyf[j];
       rhow = 0.5 * (rho[i][j] + rho[i - 1][j]);
       uw = 0.5 * (u[i][j] + u[i - 1][j]);
       Fw = rhow * uw * dyf[j];

       // north-south -- convection
       rhon = 0.5 * (rho[i][j] + rho[i][j + 1]);
       vn = 0.5 * (v[i][j] + v[i][j + 1]);
       Fn = rhon * vn * dxf[i];
       rhos = 0.0;
       vs = 0.0;
       Fs = 0.0;

#include "convective_scheme_coeffs.h"

      // Dirichlet BC on the bottom; 
       // note Fbot is added to source sometimes, but is always subtracted from aP
       // set ifac, boundary coefficient (bcS) and modify b, aP as needed -- 
       ifac = 1.0; //default
       if (convec_scheme == 1 && Fbot[i] < 0.0){ // UW
          ifac = 0.0;
       }
       bcS = (kbot[i] / dyc[j] + ifac * Fbot[i]) * dxf[i];
       b[i][j] += bcS * Tbot[i];
       aP[i][j] += bcS - Fbot[i] * dxf[i];
   }
  // ------ Step 2c :: bottom boundary done ---

  // ------ Step 2d :: top boundary ----------
  j = ny-1;
  for(i=1; i<nx-1; i++)
   {
      // east-west
       ke = 0.5 * (kdiff[i][j] + kdiff[i + 1][j]);
       De = ke * dyf[j] / dxc[i + 1];
       kw = 0.5 * (kdiff[i][j] + kdiff[i - 1][j]);
       Dw = kw * dyf[j] / dxc[i];

       // north-south
       kn = 0.0;
       Dn = 0.0;
       ks = 0.5 * (kdiff[i][j] + kdiff[i][j - 1]);
       Ds = ks * dxf[i] / dyc[j];

       // east-west -- convection
       rhoe = 0.5 * (rho[i][j] + rho[i + 1][j]);
       ue = 0.5 * (u[i][j] + u[i + 1][j]);
       Fe = rhoe * ue * dyf[j];
       rhow = 0.5 * (rho[i][j] + rho[i - 1][j]);
       uw = 0.5 * (u[i][j] + u[i - 1][j]);
       Fw = rhow * uw * dyf[j];

       // north-south -- convection
       rhon = 0.0;
       vn = 0.0;
       Fn = 0.0;
       rhos = 0.5 * (rho[i][j] + rho[i][j - 1]);
       vs = 0.5 * (v[i][j] + v[i][j - 1]);
       Fs = rhos * vs * dxf[i];

#include "convective_scheme_coeffs.h"

      // Dirichlet BC on the top; 
       // note Ftop is subtracted from source sometimes, but is always added to aP
       // set ifac, boundary coefficient (bcN) and modify b, aP as needed -- 
       ifac = 1.0; //default
       if(convec_scheme==1 && Ftop[i]>0.0) // UW
          ifac = 0.0;
       bcN = (ktop[i] / dyc[j] - ifac * Ftop[i]) * dxf[i];
       b[i][j] += bcN * Ttop[i];       aP[i][j] += bcN + Ftop[i] * dxf[i];
   }
  // ------ Step 2d :: top boundary done ---

  // ------ Step 3  :: 2-boundary points ------
  // ------ Step 3a :: top-left boundary ----------
  i = 0; j = ny-1;
      // east-west
  ke = 0.5 * (kdiff[i][j] + kdiff[i + 1][j]);
  De = ke * dyf[j] / dxc[i + 1];
  kw = 0.0;
  Dw = 0.0;

  // north-south
  kn = 0.0;
  Dn = 0.0;
  ks = 0.5 * (kdiff[i][j] + kdiff[i][j - 1]);
  Ds = ks * dxf[i] / dyc[j];

  // east-west -- convection
  rhoe = 0.5 * (rho[i][j] + rho[i + 1][j]);
  ue = 0.5 * (u[i][j] + u[i + 1][j]);
  Fe = rhoe * ue * dyf[j];
  rhow = 0.0;
  uw = 0.0;
  Fw = 0.0;

  // north-south -- convection
  rhon = 0.0;
  vn = 0.0;
  Fn = 0.0;
  rhos = 0.5 * (rho[i][j] + rho[i][j - 1]);
  vs = 0.5 * (v[i][j] + v[i][j - 1]);
  Fs = rhos * vs * dxf[i];

#include "convective_scheme_coeffs.h"

      // Dirichlet BC on the left; 
       ifac = 1.0; //default
       if(convec_scheme==1 && Fleft[j] < 0.0)
       { // UW
          ifac = 0.0;
       }
       bcW = (kleft[j] / dxc[i] + ifac * Fleft[j]) * dyf[j];
       b[i][j] += bcW * Tleft[j];
       aP[i][j] += bcW - Fleft[j] * dyf[j];

       // Dirichlet BC on the top; 
       ifac = 1.0; //default
       if (convec_scheme == 1 && Ftop[i] > 0.0)
       { // UW
          ifac = 0.0;
       }
       bcN = (ktop[i] / dyc[j] - ifac * Ftop[i]) * dxf[i];
       b[i][j] += bcN * Ttop[i];
       aP[i][j] += bcN + Ftop[i] * dxf[i];
       // ------ Step 3a :: top-left boundary done ---

       // ------ Step 3b :: top-right boundary ----------
       i = nx - 1;
       j = ny - 1;
       // east-west
       ke = 0.0;
       De = 0.0;
       kw = 0.5 * (kdiff[i][j] + kdiff[i - 1][j]);
       Dw = kw * dyf[j] / dxc[i];

       // north-south
       kn = 0.0;
       Dn = 0.0;
       ks = 0.5 * (kdiff[i][j] + kdiff[i][j - 1]);
       Ds = ks * dxf[i] / dyc[j];

       // east-west -- convection
       rhoe = 0.0;
       ue = 0.0;
       Fe = 0.0;
       rhow = 0.5 * (rho[i][j] + rho[i - 1][j]);
       uw = 0.5 * (u[i][j] + u[i - 1][j]);
       Fw = rhow * uw * dyf[j];

       // north-south -- convection
       rhon = 0.0;
       vn = 0.0;
       Fn = 0.0;
       rhos = 0.5 * (rho[i][j] + rho[i][j - 1]);
       vs = 0.5 * (v[i][j] + v[i][j - 1]);
       Fs = rhos * vs * dxf[i];

#include "convective_scheme_coeffs.h"

      // Dirichlet BC on the right; 
       ifac = 1.0; //default
       if (convec_scheme == 1 && Fright[j] > 0.0) 
       {// UW
          ifac = 0.0;
       }
       bcE = (kright[j] / dxc[i] - ifac * Fright[j]) * dyf[j];
       b[i][j] += bcE * Tright[j];
       aP[i][j] += bcE + Fright[j] * dyf[j];

       // Dirichlet BC on the top; 
       ifac = 1.0; //default
       if (convec_scheme == 1 && Ftop[i] > 0.0)
       { // UW
          ifac = 0.0;
       }
       bcN = (ktop[i] / dyc[j] - ifac * Ftop[i]) * dxf[i];
       b[i][j] += bcN * Ttop[i];
       aP[i][j] += bcN + Ftop[i] * dxf[i];
       // ------ Step 3b :: top-right boundary done ---

       // ------ Step 3c :: bottom-left boundary ----------
       i = 0;
       j = 0;
       // east-west
       ke = 0.5 * (kdiff[i][j] + kdiff[i + 1][j]);
       De = ke * dyf[j] / dxc[i + 1];
       kw = 0.0;
       Dw = 0.0;

       // north-south
       kn = 0.5 * (kdiff[i][j] + kdiff[i][j + 1]);
       Dn = kn * dxf[i] / dyc[j + 1];
       ks = 0.0;
       Ds = 0.0;

       // east-west -- convection
       rhoe = 0.5 * (rho[i][j] + rho[i + 1][j]);
       ue = 0.5 * (u[i][j] + u[i + 1][j]);
       Fe = rhoe * ue * dyf[j];
       rhow = 0.0;
       uw = 0.0;
       Fw = 0.0;

       // north-south -- convection
       rhon = 0.5 * (rho[i][j] + rho[i][j + 1]);
       vn = 0.5 * (v[i][j] + v[i][j + 1]);
       Fn = rhon * vn * dxf[i];
       rhos = 0.0;
       vs = 0.0;
       Fs = 0.0;

#include "convective_scheme_coeffs.h"

      // Dirichlet BC on the bottom; 
       ifac = 1.0; //default
       if (convec_scheme == 1 && Fbot[i] < 0.0)
       { // UW
          ifac = 0.0;
       }
       bcS = (kbot[i] / dyc[j] + ifac * Fbot[i]) * dxf[i];
       b[i][j] += bcS * Tbot[i];
       aP[i][j] += bcS - Fbot[i] * dxf[i];

       // Dirichlet BC on the left; 
       ifac = 1.0; //default
       if (convec_scheme == 1 && Fleft[j] < 0.0) // UW
          ifac = 0.0;
       
       bcW = (kleft[j] / dxc[i] + ifac * Fleft[j]) * dyf[j];
       b[i][j] += bcW * Tleft[j];
       aP[i][j] += bcW - Fleft[j] * dyf[j];

       //printf("%d %d %lf %lf %lf %lf %lf\n", i, j, kleft[j], dxc[i], Fleft[j], dyf[j], bcW);
       //printf("%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", i, j, aE[i][j], aW[i][j], aN[i][j], aS[i][j], Fe, Fw, Fn, Fs, bcW, bcS, aP[i][j], b[i][j]);

  // ------ Step 3c :: bottom-left boundary done ---

  // ------ Step 3d :: bottom-right boundary ----------
  i = nx-1;  j = 0;
      // east-west
  ke = 0.0;
  De = 0.0;
  kw = 0.5 * (kdiff[i][j] + kdiff[i - 1][j]);
  Dw = kw * dyf[j] / dxc[i];

  // north-south
  kn = 0.5 * (kdiff[i][j] + kdiff[i][j + 1]);
  Dn = kn * dxf[i] / dyc[j + 1];
  ks = 0.0;
  Ds = 0.0;

  // east-west -- convection
  rhoe = 0.0;
  ue = 0.0;
  Fe = 0.0;
  rhow = 0.5 * (rho[i][j] + rho[i - 1][j]);
  uw = 0.5 * (u[i][j] + u[i - 1][j]);
  Fw = rhow * uw * dyf[j];

  // north-south -- convection
  rhon = 0.5 * (rho[i][j] + rho[i][j + 1]);
  vn = 0.5 * (v[i][j] + v[i][j + 1]);
  Fn = rhon * vn * dxf[i];
  rhos = 0.0;
  vs = 0.0;
  Fs = 0.0;

#include "convective_scheme_coeffs.h"

      // Dirichlet BC on the bottom; 
       ifac = 1.0; //default
       if (convec_scheme == 1 && Fbot[i] < 0.0) // UW
          ifac = 0.0;
       bcS = (kbot[i] / dyc[j] + ifac * Fbot[i]) * dxf[i];
       b[i][j] += bcS * Tbot[i];
       aP[i][j] += bcS - Fbot[i] * dxf[i];

       // Dirichlet BC on the right; 
       ifac = 1.0; //default
       if (convec_scheme == 1 && Fright[j] > 0.0) // UW
          ifac = 0.0;
       bcE = (kright[j] / dxc[i] - ifac * Fright[j]) * dyf[j];
       b[i][j] += bcE * Tright[j];
       aP[i][j] += bcE + Fright[j] * dyf[j];

       // ------ Step 3d :: bottom-right boundary done ---
 
}

