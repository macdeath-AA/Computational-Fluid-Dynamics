#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void grid(int nx, double xst, double xen, double *x, double *dx){

    *dx = (xen - xst) /( nx - 1);
    
   
    int i;
    for (i = 0; i < nx; i++){
        x[i+1] = xst+ i*(*dx);
      
    }
}
void enforce_bcs(int nx, double *x, double *T) {
   // x[0] = 0;
   // x[nx-1] = 1;
    T[0] = -1.0;
    T[nx-1] = 1.0;
}

void set_initial_condition(int nx, double *x, double *T)
{ 
      
    int i;
    for ( i = 0; i < nx; i++)
    {
        //t = 0.001
        T[i] = erf((x[i] - 0.5) / 0.063245 );
       
    }
enforce_bcs(nx,x,T);
}

void get_rhs(int nx, double dx, double *x, double *T, double *rhs)
{
    int i;
    for (i = 1; i < nx-1; i++){
        rhs[i] = (T[i + 1] + T[i - 1] - 2 * T[i]) / (dx * dx);      
    }
   rhs[0] = 0;
   rhs[nx-1]=0;
}

void timestep_Euler(int nx, double dt, double dx, double *x, double *T, double *rhs)
{
       int j;

        get_rhs(nx,dx,x,T,rhs);

        for(j=0;j<nx ;j++ ){
        T[j] = dt*rhs[j] + T[j];
       
      
    }
 }
  

void output_soln(int nx, int it, double tcurr, double *x, double *T)
{
 FILE *fp;
  char filename[100];

  sprintf(filename,"/clhome/me20btech11001/cfd/q1bnx1000/T_x_%d.txt",it);
  printf("%s\n",filename);

  fp=fopen(filename,"w");

  int i;
  for(i=0;i<nx;i++){
    fprintf(fp,"%lf %lf\n",x[i],T[i]);
  }
  fclose(fp);
  printf("Done writing solution for timestep %d\n",it);

}
    

int main()
{
    int nx;
    double *x, *T, *rhs, tst, ten, xst, xen, dx, dt, tcurr;

    int i, it, num_time_steps, it_print;
    FILE *fp;

    
    fp = fopen("input.in", "r");
    fscanf(fp, "%d\n", &nx);
    fscanf(fp, "%lf %lf\n", &xst, &xen);
    fscanf(fp, "%lf %lf\n", &tst, &ten);
    fscanf(fp, "%le\n",&dt);
    fclose(fp);

    printf("nx: %d\nxo: %lf\nx1: %lf\nt0: %lf\nt1: %lf\ndt: %le\n", nx, xst, xen, tst, ten,dt);

    x = (double *)malloc(nx * sizeof(double));
    T = (double *)malloc(nx * sizeof(double));
    rhs = (double *)malloc(nx * sizeof(double));

    grid(nx, xst, xen, x, &dx);
    
    set_initial_condition(nx, x, T);
   
    
    num_time_steps = (int)((ten - tst) / dt) + 1;
    
    it_print = num_time_steps / 10;               
       
    for (it = 0; it < num_time_steps; it++)
    {
        tcurr = tst + (double)it * dt;
 
        timestep_Euler(nx,dt,dx,x,T,rhs);
        
        if(it%it_print ==0){ 
        output_soln(nx,it,tcurr,x,T);}
       
    }

  free(rhs);
  free(T);
  free(x);
return 0;

}

