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

void exact_temp(int nx, double tcurr, double *x, double *Ta){

int i;
Ta[0] = -1;
Ta[nx-1]= 1;
for (i=1; i< nx-1;i++){
 Ta[i] = erf((x[i]-0.5)/(2* sqrt(tcurr)));
}

}

void output_soln(int nx, double *Ta, double *T){

FILE *fp;
  char filename[100];

  sprintf(filename,"/clhome/me20btech11001/cfd/errors/nx%dT_fb.txt",nx);
  printf("%s\n",filename);

  fp=fopen(filename,"w");

  int i;
  for(i=0;i<nx;i++){
    fprintf(fp,"%lf %lf\n",Ta[i],T[i]);
  }
  fclose(fp);
  printf("Done writing solution for nx %d\n",nx);
}


     

int main()
{
  
    int i, it, num_time_steps, it_print;  
    int nx; 
    double xst = 0.0, xen = 1.0, tst = 0.001, ten = 0.002, dt = 3.125e-6;
    double *x, *T, *rhs,*Ta;
    double dx,tcurr;
   // printf("nx: %d\nxo: %lf\nx1: %lf\nt0: %lf\nt1: %lf\ndt: %le\n", nx, xst, xen, tst, ten,dt);

    
    for (nx= 50; nx <= 400; nx+=50){
    x = (double *)malloc(nx*sizeof(double));
    T = (double *)malloc(nx*sizeof(double));
    rhs = (double *)malloc(nx*sizeof(double));
    Ta = (double *)malloc(nx*sizeof(double));
   
    grid(nx, xst, xen, x, &dx);    
    set_initial_condition(nx, x, T);   
    
    num_time_steps = (int)((ten - tst) / dt) + 1;
    
    it_print = num_time_steps / 10;               
      
    for (it = 0; it < num_time_steps; it++)
    {
    tcurr = tst + (double)it * dt;
 
       timestep_Euler(nx,dt,dx,x,T,rhs);
       exact_temp(nx,tcurr,x,Ta); 
      
       if (it == num_time_steps -1){
       
       output_soln(nx,Ta,T);
      }
} 
     free(rhs);
  free(T);
  free(x); 
  free(Ta);    
    }
return 0;

}

