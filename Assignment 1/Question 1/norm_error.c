#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(){

int i = 0,nx,j;
double x, Ta,Tn,norm;

double count = 0;
FILE *file;
FILE *newfile;
char filename[100];

newfile = fopen("norm_nx.txt","w");

for (nx = 50 ; nx <= 500; nx +=50){


double *error = (double*)malloc(nx*sizeof(double));

  i = 0;
  count = 0;
  
  sprintf(filename,"nx%dT_fb.txt",nx);
  file = fopen(filename,"r");

  while(fscanf(file,"%lf %lf", &Ta,&Tn)==2){
   
   error[i] = fabs(Tn - Ta);
   i++ ;
  }
  for (j=0; j<nx ; j++){
  
   count +=( error[j] * error[j]);
   norm = sqrt(count);
  }
 fprintf(newfile,"%d %lf\n",nx,norm);
 free(error);
 fclose(file);
}


fclose(newfile);

return 0;
}

