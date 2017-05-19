#include <defs.h>
#include <flagcal.h>

int  group_data(float *randpar,float U[],int *ant1, int *ant2,int*isrc){
  float coeff;
  int l,m;
  
  coeff = crval4/1000.0;
  U[0]  = coeff*randpar[0]*pscal[1];
  U[1]  = coeff*randpar[1]*pscal[2];
  U[2]  = coeff*randpar[2]*pscal[3];
  
  l = (int) (randpar[3]/256.0 +0.1) ;
  m = (int) (randpar[3] - l *256.0+0.1);

  *ant1= l-1; 
  *ant2= m-1;
  *isrc = (int) randpar[6];

  if (*ant1 < 0 || *ant2 < 0){
    fprintf(stderr,"\t Error ! ant1 = %d ant2 = %d\n",*ant1,*ant2);
    return(-1);
  } 
  
  return(0);
}


