#include<stdio.h>
#include<stdlib.h>
#include<math.h>

float median (int, float[]);

int  mad(int n, float x[], float *y1,float *y2){
  int i;
  float y[n],xx[n]; 

  if(n<0) /* jnc 9/sep/11 */
    { *y1 =0.0; *y2 = 0.0;}

  if (n == 1){
    *y1 = x[0]; *y2 = 0.0;
    return(0);
  }else{
    for(i=0; i < n; i++)
      xx[i] = x[i] ;
    
    *y1 = median(n,x);
    
    for(i=0; i < n; i++)
      y[i] = fabsf(xx[i] - *y1);
    
    *y2 = median(n,y);
    return(0);
  }// for if
  
  return(0);
  
}// mad ends 


