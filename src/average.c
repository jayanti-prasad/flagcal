#include<math.h>
#include<stdio.h>
#define TINY  1.0e-8
int average(int np, float *arr, float *mean, float *rms)
{ int i;
  float sum,sumsqr,devsqr;
  
  if(np<0){
    fprintf(stderr,"ERR:average() called with negative array length!\n"); 
    *mean=0.0;*rms=0.0;return 1;
  }

  if(np==0){
    fprintf(stderr,"WRN:average() called with 0 array length!\n"); 
    *mean=0.0; *rms=0.0;return 0;
  }
  if(np==1){
    *mean=arr[0];*rms=0.0;
    return 0;
  }
  
  sum=sumsqr=0.0;
  for(i=0;i<np;i++){
    sum   += arr[i];
    sumsqr+= arr[i]*arr[i];
  }

  *mean=sum/np;
  devsqr=(sumsqr - np*(*mean)*(*mean));
  if((devsqr/(*mean))<TINY) *rms=0.0;
  else *rms = sqrt(devsqr)/(np-1.0);

  return 0;
}
