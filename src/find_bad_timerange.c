#include<stdio.h>
#include<stdlib.h>
#include<math.h>

/*
  looks for data that has a discrepant mean after averaging (using the mean, not meadian) over
  a window of length win. Data that is thrshigh*mad2 greater than mad1 or data that
  is thrslow*mad2 less than mad1 or data that is less than thrsmean*mad1 is flagged.
  The assymetric flagging is meant to account for the non gaussian statistics of the amplitude
  of the visibility for the target source.

  jnc 10/sep/11
*/
int find_bad_timerange(int np, float x[], float y[], float thrshigh, float thrslow, float thrsmean,
			float m1, float m2, int win){
  float *x1;
  int    i1,i2,winlen,n;

  
  /* don't process -- this should be caught as the whole baseline being bad*/
  if(win >np/2) return 0; 

  
  /*  fprintf(stderr,"m1=%6.2f m2=%6.2f thrshigh=%6.2f thrslow=%6.2f thrsmean=%6.2f win=%d\n",
	  m1,m2,thrshigh,thrslow,thrsmean,win);
  */
  winlen = 2*win+1; 
  if((x1=(float*)malloc(np*sizeof(float)))==NULL){
      fprintf(stderr,"find_bad_timerange() Malloc Error\n");
      return 1;
  }
  for (i1=0; i1 < np; i1++)x1[i1]=0.0;

  for(i1=0; i1 < np; i1++){
    if (i1 >=0 && i1 < win ){
      for(n=0,i2=0; i2 < winlen; i2++)
	if(y[i1+i2]>0){x1[i1] +=x[i1+i2];n++;}
    }
    if (i1 >=win && i1 < np-win){
      for(n=0,i2=-win; i2 <=win; i2++)
	if(y[i1+i2]>0){x1[i1] +=x[i1+i2];n++;}
    }
    if(i1 >= np-win && i1 < np){
      for(n=0,i2=-winlen; i2 < 0; i2++)
	if(y[i1+i2] > 0){x1[i1] +=x[i1+i2];n++;}
    }
    if(n>0)x1[i1] = x1[i1]/n;
  }// for i1
 
  if(y[0] > 0)
    x1[0]=0.5*(x[0]+x1[0]); // This is ad-hoc however I found it useful 
  if(y[np-1]>0)
    x1[np-1]=0.5*(x[np-1]+x1[np-1]);
 
 for(i1 = 0; i1 < np; i1++){
   /*
     if(x1[i1]-m1 >  thrshigh * m2){fprintf(stderr,"flag high\n");}
     if(m1-x1[i1] >  thrslow * m2){fprintf(stderr,"flag low\n");}
     if(x1[i1] < thrsmean*m1){fprintf(stderr,"flag mean\n");}
   */
   if(y[i1]<0) continue;
   if((x1[i1]-m1 >  thrshigh * m2)||(m1-x1[i1] >  thrslow * m2) || x1[i1] < thrsmean*m1)
     y[i1] = -1.0;
 }// for if  

 free(x1);
  return(0); 
}// peak_remove ends 

 
