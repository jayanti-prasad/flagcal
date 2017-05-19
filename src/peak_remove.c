#include<stdio.h>
#include<stdlib.h>
#include<math.h>

int mad(int,float[],float*,float*);

int peak_remove(int np, float x[], float y[], float thrshigh, float thrslow,int win){
  float x1[np],x2[np],m1,m2;
  int i1,i2,winlen;

  for (i1=0; i1 < np; i1++){
    x1[i1] = x[i1];
    x2[i1] = 0.0;
  }// for i1
  
  winlen = 2*win+1; 
  for(i1=0; i1 < np; i1++){
    if (i1 >=0 && i1 < win )
      for(i2=0; i2 < winlen; i2++)
	x2[i1]+=x1[i1+i2]/winlen;
    if (i1 >=win && i1 < np-win)
      for(i2=-win; i2 <=win; i2++)
	x2[i1]+=x1[i1+i2]/winlen;
    if(i1 >= np-win && i1 < np)
      for(i2=-winlen; i2 < 0; i2++)
	x2[i1]+=x1[i1+i2]/winlen;
  }// for i1
 
  x2[0]=0.5*(x1[0]+x2[0]); // This is ad-hoc however I found it useful 
  x2[np-1]=0.5*(x1[np-1]+x2[np-1]);
 
  // now subtract the smoothed one from the original
    
 for(i1 = 0; i1 < np; i1++)
   x2[i1] = x1[i1]-x2[i1];
 
 mad(np,x2,&m1,&m2); 
 
 for(i1 = 0; i1 < np; i1++){
   if((x2[i1]-m1 >  thrshigh * m2)||(m1-x2[i1] >  thrslow * m2))
     y[i1] = -1.0;
 }// for if  

  return(0); 
}// peak_remove ends 

 
