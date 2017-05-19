#include<stdio.h>
#include<stdlib.h>
#include<math.h>

void piksrt(int , float []);

float median (int n, float x[]){
  int i;
  float y[n], m;

  if(n<1) //jnc 9sep11
    return 0.0;

  for (i=0; i < n; i++)
    y[i+1] = x[i]; 

  if (n ==1){
    return(y[1]);
  }else{
    piksrt(n,y);
    if (n % 2 ==1){
      m = y[n/2];
    }else{
      m = 0.5*(y[n/2-1]+y[n/2]);
    }
    return(m);
  }// for else if
}// end program 



