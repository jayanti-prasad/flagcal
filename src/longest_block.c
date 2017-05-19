#include<stdio.h>
#include<stdlib.h>
#include<math.h>

float max(float, float);

int longest_block(int n, float T, float  x[]){
  int i, j,l;
  float tt;
  l = 0 ;
  j = 0;
  for(i = 0; i < n; i++){
    if(x[i] > T){
      j++;
      tt = max((float)l,(float)j);
      l = (int) tt;
    }else{
      j = 0;
    }
  }
  return(l);
}

