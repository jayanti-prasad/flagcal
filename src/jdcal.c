#include<stdio.h>
#include<stdlib.h>
#include<math.h>

int ftime(float time_data, float ttime[],float *tmsec){
  double xx[7],yy[7],tsec,t1,t2,t3;
  int i;
  
  for(i=0; i < 7; i++){
    xx[i]  = 0.0;
    yy[i]  = 0.0;
  }
  
  ttime[0] = floor(time_data); // day 
  
  xx[0] = time_data-ttime[0];
  yy[0] = 0.0;
  
  for(i=1; i < 7; i++){
    xx[i] = 10.0 *xx[i-1];
    yy[i] = floor(xx[i]);
    xx[i] -= yy[i];
  }// for if
  
  tsec = 0.0;
  for(i = 1; i < 7; i++)
    tsec += yy[i]  * (8640.0/pow(10.0,i-1));
  
  t1 = floor(tsec/3600.0);
  t2 = floor((tsec - t1 * 3600.00)/60.0);
  t3 = tsec -60.0*(60.0*t1 + t2);
  
  ttime[1] = t1; // hour 
  ttime[2] = t2; // min 
  ttime[3] = t3; // sec 
  
  *tmsec  = (ttime[0]*24.0 + ttime[1])*3600.0 + 60.0 *ttime[2] + ttime[3];

  return(0);

}

int fdate(double day_data, int dd[]){
  int i,j,l,n;
  
  l = day_data + 68569;
  n = ( 4 * l ) / 146097;
  l = l - ( 146097 * n + 3 ) / 4;
  i = ( 4000 * ( l + 1 ) ) / 1461001;
  l = l - ( 1461 * i ) / 4 + 31;
  j = ( 80 * l ) / 2447;
  dd[0] = l - ( 2447 * j ) / 80; // day
  l = j / 11;
  dd[1] = j + 2 - ( 12 * l ); // month 
  dd[2] = 100 * ( n - 49 ) + i + l; // year 

  return(0);

}



