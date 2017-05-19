#include<stdio.h>
#include<stdlib.h>
#include<math.h>

float gettime(float rand4, float rand5){
  float tdata,fracpart; 
  int intpart;

    intpart  = floor(rand4);

    if (intpart > 0)
       fracpart = rand4-intpart;
    else
       fracpart = rand4;

    tdata = fracpart + rand5;

    return(tdata); 

 }
 
