#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<stdbool.h>
#include<complex.h>

float  cmod(cmplx z){
   return(sqrt(z.re*z.re+z.im*z.im));
}

float cphase(cmplx z){
  float theta = atan2(z.im,z.re);
  theta *=(180.0/M_PI);
  return(theta);
}

cmplx cmulti(cmplx z1, cmplx z2){
  cmplx z;
  z.re = z1.re * z2.re - z1.im * z2.im;
  z.im = z1.re * z2.im + z1.im * z2.re;
  return(z);
}

cmplx ccmulti(cmplx z1, cmplx z2){
  cmplx z;
  z.re = z1.re * z2.re + z1.im * z2.im;
  z.im = -z1.re * z2.im+ z1.im * z2.re;
  return(z);
}

cmplx csmulti(cmplx z, float x){
  z = Cmplx(x*z.re,x*z.im); 
  return(z);
}

cmplx  cdiff(cmplx z1, cmplx z2){
  cmplx z;
  z.re = z1.re  - z2.re;
  z.im = z1.im  - z2.im;
   return(z);
}

cmplx  cadd(cmplx z1,cmplx z2){
    cmplx z;
    z.re = z1.re  + z2.re;
    z.im = z1.im  + z2.im;
    return(z);
}

cmplx cinver(cmplx z){
  float R;
  cmplx z1;
  R = z.re * z.re + z.im * z.im ;
  if(R > 0.0){
    z1.re =  z.re/R; 
    z1.im =  -z.im/R;
  }else{
    z1.re  = 0.0;
    z1.im  = 0.0;
  }
  return(z1);
 }

cmplx Cmplx(float x, float y){
    cmplx z; 
    z.re=x;
    z.im=y;
    return(z); 
}

cmplx3 Cmplx3(float x, float y, float w){
    cmplx3 z;
    z.re=x;
    z.im=y;
    z.w =w;
    return(z);
}

