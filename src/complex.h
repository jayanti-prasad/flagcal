#ifndef _COMPLEX_H
#define _COMPLEX_H

typedef struct Complex {float re,im;} cmplx;
typedef struct Complex3 {float re,im,w;} cmplx3;

float  cmod(cmplx );
float  cphase(cmplx );
cmplx cmulti(cmplx ,cmplx );
cmplx csmulti(cmplx , float );
cmplx ccmulti(cmplx ,cmplx );
cmplx cdiff(cmplx,cmplx);
cmplx cadd(cmplx,cmplx);
cmplx cinver(cmplx );
cmplx Cmplx(float,float);
cmplx3 Cmplx3(float,float,float); 

#endif /*_COMPLEX_H */
