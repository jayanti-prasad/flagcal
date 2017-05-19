#include <defs.h>
#include <flagcal.h>

int gainsol(int niter, cmplx3 v[], cmplx g[]){
  int i,l,m;
  cmplx z,z1,gg[nants]; 
  float tmp,dmax;
  
  for(l = 0; l < nants; l++){
    g[l]  = Cmplx(1.0,1.0);
    gg[l] = Cmplx(0.0,0.0);
  }// for i
  for (i = 0; i <= niter; i++){
    dmax = 0.0;
    for(l = 0; l < nants ; l++){
      z = Cmplx(0.0,0.0);
      for(tmp=0.0,m = 0; m < nants; m++){
	if(m != l) {
	  z1  = Cmplx(v[m+nants*l].re,v[m+nants*l].im);
	  z   = cadd(z,cmulti(z1,csmulti(g[m],v[m+nants*l].w))); 
	  tmp += pow(cmod(g[m]),2) * v[m+nants*l].w;
	}// for if
      } // for m
      if (tmp > 0.0)
	gg[l] = Cmplx(z.re/tmp,z.im/tmp);
    }// for l           
    for(l = 0; l < nants; l++){
      if(cmod(gg[l]) > 0.0){
	z  = cdiff(gg[l],g[l]);
	dmax = max(dmax,cmod(z)/cmod(g[l]));
	g[l] = cadd(g[l],csmulti(z,alpha)); 
      }else{
	g[l]=Cmplx(0.0,0.0); 
      }// else 
    }// for l
    //fprintf(stderr,"iter=%d dmax=%2.12f\n",i,dmax); 
    if (dmax < eps){
      //  fprintf(stderr,"convergence after %d interations: dmax = %2.8f\n",i,dmax);
      break ;
    }// fi 
    if (i == niter && dmax > eps){
      g[l] = Cmplx(0.0,0.0);
    }// for if 
  } // for i
  return(0); 
}      


