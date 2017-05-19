#include <defs.h>
#include <flagcal.h>

int gain_apply(){
  int i1,i2,i3,j,j1,k1,l,l1,l2,m,n;
  cmplx z,z1; 

  fprintf(stdout,"\t  APPLYING gain\n");
   
  for(i1 = 0; i1 < nscans; i1++){
    for(i2 = 0; i2 < nab[i1]; i2++){
      j = abid_in[i2+max_nab*i1];
      for(i3 = ab_start[j]; i3 < ab_end[j]; i3++){
	for(j1 = 0; j1 < nchans; j1++){
	  for(k1 = 0; k1 < nstokes; k1++){
	    l1 = k1 + nstokes *(j1+nchans*j);
	    l2 =  k1 + nstokes*(j1+nchans*i3);
	    for(n=0; n < nbaselines; n++){
	      l = (id[n]+1)/nants;
	      m = (id[n]+1) - nants *l;
	      z  =  ccmulti(g1[l+nants*l1],g1[m+nants*l1]);
	      z1 =  Cmplx(Visib[n+nbaselines*l2].re,Visib[n+nbaselines*l2].im);
	      if (Visib[n+nbaselines*l2].w > 0.0  &&  cmod(z) > 0.0){
                Visib[n+nbaselines*l2].re = cmulti(z1,z).re;
                Visib[n+nbaselines*l2].im = cmulti(z1,z).im;
              }else{
	       if(gs[i1] > 0) Visib[n+nbaselines*l2].w  = -1.0;
              }// for if
	    }// for n
	  }// for k1
	} // for j1
      }// for i3
    }// for i2
  } // for i1
  
  fprintf(stdout,"\t  gain applied \n");
  
  return(0);
} // gain_apply ends 
