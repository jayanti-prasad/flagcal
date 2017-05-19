#include <defs.h>
#include <flagcal.h>

int smooth_data(){
  int i1,i2,i3,j1,k1,j,n,l1,l2,intm,*ino;

  fprintf(stdout,"\t  SMOOTHING data \n"); 
  
  intm = 0;
  for (i1=0; i1 < nscans; i1++)
    intm+=nab[i1]; 
  
  ino = (int *)malloc(intm*nchans*nstokes*nbaselines*sizeof(int));
  V1 = (cmplx3 *)malloc(nablocks*nstokes*nchans*nbaselines*sizeof(cmplx3));

  for(i1 = 0; i1 < nscans; i1++){
    for(i2 = 0; i2 < nab[i1]; i2++){
      j = abid_in[i2+max_nab*i1]; 
      for(j1 = 0; j1 < nchans; j1++){
	for(k1 = 0; k1 < nstokes; k1++){
	  l1  =  k1 + nstokes * (j1+nchans*j);
	  for(n=0; n < nbaselines; n++){
	    ino[n+nbaselines*l1] = 0; 
	    V1[n+nbaselines*l1]  = Cmplx3(0.0,0.0,1.0); 
         }// for n
	} // for k1 
      } // for j1
    }// for i2
  }// for i1
  
  for(i1 = 0; i1 < nscans; i1++){
    for(i2 = 0; i2 < nab[i1]; i2++){
      j = abid_in[i2+max_nab*i1]; 
      for(i3 = ab_start[j]; i3 < ab_end[j]; i3++){
	for(j1 = 0; j1 < nchans; j1++){
	  for(k1 = 0; k1 < nstokes; k1++){
	    l1  =  k1 + nstokes * (j1+nchans*j);
	    l2  =  k1 + nstokes * (j1+nchans*i3);
	    for(n = 0; n < nbaselines; n++){
	      if (Visib[n+nbaselines*l2].w > 0.0){
		V1[n+nbaselines*l1].re += Visib[n+nbaselines*l2].re; 
		V1[n+nbaselines*l1].im += Visib[n+nbaselines*l2].im; 
		ino[n+nbaselines*l1]++;
	      }// for if
	    }// for n
	  } // for k1 
	} // for j1
      }// for i3
    }// for i2
  }// for i1
  
  for(i1 = 0; i1 < nscans; i1++){
    for(i2 = 0; i2 < nab[i1]; i2++){
      j = abid_in[i2+max_nab*i1]; 
      for(j1 = 0; j1 < nchans; j1++){
	for(k1 = 0; k1 < nstokes; k1++){
	  l1  =  k1 + nstokes * (j1+nchans*j);
	  for(n = 0; n < nbaselines; n++){
	    if (ino[n+nbaselines*l1] > 0){     
	      V1[n+nbaselines*l1].re /= ino[n+nbaselines*l1];  
	      V1[n+nbaselines*l1].im /= ino[n+nbaselines*l1];
	      V1[n+nbaselines*l1].w   = 1.0;
	    }else{
	      V1[n+nbaselines*l1].w   = 0.0;
	    }// for if
	  }// for n
	} // for k1 
      } // for j1
    }// for i2
  }// for i1
  free(ino); 
  free(vsr);
  free(dant);
  free(bad_base);
  free(bad_chans);
  
  return(0);
  
}// smooth_data ends 
