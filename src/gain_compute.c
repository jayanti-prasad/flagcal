#include <defs.h>
#include <flagcal.h>
/*------------------------------------------------------------------------------
- This program computes antenna based complex gains on the basis of visibilties.
- Gain computation for different channels is done in parallel using OpenMP
- On some computers this program reports memory errors which I do not understand.
- It should be possible to check where things are going.
- If you know parallel programming you can use something also in place of OpenMP
                                                  ----- Jayanti Prasad          
                                                      Mon Jun 28 10:20:32 IST 2010
 -------------------------------------------------------------------------------*/
int phase_calib(cmplx []);
int gainsol(int ,cmplx3 [], cmplx []);

int gain_compute(){
  cmplx  g[nants]; 
  cmplx3 V[nants*nants]; 
  int i1,j1,k1,i2,j,l,l1,m,n;
   
  fprintf(stdout,"\t  COMPUTING GAIN  \n");
   
  for(i1 = 0; i1 < nscans; i1++){
    for(i2 = 0; i2 < nab[i1] ; i2++){
      j = abid_in[i2+max_nab*i1];
      for(j1 = 0; j1 < nchans; j1++){
	for(k1 = 0; k1 < nstokes; k1++){
	  l1  =  k1 + nstokes * (j1+nchans*j);
	  for(l = 0; l < nants; l++)
	    g1[l+nants*l1] = Cmplx(0.0,0.0);
	} // for k1
      } // for j1
    } // for i1
  }// for i1 
  
  omp_set_num_threads(nthreads);
  
  for(i1 = 0; i1 < nscans; i1++){
    if (src[idmap[sid[i1]]].type != 2  && bs[i1] > 0){ // source is not the science source (is a calib) 
      for(i2 = 0; i2 < nab[i1]; i2++){ //  and the scan is not the flagged one 
  	j = abid_in[i2+max_nab*i1]; 
#pragma omp parallel for private(j1,k1,l,m,n,l1,V,g) shared(i1,i2,j,g1)
	for(j1 = 0; j1 < nchans; j1++){
	  for(k1=0; k1 < nstokes; k1++){
	    l1 = k1 + nstokes *(j1+nchans*j);
	    for(l=0; l < nants; l++)
	      g[l]  = Cmplx(0.0,0.0);
	    for(n=0; n < nbaselines; n++){
	      l = (id[n]+1)/nants;
	      m = (id[n]+1) - nants *l;
	      if(l==m)
		V[m+nants*l]  = Cmplx3(0.0,0.0,0.0);
	      else
		V[m+nants*l]  =  V1[n+nbaselines*l1];
	      V[l+nants*m]    = Cmplx3(V[m+nants*l].re,-V[m+nants*l].im,V[m+nants*l].w);
	    } // for n 
	    gainsol(niter,V,g);
	    // compute the gains for each time,stokes & chan 
	    for(l=0; l < nants; l++) 
	      g[l] = cinver(g[l]);
	    phase_calib(g);
	    // set the phase with respect to a ref antenna
	    for(l=0; l < nants; l++)
	      g1[l+nants*l1] = g[l]; 
	  } // for k1
	} // for j1
	// let us wait for everbody
#pragma omp barrier
	// paralle section ends 
      } // for i2
    } // for if 
  } // for i1 
  free(V1);
  fprintf(stdout,"\t gain computed successfully \n");
  
  return(0);
  
}// end gain_compute 


int phase_calib(cmplx g[]){
  /* This normalizes the gain phases with respect to a ref ant */
  int l;
  float theta,R, theta_ref; 
  
  theta_ref = cphase(g[ref_ant]);
  
  for(l = 0; l < nants; l++){
    theta = cphase(g[l]);
    R = cmod(g[l]);
    theta = (M_PI/180.0)*(theta-theta_ref);
    g[l]  = Cmplx(R*cos(theta),R*sin(theta));
  }
  return(0);
}// end phase_calib


