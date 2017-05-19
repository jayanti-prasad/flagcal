#include <defs.h>
#include <flagcal.h> 
/*----------------------------------------------------------------------------------
- This program uses OpenMP 
- This program and the associated peak_remove have been sufficiently modified. Now the 
  size of the smoothing windows has been frozen. However, in order to make the 
  program computationally cheap in the n+1 th iteration the input data for smoothing is 
  taken from n th iteration. The module has been tested and found working. Now the 
  only free paramedters for the smoothing & thresholding are the initial threshold, 
  lower factor of threshold,  the number of iterations and the initial smoothing
  only.

                      ---------   Jayanti Prasad,  Fri Jun 18 11:38:50 IST 2010

---------------------------------------------------------------------------------------*/

int peak_remove(int , float [], float [], float , float,int );

int remove_peaks_time(){
  int i1,j1,k1,n,n1,n2,i2,l1,i,win;
  float x[max_times],thrshigh,thrslow,y[max_times],y1[max_times]; 
  cmplx z; 
  
  if (!rfi_mode_t) return(0);
  
  omp_set_num_threads(nthreads);
  
  fprintf(stdout,"\t removing RFI peaks from time \n");
  
  for(n1=0,n2=0,i1 = 0; i1 < nscans; i1++){
    if (gs[i1] > 0){
      for(j1 = 0; j1 < nchans; j1++){
	for(k1 = 0;  k1 < nstokes; k1++){
#pragma omp parallel for shared(i1,j1,k1,n1,n2) private(n,i2,l1,x,y,y1,z,thrshigh,thrslow,win,i)     
	  for (n = 0; n < nbaselines; n++) {
	    for(i2 = tid[i1]; i2 < tid[i1+1]; i2++){
	      l1   = k1 + nstokes *(j1+nchans*i2);
	      z    = Cmplx(Visib[n+nbaselines*l1].re,Visib[n+nbaselines*l1].im);
	      x[i2-tid[i1]]  = cmod(z);
	      y[i2-tid[i1]]  = Visib[n+nbaselines*l1].w;
	      y1[i2-tid[i1]] = y[i2-tid[i1]];
	    }// for i2 
	    // the array has been populated 
	    if(src[idmap[sid[i1]]].type==2){
	      thrshigh = thrs_peak_t;
	      thrslow  = thrs_peak_tlow;
	    }else{
	      thrshigh = thrs_peak_t;
	      thrslow  = thrs_peak_t;
	    }
	    win  = win_min_t;
	    for(i = 0; i < nriter_t; i++){
	      peak_remove(tpnt[i1],x,y,thrshigh,thrslow,win);
	      win      *= dwin_t;
	      thrshigh *= dthrs_peak_t; 	    
	      thrslow  *= dthrs_peak_t; 	    
	    }// for i
	    for(i2 = tid[i1]; i2 < tid[i1+1]; i2++){
	      l1   = k1 + nstokes *(j1+nchans*i2);
	      Visib[n+nbaselines*l1].w = y[i2-tid[i1]];
	      if (y1[i2-tid[i1]] > 0.0  &&  y[i2-tid[i1]] < 0.0)
	      n2++;
	      n1++;
	    }// for i2 
	  }// for n
	  // let us wait for everbody
#pragma omp barrier
	  // end of parallel section 
	}// for k1
      }// for j1
    }
  }// for i1
  
  fprintf(stdout,"\t Flagging by removing peaks  [ALL]: %2.2f \n",(float)n2/n1); 
  return(0);
  
}// for remove peaks 
 
int remove_peaks_freq(){
  int i,i1,i2,j1,k1,n,l1,win,n1,n2;
  float x[nchans],thrshigh,thrslow,y[nchans],y1[nchans];
  cmplx z; 


// the j loop has been changed from 1,nchans-2, to 0 nchans-1//jayanti/30-03-2011
  
  if (!rfi_mode_f) return(0);
  
  omp_set_num_threads(nthreads);
  
  fprintf(stdout,"\t removing RFI peaks from frequency \n");
  
  for(n1=0,n2=0,i1 = 0; i1 < nscans; i1++){
    if (gs[i1] > 0){
      for(i2 = tid[i1]; i2 < tid[i1+1]; i2++){
	for(k1 = 0;  k1 < nstokes; k1++){
#pragma omp parallel for  shared(i1,i2,k1,n1,n2) private(i,n,j1,l1,x,y,y1,z,win,thrshigh,thrslow) 	    
	  for (n = 0; n < nbaselines; n++) {
	    for(j1 = 0; j1 < nchans; j1++){
	      l1   = k1 + nstokes *(j1+nchans*i2);
	      z    = Cmplx(Visib[n+nbaselines*l1].re,Visib[n+nbaselines*l1].im);
	      x[j1]= cmod(z);
	      y[j1]= Visib[n+nbaselines*l1].w; 
	      y1[j1]= y[j1]; 
	    }//for  j1
	    // the array has been populated 
	    if(src[idmap[sid[i1]]].type==2){
	      thrshigh = thrs_peak_f;
	      thrslow  = thrs_peak_flow;
	    }else{
	      thrshigh = thrs_peak_f;
	      thrslow  = thrs_peak_f;
	    }
	    win  = win_min_f;
	    for(i = 0; i < nriter_f; i++){
	      peak_remove(nchans,x,y,thrshigh,thrslow,win);
	      win      *= dwin_f;
	      thrshigh *= dthrs_peak_f; 	    
	      thrslow  *= dthrs_peak_f; 	    
	    }// for i
	    for(j1= 0; j1 < nchans; j1++){
	      l1   = k1 + nstokes *(j1+nchans*i2);
	      Visib[n+nbaselines*l1].w = y[j1];
	      if (y1[j1] > 0.0 && y[j1] < 0.0)
		n2++;
	      n1++;
	    }// for j1 
	  }// for n
#pragma omp barrier
	  // let us wait for everbody
	}// for k1
      }// for i2
    }//for if
  } // for i1
  
  fprintf(stdout,"\t Flagging by removing peaks  [ALL]: %2.2f \n",(float)n2/n1); 
  
  return(0);
  
}// for remove peaks 

 
