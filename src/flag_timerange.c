#include <defs.h>
#include <flagcal.h> 
/*----------------------------------------------------------------------------------
- This program uses OpenMP 
- Flags time ranges whose median value is either 

  > scan_median + thrshigh*scan_mad
  < scan_median - thrslow*scan_mad
  < thrsmean*scan_median

  The assymetric flagging about the mean is done only for the target source,
  where for weak sources the statistics is non gaussian. For the calibrator
  sources, thrshigh is used to do symmetric flagging about the scan_median.

  jnc sep 2011

---------------------------------------------------------------------------------------*/

int  find_bad_timerange(int npt,float *x, float *y,float thrshigh,float thrslow, float thrsmean,
			float m1,float m2, int win);

int flag_timerange( float *mad1, float *mad2){
  int i1,j1,k1,n,n1,n2,np,i2,l1,i,win;
  float x[max_times],thrshigh,thrslow,thrsmean,y[max_times],y1[max_times]; 
  float dthrs;
  int   dwin;
  float m1,m2;
  cmplx z; 
  
  if (!range_t_mode) return(0);
  
  omp_set_num_threads(nthreads);
  
  fprintf(stdout,"\t flagging bad time ranges \n");
  
  for(n1=0,n2=0,i1 = 0; i1 < nscans; i1++){
    if (gs[i1] > 0){
      for(j1 = 0; j1 < nchans; j1++){
	for(k1 = 0;  k1 < nstokes; k1++){
#pragma omp parallel for shared(i1,j1,k1,n1,n2) private(n,i2,l1,x,y,y1,z,i,thrshigh,thrslow,thrsmean,m1,m2,win,dwin,dthrs,np)
	  for (n = 11; n < nbaselines; n++) { 
	    m1=mad1[i1]; 
	    m2=mad2[i1];
	    if(src[idmap[sid[i1]]].type==2){
	      thrshigh = thrs_range_t;
	      thrslow  = thrs_low_range_t;
	      thrsmean  = thrs_mean_range_t;
	    }else{
	      thrshigh = thrs_range_t;
	      thrslow  = thrs_range_t;
	      thrsmean  = thrs_mean_range_t;
	    }
	    win   = win_min_rt;
	    dwin  = dwin_rt;
	    dthrs = dthrs_rt;
	    for(np=0,i2 = tid[i1]; i2 < tid[i1+1]; i2++){
	      l1   = k1 + nstokes *(j1+nchans*i2);
	      if(Visib[n+nbaselines*l1].w<0) continue;
	      z    = Cmplx(Visib[n+nbaselines*l1].re,Visib[n+nbaselines*l1].im);
	      x[np]  = cmod(z);
	      y[np]  = Visib[n+nbaselines*l1].w;
	      y1[np] = y[i2-tid[i1]];
	      np++;
	    }// for i2 
	    if(np<1) continue;
	    for(i = 0; i < nriter_t; i++){
	      find_bad_timerange(np,x,y,thrshigh,thrslow,thrsmean,m1,m2,win);
	      win      *= dwin;
	      thrshigh *= dthrs; 	    
	      thrslow  *= dthrs; 	    
	    }// for i
	    for(np=0,i2 = tid[i1]; i2 < tid[i1+1]; i2++){
	      l1   = k1 + nstokes *(j1+nchans*i2);
	      if(Visib[n+nbaselines*l1].w <0) continue;
	      if(y[np] < 0)
		Visib[n+nbaselines*l1].w = -1.0;
	      if (y1[np] > 0.0  &&  y[np] < 0.0) n2++;
	      np++;
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
  
    fprintf(stdout,"\t Flagging of bad timerange  [ALL]: %2.2f \n",(float)n2/n1); 
  return(0);
  
}// for remove flag_timerange
 
