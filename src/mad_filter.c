#include<defs.h>
#include<flagcal.h>
/*-----------------------------------------------------------------------------
  - This program contains modules which are used to filter very low/high data 
   points.
  - Some of the modules are called before the calibration and others after the
  calibration.
  - The main modules are as follws:
  (a)  global_mad_filter : computes a set of numbers (m1,m2) for every source, and 
  on the basis of that the data is filtered. At present this modules is called just 
  after reading and indexing the data.
  (b)  pre_mad_fileter   : this filters the data on the basis of a set of representative 
  numbers (m1,m2) for every scan.
  (c)  post_milter       : this filters data after calibration and is expected to 
  compare the visibility amplitude for different baselines and flag the baselines 
  for which the values are two high or too low. 
  -- Note that all of the module use OpenMP and carry out computation in parallel.
  -- Jayanti Prasad (August 14, 2010)
  ---------------------------------------------------------------------------*/

int flag_timerange( float *mad1, float *mad2);

int   mad(int , float [], float *,float *);
int   compute_mad(float[],float[]);
float median (int, float[]);

int mad_filter_global(){
  int i1,i2,j1,k1,l1,n,n1,n2,np,np1[nsrc+1];
  float *m1,*m2,*mm1,*mm2,*y,*y1,*y2,*mad1,*mad2,tmp,x[nbaselines];
  cmplx z; 
  
  if (!gmad_mode) return(0); 
  
  fprintf(stdout,"\t GLOBAL MAD FILTETING \n"); 
  
  sprintf(err_msg,"ERROR ! memory allocation in MAD_FILTER_GLOBAL failed for");
  
  if((mad1 = (float *)malloc(ntimes*nchans*nstokes*sizeof(float)))==NULL)
    {fprintf(stdout,"\t %s mad1\n",err_msg);return(-1);}
  if((mad2 = (float *)malloc(ntimes*nchans*nstokes*sizeof(float)))==NULL)
    {fprintf(stdout,"\t %s mad2\n",err_msg);return(-1);}
  if((m1 = (float *)malloc(ntimes*sizeof(float)))==NULL)
    {fprintf(stdout,"\t %s m1\n",err_msg);return(-1);}
  if((m2 = (float *)malloc(ntimes*sizeof(float)))==NULL)
    {fprintf(stdout,"\t %s m2\n",err_msg);return(-1);}
  if((mm1 = (float *)malloc((nsrc+1)*sizeof(float)))==NULL)
    {fprintf(stdout,"\t %s mm1\n",err_msg);return(-1);}
  if((mm2 = (float *)malloc((nsrc+1)*sizeof(float)))==NULL)
    {fprintf(stdout,"\t %s mm2\n",err_msg);return(-1);}
  
  omp_set_num_threads(nthreads);

  // initilization 
  
  for(i1=0; i1 <=nsrc; i1++){
    mm1[i1] = 0.0;
    mm2[i1] = 0.0;
  }

  for(i1 = 0; i1 < ntimes; i1++){
    m1[i1] = 0.0;
    m2[i1] = 0.0;
    for(j1 = 0; j1 < nchans; j1++){
      for(k1 = 0; k1 < nstokes; k1++){
	l1 = k1+nstokes*(j1+nchans*i1);
	mad1[l1] = 0.0;
	mad2[l1] = 0.0;
      }// for k1
    }// for j1
  }// for i1
  
  // this is the main loop 
  
  for(i1 = 0; i1 < ntimes; i1++){
#pragma omp  parallel for shared(i1) private(j1,k1,n,np,l1,z,x) 
    for(j1 = 0; j1 < nchans; j1++){
      for(k1 = 0; k1 < nstokes; k1++){
	l1 = k1+nstokes*(j1+nchans*i1);
	np = 0;
        for(n = 0; n < nbaselines; n++){
	  z   = Cmplx(Visib[n+nbaselines*l1].re,Visib[n+nbaselines*l1].im);
	  if (Visib[n+nbaselines*l1].w < 0.0) continue;
	  x[np] = cmod(z);
	  np++;
	}// for n
	mad(np,x,&mad1[l1],&mad2[l1]); // this is about baselines 
      }// for k1
    }// for j1
#pragma omp barrier   
  }// for i1    
  
  
 if((y =(float *)malloc(ntimes*(nsrc+1)*sizeof(float)))==NULL)
   {fprintf(stdout,"\t %s y1\n",err_msg);return(-1);}
 if((y1=(float *)malloc(nchans*nstokes*sizeof(float)))==NULL)
   {fprintf(stdout,"\t %s y1\n",err_msg);return(-1);}
 if((y2=(float *)malloc(nchans*nstokes*sizeof(float)))==NULL)
   {fprintf(stdout,"\t %s y1\n",err_msg);return(-1);}
 
 for(i1 = 0; i1 < ntimes; i1++){
   for(np=0,j1 = 0; j1 < nchans; j1++){
     for(k1 = 0; k1 < nstokes; k1++){
       l1 = k1+nstokes*(j1+nchans*i1);
       if (mad1[l1] > 0.0 && mad2[l1] > 0.0){
	 y1[np] = mad1[l1]; 
	 y2[np] = mad2[l1]; 
	 np++;
       }// for if
     }// for k1
   }// for j1 // average over channles & stokes 
   m1[i1] = median(np,y1);
   m2[i1] = median(np,y2);
 }// for i1
 
 if((y1=realloc(y1,ntimes*(nsrc+1)*sizeof(float)))==NULL)
   {fprintf(stdout,"\t %s y1\n",err_msg);return(-1);}
 if((y2=realloc(y2,ntimes*(nsrc+1)*sizeof(float)))==NULL)
   {fprintf(stdout,"\t %s y2\n",err_msg);return(-1);}
 
  for(i1=0; i1 <= nsrc; i1++){
    for(i2 =0; i2 < ntimes; i2++){
      y1[i2+ntimes*i1] = 0.0;
      y2[i2+ntimes*i1] = 0.0;
    }// for i2
    np1[i1]=0;
    mm1[i1]=0.0;
    mm2[i1]=0.0;
  }// for i1
  
  for(i1=0; i1 < nscans; i1++){
    for(i2 = tid[i1]; i2 < tid[i1+1]; i2++){ 
      l1 = np1[idmap[sid[i1]]]+ntimes*idmap[sid[i1]]; 
      y1[l1] = m1[i2];
      y2[l1] = m2[i2];
      np1[idmap[sid[i1]]]++; 
    }// for i2
  }// for i2
  
  for(i1=1; i1 <= nsrc; i1++){
    for(i2 = 0; i2 < np1[i1]; i2++)
      y[i2] = y1[i2+ntimes*i1];
    mm1[i1] = median(np1[i1],y);
    for(i2 = 0; i2 < np1[i1]; i2++)
      y[i2] = y2[i2+ntimes*i1];
    mm2[i1] = median(np1[i1],y);
    fprintf(stdout,"\t Src = %12s:  m1: %6.3f m2: %6.3f\n",src[i1].name,mm1[i1],mm2[i1]);  
  }// for i1
  
  for(n1=0,n2=0,i1=0; i1 < nscans; i1++){
    if (gs[i1] > 0){
      for(i2=tid[i1]; i2 < tid[i1+1]; i2++){
	for(j1=0; j1 < nchans; j1++){
	  for(k1=0; k1 < nstokes; k1++){
	    l1=k1+nstokes*(j1+nchans*i2);
	    for(n=0; n < nbaselines; n++){
	      z   = Cmplx(Visib[n+nbaselines*l1].re,Visib[n+nbaselines*l1].im);
	      tmp = fabs(cmod(z)-mm1[idmap[sid[i1]]])/mm2[idmap[sid[i1]]];
	      if (Visib[n+nbaselines*l1].w > 0.0  && (tmp> thrs_gmad)){
		Visib[n+nbaselines*l1].w = -1.0;
		n2++;
	      }// for if
	      n1++; 
	    }// for n
	  }// for k1
	}// for j1
      }// for i2
    }// for if
  }// for i1
  
  free(mad1); 
  free(mad2);
  
  fprintf(stdout,"\t data filtered by GLOBAL mad filter = %2.6f\n",(float)n2/n1); 
  
  free(y);
  free(y1);
  free(y2);
  free(m1);
  free(m2);
  free(mm1);
  free(mm2);
  
  return(0); 
}// end mad_filter_global 

int pre_mad_filter(){
  int i1, i2,j1, k1,n,n1,n2,l1,l2;
  float tmp,*mad1,*mad2;
  cmplx z; 
  
  if (!pmad_mode) return(0);  
  
  fprintf(stdout,"\t  pre MAD FILTERRING input data\n");
  
  sprintf(err_msg,"\t ERROR ! Memory allocation failed in pre_mad_filter for");
  
  if((mad1 = (float *)malloc(nscans*nchans*nstokes*nbaselines*sizeof(float)))==NULL)
    {fprintf(stdout,"\t %s mad1 !\n",err_msg);return(-1);}
  if((mad2 = (float *)malloc(nscans*nchans*nstokes*nbaselines*sizeof(float)))==NULL)
    {fprintf(stdout,"\t %s mad2 !\n",err_msg);return(-1);}
  
  compute_mad(mad1,mad2); // this is over all the time samples of a scan 
  
  for(n1=0,n2=0,i1 = 0; i1 < nscans; i1++){ 
    if (gs[i1] > 0){
      if (src[idmap[sid[i1]]].type == 0 || src[idmap[sid[i1]]].type == 1){
	for(i2 = tid[i1]; i2 < tid[i1+1]; i2++){ 
	  for(j1 = 0; j1 < nchans; j1++){
	    for(k1= 0; k1 < nstokes; k1++){
	      l1 = k1 +  nstokes*(j1+nchans*i1);
	      l2 = k1 +  nstokes*(j1+nchans*i2); 
	      for(n = 0; n < nbaselines; n++){
		z   = Cmplx(Visib[n+nbaselines*l2].re,Visib[n+nbaselines*l2].im);
		tmp = fabsf(cmod(z)- mad1[n+nbaselines*l1])/mad2[n+nbaselines*l1];
		if (Visib[n+nbaselines*l2].w > 0.0  && (tmp > thrs_pmad) ){
		  Visib[n+nbaselines*l2].w = -1.0;
		  n2++; 
		}// for if
		n1++;
	      }// for n 
	    } // for k1
	  }// for j1
	}// for i2
      }// for if
    }// for if
  }// for i1
  
  fprintf(stdout,"\t data filtered by pre-MAD filter = %2.6f \n",(float)n2/n1);
  
  free(mad1);
  free(mad2);
  
  return(0);
}

int post_mad_filter(){
  int i1,i2,j1, k1,n,np,n1,n2,l1,l2;
  float tmp1,tmp2,tmp3,*mm1,*mm2,*mad1,*mad2,*y1,*y2;
  float y[nbaselines],m1[nscans],m2[nscans]; 
  float thrshigh,thrslow;
  cmplx z;
 
 if (!qmad_mode) return(0);
  
  fprintf(stdout,"\t  post MAD FILTERRING all the data \n");
  sprintf(err_msg,"\t ERROR ! Memory allocation failed in post_mad_filter for");
  
  if((mad1 = (float *)malloc(nscans*nchans*nstokes*nbaselines*sizeof(float)))==NULL)
    {fprintf(stdout,"\t %s  mad1 !\n",err_msg); return(-1);}
  if((mad2 = (float *)malloc(nscans*nchans*nstokes*nbaselines*sizeof(float)))==NULL)
    {fprintf(stdout,"\t %s  mad2 !\n",err_msg); return(-1);}
  if((mm1  = (float *)malloc(nscans*nchans*nstokes*sizeof(float)))==NULL)
    {fprintf(stdout,"\t %s  mm1!\n",err_msg); return(-1);}
  if((mm2  = (float *)malloc(nscans*nchans*nstokes*sizeof(float)))==NULL)
    {fprintf(stdout,"\t %s  mm1!\n",err_msg); return(-1);}
  
  for(i1 = 0; i1 < nscans; i1++){ 
    for(j1 = 0; j1 < nchans; j1++){
      for(k1= 0; k1 < nstokes; k1++){
	l1 = k1 +  nstokes*(j1+nchans*i1);
	mm1[l1] = 0.0;
	mm2[l1] = 0.0;
	for(n  = 0; n < nbaselines; n++)  {
	  mad1[n+nbaselines*l1] = 0.0;
	  mad2[n+nbaselines*l1] = 0.0;
	}// for n
      }// for k1
    }// for j1
  }// for i1
  
  if(compute_mad(mad1,mad2))
    return(0);
 
  for(i1 = 0; i1 < nscans; i1++){ 
    for(j1 = 0; j1 < nchans; j1++){
      for(k1= 0; k1 < nstokes; k1++){
	l1 = k1 +  nstokes*(j1+nchans*i1);
	for(np=0, n= 0; n < nbaselines; n++){  
	  if (mad1[n+nbaselines*l1] > 0.0){
	    y[np] = mad1[n+nbaselines*l1];
	    np++;
	  }// for if
	}// for n
	mm1[l1] = median(np,y);
	for(np=0,n  = 0; n < nbaselines; n++){  
	  if (mad2[n+nbaselines*l1] > 0.0){
	    y[np] = mad2[n+nbaselines*l1];
	    np++;
	  }// for if
	}// for n
	mm2[l1] = median(np,y);
      }// for k1
    }// for j1
  }// for i1

  if((y1 =(float *)malloc(nchans*nstokes*sizeof(float)))==NULL)
    {fprintf(stdout,"\t %s  y1 !\n",err_msg); return(-1);}
  if((y2 =(float *)malloc(nchans*nstokes*sizeof(float)))==NULL)
    {fprintf(stdout,"\t %s  y2 !\n",err_msg); return(-1);}
  
  for(i1=0; i1 < nscans; i1++){
    m1[i1] = 0.0;
    m2[i1] = 0.0;
    for(np=0,j1  = 0; j1 < nchans; j1++){
      for(k1= 0; k1 < nstokes; k1++){
	l1  = k1 +  nstokes*(j1+nchans*i1);
	y1[np] = mm1[l1];
	y2[np] = mm2[l1];
	np++;
      }// for k1
    }// for j1
    m1[i1] = median(np,y1);
    m2[i1] = median(np,y2);
    if(gs[i1] > 0) fprintf(stdout,"\tSrc:%12s  m1:%6.3f  m2:%6.3f\n",src[idmap[sid[i1]]].name,m1[i1],m2[i1]); 
  }// for i1


  /* flag baselines if their median value is too low. */
  for(n1=0,n2=0,i1 = 0; i1 < nscans; i1++){ 
    if(src[idmap[sid[i1]]].type==2){
      thrshigh = thrs_qmadb;
      thrslow  = thrs_qmadblow;
    }else{
      thrshigh = thrs_qmadb;
      thrslow  = thrs_qmadb;
    }
    for(j1 = 0; j1 < nchans; j1++){
      for(k1= 0; k1 < nstokes; k1++){
	l1 = k1 +  nstokes*(j1+nchans*i1);
	for(np=0, n= 0; n < nbaselines; n++){  
	  tmp1= mad1[n+nbaselines*l1];
	  tmp2=(tmp1 -m1[i1])/m2[i1];
	  tmp3=(tmp1 - mad1[n+nbaselines*l1])/m2[i1];
	  if((tmp1/m1[i1] < thrs_qmadbmean)|| tmp2 > thrshigh || tmp3 > thrslow){
	    for(i2 = tid[i1]; i2 < tid[i1+1]; i2++){ 
	      l2 = k1 +  nstokes*(j1+nchans*i2); 
	      if(Visib[n+nbaselines*l2].w > 0.0){
		Visib[n+nbaselines*l2].w = -1.0;
		n2++;
	      }
	    }
	  }
	  n1 += tid[i1+1]-tid[i1];
	}
      }
    }
  }
  fprintf(stderr,"\t Data filtered because of bad baseline %6.2f\n", (float)n2/n1);
  
  if(flag_timerange(m1,m2))         /* Identify and flag band time ranges*/
    return(0);

  for(n1=0,n2=0,i1 = 0; i1 < nscans; i1++){ 
    if (gs[i1] > 0){
      if(src[idmap[sid[i1]]].type==2){
	thrshigh = thrs_qmad;
	thrslow  = thrs_qmadlow;
      }else{
	thrshigh = thrs_qmad;
	thrslow  = thrs_qmad;
      }
      for(i2 = tid[i1]; i2 < tid[i1+1]; i2++){ 
	for(j1 = 0; j1 < nchans; j1++){
	  for(k1= 0; k1 < nstokes; k1++){
	    l1 = k1 +  nstokes*(j1+nchans*i1);
	    l2 = k1 +  nstokes*(j1+nchans*i2); 
	    for(n  = 0; n < nbaselines; n++)  {
	      z    = Cmplx(Visib[n+nbaselines*l2].re,Visib[n+nbaselines*l2].im);
	      tmp1 = cmod(z);
	      tmp2 = (tmp1-m1[i1])/m2[i1];
	      tmp3 = (m1[i1]-tmp1)/m2[i1];
	      if (Visib[n+nbaselines*l2].w > 0.0  && 
		  ((tmp2 > thrshigh)||(tmp3 > thrslow)||(tmp1/m1[i1]<thrs_qmadmean))){
		Visib[n+nbaselines*l2].w = -1.0;
		n2++;
	      }// for if 
	      n1++;
	    }// for i2 
	  } // for n
	}// for k1
      }// for j1
    }// for if 
  }// for i1

  fprintf(stdout,"\t data filtered by post-MAD filter = %6.2f \n",(float)n2/n1);
  
  /* if the data for a particular channel/stokes is quite off from others 
     that can be flagged here   
     
  */

  
  free(y1);
  free(y2); 
  free(mm1);
  free(mm2); 
  free(mad1);
  free(mad2);
  
  
  return(0);
}

int compute_mad(float mad1[],float mad2[]){
  int i1,i2,i3,j1,k1,n,np,l1,l2;
  float y[max_times]; 
  cmplx z;
 
  omp_set_num_threads(nthreads);
  
  for(i1 = 0; i1 < nscans; i1++){
    for(j1 = 0; j1 < nchans; j1++){
      for(k1 = 0; k1 < nstokes; k1++){
	l2 = k1 + nstokes*(j1+nchans*i1);
	for(n=0; n < nbaselines; n++){
	  mad1[n+nbaselines*l2] = 0.0;
	  mad2[n+nbaselines*l2] = 0.0;
	} // for n
      } // for k1
    }// for j1
  }// for i1
  
  for(i1 = 0; i1 < nscans; i1++){
#pragma omp parallel for private(j1,k1,i2,i3,l1,l2,n,np,z,y)shared(i1,mad1,mad2)
    for(j1 = 0; j1 < nchans; j1++){
      for(i2  = 0; i2 < tpnt[i1]; i2++)
	y[i2] = 0.0; 
      for(k1= 0; k1 < nstokes; k1++){
	l1 = k1 + nstokes *(j1+nchans*i1);
	for(n  = 0; n < nbaselines; n++){
	  for(i2=0; i2 < tpnt[i1]; i2++)
	    y[i2] = 0.0;
	  for(np=0,i2 = tid[i1]; i2 <  tid[i1+1]; i2++){
	    l2 = k1+ nstokes *(j1+nchans*i2); 
	    i3 = i2 - tid[i1];
	    if (Visib[n+nbaselines*l2].w < 0.0) continue;
	    z     = Cmplx(Visib[n+nbaselines*l2].re,Visib[n+nbaselines*l2].im);
	    y[i3] = cmod(z);
	    np++; 
	  }// for i2
	  mad(np,y,&mad1[n+nbaselines*l1],&mad2[n+nbaselines*l1]); //jnc 9/sep/11 removed np>0 check
	} // for n
      }// for k1
    }// for j1
    // let us wait for everbody 
#pragma omp barrier
    // parallel section ends 
  }// for i1
  return(0);
}
