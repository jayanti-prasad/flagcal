#include <defs.h>
#include <flagcal.h>

float calflux(struct source, float);
  
int setjy(){
  int i;

  fprintf(stdout,"\t computing the flux for the fluxcal \n");

  for(i=1; i <= nsrc; i++){
    src[i].flux = 1.0;
  }
  
  for(i=1; i <= nsrc; i++){
    switch(src[i].type){
    case 0:
      src[i].flux = calflux(src[i],freq); 
      printf("\t Src: %12s Flx: %2.6f\n",src[i].name,src[i].flux); 
      break; 
    case 1:
      src[i].flux = 0.0;
      break;
    case 2:
      src[i].flux = 1.0;
      break;
    }
  }// for i2
  
  return(0);
}// end setjy 

int getjy(){
  int i1,i2,k1,l,l1,l2,m1,j,j1,np1,pidx,fidx,ipcal;
  int *flx_scan,flx_scans,*phs_scan,phs_scans;
  float *gain_ratio,*flx_vals,*scn_gain,*ratio_arr,*tmp_arr,*tmp_farr;
  float median (int n, float x[]);
  int   average(int np,float *x,float *mean,float *rms);
  int   debug=0;

  ipcal=0;
  for(i1=0; i1 <=nsrc; i1++){
    if (src[i1].type==1)
    ipcal = 1;
  }
  
  if(!ipcal) return(0); // if there is no phasecal return   

  fprintf(stdout,"\t Computing FLUX for phase cal \n"); 

  /* only for debug purposes */
  if(debug){
    if((g1=(cmplx  *)malloc(nablocks*nstokes*nchans*nants*sizeof(cmplx)) ) == NULL){
      fprintf(stdout,"\t getjy()memory allocation for  Visib");
      return(-1);
    }
    /* put in fake gain values */
    for(i1 = 0; i1 < nscans ; i1++){
      if ((src[idmap[sid[i1]]].type != 0) && (src[idmap[sid[i1]]].type != 1))continue;
      for(j1 = 1; j1 < nchans; j1++){
	for(k1 = 0; k1 < nstokes_half; k1++){
	  l2 = k1+nstokes_half *(j1+nchans*i1);
	  for(l = 0; l < nants; l++){
	    for(i2 = 0; i2 < nab[i1]; i2++){
	      j = abid_in[i2+max_nab*i1]; 
	      l1 = k1+nstokes * (j1+nchans*j);
	      g1[l+nants*l1].im=0.0;
	      if(src[idmap[sid[i1]]].type==0)
		g1[l+nants*l1].re=sqrt(src[idmap[sid[i1]]].flux);
	      else
		g1[l+nants*l1].re=1.0;
	    }// for i2
	  } // for l
	}// for k1
      } // for j1
    } //for i1 
  }// debug

  if((flx_scan=(int *)malloc(sizeof(int)*nscans))==NULL)
  { fprintf(stdout,"\t getjy()Malloc failed for flx_scan\n"); return 1;}
  if((phs_scan=(int *)malloc(sizeof(int)*nscans))==NULL)
  {fprintf(stdout,"\t getjy()Malloc failed for phs_id\n"); return 1;}
  
  flx_scans=phs_scans=0.0;
  for(i1=0; i1 < nscans; i1++)
  { if (src[idmap[sid[i1]]].type==0){flx_scan[flx_scans]=i1;flx_scans++;}
    if (src[idmap[sid[i1]]].type==1){phs_scan[phs_scans]=i1;phs_scans++;}
  }

  if(flx_scans==0)
  { fprintf(stdout,"\t getjy() Can't calibrate: No Flux calibrator scans\n"); return 1;}
  if(phs_scans==0)
  { fprintf(stdout,"\t getjy() Can't calibrate: No Phase calibrator scans\n"); return 1;}

  if((gain_ratio=(float *)calloc(flx_scans*(flx_scans+phs_scans),sizeof(float)))==NULL)
  {fprintf(stdout,"\t getjy()Malloc failed for gain_ratio\n"); return 1;}
  if((flx_vals=(float *)calloc(flx_scans*(flx_scans+phs_scans),sizeof(float)))==NULL)
  {fprintf(stdout,"\t getjy()Malloc failed for flx_vals\n"); return 1;}
  if((tmp_farr=(float *)calloc(flx_scans*(flx_scans+phs_scans),sizeof(float)))==NULL)
  {fprintf(stdout,"\t getjy()Malloc failed for tmp_farr\n"); return 1;}
  if((scn_gain=(float *)calloc(nscans*nants*nchans*nstokes_half,sizeof(float)))==NULL)
  {fprintf(stdout,"\t getjy()Malloc failed for scn_gain\n"); return 1;}
  if((ratio_arr=(float *)calloc(nscans*nants*nchans*nstokes_half,sizeof(float)))==NULL)
  {fprintf(stdout,"\t getjy()Malloc failed for ratio_array\n"); return 1;}
  if((tmp_arr=(float *)calloc(max_nab,sizeof(float)))==NULL)
  {fprintf(stdout,"\t getjy()Malloc failed for tmp_arr\n"); return 1;}


  /* get median gain per scan for flx,phs cal */
  for(i1 = 0; i1 < nscans ; i1++){
    if ((src[idmap[sid[i1]]].type != 0) && (src[idmap[sid[i1]]].type != 1))continue;
    for(j1 = 1; j1 < nchans; j1++){
      for(k1 = 0; k1 < nstokes_half; k1++){
	l2 = k1+nstokes_half *(j1+nchans*i1);
	for(l = 0; l < nants; l++){
	  int np=0;
	  for(i2 = 0; i2 < nab[i1]; i2++){
	    float g;
	    j = abid_in[i2+max_nab*i1]; 
	    l1 = k1+nstokes * (j1+nchans*j);
	    if(!debug && bad_ant[l+nants*i1] > thrs_ant)continue;
	    if((g=cmod(g1[l+nants*l1]))<= 0.0) continue;
	    tmp_arr[np] =g ;np++;
	  }// for i2
	  if(np>1)scn_gain[l+nants*l2]=median(np,tmp_arr); /* initialized to 0.0*/
	} // for l
      }// for k1
    } // for j1
  } //for i1 


  /* compute the gain ratios for all flx, phs cal scans */

  for(np1=0,m1=0;m1<flx_scans;m1++){
    for(i1 = 0; i1 < nscans ; i1++){
      int np=0;
      if ((src[idmap[sid[i1]]].type != 0) && (src[idmap[sid[i1]]].type != 1))continue;
      for(j1 = 1; j1 < nchans; j1++){
	for(k1 = 0; k1 < nstokes_half; k1++){
	  l1 = k1+nstokes_half *(j1+nchans*i1);
	  l2 = k1+nstokes_half *(j1+nchans*flx_scan[m1]);
	  for(l = 0; l < nants; l++){
	    pidx=l+nants*l1;
	    fidx=l+nants*l2;
	    if(scn_gain[pidx]>0.0 && scn_gain[fidx] > 0.0)
	    { ratio_arr[np] = scn_gain[fidx]/scn_gain[pidx];np++;}
	  }// for l
	} // for k1
      }// for j1
      if(np>0)gain_ratio[np1]=median(np,ratio_arr);
      np1++;
    } // for i1
  }// for m1


  /* compute the cross calibrated flux per scan,FlxCal combination */
  for(np1=0,m1=0;m1<flx_scans;m1++){
    fidx=idmap[sid[flx_scan[m1]]];
    printf("\t [FlxCal:%s Scn: %4d Flx %10.3f]\n",src[fidx].name,flx_scan[m1],src[fidx].flux);
    for(i1=0;i1<nscans;i1++){
      if ((src[idmap[sid[i1]]].type != 0) && (src[idmap[sid[i1]]].type != 1))continue;
      printf("\t Scan: %4d Src: %16s",i1,src[idmap[sid[i1]]].name);
      flx_vals[np1] = (gain_ratio[np1]*gain_ratio[np1])*src[fidx].flux;
      printf(" %10.3f ",flx_vals[np1]);
      if(flx_vals[np1]>0.0) printf("\n");
      else printf("\t [ignored]\n");
      np1++;
    }
    printf("\n");
  }

  /* get the average flux for each PhsCal */
  for(j1=1;j1<=nsrc;j1++){
    int np=0;
    float mean,rms;
    if(src[j1].type !=1) continue;
    for(np1=0,m1=0;m1<flx_scans;m1++){
      for(i1=0;i1<nscans;i1++){
	if ((src[idmap[sid[i1]]].type != 0) && (src[idmap[sid[i1]]].type != 1))continue;
	if(idmap[sid[i1]]==j1 && flx_vals[np1]>0.0) /* flux is 0 if no data in scan */
	  {tmp_farr[np]=flx_vals[np1];np++;}
	np1++;
      }
    }
    if(average(np,tmp_farr,&mean,&rms)) return 1;
    printf("\t Src: %s %10.3f +- %10.3f\n",src[j1].name,mean,rms);
    src[j1].flux=mean;
  }

  if(debug) free(g1);
  free(flx_scan); free(phs_scan); free(gain_ratio);   
  free(flx_vals); free(tmp_farr);  free(scn_gain);
  free(ratio_arr);  free(tmp_arr);

  return 0;
} // flux_cal ends 


int gain_correct(){
  int i1, i2,j1, k1, j,l,l1; 

  for(i1=0; i1 < nscans; i1++){
    if (src[idmap[sid[i1]]].type != 2){
      for(i2 = 0; i2 < nab[i1]; i2++){
	j = abid_in[i2+max_nab*i1]; 
	for(j1=0; j1 < nchans; j1++){
	  for(k1=0; k1 < nstokes; k1++){
	    l1 = k1+nstokes*(j1+nchans*j);
	    for(l=0; l < nants; l++){
	      g1[l+nants*l1].re *= sqrt(src[idmap[sid[i1]]].flux);
	      g1[l+nants*l1].im *= sqrt(src[idmap[sid[i1]]].flux);				
	    } // for l
	  }// for k1
	} // for j1
      }// for i2
    } // for if 
  } // for i1
  return(0);
}

float calflux(struct source fcal, float f){
  int i, nfcal = 6, found=0;
  struct source *stdfcal;
  float A[6], B[6],C[6];

  if((stdfcal = (struct source *)malloc(nfcal*sizeof(struct source)))==NULL)
  {fprintf(stdout,"calflux()failled to malloc for stdfcal\n");return 1;}
  
  strcpy(stdfcal[0].name,"3C286");
  A[0] =  1.480000;
  B[0] =  0.292000;
  C[0] = -0.124000;

  strcpy(stdfcal[1].name,"3C48");
  A[1]= 2.345000;
  B[1]= 0.071000;
  C[1]=-0.138000; 

  strcpy(stdfcal[2].name,"3C147");
  A[2] =  1.766000;
  B[2] =  0.447000;
  C[2] = -0.184000; 
  
  strcpy(stdfcal[3].name,"3C138");
  A[3] =  2.00900;
  B[3] = -0.07176;
  C[3] = -0.08620; 

  strcpy(stdfcal[4].name,"1934-638");
  A[4] = -23.839000;
  B[4] =  19.569000;
  C[4] = -4.8168000; 

  strcpy(stdfcal[5].name,"3C295");
  A[5] =  1.485000;
  B[5] =  0.759000;
  C[5] = -0.255000; 

  
  for(i=0; i < nfcal; i++){
    stdfcal[i].flux = A[i] + B[i] * log10(f)+ C[i] * log10(f)*log10(f);
    stdfcal[i].flux = pow(10.0,stdfcal[i].flux);
  } 
  
  for(i=0; i < nfcal; i++){
    if(strcmp(fcal.name,stdfcal[i].name)==0){
      fcal.flux = stdfcal[i].flux;
      found = 1;
    }
  }

  if(found == 0){
    fprintf(stdout,"\t your flux cal %s is not found please give its flux in Jy\n",fcal.name);
    scanf("%f",&fcal.flux);
    fprintf(stdout,"\t you have given = %2.6f\n",fcal.flux);
  }

  return(fcal.flux);
}// end calflux 


