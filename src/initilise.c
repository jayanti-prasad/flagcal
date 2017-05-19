#include <defs.h>
#include <flagcal.h>
/*------------------------------------------------------------------------------------
  -   This program sets the properties of the source (type etc.) and arrange the 
      data in the form of flagging (sblocks) and solution blocks (ablocks).
  -   Memory allocation for some global varibles is also done here. 
  
  Jayanti Prasad
  Sat Jul 24 09:58:42 IST 2010
  
  --------------------------------------------------------------------------------------*/
int setup_idmap(){
  /* JNC 5/Sep/11 Maps from the source id as listed in the random parameters
     of the UV data to the index in the source table src
  */
  int i,j;
  for(i=0;i<nscans;i++){
    for(j=1;j<=nsrc;j++){
      if(sid[i]==src[j].id){idmap[sid[i]]=j;break;}
      if(j==nsrc){
	fprintf(stderr,"Cannot find a match for Src_ID %ld in ScanNo %d\n",
		sid[i],i);
	return(-1);
      }
    }
  }
  return(0);
}
int initilise(){
  int i,j,k,i1,i2,block_size;
  
  /* ask for the identification of the souces */ 
  
  for(i=0; i <3; i++){
    switch(i){
    case 0:
      strcpy(srtype[i].val,"FCAL");
      break;
    case 1:
      strcpy(srtype[i].val,"PCAL");
      break;
    case 2:
      strcpy(srtype[i].val,"TSRC");
      break;  
    }// for switch 
  }// for i 

   for(i = 1; i <= nsrc; i++){
     //      src[i].id = i-1; ID now read from the table JNC 5Sep11
      src[i].flux = 1.0;
      src[i].type = -1;
   }// for i
   if(setup_idmap()) return -1;

  if(access(srcfile, F_OK) !=-1)
    read_infile(0); 
  
  fprintf(stderr,"\t ----------------source info ----------------\n");
  
  for(i = 1; i <= nsrc; i++){
    if (src[i].type < 0){
      fprintf(stdout,"\n \t What is  source %s ? \n",src[i].name);
      fprintf(stdout,"\t 0 : FCAL  \n");
      fprintf(stdout,"\t 1 : PCAL  \n");
      fprintf(stdout,"\t 2 : TSRC  \n");
      scanf("%d",&src[i].type);
    }// for if
  }// for i2
 
  
  for(i=0; i < nscans; i++){
    if (src[idmap[sid[i]]].type < 0 || src[idmap[sid[i]]].type > 2){
      fprintf(stdout,"\t  Scn:%2d   Src <id:%2ld  type: UNKNOWN  Nm: %s>\n",
	      i,sid[i],src[idmap[sid[i]]].name);
      return(-1);  
    }// for if
   fprintf(stdout,"\t  Scn:%2d   Src <id:%2ld  type: %s  Nm: %s>\n",
	   i,sid[i],srtype[src[idmap[sid[i]]].type].val,src[idmap[sid[i]]].name);
  }// for i 
  fprintf(stdout,"\t ---------------------------------------------------\n");
  
  /* indexing for the flagging/averaging blocks */
  
  for(i=0; i < nscans; i++)
    tpnt[i] = tid[i+1] - tid[i];
  
  max_times = 0;
  for (i1=0; i1 < nscans; i1++){
    if(tpnt[i1] > max_times)
      max_times = tpnt[i1];
  }
  
  max_nsb  =  0;     
  max_nab  =  0;     
  nsb      =   (int   *)malloc(nscans*sizeof(int));
  last_sbl =   (int   *)malloc(nscans*sizeof(int));
  nab      =   (int   *)malloc(nscans*sizeof(int));
  last_abl =   (int   *)malloc(nscans*sizeof(int));
  
  nsblocks = 0;
  nablocks = 0;
  for(i1 = 0; i1 < nscans; i1++){
    nsb[i1] = tpnt[i1]/nflag; 
    nab[i1] = tpnt[i1]/nsol; 
    last_sbl[i1] =  nflag + (tpnt[i1]-nsb[i1] * nflag);
    last_abl[i1] =  nsol  + (tpnt[i1]-nab[i1] * nsol);
    if (nsb[i1] > max_nsb)
      max_nsb = nsb[i1];
    if (nab[i1] > max_nab) 
      max_nab = nab[i1];
    nsblocks +=nsb[i1]; 
    nablocks +=nab[i1]; 
  }//for i1
  nsblocks+=1;
  nablocks+=1;
  
  sbid_in  = (int *)malloc(nscans*max_nsb*sizeof(int));
  sb_start = (int *)malloc(nsblocks*sizeof(int)); 
  sb_end   = (int *)malloc(nsblocks*sizeof(int));
  abid_in  = (int *)malloc(nscans*max_nab*sizeof(int));
  ab_start = (int *)malloc(nablocks*sizeof(int)); 
  ab_end   = (int *)malloc(nablocks*sizeof(int));
  abmt     = (float *)malloc(nablocks*sizeof(float));
 
  if((g1=(cmplx  *)malloc(nablocks*nstokes*nchans*nants*sizeof(cmplx)) ) == NULL){
    perror("memory allocation for  Visib");
    return(-1);
   }// for if
  
  for(i1 = 0; i1 < nscans; i1++){
    for(i2 = 0; i2 < max_nsb; i2++){
      sbid_in[i2+max_nsb*i1] = 0; 
    }// for i2
    for(i2 = 0; i2 < max_nab; i2++){
      abid_in[i2+max_nab*i1] = 0; 
    }// for i2
  }// for i1
  
  j = 0; 
  k = 0;
  for(i1 = 0; i1 < nscans; i1++){
    for(i2 = 0; i2 < nsb[i1]; i2++){
      i = i2 + max_nsb * i1; 
      sbid_in[i] = j;
      j++;
    }// for i2
    for(i2 = 0; i2 < nab[i1]; i2++){
      i = i2 + max_nab * i1; 
      abid_in[i] = k;
      k++;
    }// for i2
  }// for i1
  
  for(i1 = 0; i1 < nscans; i1++){
    for(i2 = 0; i2 < nsb[i1]; i2++){	  
      j = sbid_in[i2+max_nsb*i1];
      if (i2  == (nsb[i1]-1)){
	block_size = last_sbl[i1];
      }else{
	block_size = nflag;
      }
      sb_start[j]  =  tid[i1]+nflag*i2;
      sb_end[j]    =  sb_start[j] + block_size;
    }// for i2
    for(i2 = 0; i2 < nab[i1]; i2++){	  
      k = abid_in[i2+max_nab*i1];
      if (i2  == (nab[i1]-1)){
	block_size = last_abl[i1];
      }else{
	block_size = nsol;
      }
      ab_start[k]  =  tid[i1]+nsol*i2;
      ab_end[k]    =  ab_start[k] + block_size;
    }// for i2
  }// for i1
  
  for(i1 = 0; i1 < nscans; i1++){
    for(i2 = 0; i2 < nab[i1]; i2++){
      j = abid_in[i2+max_nab*i1]; 
      abmt[j] = (time_data[ab_start[j]] + time_data[ab_end[j]])/2.0;
    } // for i2 
  } // for i1 
  
  write_index();
  
  return(0);
  
} // index_data ends 


int clip(){
  int i1, j1, k1, n, l1,n1=0,n2=0;
  cmplx z; 
  
  if (!clip_mode) return(0); 
  
  for(i1 = 0; i1 < ntimes; i1++){
    for(j1 = 0; j1 < nchans; j1++){
      for(k1 = 0; k1 < nstokes; k1++){
        l1=k1+nstokes*(j1+nchans*i1); 
	for(n = 0; n < nbaselines; n++){
	  z.re = Visib[n+nbaselines*l1].re; 
	  z.im = Visib[n+nbaselines*l1].im;
	  if(cmod(z) > thrs_clip){
	    Visib[n+nbaselines*l1].w = -1.0;
	    n2++;
	  }// for if
          n1++;
        }// for n
      }// for k1
    }// for j1
  }// for i1
  
  fprintf(stdout,"\t data clipped = %6.3f\n",(float)n2/n1); 
  
  return(0); 
}

int init_flag(){
  /*  This is the place where you can apply flagging on the basis of your
   *  prior information about the bad antennas, basebaselines & bad channels */
  int i,l,m,n,n1,n2,i1,j1,k1,l1,i2;
 
  ba=(int *)malloc(nants*sizeof(int));
  bb=(int *)malloc(nbaselines*sizeof(int));
  bc=(int *)malloc(nchans*sizeof(int));
  gs=(int *)malloc(nchans*sizeof(int));
  bs=(int *)malloc(nscans*sizeof(int));
  rc=(int *)malloc(nscans*sizeof(int));
  lc=(int *)malloc(nscans*sizeof(int));
  
  for(i = 0; i < nants; i++)
    ba[i] = 1;
  for(i=0; i < nbaselines; i++)
    bb[i] = 1;
  for(i=0; i < nchans; i++)
    bc[i] = 1;
  for(i=0; i < nscans; i++)
    bs[i] = 1;
  for(i=0; i < nscans; i++)
    gs[i] = 1;	  
  
  for(i=0; i < nscans; i++){
    lc[i]=0;
    rc[i]=0;
  }
  
  // now a calib file(optional) also can be passed to overrule the 
  // gain intrapolation.
  
  if(access(calfile, F_OK) !=-1){
    read_infile(1); // read cal file 
    for(i=0; i < nscans; i++)
     printf("scan=%d lc=%d rc=%d\n",i,lc[i],rc[i]);	   
  }
  
  nbadtimes=0;
  if(access(flgfile, F_OK) !=-1){
    read_infile(2); // read flag file
    
    for(n1=0,n2=0,i1=0; i1 < nscans; i1++){
      for(i2 = tid[i1]; i2 < tid[i1+1]; i2++){
	for(j1=0; j1 < nchans; j1++){
	  for(k1=0; k1 < nstokes; k1++){
	    l1 = k1+nstokes*(j1+nchans*i2);
	    for(n = 0; n < nbaselines; n++){
	      l = (id[n]+1)/nants;
	      m = (id[n]+1) - nants *l;
	      if(ba[l] < 0 || ba[m] < 0 || bb[n] < 0  || bc[j1] < 0 || bs[i1] < 0){
		Visib[n+nbaselines*l1].w = -1.0;
		n2++;
	      }// for if
	      n1++;
	    }// for n
	  }// for k1
	}// for j1
      }// for i2
    }// for i1
    fprintf(stdout,"\t Initial flagging [INPUT] = %2.6f\n",(float)n2/n1);
  } 
  
 return(0); 
}


