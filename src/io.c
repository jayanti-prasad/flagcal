#include <defs.h>
#include<flagcal.h>
#include <sys/unistd.h>

/*-------------------------------------------------------------------------
  This file contain the functions which are used to read command line options,
  in files. It also contains functions which are used to read/write 
  auxiliary file, like flag files, gain file and the index files.
                             JAAYNTI PRASAD, Fri Apr 22 11:54:37 IST 2011
http://www.iucaa.ernet.in/~jayanti/flagcal.html
  --------------------------------------------------------------------------*/

int get_timerange(char *string1){
  char *q;
  int   d1,h1,m1,s1;
  int   d2,h2,m2,s2;

  if((q=strchr(string1,'='))==NULL){
    fprintf(stderr,"Missing '=' in %s\n",string1);
    return 1;
  }
  q++;
  if(sscanf(q,"%d,%d,%d,%d,%d,%d,%d,%d",&d1,&h1,&m1,&s1,
	    &d2,&h2,&m2,&s2)!=8){
    fprintf(stderr,"%s has the wrong format\n",string1);
    return 1;
  }

  badtimes[nbadtimes].start = d1+((float)h1)/24.0+((float)m1)/(60.0*24.0)
                                +((float)s1)/(24.0*3600.0);
  badtimes[nbadtimes].end   = d2+((float)h2)/24.0+((float)m2)/(60.0*24.0)
                                +((float)s2)/(24.0*3600.0);
  nbadtimes++;

  if(nbadtimes==MAX_BAD_TIMES){
    fprintf(stderr,"MAX_BAD_TIMES %d exceeded!\n",MAX_BAD_TIMES);
    return 1;
  }
  return 0;
}
int read_input(int nargc, char *nargv[]){
  int i;
  bchan = 0;
  echan = 0;
  src_out  = 0;  // all source 
  nthreads = 1;  
  imode    = 0;
 
   fprintf(stdout,"reading command line options\n");
  
  // this is the default parameter file 
  sprintf(parafile,"parameters.in");
  
  if(nargc < 7){
    print_message(); 
    return(-1);
  }else{
    fprintf(stdout,"\t -------- command line options ------------- \n");
    for (i=0; i < nargc; i++){
      if (!strcmp(nargv[i],"-i")){	      
	sprintf(infits,"%s",nargv[++i]);
	fprintf(stdout,"\t input file        = %s\n",infits);
      }if (!strcmp(nargv[i],"-o")){
	sprintf(outfits,"%s",nargv[++i]);
	fprintf(stdout,"\t output file       = %s\n",outfits);
      }if(!strcmp(nargv[i],"-p")){
	sprintf(parafile,"%s",nargv[++i]);
	fprintf(stdout,"\t parameter file    = %s\n",parafile);
      } if(!strcmp(nargv[i],"-s")){
        sprintf(srcfile,"%s",nargv[++i]);
        fprintf(stdout,"\t source file       = %s\n",srcfile);
      }if(!strcmp(nargv[i],"-f")){
        sprintf(flgfile,"%s",nargv[++i]);
        fprintf(stdout,"\t flag file         = %s\n",flgfile);
      }if(!strcmp(nargv[i],"-c")){
        sprintf(calfile,"%s",nargv[++i]);
        fprintf(stdout,"\t calibration  filr = %s\n",calfile);
      }if(!strcmp(nargv[i],"-g")){
        sprintf(gainfile,"%s",nargv[++i]);
        fprintf(stdout,"\t gain filr         = %s\n",calfile);
      } if(!strcmp(nargv[i],"-n") ){
	nthreads = atoi(nargv[++i]);
	fprintf(stdout,"\t # threads         = %d\n",nthreads);
      }if(!strcmp(nargv[i],"-bchan") ){
	bchan = atoi(nargv[++i]);
	fprintf(stdout,"\t bchan             = %d\n",bchan);  
      }if(!strcmp(nargv[i],"-echan")){
	echan = atoi(nargv[++i]);
	fprintf(stdout,"\t echan             = %d\n",echan);
      }if(!strcmp(nargv[i],"-osrc")){
	src_out= atoi(nargv[++i]);
	fprintf(stdout,"\t output source     = %d\n",src_out);
      }if(!strcmp(nargv[i],"-imode")){
	imode = atoi(nargv[++i]);
	fprintf(stdout,"\t imode             = %d\n",imode);
      }
    }// for
    if(strlen(infits)==0){
      fprintf(stdout,"\t ERROR ! You have not given a valid INPUT file \n");
      return(-1);
    }if(strlen(outfits)==0){
      fprintf(stdout,"\t ERROR ! You have not given a valid OUTPUT file \n");
      return(-1);
    }if(strlen(parafile)==0){
      fprintf(stdout,"\t ERROR ! You have not given a valid PARAMETER file \n");
      return(-1);
    }if (echan !=0 && bchan !=0  &&  bchan > echan ){
      fprintf(stdout,"\t ERROR ! echan < bchan !\n");
      return(-1);
    } if(access(outfits, F_OK)!= -1){
      printf("\t ERROR ! file %s exists, delete that !\n",outfits);
        return(-1);
    }
    return(0);
  }// for if
}

int read_parameters(char parafile[]){
  int i,n1=23,n2=28;
  FILE *inpfile;
  char string1[max_len],c[max_len], *p,*p1; 
  int  iry[n1];
  float fry[n2];
  //struct sarray  ikwrd[n1],fkwrd[n2], switch_val[2];
  struct sarray  switch_val[2];
  
  sprintf(switch_val[0].val,"OFF");
  sprintf(switch_val[1].val,"ON");

 
  char *ikwrd[23]={"nsol","niter","ref_ant","clip_mode","gmad_mode","pmad_mode","nflag","aflag_mode","bflag_mode","cflag_mode","qmad_mode","rfi_mode_t","nriter_t","win_min_t","dwin_t","rfi_mode_f","nriter_f","win_min_f","dwin_f","range_t_mode","win_min_rt","dwin_rt","nriter_rt"};
  char *fkwrd[28]={"alpha","eps","thrs_clip","thrs_gmad","thrs_pmad","thrs_vsr","thrs_ngood","thrs_ant","thrs_base","thrs_chan","thrs_qmad","thrs_low_qmad","thrs_mean_qmad","thrs_b_qmad","thrs_b_low_qmad","thrs_b_mean_qmad","thrs_comp","thrs1_comp","thrs_peak_t","thrs_low_peak_t","dpt","thrs_peak_f","thrs_low_peak_f","dpf","thrs_range_t","thrs_low_range_t","thrs_mean_range_t","dprt"};


  /*    
  inpfile=fopen("keywords.inc","r");
  for(i =0; i < n1; i++){
    fgets(c,max_len, inpfile);  
    p=c+(strlen(c)-1);
    while(isspace(*p)){*p='\0';p--;}
    if(strlen(c)>31) // length defined in sarray
      { fprintf(stderr,"keyword %s too long\n",c);  return 1; }
    strcpy(ikwrd[i].val,c);
  } 
  for(i =0; i < n2; i++){
    fgets(c,max_len,inpfile);
    p=c+(strlen(c)-1);
    while(isspace(*p)){*p='\0';p--;}
    if(strlen(c)>31) // length defined in sarray
      { fprintf(stderr,"keyword %s too long\n",c);  return 1; }
    strcpy(fkwrd[i].val,c);
  } 
  fclose(inpfile); 
  */

  if(access(parafile, F_OK) != -1){
    fprintf(stdout,"\t Reading parameter file %s \n",parafile);
    inpfile=fopen(parafile,"r");
    fprintf(stdout,"\t---------------parameters--------------------\n");
    while(fgets(c,max_len, inpfile)!=NULL) { 
      sprintf(string1,"%s",c);
      //fprintf(stdout,"read %s\n",string1);
      p=string1;
      while(isspace(*p))p++;
      if(!strlen(p)) continue;
      p = strchr(string1,'!');
      if (p == NULL){
	for(i=0; i < n1; i++){
          if((p1 = strchr(ikwrd[i], '\n')) != NULL)
            *p1 = '\0';
	  if (strstr(string1,ikwrd[i]) != NULL){
	    p = strchr(c,'=');
	    iry[i]= atoi(p+1);
	    //fprintf(stdout,"\t %s         = %d\n",ikwrd[i].val,iry[i]); 
	  } 
	}// for i
	for(i=0; i < n2 ; i++){
	  if((p1 = strchr(fkwrd[i], '\n')) != NULL)
	    *p1 = '\0';
	  if (strstr(string1,fkwrd[i]) != NULL){
	    p = strchr(c,'=');
	    fry[i]= atof(p+1);
	    //fprintf(stdout,"\t %s         = %2.6f\n",fkwrd[i].val,fry[i]); 
	  } 
	}// for 
      }// for if
    }// for while
  }// for if	
  
  nsol        = iry[0];
  niter       = iry[1];
  ref_ant     = iry[2];
  clip_mode   = iry[3];
  gmad_mode   = iry[4];
  pmad_mode   = iry[5];
  nflag       = iry[6];
  aflag_mode  = iry[7];
  bflag_mode  = iry[8];
  cflag_mode  = iry[9];
  qmad_mode   = iry[10];
  rfi_mode_t  = iry[11];
  nriter_t    = iry[12];
  win_min_t   = iry[13];
  dwin_t      = iry[14];
  rfi_mode_f  = iry[15];
  nriter_f    = iry[16];
  win_min_f   = iry[17];
  dwin_f      = iry[18];
  range_t_mode= iry[19];
  win_min_rt  = iry[20];
  dwin_rt     = iry[21];
  nriter_rt   = iry[22];

  
  alpha             =  fry[0];
  eps               =  fry[1];
  thrs_clip         =  fry[2];
  thrs_gmad         =  fry[3];
  thrs_pmad         =  fry[4];
  thrs_vsr          =  fry[5];
  thrs_ngood        =  fry[6];
  thrs_ant          =  fry[7];
  thrs_base         =  fry[8];
  thrs_chan         =  fry[9];
  thrs_qmad         =  fry[10];
  thrs_qmadlow      =  fry[11];  
  thrs_qmadmean     =  fry[12];  
  thrs_qmadb        =  fry[13];
  thrs_qmadblow     =  fry[14];  
  thrs_qmadbmean    =  fry[15];  
  thrs_comp         =  fry[16];
  thrs1_comp        =  fry[17];
  thrs_peak_t       =  fry[18];
  thrs_peak_tlow    =  fry[19];
  dthrs_peak_t      =  fry[20];
  thrs_peak_f       =  fry[21];
  thrs_peak_flow    =  fry[22];
  dthrs_peak_f      =  fry[23];
  thrs_range_t      =  fry[24];
  thrs_low_range_t  =  fry[25];
  thrs_mean_range_t =  fry[26];
  dthrs_rt          =  fry[27];
  
  printf("\t (1)---calibration parameters-----------------\n");
  printf("\t nsol = %d, niter = %d, ref_an = %d, nflag = %d\n",nsol,niter,ref_ant,nflag);
  printf("\t eps = %6.2e, alpha = %2.6f\n",eps,alpha); 
  printf("\t    --- flagging parameters-----------------\n");
  printf("\t (2) clip : %s, clip threshold = %6.2e\n",switch_val[clip_mode].val,thrs_clip);
  printf("\t (3) gloabl mad filtering: %s, threshold = %6.2f\n",switch_val[gmad_mode].val,thrs_gmad);
  printf("\t (4) pre-mad filtering : %s, threshold = %6.2f\n",switch_val[pmad_mode].val,thrs_pmad);
  printf("\t (5) vsr flagging [always ON]:threshold = %6.2f, thrs_ngood = %6.2f\n",thrs_vsr,thrs_ngood);
  printf("\t (6) antena flagging : %s, threshold = %6.2f\n",switch_val[aflag_mode].val,thrs_ant);
  printf("\t (7) baseline flagging : %s, threshold = %6.2f\n",switch_val[bflag_mode].val,thrs_base);
  printf("\t (8) channel flagging : %s, threshold = %6.2f\n",switch_val[cflag_mode].val,thrs_chan);
  printf("\t (9) post mad filtering : %s\n",switch_val[qmad_mode].val);
  printf("\t\t thrs_qmad=%6.2f thrs_qmadlow=%6.2f, thrs_qmadmean = %6.2f\n",thrs_qmad,thrs_qmadlow,
	 thrs_qmadmean);
  printf("\t\t thrs_qmadb=%6.2f thrs_qmadblow=%6.2f, thrs_qmadbmean = %6.2f\n",thrs_qmadb,thrs_qmadblow,
	 thrs_qmadbmean);
  printf("\t\t thrs_comp=%2.6f, thrs1_comp = %6.2f\n",thrs_comp,thrs1_comp);
  printf("\t (10) remove RFI peak in time   : %s\n",switch_val[rfi_mode_t].val);
  if (rfi_mode_t==1){
    printf("\t              iter            = %d\n",nriter_t);
    printf("\t              win_start       = %d\n",win_min_t);
    printf("\t              dwin            = %d\n",dwin_t);
    printf("\t              thrs_peak_t     = %2.2f\n",thrs_peak_t);
    printf("\t              thrs_peak_tlow  = %2.2f\n",thrs_peak_tlow);
    printf("\t              dpt             = %2.2f\n",dthrs_peak_t);
  }
  printf("\t (11) remove RFI peak in frequency    : %s\n",switch_val[rfi_mode_f].val);
  if (rfi_mode_f==1){
    printf("\t               iter           = %d\n",nriter_f);
    printf("\t               win_start      = %d\n",win_min_f);
    printf("\t               dwin           = %d\n",dwin_f);
    printf("\t               thrs_peak_f    = %2.2f\n",thrs_peak_t);
    printf("\t              thrs_low_peak_f = %2.2f\n",thrs_peak_t);
    printf("\t               dpf            = %2.2f\n",  dthrs_peak_f);
  }
  printf("\t (9) time range filtering : %s\n",switch_val[range_t_mode].val);
  printf("\t\t thrs_range_t=%6.2f thrs_low_range_t=%6.2f, thrs_mean_range_t = %6.2f\n",
	 thrs_range_t,thrs_low_range_t, thrs_mean_range_t);
  printf("\t\t dthrs_rt=%6.2f win_min_rt=%6d, dwin_rt = %6d nriter_rt=%6d\n",
	 dthrs_rt,win_min_rt,dwin_rt,nriter_rt);
  return(0);
}

int read_infile(int icode){
  int i,i1,j1,k1,a1,a2;
  FILE *inpfile;
  char string1[max_len],c[max_len], *p;
  
  switch (icode){
  case 0:
    printf("\t scanning  source list from file: %s \n",srcfile);
    inpfile=fopen(srcfile,"r");
    while(fgets(c, 100, inpfile)!=NULL) {
      sprintf(string1,"%s",c);
      p = strchr(string1,'!');
      if (p == NULL){
	for(i=1; i <=nsrc; i++){
	  if (strstr(string1,src[i].name) != NULL){
	    p = strchr(c, '=');
	    src[i].type = atoi(p+1);
	    switch(src[i].type){
	    case 0:
	      fprintf(stdout,"\t Src: %12s Type: FCAL \n",src[i].name);
	      break;
	    case 1:
	      fprintf(stdout,"\t Src: %12s Type: PCAL \n",src[i].name);
	      break;
	    case 2:
	      fprintf(stdout,"\t Src: %12s Type: TSRC \n",src[i].name);
	      break;
	    }// for switch 
	  }// for if
	}// for i 
      }// for if
    }// for while
    fclose(inpfile);
    return(0);
    break;
  case 1:
    printf("\t Reading CALIB info from the input  file: %s \n",calfile);
    inpfile = fopen(calfile,"r");
      while(!feof(inpfile)){
	fscanf(inpfile,"%d %d %d\n",&i1,&j1,&k1);
	if(i1 < 0 || i1 > nscans-1 || src[idmap[sid[i1]]].type !=2){
	  fprintf(stderr,"Error: first entry of your calib file wrong !\n");
	  return(-1); 
	}if(j1 < 0 || j1 > nscans-1 || src[idmap[sid[j1]]].type < 0  || src[idmap[sid[j1]]].type > 1){
	  fprintf(stderr,"Error: second entry of your calib file wrong !\n");
	  return(-1);
	} if(k1 < 0 || k1 > nscans-1 || src[idmap[sid[j1]]].type < 0  || src[idmap[sid[j1]]].type > 1){
	  fprintf(stderr,"Error: third entry of your calib file wrong !\n");
	  return(-1);
	}
	lc[i1]=j1;      
	rc[i1]=k1;	   
	printf("\t scan:%d lc=%d rc=%d\n",i1,lc[i1],rc[i1]); 
      }
      fclose(inpfile); 
      return(0);
    break;
  case 2: 
    printf("\t FLAGGING from the input flag file : %s \n",flgfile);
    inpfile = fopen(flgfile,"r");
    if(check_badtimes_only){
      printf("\t CHECKING ONLY BADTIMES\n");
      while(fgets(c, 100, inpfile)!=NULL) {
	sprintf(string1,"%s",c);
	p = strchr(string1,'!');
	if (p == NULL) {
	  if (strstr(string1,"BADTIME") != NULL){
	    if(get_timerange(string1)){
	      fprintf(stderr,"can't parse %s\n",string1);
	      return 1;
	    }
	    fprintf(stdout,"\t BADTIME start= %f end =%f (days)\n",
		    badtimes[nbadtimes-1].start,badtimes[nbadtimes-1].end);
	  }
	}
      }
      fclose(inpfile);
      return 0;
    }
    while(fgets(c, 100, inpfile)!=NULL) {
      sprintf(string1,"%s",c);
      p = strchr(string1,'!');
      if (p == NULL) {
	if (strstr(string1,"ANT") != NULL){
	  string_decode(string1,&a1,&a2);
	  for(i=a1; i <=a2; i++){
	    if ( i >=0 && i < nants){
	      ba[i] = -1;
	      fprintf(stdout,"\t BAD ANT     = %d\n",i);
	    }// for if
	  }// for i
	}if (strstr(string1,"BASE") != NULL){
	  string_decode(string1,&a1,&a2);
	  for(i=a1; i <=a2; i++){
	    if ( i >=0 && i < nbaselines){
	      bb[i] = -1;
	      fprintf(stdout,"\t BAD BASELINE = %d\n",i);
	    }// for if
	  }// for i
	}if (strstr(string1,"CHAN") != NULL){
	  string_decode(string1,&a1,&a2);
	  for(i=a1; i <=a2; i++){
	    if ( i >=0 && i < nchans){
	      bc[i] = -1;
	      fprintf(stdout,"\t BAD CHANNEL = %d\n",i);
	    }// for if
	  }// for for
	}// for if
	if (strstr(string1,"BSCAN") != NULL){
	  string_decode(string1,&a1,&a2);
	  for(i=a1; i <=a2; i++){
	    if (i >=0 && i < nscans){
	      bs[i] = -1;
	      fprintf(stdout,"\t BAD SCAN = %d\n",i);
	    }// for if
	  }// for for
	}// for if
	if (strstr(string1,"GSCAN") != NULL){
	  string_decode(string1,&a1,&a2);
	  for(i=a1; i <=a2; i++){
	    if (i >=0 && i < nscans){
	      gs[i] = -1;
	      fprintf(stdout,"\t GOOD SCAN = %d\n",i);
	    }// for if
	  }// for for
	}// for i
      }// for if
    }// for while 
    fclose(inpfile);
    return(0);
    break;
  }
  return(0); 
}

int read_gains(){
  /* this program is optional */
  int j,j1,k1,l,l1,ntm,nchans2,nstokes2;
  FILE *outf;
  cmplx z; 
  
  outf   = fopen(gainfile,"rb");
  
  fread(&ntm,sizeof(int),1,outf);
  fread(&nchans2,sizeof(long),1,outf);
  fread(&nstokes2,sizeof(long),1,outf);
  
  if (nchans2==nchans && nstokes2==nstokes){
    if (ntm > nablocks)
      g1 = realloc(g1,ntm*nstokes*nchans*nants*sizeof(cmplx));
  }else{
    perror(" # channels/stokes in FITS file and gain file do not match !");
    return(-1);
  }
  for(j=0; j < ntm ; j++){
    fread(&abmt[j],sizeof(float),1,outf);
    for(j1 = 0; j1 < nchans; j1++){
      for(k1 = 0; k1 < nstokes; k1++){
	l1 = k1 + nstokes *(j1+nchans*j);
	for(l = 0; l < nants; l++){
	  fread(&z,sizeof(cmplx),1,outf);
	  g1[l+nants*l1]=Cmplx(z.re,-z.im);
	}// for l
      }// for k1
    }// for j1
  } // for j
  
  fclose(outf);
  free(g1); 
  fprintf(stdout,"\t gains have been computed and  written successfully !\n");
  
  return(0);
}// gain_write ends


int write_index(){
  FILE *outf,*outf1;
  int i1,i2,j,l,m,n;
  
  outf  = fopen("sindex.dat","w");
  outf1 = fopen("aindex.dat","w");
  for(i1=0; i1 < nscans; i1++){
    for(i2=0; i2 < nsb[i1]; i2++){
      j = sbid_in[i2+max_nsb*i1];
      fprintf(outf,"%d %d %d %d %d\n",i1,i2,j,sb_start[j],sb_end[j]); 
    }// for i2
    for(i2=0; i2 < nab[i1]; i2++){
      j = abid_in[i2+max_nab*i1];
      fprintf(outf1,"%d %d %d %d %d\n",i1,i2,j,ab_start[j],ab_end[j]); 
    }// for i2
  }// for i1
  fclose(outf); 
  fclose(outf1);
  
  outf=fopen("baselines.dat","w");
  for(n=0; n < nbaselines; n++){
    l = (id[n]+1)/nants;
    m = (id[n]+1) - nants *l;
    fprintf(outf,"%d %d %d\n",n,l,m); 
  }
  fclose(outf);  
  return(0);
}


int write_gains(){
  /* this program is optional */
  int i1,i2,j,j1,k1,l,l1,ntm;
  FILE *outf;
  cmplx z;
  
  outf   = fopen("gain.dat","wb");
  ntm    = 0;
  for(i1 = 0; i1 < nscans; i1++)
    ntm += nab[i1];

  fwrite(&ntm,sizeof(int),1,outf);
  fwrite(&nchans,sizeof(long),1,outf);
  fwrite(&nstokes,sizeof(long),1,outf);
  
  for(i1 = 0; i1 < nscans; i1++){
    for(i2 = 0; i2 < nab[i1]; i2++){
      j = abid_in[i2+max_nab*i1];
      fwrite(&abmt[j],sizeof(float),1,outf);
      for(j1 = 0; j1 < nchans; j1++){
        for(k1 = 0; k1 < nstokes; k1++){
          l1 = k1 + nstokes *(j1+nchans*j);
          for(l = 0; l < nants; l++){
            z = Cmplx(g1[l+nants*l1].re,-g1[l+nants*l1].im); // just for AIPS convention 
             fwrite(&z,sizeof(cmplx),1,outf);
          }// for l
        }// for k1
      } // for j1
    }// for i2
  } // for i1
  fclose(outf);
  free(g1);
  fprintf(stdout,"\t gains have been computed and  written successfully !\n");
  
  return(0);
}// gain_write ends
  
int write_flag_info(){
  int i1,i2,j,j1,k1,l,l1,m,n,ntm,n1; 
  FILE *outf1,*outf2,*outf3,*outf4; 
  float  ttm;
 
  fprintf(stdout,"\t Writing FLAG files\n"); 
  
  outf1 = fopen("vsr.dat","wb");
  outf2 = fopen("bad_ant.dat","w");
  outf3 = fopen("bad_bas.dat","w");
  outf4 = fopen("bad_chn.dat","w");
  
  ntm    = 0;
  for(i1 = 0; i1 < nscans; i1++){
    if (src[idmap[sid[i1]]].type != 2)
      ntm += nsb[i1];
  }// for i1

  fwrite(&ntm,sizeof(int),1,outf1);
  fwrite(&nchans,sizeof(long),1,outf1);
  fwrite(&nstokes_half,sizeof(long),1,outf1);
  fwrite(&thrs_vsr,sizeof(float),1,outf1);
  
  n1 = 0;
  for(i1 = 0; i1 < nscans; i1++){
    if (src[idmap[sid[i1]]].type != 2){
      for(i2 = 0; i2 < nsb[i1]; i2++){
	j = sbid_in[i2+max_nsb*i1];
	ttm = (time_data[sb_start[j]] + time_data[sb_end[j]])/2.0;
	fwrite(&ttm,sizeof(float),1,outf1);
	for(j1 = 0; j1 < nchans; j1++){
	  for(k1 = 0; k1 < nstokes_half; k1++){
	    l1 = k1 + nstokes *(j1+nchans*j);
	    for(n = 0; n < nbaselines;n++){
	      fwrite(&vsr[n+nbaselines*l1],sizeof(float),1,outf1);
	    }// for n
	  }// for k1
	} // for j1
	n1++;
      }// for i2
      for(l = 0;  l < nants; l++)
	fprintf(outf2,"%d %d %2.6f\n",i1,l+1,bad_ant[l+nants*i1]);
      for(n = 0;  n < nbaselines; n++){
	l = (id[n]+1)/nants;
	m = (id[n]+1) - nants *l;
	if ( bad_ant[l+nants*i1] > thrs_ant || bad_ant[m+nants*i1] > thrs_ant)
	  ttm  = 0.0;
	else
	  ttm  =  bad_base[n+nbaselines*i1]; 
	fprintf(outf3,"%d %d  %2.6f\n",i1,n,ttm);
      }// for n
      for(j1=0;  j1 < nchans; j1++)
	fprintf(outf4,"%d %d %2.6f\n",i1,j1,bad_chans[j1+nchans*i1]);
    }// for if
  } // for i1

  fclose(outf1);
  fclose(outf2);
  fclose(outf3);
  fclose(outf4);
  fprintf(stderr,"\t Flag files have been written successfully !\n");
  
  return(0);
}

int print_message(){
  
  fprintf(stdout," ---------------------------- FLAGCAL Version-1.0 (Nov  2010 )-------------------------------\n");
  fprintf(stdout,"\t - FLGACAL  flags and calibrates a multi-source, multi-channel data file\n");
  fprintf(stdout,"\t   which should be in the FITS format only. \n");
  fprintf(stdout,"\t - The effectiveness of the program depends on the set of parameters which \n");
  fprintf(stdout,"\t   are given in the parameter file (parameters.in).  \n");
  fprintf(stdout,"\t - Before running the full pipeline, run the pipeline with imode 1 \n");
  fprintf(stdout,"\t   and look at the information about the bad antennas, baselines & channels and \n");
  fprintf(stdout,"\t   set various parametesr accordingly. However the program may work fine with \n");
  fprintf(stdout,"\t   the default values. For detail check README and defs.h \n");
  fprintf(stdout,"\t - If you need more help check README, parameters.txt, defs.h & parameters.in \n");
  fprintf(stdout,"\t - You can always disable various modules by commentig their calls in the main program\n");
  fprintf(stdout,"\t   or by changing the valued of [a-h]modes in the parameter file. \n");
  fprintf(stdout,"\t -  If you jump from read_fits to write_fits (by commenting everything between in main program)\n");
  fprintf(stdout,"\t   you can use the flagcal for splitting files (single-channel/single -source etc.) also  \n");
  fprintf(stdout,"\t - You can feed a FITS file generated by flagcal to itself (convergence of this is not checked) \n");
  fprintf(stdout,"\t - Use the following options to run the pipeline: \n");
  fprintf(stdout,"-------------------------------------------------------------------------------------\n");
  fprintf(stdout,"\t  -i     : Input FITS file name (no default !) \n");
  fprintf(stdout,"\t  -o     : Output Fits File name (no default !)\n");
  fprintf(stdout,"\t  -p     : Input parameter file name  <default=parameters.in>\n");
  fprintf(stdout,"\t  -s     : Source file name <optional, no default>)\n");
  fprintf(stdout,"\t  -f     : FLAG  file name  <optiona,  no default>)\n");
  fprintf(stdout,"\t  -c     : CALIB  file name <optional, no default>)\n");
  fprintf(stdout,"\t  -g     : GAIN   file name <optional, no default>)\n");
  fprintf(stdout,"\t  -n     : # of openMP threads <default=1> \n");
  fprintf(stdout,"\t  -bchan : starting channel for output <default=0>\n");
  fprintf(stdout,"\t  -echan : ending channel for output <default=nchans>\n");
  fprintf(stdout,"\t  -osrc  : source [id] for output <default=0: all>\n");
  fprintf(stdout,"\t  -imode : mode for running  <default=0>\n");
  fprintf(stdout,"\t         :  0  => run full pipeline   \n");
  fprintf(stdout,"\t         :  1  => run only up to flagging \n");
  fprintf(stdout,"----------------------------------------------------------------------------------------\n");
  fprintf(stdout,"./flagcal -i ZCOSMOSA.FITS -o TEMP.FITS -p parameters.in -s sources.in -n 4\n");
  fprintf(stdout,"----------------------------------------------------------------------------------------\n");
  fprintf(stdout,"                                        Jayanti Prasad,   Sun Nov  7 17:19:33 IST 2010 \n");
  fprintf(stdout," ======================================================================================\n");
 
  return(0);   

}


int string_decode(char string1[], int *i1, int *i2){
  char  *p,*p1;
  p  = strchr(string1,'=');
  *i1 = atoi(p+1);
  if (strstr(string1,":") != NULL){
    p1 = strchr(string1,':');
    *i2 = atoi(p1+1);
  }else{
    *i2 = *i1;
  }
  return(0);
}
