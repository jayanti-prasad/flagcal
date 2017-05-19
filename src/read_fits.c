#include <defs.h>
#include <flagcal.h>
#include<fitsio.h>
#include<string.h>

/*-------------------------------<read_fits.c>-------------------------------------------------
 * This module is a part of flagcal v(0.0) global variables are defined in flagcal/defs.h * 
 *    * Report suggestions & bugs to prasad.jayanti@gmail.com 
 *  -----------------------------------------------------------------------------------
 *       
  This module is written to read a Multi-source-Multi-channel radio-interferrometric 
observation FITS file. At present the module is very conservative and make the 
following assumtions: 
- Header in the FITS file contains keywords GCOUNT, PCOUNT
- The file contains > four HDUs,  one primary HDU (containing visibilities) and 
  other containg AN, FQ and SU tables.
- Visibility data MUST be in time order ! 
- In the SU table the indexing of the sources must start from 1
- I do not check the random group paramters and assume the following:
  rand(0,1,2) = > U,V,W
  rand(3)     = > baseline 
  rand(4,5)   = > time 
  rand(6)     = > source id  
 . integer part of rand(4) is assumed to contain day information and fractional 
   part hh/mm/ss 
 . rand(5) is assumed to contain hh/mm/ss information  
  
- The data is indexed into scans and the logic is as follows:
  when time and source id (rand(6)) changes it is assumed that a new scan 
  has started.
- Every scan is assigned the following variables:
  starting-ending time, source id, starting and ending group number. 

WARNING ! It is dangerous to open a FITS file more than one times 
          in a program I learned it after spending months !
  
                                                             -- Jayanti Prasad
                                                             Wed Nov 10 12:42:52 IST 2010
--------------------------------------------------------------------------------*/
int read_fits(char filename[]){
  fitsfile *fptr;
  char key_sutabl[FLEN_KEYWORD]="EXTNAME",key_simple[FLEN_KEYWORD]="SIMPLE";
  char key_bitpix[FLEN_KEYWORD]="BITPIX",key_naxis[FLEN_KEYWORD] ="NAXIS";
  char key_naxis2[FLEN_KEYWORD]="NAXIS2",key_naxis3[FLEN_KEYWORD]="NAXIS3";
  char key_naxis4[FLEN_KEYWORD]="NAXIS4",key_gcount[FLEN_KEYWORD]="GCOUNT"; 
  char key_pcount[FLEN_KEYWORD]="PCOUNT",key_crval4[FLEN_KEYWORD]="CRVAL4";
  char key_ra[FLEN_KEYWORD]="CRVAL6"; 
  char key_dec[FLEN_KEYWORD]="CRVAL7",keyvalue[FLEN_VALUE];
  char comment[FLEN_COMMENT], *strval[max_len],this_key[FLEN_KEYWORD];;
  char *val, value[2000], nullstr[]="*";
  float *randpar,*data,U[3];
  float t1,t2,nulval=0.0;
  long gg[MAX_SCANS],el1=1,group,l1,nel, jj, nrows;
  int i1,i,j,l,m,n,j1,k1,iday,iday1,nday,isrc,isrc1;
  int nsu,ncols,ii,dispwidth[2000];
  int firstcol, lastcol = 0, linewidth,status  = 0,anynul = 0,hdutype;
  int bad_time_flag;
  if (fits_open_file(&fptr,filename,READONLY,&status))
    printerror( status );
  
  fits_read_keyword(fptr,key_simple,keyvalue,comment,&status);
  if(*keyvalue!='T'){ fprintf(stdout,"Input File is NOT a SIMPLE FITS file\n");}
  
  fits_read_key_lng(fptr,key_bitpix,&nbit,comment,&status); 
  fits_read_key_lng(fptr,key_naxis, &naxis,comment,&status);
  fits_read_key_lng(fptr,key_naxis2,&ncmplx,comment,&status);
  fits_read_key_lng(fptr,key_naxis3,&nstokes,comment,&status);
  fits_read_key_lng(fptr,key_naxis4,&nchans,comment,&status);
  fits_read_key_lng(fptr,key_gcount,&gcount,comment,&status);
  fits_read_key_lng(fptr,key_pcount,&pcount,comment,&status);
  fits_read_key_flt(fptr,key_crval4,&crval4,comment,&status);
  fits_read_key(fptr,TFLOAT,key_ra,&ra, comment,&status);
  fits_read_key(fptr,TFLOAT,key_dec,&dec, comment,&status); 
  fits_get_num_hdus(fptr,&nhdu,&status);
  freq  = crval4 /1.0e6; // in MHz
  
  fprintf(stdout,"\t ------------------- KEYWORDS--------------------------\n");
  fprintf(stdout,"\t nbit      : %12ld\n",nbit);
  fprintf(stdout,"\t nHDUs     : %12d\n",nhdu);
  fprintf(stdout,"\t naxis     : %12ld\n",naxis);
  fprintf(stdout,"\t ncmplx    : %12ld\n",ncmplx);
  fprintf(stdout,"\t nstokes   : %12ld\n",nstokes);
  fprintf(stdout,"\t nchans    : %12ld\n",nchans);
  fprintf(stdout,"\t gcount    : %12ld\n",gcount);
  fprintf(stdout,"\t pcount    : %12ld\n",pcount);
  fprintf(stdout,"\t freq[Mhz] : %12.2f\n",freq);
  
  tgrp=(int *)malloc((gcount+1)*sizeof(int));
  bgrp=(int *)malloc((gcount+1)*sizeof(int));
  randpar=(float *)malloc(pcount*sizeof(float));
  pscal=(float *)malloc(pcount*sizeof(float));
  pzero=(float *)malloc(pcount*sizeof(float));
  data=(float *)malloc(ncmplx*nstokes*nchans*sizeof(float));
  
  for(i=1; i <=pcount; i++){
    sprintf(this_key,"PSCAL%d",i);
    fits_read_key_flt(fptr,this_key,&pscal[i],comment,&status);
  //  fprintf(stdout,"\t pscal[%d]  :%12.2e\n",i,pscal[i]);
    sprintf(this_key,"PZERO%d",i);
    fits_read_key_flt(fptr,this_key,&pzero[i],comment,&status);
   // fprintf(stdout,"\t pzero[%d]  :%12.2e\n",i,pzero[i]);
  }// for if
  
  /* This is for the maping for baselines */
  
  for(j=0,l=0; l <nants ; l++){
    for(m=l+1;  m <nants; m++){
      id[j] = m + nants *l-1;
      dc[id[j]] = j;
      j++;
    }// for m
  }// for l

  /*If there are more than 2 stokes then RL & LR ignored for flagging */
  
  if (nstokes > 2)
    nstokes_half = nstokes/2;
  else
    nstokes_half = nstokes;
  
  /* Initilize the variables */
  
  t1=0.0;
  t2=0.0;
  ntimes=-1;
  nscans=-1;
  isrc=-1;
  isrc1=-1;
  nday = 0;
  
  /* This is the first loop which will find out how many scans, time samples etc. are there */
  
  for(group = 1; group <= gcount; group++){
    if (fits_read_grppar_flt(fptr,group,el1,pcount,randpar,&status))
      printerror(status);
    
    // decode u,v,w,l.m,source 
    group_data(randpar,U,&l,&m,&isrc);
    
    if (l<nants && m < nants)    // for software backend data l & m be greater than 30
      n = dc[m+nants*l-1];       // we ignore all the data for l > 30 and m > 30  
    
    // in general the fourth random group parameter is day and does not have 
    // any fractional part 
  
    iday = (int) randpar[4];
    
    if(iday > iday1)
      nday++; 
    
    // in some cases it was found that people are packing all the info regarding
    // date and time in the fourth parameter and the fifth one is zero 
    // in that case random group parameter has a fractional part also which 
    // represent  time 

     t2 = (float) nday + gettime(randpar[4],randpar[5]) ;

    // in any case t2 must be monotonically increasing. if that is not the case 
    // the whole thing will break dowon. 
    
    if(t2 > t1){
      ntimes++;
      time_data[ntimes]=t2;
      if (isrc != isrc1){
	nscans++;
	tid[nscans]=ntimes;
	gid[nscans]=group;
	sid[nscans] =isrc;
      }
    }// for if
    tgrp[group] = ntimes; 
    bgrp[group]=n;
    iday1 = iday;  
    isrc1 = isrc; 
    t1    = t2; 
  } // for group  
  ntimes++;
  nscans++;
  tid[nscans] = ntimes;
  gid[nscans] = gcount;
  time_data[ntimes] = time_data[ntimes-1];

  fprintf(stdout,"\t -----------------------------------------------------\n");  
  fprintf(stdout,"\t nscans      : %d \n",nscans);
  fprintf(stdout,"\t ntimes      : %d \n",ntimes);
  fprintf(stdout,"\t tgrp[gcount : %d\n",tgrp[gcount]);
  fprintf(stdout,"\t bgrp[gcount : %d\n",bgrp[gcount]);

  fprintf(stdout,"\t ------------------------------Index-----------------------\n");
  for(i=0; i < nscans; i++)
    fprintf(stdout,"\t Scan: %6d SrcId:%6ld t_begin :%6ld t_end :%6ld\n",i,sid[i],tid[i],tid[i+1]);
  fprintf(stdout,"\t ----------------reading visibilities-------------------------------------\n");
  
  Visib=(cmplx3 *)malloc(ntimes*nbaselines*nstokes*nchans*sizeof(cmplx3));
  for(i1=0; i1 < ntimes; i1++){
    for(j1=0; j1 < nchans; j1++){
      for(k1=0; k1 < nstokes; k1++){
	l1=k1+nstokes*(j1+nchans*i1);
	for(n=0; n < nbaselines; n++)
	  Visib[n+nbaselines*l1]   = Cmplx3(0.0,0.0,1.0);
      }// for k1
    }// for j1
  }// for i1
  
  /*  This is the second and the main  loop which which reads Visibility data 
      Timeranges indicated as bad are flagged. The day count is forced to
      zero for the first record
  */
    
  nel     = ncmplx * nstokes * nchans;
  for(nday=0,group = 1; group <= gcount; group++){
    double t1,t2;
    long   d1,d2;
    if (fits_read_grppar_flt(fptr,group,el1,pcount,randpar,&status))
      printerror(status);
    t2=randpar[4];
    d2=(long)floor(randpar[4]);
    t2=t2-d2;
    t2= t2+randpar[5];
    if(group==1){t1=t2;d1=d2;}
    t2=t2-t1+(d2-d1);
    bad_time_flag=0;
    for(i=0;i<nbadtimes;i++){
      if(t2 >= badtimes[i].start && t2 <= badtimes[i].end){
	bad_time_flag=1;
	break;
      }
    }
    if(fits_read_img_flt(fptr,group,el1,nel,nulval,data,&anynul,&status))
     printerror(status );
   
    
    if (bgrp[group] > 0 && bgrp[group]  <nbaselines){
      for(j1=0; j1 < nchans; j1++){
	for(k1 = 0; k1< nstokes; k1++){
	  l1 = k1+nstokes*(j1+nchans*tgrp[group]);   
	  Visib[bgrp[group]+nbaselines*l1]  = 
            Cmplx3(data[0+ncmplx*(k1+nstokes*j1)],data[1+ncmplx*(k1+nstokes*j1)],data[2+ncmplx*(k1+nstokes*j1)]);
	  if(bad_time_flag)
	    Visib[bgrp[group]+nbaselines*l1].w=-1.0;
	} // for k1
      } // for j1
    } // for if
  } // for group  
  
  fprintf(stdout,"\t read visibility data \n\t reading SU table \n");
  
  /* Now find out what is the order of the SU table */
  
  for(nsu=0,i=2; i <= nhdu; i++){
    if (fits_movabs_hdu(fptr,i,&hdutype, &status) )
      printerror(status );
    fits_read_key_longstr(fptr,key_sutabl,strval,comment,&status);
    if (strstr(*strval,"SU") != NULL)
      nsu = i;
  }// for i
  
  fprintf(stdout,"\t SU table order : %6d\n",nsu);
  
  if (fits_movabs_hdu(fptr,nsu,&hdutype, &status) )
    printerror(status );
  
  val = value;
  fits_get_num_rows(fptr, &nrows, &status);
  fits_get_num_cols(fptr, &ncols, &status);
  nsrc = nrows;
  
  while(lastcol < ncols) {
    linewidth = 0;
    firstcol = lastcol+1;
    for (lastcol = firstcol; lastcol <= ncols; lastcol++) {
      fits_get_col_display_width(fptr, lastcol, &dispwidth[lastcol], &status);
      linewidth += dispwidth[lastcol] + 1;
      if (linewidth > 80)break;  
    }// for 
    
    if (lastcol > firstcol)lastcol--;  /* the last col didn't fit */
    
    /* read each column, row by row (there are faster ways to do this) */
    for (jj = 1; jj <= nrows && !status; jj++) {
      for (ii = firstcol; ii <= lastcol; ii++){
	fits_read_col_str(fptr,ii,jj, 1, 1, nullstr, &val, &anynul, &status);
	if (ii==1)  src[jj].id=atoi(value);
	if (ii==2)strcpy(src[jj].name,value);
	if (ii==11) src[jj].ra = atof(value);
	if (ii==12) src[jj].dec= atof(value);
      }// for ii
    }// for jj
  }// for while 
  
  if (fits_close_file(fptr,&status))
    printerror( status );

  for(i=0; i < nscans; i++)
    gg[i] = gid[i+1] - gid[i];
  
  for(i=0; i <= nsrc; i++)
    ggcount[i] = 0;
  
  for(i=0; i < nscans; i++)
    ggcount[sid[i]]+=gg[i]; 
  
  for(i=1 ;i <=nsrc; i++)
    fprintf(stdout,"\t src=%18s gcount=%12ld\n",src[i].name,ggcount[i]); 

  fprintf(stdout,"\t FITS file has been read successfully !\n"); 
  
  return(0);
}

