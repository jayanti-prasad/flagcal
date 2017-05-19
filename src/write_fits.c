#include <defs.h>
#include <flagcal.h>
#include<fitsio.h>

/*-----------------------------<write_fits.c>-----------------------------------------------
   * This module is a part of flagcal v(0.0) global variables are defined in flagcal/defs.h * 
   * Report suggestions & bugs to prasad.jayanti@gmail.com 
   ------------------------------------------------------------------------------------------
  -  This module writes a single/multi channel single/multi source FITS file.
-  If no bchan, ehchan, isrc etc., are give the structure of the input and 
   output file is exactly the same.
-  There is a lot of reshuffling done with random group parameters which is 
   prone to break  down ! 
-  The logic followed is identical to that in read_fits.
                                                       --- Jayanti Prasad
                                                           Wed Nov 10 12:42:52 IST 2010
 ------------------------------------------------------------------------------*/   
int write_fits(char infits[], char outfits[]){
  char comment1[]="Created by FLGCAL ---- Jayanti Prasad";
  fitsfile *infptr,*outfptr;
  char  *strval[max_len],key_sutabl[FLEN_KEYWORD]="EXTNAME",comment[FLEN_COMMENT];
  int i1,j1,k1,l,l1,l2,m,isrc,osrc[nsrc+1],status=0,hdutype,morekeys=0; 
  long nchans1,pcount1,gcount1,group,group1,el1 = 1,nel; 
  float randpar[pcount],randpar1[pcount],*data,U[3];
  
  fprintf(stdout,"\t  WRITTING output FITS file\n");
  
  for(i1=0; i1 <= nsrc; i1++)
    osrc[i1] = 0;
  
  if(fits_open_file(&infptr,infits,READONLY,&status))
    printerror(status);    
  
  if (fits_create_file(&outfptr,outfits, &status)) 
    printerror(status);           
  
  if (fits_copy_header(infptr,outfptr, &status)) 
    printerror(status);           
  
  if (fits_write_history(outfptr,comment1,&status))
    printerror(status);
  
  if (fits_write_date(outfptr,&status))
    printerror(status);
  
  if(echan ==0)
    echan = nchans; 
  
  if(echan != bchan )
    nchans1 = echan-bchan;
  else
    nchans1 = 1;
  
  if(src_out != 0 ){
    osrc[src_out] = 1;
    gcount1  = ggcount[src_out];
    if (fits_update_key(outfptr,TSTRING,"OBJECT",&src[src_out].name,"",&status))
      printerror(status);
    if (fits_update_key(outfptr,TFLOAT,"CRVAL6",&src[src_out].ra,"",&status))
      printerror(status);
    if (fits_update_key(outfptr,TFLOAT,"CRVAL7",&src[src_out].dec,"",&status))
      printerror(status);
    if (fits_delete_key(outfptr,"PTYPE7",&status))
      printerror( status );
    if (fits_delete_key(outfptr,"PSCAL7",&status))
      printerror( status );
    if (fits_delete_key(outfptr,"PZERO7",&status))
      printerror( status );
    if(nchans1 == 1){
      pcount1  = pcount-2; 
      fprintf(stdout,"\t single source - single channel\n");
      if (fits_delete_key(outfptr,"PTYPE8",&status))
	printerror( status );
      if (fits_delete_key(outfptr,"PSCAL8",&status))
	printerror( status );
      if (fits_delete_key(outfptr,"PZERO8",&status))
	printerror( status );
    }else{
      pcount1  = pcount-1; 
      fprintf(stdout,"\t single source - multi  channel\n");
      if (fits_modify_name(outfptr,"PTYPE8","PTYPE7",&status))
	printerror( status );
      if (fits_modify_name(outfptr,"PSCAL8","PSCAL7",&status))
	printerror( status );
      if (fits_modify_name(outfptr,"PZERO8","PZERO7",&status))
	printerror( status );
    }// for if
  }else{
    for(i1=1; i1 <= nsrc; i1++)
      osrc[i1] = 1;
    gcount1 = gcount;
    if (nchans1==1){
      pcount1=pcount-1;
      fprintf(stdout,"\t multi source - single channel\n");
      if (fits_delete_key(outfptr,"PTYPE8",&status))
	printerror( status );
      if (fits_delete_key(outfptr,"PSCAL8",&status))
	printerror( status );
      if (fits_delete_key(outfptr,"PZERO8",&status))
	printerror( status );
    }else{
      pcount1=pcount;
      fprintf(stdout,"\t multi source - multi  channel\n");
    }// for else
  }// for else 
   
  fprintf(stdout,"\t output FITS file parameters \n");
  fprintf(stdout,"\t # output chans     = %ld\n",nchans1);
  fprintf(stdout,"\t # output gcount    = %ld\n",gcount1);
  fprintf(stdout,"\t # output pcount    = %ld\n",pcount1);
  
  for(i1=1; i1 <= nsrc; i1++)
    fprintf(stdout,"\t source = %d out = %d\n",i1,osrc[i1]);
  
  if (fits_update_key(outfptr,TLONG,"GCOUNT",&gcount1,"GCOUNT = ",&status))
    printerror(status);
  if (fits_update_key(outfptr,TLONG,"PCOUNT",&pcount1,"PCOUNT = ",&status))
    printerror(status);
  if (fits_update_key(outfptr,TLONG,"NAXIS4",&nchans1,"NAXIS4 = ",&status))
    printerror(status);
  
   
  data = (float *)malloc(ncmplx*nstokes*nchans1*sizeof(float));
  nel      = ncmplx * nstokes * nchans1; 

  for(l2=0,group1=1,group = 1; group <= gcount; group++){  
    
    if(fits_read_grppar_flt(infptr,group,el1,pcount,randpar,&status))
      printerror( status );
    
    group_data(randpar,U,&l,&m,&isrc);
    
    if(osrc[idmap[isrc]] == 1 && bgrp[group] < nbaselines ){
      for(k1 = 0; k1 < pcount1; k1++)
	randpar1[k1] = randpar[k1]; 
      
      if(src_out !=0 && nchans1 > 1)
	randpar1[pcount1-1] = randpar[pcount-1]; 
      
      if(fits_write_grppar_flt(outfptr,group1,el1,pcount1,randpar1,&status))
      	printerror(status );
      
      for(j1 = bchan; j1 < echan; j1++){
	for(k1 = 0; k1 < nstokes; k1++){
	  l1  =  k1 + nstokes*(j1+nchans*tgrp[group]);
	  data[0+ncmplx*(k1+nstokes*(j1-bchan))] = Visib[bgrp[group]+nbaselines*l1].re;
	  data[1+ncmplx*(k1+nstokes*(j1-bchan))] = Visib[bgrp[group]+nbaselines*l1].im;
	  data[2+ncmplx*(k1+nstokes*(j1-bchan))] = Visib[bgrp[group]+nbaselines*l1].w;
	}// for k1
      }// for j1

      if(fits_write_img_flt(outfptr,group1,el1,nel,data,&status))
    	printerror(status);
      group1++;
      l2++;
    } // for if source 
  } // for group 
  
  fprintf(stderr,"\t *** wrote visibility data, writing extra HDUs ***\n");
  
  for(i1=2; i1 <=nhdu; i1++){
    if (fits_movabs_hdu(infptr,i1,&hdutype, &status) )
      printerror( status );
    if(fits_read_key_longstr(infptr,key_sutabl,strval,comment,&status))
      printerror( status );
    if (strstr(*strval,"SU") != NULL){
      if (src_out == 0)
        if(fits_copy_hdu(infptr, outfptr, morekeys, &status) )
          printerror( status );
    }else{
      if(fits_copy_hdu(infptr, outfptr, morekeys, &status) )
        printerror( status );
    }// for if
  }// for i1
  
  fprintf(stderr,"\t *** wrote  extra HDUs ***\n");
  
  if (fits_close_file(outfptr,&status))
    printerror( status );
  
  fprintf(stdout,"\t Every thing written \n");
  
  return(0);
}// write_output ends 

