#include<fitsio.h>
#include <defs.h>
#include <flagcal.h> 

int main(int argc, char *argv[]){
  fitsfile *fptr;
  char key_simple[FLEN_KEYWORD]="SIMPLE";
  char key_naxis2[FLEN_KEYWORD]="NAXIS2";
  char key_naxis3[FLEN_KEYWORD]="NAXIS3";
  char key_naxis4[FLEN_KEYWORD]="NAXIS4";
  char key_gcount[FLEN_KEYWORD]="GCOUNT";  // # of visibilities
  char key_pcount[FLEN_KEYWORD]="PCOUNT";  // # of random parameters
  char key_crval4[FLEN_KEYWORD]="CRVAL4";  // frequency 
  char keyvalue[FLEN_VALUE],comment[FLEN_COMMENT],this_key[FLEN_KEYWORD];;
  int i,l,m,j1,k1,isrc,status=0,anynul=0,tday,tday1,iday;
  long group,el1=1,nel;
  float *randpar,*data,tdata,U[3],w,nulval = 0.0;
  cmplx z; 
  
  if(argc < 4){
    fprintf(stdout,"\n\t - This program reads single channel single stokes data from\n");
    fprintf(stdout,"\t   a FITS file which can be directed & then plotted using the ivpot.   \n");
    fprintf(stdout,"\t - Use the following options: \n\n");
    fprintf(stderr,"\t ./readfits <inputfile> <chanel> <stokes> \n\n");
    fprintf(stdout,"\t  - Last modified on Thu Jul  8 12:03:54 IST 2010 --Jayanti Prasad \n");
    return(-1);
  }// for if
  
  if (fits_open_file(&fptr,argv[1],READONLY,&status))
    printerror(status);
  
  fits_read_keyword(fptr,key_simple,keyvalue,comment,&status);
  if(*keyvalue!='T'){fprintf(stderr,"Input File is NOT a SIMPLE FITS file\n");}
  
  fits_read_key_lng(fptr,key_naxis2,&ncmplx,comment,&status);
  fits_read_key_lng(fptr,key_naxis3,&nstokes,comment,&status);
  fits_read_key_lng(fptr,key_naxis4,&nchans,comment,&status);
  fits_read_key_lng(fptr,key_gcount,&gcount,comment,&status);
  fits_read_key_lng(fptr,key_pcount,&pcount,comment,&status);
  fits_read_key_flt(fptr,key_crval4,&crval4,comment,&status);

  pscal=(float *)malloc(pcount*sizeof(float));
  pzero=(float *)malloc(pcount*sizeof(float));
 
  fprintf(stderr,"\t ncmplx      = %ld\n",ncmplx);
  fprintf(stderr,"\t nstokes     = %ld\n",nstokes);
  fprintf(stderr,"\t nchans      = %ld\n",nchans);
  fprintf(stderr,"\t gcount      = %ld\n",gcount);
 
  nel     = ncmplx * nstokes * nchans; 
  data    = (float *)malloc(sizeof(float)*nel);
  randpar = (float *)malloc(pcount*sizeof(float)); 
  for(i=1; i <=pcount; i++){
    sprintf(this_key,"PSCAL%d",i);
    fits_read_key_flt(fptr,this_key,&pscal[i],comment,&status);
  //  fprintf(stderr,"\t pscal[%d]  :%12.2e\n",i,pscal[i]);
    sprintf(this_key,"PZERO%d",i);
    fits_read_key_flt(fptr,this_key,&pzero[i],comment,&status);
  //  fprintf(stderr,"\t pzero[%d]  :%12.2e\n",i,pzero[i]);
   }// for if
 
  iday  = -1;
  tday  = -1;
  tday1 = -1;
  tdata = 0.0;

  for(group = 1; group <= gcount; group++){
    if (fits_read_grppar_flt(fptr,group,el1,pcount,randpar,&status))
      printerror( status );
    if(group_data(randpar,U,&l,&m,&isrc))
      return(-1);

    // in some cases it was found that people are packing all the info regarding
    // date and time in the fourth parameter and the fifth one is zero 
    // in that case random group parameter has a fractional part also which 
    // represent  time 

    tday = (int) randpar[4]; 

    if (tday > tday1)
      iday++; 

    tdata = (float) iday + gettime(randpar[4],randpar[5]) ; 

//    printf("tdata=%.6f\n",tdata); 


    if (pcount > 7) 
      isrc = (int) randpar[6];
    else
      isrc = 0;

    if (l < nants && m < nants){ 
      if(fits_read_img_flt(fptr,group,el1,nel,nulval,data,&anynul,&status))
	printerror(status );
      for(j1= 0; j1 < nchans; j1++){
	for(k1 = 0; k1< nstokes; k1++){
    if(j1 == atoi(argv[2])  && k1== atoi(argv[3]) ){
//	  if(l == atoi(argv[2])  && m== atoi(argv[3]) && k1==0){
	    z.re = data[0+ncmplx*(k1+nstokes*j1)];
	    z.im = data[1+ncmplx*(k1+nstokes*j1)];
	    w    = data[2+ncmplx*(k1+nstokes*j1)];
	    printf("%2.6f %d %d %d %6.2f %6.2f %6.2f %6.3f %6.0f %2.4f\n",
		   tdata,isrc,l,m,U[0],U[1],U[3],cmod(z),cphase(z),w);

//            printf("%2.6f %2.6f %2.6f %2.6f\n",tdata,(float)j1,cmod(z),w);

	  }// for if j1 & k1 
	} // for k1
      } // for j1
    }// for if 
    tday1=tday; 
  }// for group  
  if(fits_close_file(fptr,&status))
    printerror(status );
  return(0); 
}


