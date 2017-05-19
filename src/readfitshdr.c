#include<stdio.h>
#include<string.h>
#include<fitsio.h>

int main(int argc, char *argv[]){

  fitsfile *fptr;
  char card[FLEN_CARD];
  int status = 0,  nkeys, ii,l, numhdu,hdunum,hdutype;

  if(argc < 2){
   fprintf(stderr,"./readfitshdr  <FITS FILE>\n");
   return(-1);
  }


  fits_open_file(&fptr,argv[1],READONLY,&status);
  
  fits_get_hdrspace(fptr,&nkeys, NULL, &status);

  fits_get_num_hdus(fptr,&numhdu,&status);

  printf("#  HDUs=%d\n",numhdu);

  for(l=1; l <= numhdu; l++){
    fits_movabs_hdu(fptr,l,&hdutype,&status);
    fits_get_hdu_num(fptr,&hdunum);
    
    
    printf("HDU = %d, HDU type=%d\n",hdunum,hdutype); 
    
    fits_get_hdrspace(fptr,&nkeys, NULL, &status);
    
    for(ii=1; ii <=nkeys; ii++){
      fits_read_record(fptr,ii,card,&status);
      printf("%s\n",card);
    } 
    printf("END\n\n");
    
  }
  
  fits_close_file(fptr,&status);
  
  if(status)
    fits_report_error(stderr,status);
  return(status);
}
