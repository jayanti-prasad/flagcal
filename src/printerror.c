#include<stdio.h>
#include<fitsio.h>
#include<stdlib.h>

void printerror( int status){
  if (status){
    fits_report_error(stderr, status);
    exit(status ); 
  }
  return;
}
