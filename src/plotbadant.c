#include <defs.h>
#include<cpgplot.h>
int main(int argc, char*argv[]){
  int i,j,k,nall,nscans,*iscan,*x,*y;
  float *z,*xx,*yy,x1[2],y1[2]; 
  FILE *inpf;
  char text[100];

  if (argc < 2){
     fprintf(stderr,"\t -------- PLOTBADANT-----------------------\n");
     fprintf(stderr,"\t - This program plots the performance of antennas using flag file (bad_ant.dat).  \n");
     fprintf(stderr,"\t   produced by flagcal.  \n");
     fprintf(stderr,"\t - The value is higher (btween 0 & 1) for bad antennas then for good antennas.  \n");
     fprintf(stderr,"\t - This plot can be used to set the threshold for thrs_ante (see paramrters.in file).  \n");
     fprintf(stderr,"\t  ./plotbadant bad_ant.dat  \n\n");
     return(-1); 
  }
  
  x  =  (int *)malloc(sizeof(int));
  y  =  (int *)malloc(sizeof(int));
  z  =  (float *)malloc(nants*sizeof(float));
  xx =  (float *)malloc(nants*sizeof(float));
  yy =  (float *)malloc(nants*sizeof(float));
  iscan =  (int *)malloc(sizeof(int));
  
  inpf = fopen(argv[1],"r"); 

  i = 0;
  j = -1;
  nscans =0;
  while(!feof(inpf)){
    fscanf(inpf,"%d %d %f\n",&x[i],&y[i],&z[i]);
    if (x[i] !=j){
      iscan[nscans] = x[i];
      nscans++;
      iscan = realloc(iscan,(nscans+1)*sizeof(int));
      j = x[i];
    }
    i++;
    x = realloc(x,(i+1)*sizeof(int));
    y = realloc(y,(i+1)*sizeof(int));
    z = realloc(z,(i+1)*sizeof(float));
  }//for while 
  nall = i;
  
  cpgbeg(1,"?",1,1);
  cpgpap(10.0,0.75);
  cpgsvp(0.0,1.0,0.0,1.0);
  cpgsch(0.75);
  cpgenv(1.0,(float)nants,0.0,1.1,2,1);
  cpglab("ANT","MEASURE OF BADNESS","ANTENNA PROFILE"); 
  //cpgscr(0,0.0,0.0,0.0);
 // cpgscr(1,1.0,1.0,1.0);
  cpgslw(2);
  cpgtext(0.923*nants,1.04,"scan");
 // for(i = 0; i < 4; i++){
   for(i = 0; i < nscans; i++){
    if (iscan[i] < 10)
      sprintf(text,"0%d",iscan[i]); 
    else
      sprintf(text,"%d",iscan[i]);
    for(j=0; j < nants; j++){
      k = j + nants  * i;   
      xx[j] =  (float ) y[k];
      yy[j] =  z[k];  
    }
    x1[0] = 0.88 * nants;
    x1[1] = 0.92 * nants;
    y1[0] = 1.0 - 0.35 *(float)i/nscans; 
    y1[1] = y1[0];  
    cpgsci(i+2);
  //  cpgsls(i+1);
    cpgline(2,x1,y1); 
    cpgtext(0.94*nants,y1[1],text);
    cpgline(nants,xx,yy);
  }// for i
  cpgend();
  return(0);
}

