#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<cpgplot.h>

float max(float, float);
float min(float, float);

int main(int argc, char *argv[]){
  int i,j,k,l,m,nr,nc,nall,nscans,nbin,n1,*sci,*iscan;
  float *x,*y,x1=0.0,x2=1.0;
  float *aa,*bb,*cc,*dd; 
  FILE *inpf;
  char label[100];
  
  if (argc < 2){
   fprintf(stderr,"./plothist <inputfile> <#bins> <#rows> <# cols> \n");
   fprintf(stderr," ./plothist vsr.dat 100 3 3 \n");
   return(-1); 
  }
  
  inpf = fopen(argv[1],"r"); 
  nbin = atoi(argv[2]);
  nr = atoi(argv[3]);
  nc = atoi(argv[4]);

  x = (float *)malloc(sizeof(float));
  y = (float *)malloc(sizeof(float));
  sci = (int *)malloc(sizeof(int));
  iscan = (int *)malloc(sizeof(int));
  
  i =0;
  nscans=0;
  j = -1;
  iscan[0] = 0;
  while(!feof(inpf)){
    fscanf(inpf,"%d %f\n",&sci[i],&x[i]);
    if (sci[i] !=j){
      iscan[nscans]=sci[i]; 
      nscans++;
      iscan = realloc(iscan,(nscans+1)*sizeof(int));
    }   
    j=sci[i];
    i++;
    x = realloc(x,(i+1)*sizeof(float));
    sci = realloc(sci,(i+1)*sizeof(int));
  }
  nall = i;
  
  aa = (float *)malloc(nscans*sizeof(float));
  bb = (float *)malloc(nscans*sizeof(float));
  cc = (float *)malloc(nscans*sizeof(float));
  dd = (float *)malloc(nscans*sizeof(float));
  
  for(i=0; i < nr; i++){
    for(j=0; j < nc; j++){
      k = j + nc *i;
      aa[k] = (1.0/nc) * i;
      bb[k] = (1.0/nc) * (i+1);
      cc[k] = (1.0/nr) * j;
      dd[k] = (1.0/nr) * (j+1);
    }
  }
  
  cpgbeg(0,"?",nr,nc);
  cpgpap(12.0,0.8);
  
  for(i=0; i < nr; i++){
    for(j=0; j < nc; j++){
      
      k = j + nc *i;
      
      sprintf(label,"scan %d",iscan[k]);
      
      
      l = 0;
      for(m=0; m < nall; m++){
	if(sci[m] == iscan[k]){
	  y[l] = x[m];
	  l++;
	  y = realloc(y,(l+1)*sizeof(float));
	}// for if
      }// for 
    
      cpgsvp(aa[k],bb[k],cc[k],dd[k]);
      cpgenv(x1,x2,0.0,(float)l,2,1);
      cpglab("VSR", "N",label);
      cpghist(l,y,x1,x2,nbin,1);

    }// for j 
  }// for i
  cpgiden();

  cpgend();

  return(0); 

}

