#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<cpgplot.h>

#include "defs.h"

float max(float, float);
float min(float, float);

int main(int argc, char*argv[]){
  int nall,nscans,*iscan;
  int  i,j,k,justdoit=1;
  int *x, *y;
  float *z,*xx,*yy; 
  FILE *inpf;
  char text[100];
  float x1,x2,y1,y2,tx1,tx2,ty1,ty2,tty;
  char ch='a';

  x   =  (int *)malloc(sizeof(int));
  y   =  (int *)malloc(sizeof(int));
  z   =  (float *)malloc(nbaselines*sizeof(float));
  xx  =  (float *)malloc(nbaselines*sizeof(float));
  yy  =  (float *)malloc(nbaselines*sizeof(float));
  
  iscan =  (int *)malloc(sizeof(int));
  
  if (argc < 2){
    fprintf(stderr,"\t ./plotbadbase  <input file>  \n");
    fprintf(stderr,"\t ./plotbadbase  bad_base  \n");
    return(-1);
  }
  
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
      j=x[i];
    }
    i++;
    x = realloc(x,(i+1)*sizeof(int));
    y = realloc(y,(i+1)*sizeof(int));
    z = realloc(z,(i+1)*sizeof(float));
  }
  nall = i;
  
 loop1:
  x1 = 0.0;
  x2 = (float) nbaselines;
  y1 = 0.0;
  y2 = 1.1;
  ch = 'a';
  while (justdoit == 1){
    cpgbeg(1,"/xs",1,1);
    cpgpap(10.0,0.75);
    cpgsch(0.75);
    cpgsvp(0.0,1.0,0.0,1.0);
    
    cpgenv(x1,x2,y1,y2,2,1);
    cpglab("BASELINE","MEASURE OF BADNESS","BASELINE PROFILE"); 
    
    cpgscr(0,0.0,0.0,0.0);
    cpgscr(1,1.0,1.0,1.0);
    
    for(i = 0; i < nscans; i++){
      if (iscan[i] < 10)
	sprintf(text,"0%d",iscan[i]); 
      else
      sprintf(text,"%d",iscan[i]);
      for(j = 0; j < nbaselines; j++){
	k = j + nbaselines * i;   
	xx[j] = (float ) y[k];
	yy[j] = z[k];  
      }
      tty = 1.0 - 0.35 *(float)i/nscans; 
      cpgsci(i+2);
      cpgtext(x1+0.95*(x2-x1),tty,text);
      cpgline(nbaselines,xx,yy);
    }
  
    cpgsci(2);
    cpgcurs(&tx1,&ty1,&ch);
    if (ch=='a' || ch=='A'){
      cpgband(2,1,tx1,ty1,&tx2,&ty2,&ch);
      x1 = min(tx1,tx2);
      x2 = max(tx1,tx2);
    }
    switch (ch){
    case 'b':
      justdoit = 1;
      goto loop1;
      break;
    case 'd':
      cpgend();
      return(0);
      break;
    }//for switch
  }
  cpgend();
  return(0);
}
