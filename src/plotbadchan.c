#include <defs.h>
#include <flagcal.h>
#include<cpgplot.h>

int main(int argc, char*argv[]){
  int nall,*iscan;
  int  i,j,k,justdoit=1;
  int *x, *y;
  float *z,*xx,*yy; 
  FILE *inpf;
  char text[100];
  float x1,x2,y1,y2,tx1,tx2,ty1,ty2,tty;
  char ch='a';

  if (argc < 3){
    fprintf(stderr,"\t -------- PLOTBADCHAN-----------------------\n");
    fprintf(stderr,"\t - This program plots the performance of channels  using flag file (bad_chn.dat).  \n");
    fprintf(stderr,"\t   produced by flagcal.  \n");
    fprintf(stderr,"\t - The value is higher (btween 0 & 1) for bad channels  than for good channels .  \n");
    fprintf(stderr,"\t - This plot can be used to set the threshold for thrs_chan  (see paramrters.in file).  \n");
    fprintf(stderr,"\t  ./plotbadchan bad_chn.dat  <#chans>  \n\n");
    return(-1);
  }
  
  inpf = fopen(argv[1],"r");
  nchans = atoi(argv[2]);
  
  x   =  (int *)malloc(sizeof(int));
  y   =  (int *)malloc(sizeof(int));
  z   =  (float *)malloc(nchans*sizeof(float));
  
  xx  =  (float *)malloc(nchans*sizeof(float));
  yy  =  (float *)malloc(nchans*sizeof(float));
  
  iscan =  (int *)malloc(sizeof(int));
  
  i = 0;
  j = -1;
  nscans =0;
  while(!feof(inpf)){
    fscanf(inpf,"%d %d %f\n",&x[i],&y[i],&z[i]);
    if (x[i]!=j){
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
  x2 = (float) nchans;
  y1 = 0.0;
  y2 = 1.1;
  ch = 'a';
  while (justdoit == 1){
    cpgbeg(1,"/xs",1,1);
    cpgpap(10.0,0.75);
    cpgsch(0.75);
    cpgsvp(0.0,1.0,0.0,1.0);
    
    cpgenv(x1,x2,y1,y2,2,1);
    cpglab("CHAN","MEASURE OF BADNESS","CHANNEL PROFILE"); 
    
    cpgscr(0,0.0,0.0,0.0);
    cpgscr(1,1.0,1.0,1.0);
    
    for(i = 0; i < nscans; i++){
      if (iscan[i] < 10)
	sprintf(text,"0%d",iscan[i]); 
      else
      sprintf(text,"%d",iscan[i]);
      for(j = 0; j < nchans; j++){
	k = j + nchans * i;   
	xx[j] = (float ) y[k];
	yy[j] = z[k];  
      }
      tty = 1.0 - 0.35 *(float)i/nscans; 
      cpgsci(i+2);
      cpgtext(x1+0.95*(x2-x1),tty,text);
      cpgline(nchans,xx,yy);
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
