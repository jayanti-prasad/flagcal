#include <defs.h>
#include <flagcal.h>
#include <cpgplot.h>

#define  nc  256  
#define  eps 0.0001;

int   read_fits(char []);
float min(float,float);
float max(float,float);
int   ftime(float, float[],float*);

int bresolve(int , char *[],int *, float []);
int pmod = 0,ihelp=0; 

int read_input_opt(int , char *[],int *);

double prd(double); 

int flip(int );


int main(int argc, char *argv[]) {
  int ntimes1,istck,nbl,i,j,l,m,n,ii,jj,i1,i2,j1,l1,fmod=0; 
  int nxsub=0,nysub=0;;
  float *data,tr[6],*xt,base[nbaselines],tmm[4],scl;
  float x1,x2,y1,y2,tx1,tx2,ty1,ty2,d1,d2,r;
  float xtick = 0.0,ytick=0.0;
  cmplx z;
  char ch='a',side[]="RI",label[100],xopt[]= "BSTNZH",yopt[]="BSTN",value[]="";
  char extra_text[max_len],ftext[max_len],ltext[max_len],help_text[2*max_len]; 
 
  if (read_input_opt(argc,argv,&istck))
    return(0);
  //read the options 
  
  if (bresolve(argc,argv,&nbl,base))
    return(0);
  // find out baselines to be plotted
  
  read_fits(argv[1]);
  // read the input fits file 
  
  fprintf(stdout,"\t ntimes=%d nchans=%ld nstokes=%ld nbaselines=%d end_time=%2.6f\n",
	  ntimes,nchans,nstokes,nbaselines,time_data[ntimes]);
//  sprintf(extra_text,"For help press h ");
  sprintf(help_text,"<n:next><p:previous><t:lin/log><b:reset><d:quit><f:flag on/off>");
 
 
  xt   =(float *)malloc(ntimes*sizeof(float));
  
  for(i1 = 0; i1 < ntimes; i1++){
    ftime(time_data[i1],tmm,&xt[i1]);
  }// for i1
  

  ntimes1 = (int)((xt[ntimes-1]-xt[0])/(xt[1]-xt[0]));
  data =(float *)malloc(ntimes1*nchans*sizeof(float));
  
  tr[1] = xt[1]-xt[0];
  tr[2] = 0.0;
  tr[4] = 0.0;
  tr[5] = 1.0;
  scl   = 1.0;

  i = 1;

  
  x1 = xt[0];
  x2 = xt[ntimes-1];
  
  y1 = 0.0;
  y2 = (float) nchans; 
  
  tr[0] = x1;
  tr[3] = y1;
  
  l = 5;
  m = 8;

  n = dc[m+nants*l-1];

  sprintf(label,"%d - %d",l,m);
  
  //antenna pair to be plotted 
  
  for(i1 = 0; i1 < ntimes1; i1++){
    for(j1 = 0; j1 < nchans; j1++){
      data[i1+ntimes1*j1] = 0.0;
    }// for j1
  }// for i1
  
  for(i1 = 0; i1 < ntimes; i1++){
    i2 = (int) ((xt[i1]-xt[0])/tr[1]); 
    for(j1 = 0; j1 < nchans; j1++){
      l1 = istck  + nstokes *(j1+nchans*i1);
      z.re = Visib[n+nbaselines*l1].re;
      z.im = Visib[n+nbaselines*l1].im;
      
      if (Visib[n+nbaselines*l1].w > 0.0)
	data[i2+ntimes1*j1] = log10(cmod(z));
      else
	data[i2+ntimes1*j1] = 0.0;
      
	
    }// for j1
  }// for i1
       
 
  d1 = 0.0;
  d2 = 0.0;
  for(i1=0; i1 < ntimes1; i1++){
    for(j1 = 0; j1 < nchans; j1++){
      d1 = min(d1,data[i1+ntimes1*j1]);
      d2 = max(d2,data[i1+ntimes1*j1]);
    }
  }
  
  if (d1 == d2)
    d2 = d2+1.0;    
  
    cpgbeg(0,"?",1,1);
  
    cpgscr(0,0.0,0.0,0.0);
    cpgscr(1,1.0,1.0,1.0);


    for(j = 1; j < nc; j++){
      r = (float)(j-1)/(nc-1); 
      cpgscr(j,r,r,r);
    }
    
    cpgsci(1);

    cpgpap(10.0,0.75);
    cpgsch(0.75);
    cpgsvp(0.0,1.0,0.0,1.0); 
    cpgenv(x1,1.2*x2,y1,y2,2,-1);
    cpgswin(x1,x2,y1,y2);
    cpgtbox(xopt,xtick,nxsub,yopt,ytick,nysub);
    cpglab("time","chan",label);
    
    cpgscir(5,nc-5); 
    cpgimag(data,ntimes1,nchans,1,ntimes1,1,nchans,d1,d2,tr);
    
    //cpggray(data,ntimes1,nchans,1,ntimes1,1,nchans,d1,d2,tr);
    cpgsch(0.75);

    cpgsci(4);
    cpgwedg(side,0.2,3.0,d1, d2,value); 
      
    cpgend();
  
    return(0);
}//main


int read_input_opt(int nargc, char *nargv[],int *istk){
  int i;
  if (nargc < 2){
    fprintf(stderr,"\t use the following options \n");
    fprintf(stderr,"\t ./igrplot <FITS FILE> <pmod>  -stk <stoke> -ant <ant1> -base <ant2> \n");
    fprintf(stderr,"\t FITS FILE : input FITS File \n");
    fprintf(stderr,"\t pmod      : l for log  < default=non log> \n");
    fprintf(stderr,"\t stoke     : stoke parameters <default=0>\n");
    fprintf(stderr,"\t ant1      : First antenna    <default=all>\n");
    fprintf(stderr,"\t ant2      : Second antennas  <default=all> \n");
    return(-1);
  }
  
  for (i = 0; i < nargc; i++){
    if (!strcmp(nargv[i],"-stk")){
      *istk = atoi(nargv[++i]);
    }else{
      *istk = 0;
    }
    if (!strcmp(nargv[i],"-l"))
      pmod = 1;
  }// for i 
  return(0);
}

int flip(int i){
  int j ;
  if (i==0) 
    j=1;
  else
    j=0;
  return(j);
}// for flip 


double prd(double x){

 double day,y; 
  
 day = 24.0 * 3600.00;
 
 if (x < 0)
   y=x+day;
 else
   y=x; 
 
 return(y); 

}
