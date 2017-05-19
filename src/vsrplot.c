#include<defs.h>
#include<flagcal.h>
#include<cpgplot.h>


float  *tdata,*vsr;    

int ichan,istock=0,nplots=4;

int nchans1, nstokes1;
int   ftime(float, float[],float *);
int   bresolve(int , char *[],int *, float []);
int   read_data(char []);

int main(int argc, char*argv[]){
  char ch='A', label[max_len],xopt[]="BSTNZH",yopt[]="BSTN";
  int i,j,l,m,n,i1,l1;
  int nxsub = 10,nysub = 10,*isymb1,*isymb2,p1=2,p2=5;
  float *x,*y,tmm[4],xtick=0.0,ytick=0.0,scl;
  float x1,x2,y1=0.0,y2=1.1,t_vsr=0.95;
  float base[nbaselines];
  int a1=4,b1=10,a2=12,b2=16,a3=12,b3=17;

  isymb1 =&p1;
  isymb2 =&p2;

  j   = 0;
  for(l=0; l <nants ; l++){
    for(m=l+1;  m <nants; m++){
      id[j] = m + nants *l-1;
      dc[id[j]] = j;
      j++;
    }// for m
  }// for l
  
  if(read_data(argv[1]) !=0)
    return(-1);
  
  fprintf(stdout,"\t nchans = %ld nstokes = %ld ntimes = %d\n",nchans,nstokes,ntimes);
  
  x = (float *)malloc(ntimes*sizeof(float));
  y = (float *)malloc(ntimes*sizeof(float));
  
  ftime(tdata[0],tmm,&x1); 
  ftime(tdata[ntimes-1],tmm,&x2); 
  
  printf("%2.6f %2.6f\n",tdata[0],tdata[ntimes-1]);
 
  istock = 0;
  ichan  = 12;
  
 
   
  cpgbeg(0,"?",1,3);
  cpgpap(8.0,0.75);
  cpgsch(2.0);

  // first panel 
  
  for(i1 = 0; i1 < ntimes; i1++){
    ftime(tdata[i1],tmm,&x[i1]); 
    l1 = istock + nstokes1 *(ichan+nchans1*i1);
    n  = dc[b1+nants*a1-1];
    y[i1] = vsr[n+nbaselines*l1]; 
  }// for i1
  
  sprintf(label,"%2d-%2d",a1,b1); 

  cpgsvp(0.0,1.0,0.0,0.3);
  cpgenv(x1,x2,y1,y2,2,-1);
  
  cpgswin(x1,x2,y1,y2);
  cpgtbox(xopt,xtick,nxsub,yopt,ytick,nysub);
  cpglab("", "vsr",label);
  for(i1 = 0; i1 < ntimes; i1++){
    if(y[i1] > t_vsr) 
      cpgpnts(1,&x[i1],&y[i1],isymb1,1);
    else
      cpgpnts(1,&x[i1],&y[i1],isymb2,1);
  }// for i1
   
  for(i1 = 0; i1 < ntimes; i1++){
    ftime(tdata[i1],tmm,&x[i1]); 
    l1 = istock + nstokes1 *(ichan+nchans1*i1);
    n  = dc[b2+nants*a2-1];
    y[i1] = vsr[n+nbaselines*l1]; 
  }// for i1
 
  sprintf(label,"%2d-%2d",a2,b2); 

  cpgsvp(0.0,1.0,0.3,0.6);
  cpgenv(x1,x2,y1,y2,2,-1);
  cpgswin(x1,x2,y1,y2);
  cpgtbox(xopt,xtick,nxsub,yopt,ytick,nysub);
  cpglab("", "vsr",label);
  for(i1 = 0; i1 < ntimes; i1++){
    
    if(y[i1] > t_vsr) 
      cpgpnts(1,&x[i1],&y[i1],isymb1,1);
    else
      cpgpnts(1,&x[i1],&y[i1],isymb2,1);

  }// for i1

  for(i1 = 0; i1 < ntimes; i1++){
    ftime(tdata[i1],tmm,&x[i1]); 
    l1 = istock + nstokes1 *(ichan+nchans1*i1);
    n  = dc[b3+nants*a3-1];
    y[i1] = vsr[n+nbaselines*l1]; 
  }// for i1

  sprintf(label,"%2d-%2d",a3,b3);
  cpgsvp(0.0,1.0,0.6,0.9);
  cpgenv(x1,x2,y1,y2,2,-1);
  cpgswin(x1,x2,y1,y2);
  cpgtbox(xopt,xtick,nxsub,yopt,ytick,nysub);
  cpglab("", "vsr",label);
  for(i1 = 0; i1 < ntimes; i1++){
    if(y[i1] > t_vsr) 
      cpgpnts(1,&x[i1],&y[i1],isymb1,1);
    else
      cpgpnts(1,&x[i1],&y[i1],isymb2,1);
  }// for i1
  
  cpgend();
  
  i++;
  
  return(0);
}

int read_data(char gfile[]){
  FILE *inpf;
  int l,i1,j1,k1;
  long l1;
  
  inpf = fopen(gfile,"rb");
  fread(&ntimes,sizeof(int),1,inpf);
  fread(&nchans,sizeof(long),1,inpf);
  fread(&nstokes,sizeof(long),1,inpf);
  fread(&thrs_vsr,sizeof(float),1,inpf);
 
  nstokes1 = (int) nstokes; 
  nchans1  = (int) nchans; 
  vsr      = (float *)malloc(ntimes*nchans1*nstokes1*nbaselines*sizeof(float));
  tdata    = (float *)malloc(ntimes*sizeof(float));
  
  for(i1 = 0; i1 < ntimes; i1++){
    tdata[i1] = 0.0;
    for(j1 = 0; j1 < nchans1; j1++){
      for(k1 = 0; k1 < nstokes1; k1++){
	l1 = k1 + nstokes1 *(j1+nchans1*i1);
	for(l = 0; l < nbaselines; l++){
	  vsr[l+nbaselines*l1] = 0.0;
	}// for l
      }// for k1
    }// got j1
  }// for i1
  
  for(i1 = 0; i1 < ntimes; i1++){
    fread(&tdata[i1],sizeof(float),1,inpf);
    for(j1 = 0; j1 < nchans1; j1++){
      for(k1 = 0; k1 < nstokes1; k1++){
        l1 = k1 + nstokes1 *(j1+nchans1*i1);
        for(l = 0; l < nbaselines; l++){
          fread(&vsr[l+nbaselines*l1],sizeof(float),1,inpf);
	}// for l
      }// for k1
    } // for j1
  }// for i1
  fclose(inpf);
  return(0);
}//read gain 

