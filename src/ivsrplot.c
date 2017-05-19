#include <defs.h>
#include <flagcal.h>
#include<cpgplot.h>

float *tdata; 
float  *vsr;    
int ichan,istock=0,nplots=4;

int nchans1, nstokes1;
int  ftime(float, float[],float *);
int  bresolve(int , char *[],int *, float []);
int  read_data(char []);
int read_input_option(int , char *[]);
int ihelp=0; 
int flip(int );

int main(int argc, char*argv[]){
  int i,j,m,n,i1,l,l1,npage,nbl;
  int p = 5,justdoit = 1,nxsub = 10,nysub = 10,*isymb;
  float *x,*y,tmm[4],xtick=0.0,ytick=0.0,scl;
  float x1,x2,y1,y2,tx1,tx2,ty1,ty2,a1,b1; 
  char ch='A', label[max_len],xopt[]="BSTNZH",yopt[]="BSTN",help_text[max_len];
  char extra_text[max_len];
  float base[nbaselines]; 

  isymb =&p;
 
  j   = 0;
  for(l=0; l <nants ; l++){
    for(m=l+1;  m <nants; m++){
      id[j] = m + nants *l-1;
      dc[id[j]] = j;
      j++;
    }// for m
  }// for l

  sprintf(help_text,"<n:next><p:previous><=:zoom in><-:zoom out><b:reset><d:quit><h:help>");
  sprintf(extra_text,"For help press h ");


  if (read_input_option(argc,argv) != 0)
    return(-1);
  
  fprintf(stdout,"\t chan=%d stock=%d nplots=%d\n",ichan,istock,nplots); 
  
  if(bresolve(argc,argv,&nbl,base) !=0) 
    return(-1);
  
  npage  = nbl/nplots;
  if (nbl > npage * nplots)
    npage++; 
  
  fprintf(stdout,"\t nbl=%d npage=%d\n",nbl,npage); 
  
  if(read_data(argv[1]) !=0)
    return(-1);
  
  fprintf(stdout,"\t nchans = %ld nstokes = %ld ntimes = %d\n",nchans,nstokes,ntimes);
  
  x = (float *)malloc(ntimes*sizeof(float));
  y = (float *)malloc(ntimes*sizeof(float));
  i = 0;   

  scl = 1.0;

  while (i < npage ){
  loop1:

    ftime(tdata[0],tmm,&x1); 
    ftime(tdata[ntimes-1],tmm,&x2); 

    //initlization for every page
    ch='A';
    justdoit = 1; 
    while(justdoit==1 && ch=='A' && x1 !=x2){
      cpgbeg(0,"/xs",1,nplots);
      cpgpap(scl*10.0,0.8);
      if (nplots < 4)
        cpgsch(1.0); 
      else
        cpgsch(2.0);
     
      cpgscr(0,0.0,0.0,0.0);//black
      cpgscr(1,1.0,1.0,1.0);//white
      cpgscr(2,1.0,0.0,0.0);//red
      cpgscr(3,0.0,1.0,0.0);//green
      cpgscr(4,1.0,1.0,0.0);//yellow
      
      for(j = 0; j < nplots; j++){
	if (j+nplots*i < nbl){
	  y1 = 0.0;
	  y2 = 0.0;
	  l  = base[j+nplots*i]/(nants+1) ;
	  m  = base[j+nplots*i] - l * nants;
	  n = dc[m+nants*l-1]; 
	  sprintf(label,"%d-%d",l,m);
	  for(i1 = 0; i1 < ntimes; i1++){
	    ftime(tdata[i1],tmm,&x[i1]); 
	    l1 = istock + nstokes1 *(ichan+nchans1*i1);
	    y[i1] = vsr[n+nbaselines*l1]; 
	    y1=min(y1,y[i1]);
	    y2=max(y2,y[i1]);
	  }// for i1
	  a1 = (float) j/nplots;
	  b1 = (float) (j+1)/nplots;
	  cpgsvp(0.0,1.0,a1,b1);
	  cpgsci(4); 
	  if(y1==y2)
	    y2 = 1.0;
	  cpgenv(x1,x2,y1,1.1*y2,2,-1);
	  cpgswin(x1,x2,y1,1.1*y2);
	  cpgtbox(xopt,xtick,nxsub,yopt,ytick,nysub);
          cpglab("", "vsr",label);
	  if (j==0){
            if (nplots >=3)
              cpgsch(2);
            else
              cpgsch(1);
            cpgsci(2);
            if(ihelp==0)
              cpgmtxt("T",1,0.1,0.5,extra_text);
          }// for if
	  cpgsci(3); 
	  for(i1 = 0; i1 < ntimes; i1++){
	    if (y[i1] < thrs_vsr)
	      cpgsci(2);
	    else
	      cpgsci(3);
	    cpgpnts(1,&x[i1],&y[i1],isymb,1);
	  }// for i1
         cpgsci(2);
        if(j==(nplots-1) && ihelp==1){
           if (nplots >=3)
             cpgsch(2);
            else
              cpgsch(1);
           cpgmtxt("B",3,0.5,0.5,help_text);
        }// for of

	}// for if
      }// for j 
      cpgsci(1);
      cpgcurs(&tx1,&ty1,&ch);
      if (ch=='a' || ch=='A'){
	cpgband(2,1,tx1,ty1,&tx2,&ty2,&ch);
        x1=min(tx1,tx2);
	x2=max(tx1,tx2);
      }
      switch (ch){
      case 'b':
        justdoit=0;
        goto loop1;
        break;
      case 'p':
        i--;
	if(i <0 ) i=npage-1;
        goto loop1;
        break;
      case 'n':
        i++;
	if (i ==npage) i=0;
        justdoit=0;
        goto loop1;
	break;
      case '-':
	scl *=0.9;
	goto loop1;
	break;
      case '=':
        scl *=1.1;
        goto loop1;
        break;
      case 'd':
        cpgend();
        return(0);
	break;
        case 'h':
        ihelp=flip(ihelp);
        goto loop1;
        break;
      }//for switch
    }// for while 
    cpgend();
    i++;
  }// for i
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


int read_input_option(int nargc, char *nargv[]){
  int i;
  if (nargc < 4){
    fprintf(stderr,"\t use the following options \n");
    fprintf(stderr,"\t ./ivsrplot <inputfile> -c <chan> -s <stoke> -n <nplots> -ant <> -base <> \n");
    fprintf(stderr,"\t input file: file created by flagcal (vsr.dat)\n");
    fprintf(stderr,"\t chan      : channel for which gains are plotted\n");
    fprintf(stderr,"\t stoke     : stoke for which gains are plotted <default=0>\n");
    fprintf(stderr,"\t nplots    : number of plots per page <default=4>\n"); 
    return(-1);
  }// for if
  
  for(i = 0; i < nargc; i++){
    if (!strcmp(nargv[i],"-c")){
      ichan = atoi(nargv[++i]);
    }if (!strcmp(nargv[i],"-s")){
      istock = atoi(nargv[++i]);
    }if (!strcmp(nargv[i],"-n")){
      nplots = atoi(nargv[++i]);
    }
  } // for i

  return(0);
}

int flip(int i){
  int j ;
  if (i==0)
    j=1;
  else
    j=0;
  return(j);
}






