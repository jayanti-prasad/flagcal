#include <defs.h>
#include <flagcal.h>
#include<cpgplot.h>

int nplots = 4, imod = 0, isrc = 0, ihelp=0; // default values

int   bresolve(int , char *[],int *, float []);
//int   strcmp ( const char * , const char *  );
int   read_input_opt(int, char *[]);
int   flip(int); 

int main(int argc, char *argv[]){
  FILE *inpf;
  int i,j,k,l,n,l1,m1,iter,*src,*a1,*a2,*sr;
  int p=5,nxsub=10,nysub=10;
  float *uu,*vv,*ww,*a,*t,*w,*fdata,tmm[4],base[nbaselines];
  float ttt,x1,x2,y1,y2,tx1,tx2,ty1,ty2,*x,*y,*z,cc,dd,scl; 
  float xtick = 0.0,ytick = 0.0; 
  char ch, label[max_len],extra_text[max_len];
  int ntimes,npage,*isymb,nbl;
  char xopt[]="BSTNZH",yopt[]="BSTN",help_text[max_len];

  isymb=&p;
  
  if (read_input_opt(argc,argv) < 0)
    return(0);
  if(bresolve(argc,argv,&nbl,base) !=0)
    return(0);
  if (nbl < nplots)
    nplots = nbl;
  
  fprintf(stdout,"\t  # of baselines to be plotted =%d \n",nbl);
   sprintf(extra_text,"For help press h ");
 
 
  a1    = (int *)malloc(sizeof(int));
  a2    = (int *)malloc(sizeof(int));
  src   = (int *)malloc(sizeof(int));
  uu    = (float *)malloc(sizeof(float));
  vv    = (float *)malloc(sizeof(float));
  ww    = (float *)malloc(sizeof(float));
  w     = (float *)malloc(sizeof(float));
  a     = (float *)malloc(sizeof(float));
  t     = (float *)malloc(sizeof(float));
  fdata = (float *)malloc(sizeof(float));
  
  k = 0;
  i = 0;
  inpf = fopen(argv[1],"r");
  while(!feof(inpf)){
    fscanf(inpf,"%f %d %d %d %f %f %f %f %f %f\n",
	   &ttt,&src[i],&a1[i],&a2[i],&uu[i],&vv[i],&ww[i],&a[i],&t[i],&w[i]);
    ftime(ttt,tmm,&fdata[i]);
    i++;
    a1=realloc(a1,(i+1)*sizeof(int));
    a2=realloc(a2,(i+1)*sizeof(int));
    src=realloc(src,(i+1)*sizeof(int));
    uu=realloc(uu,(i+1)*sizeof(float));
    vv=realloc(vv,(i+1)*sizeof(float));
    ww=realloc(ww,(i+1)*sizeof(float));
    a=realloc(a,(i+1)*sizeof(float));
    t=realloc(t,(i+1)*sizeof(float));
    fdata=realloc(fdata,(i+1)*sizeof(float));
    w=realloc(w,(i+1)*sizeof(float));
  }// for while 
  n = i;

  ntimes = n/nbaselines; 
  
  sprintf(help_text,"<n:next><p:previous><=:zoom in><-:zoom out><t:amp/phase><b:reset><d:quit><0/1/2:all/1/2 src><h:help bar>");

  npage  = nbl/nplots;
  
  fprintf(stdout,"\t imod = %d nplots = %d isrc=%d npage=%d\n",imod,nplots,isrc,npage);
 
 
  if (nbl - nplots * npage > 0)
    npage++; 
  
  x  = (float *)malloc(ntimes*sizeof(float));
  y  = (float *)malloc(ntimes*sizeof(float));
  z  = (float *)malloc(ntimes*sizeof(float));
  sr = (int  *)malloc(ntimes*sizeof(int));
 
  i = 0;
  scl = 1.0;
  while (i < npage){
  loop1:
    x1=fdata[0];
    x2=fdata[n-1]; 
    iter=0;
    ch='A';
    while(ch=='A'){
      cpgbeg(0,"/xs",1,nplots);
      cpgpap(scl*10.0,0.75);
      if (nplots >=4) 
	cpgsch(2.0); 
      else
	cpgsch(1.0);
      cpgscr(0,0.0,0.0,0.0);
      cpgscr(1,1.0,1.0,1.0);
      cpgscr(2,1.0,0.0,0.0);//red
      cpgscr(3,0.0,1.0,0.0);//green 
      cpgscr(4,0.0,0.0,1.0);//blue
      cpgscr(5,1.0,1.0,0.0);//yellow 
      for(j = 0; j < nplots; j++){
	if (j+nplots*i < nbl){
	  l1  = base[j+nplots*i]/(nants+1) ;
	  m1  = base[j+nplots*i] - l1 * nants; 
	  sprintf(label,"%d-%d",l1,m1);
	  k  = 0;
	  y1 = 0.0;
	  y2 = 0.0;
	  for(l=0; l < n ; l++){
	    if(a1[l]==l1 && a2[l]==m1) {
	      x[k] = fdata[l];
	      if(imod==0) 
		y[k] = a[l];
	      else
	      y[k] = t[l];
	      z[k] = w[l]; 
	      sr[k] = src[l]; 
	      if(z[k] > 0.0){
		if (isrc == 0){
		  y1 = min(y1,y[k]);
		  y2 = max(y2,y[k]);
		}else{
		if(sr[k] == isrc){
		  y1 = min(y1,y[k]);
		  y2 = max(y2,y[k]);
		}// for if
		}// for else
	      }// for if 
	      k++;
	    }// for if
          }// for l
	  if (y1 == y2)
	    y2+=1.0;
	  cc = (float) j/nplots;
	  dd = (float) (j+1)/nplots;
	  cpgsvp(0.0,1.0,cc,dd);
	  cpgsci(5); 
	  cpgenv(x1,x2,y1,1.1*y2,2,-1);
	  cpgswin(x1,x2,y1,1.1*y2);
	  cpgtbox(xopt,xtick,nxsub,yopt,ytick,nysub); 
	  if(imod==0)
	    cpglab("","Jy",label);
	  else
	    cpglab("","Deg",label);
	  cpgsci(0);
	  cpgline(ntimes,x,y);
	  for(k=0; k < ntimes; k++){
	    if (z[k] < 0.0)
	      cpgsci(2);
	    else
	      cpgsci(3);
	    if (isrc == 0)
	      cpgpnts(1,&x[k],&y[k],isymb,1);
	    else{
	      if(sr[k] == isrc)
		cpgpnts(1,&x[k],&y[k],isymb,1);
	    }
	  }//for k 
	}// for if 
	 cpgsci(1); 
	
	 cpgsci(2); 
	 if (nplots >=3)
	   cpgsch(2);
	 else
	   cpgsch(1);
	 
	 if(j==0 && ihelp==0)
	   cpgmtxt("T",1,0.1,0.5,extra_text);
	 if(j==(nplots-1) && ihelp==1)
	   cpgmtxt("B",3,0.5,0.5,help_text);
	 
      }// for j 
     cpgsci(1); 
     cpgcurs(&tx1,&ty1,&ch);
      if (ch=='a' || ch =='A'){
	cpgband(2,1,tx1,ty1,&tx2,&ty2,&ch);
	x1=min(tx1,tx2);
	x2=max(tx1,tx2);
      }
      switch (ch){
      case 'b': 
       goto loop1; 
       break;
      case 'p':
	i--;
	if (i < 0)
	  i=npage-1;
        goto loop1;
	break;
      case 'n':
	i++;
	if (i==npage)
	  i=0;
      	goto loop1; 
	break;
      case 'd':
	cpgend();
	return(0);
      case '-':
	scl *=0.9;
	goto loop1; 
	break;
      case '=':
	scl *=1.1;
	goto loop1; 
	break;
      case '0':
	isrc = 0;
	goto loop1; 
	break;
      case '1':
	isrc = 1;
	goto loop1; 
	break;
      case '2':
	isrc = 2;
	goto loop1; 
	break;
      case '3':
	isrc = 3;
	goto loop1; 
	break;
      case 't':
	imod = flip(imod);
	goto loop1;
	break;
      case 'h':
	ihelp=flip(ihelp);
        goto loop1;
	break;  
      }// for switch 
      iter++;
    }// for while 
     cpgend();
  }// for i
  return(0);
}

int read_input_opt(int nargc, char *nargv[]){
  int i;
  if (nargc < 2){
    fprintf(stdout,"\t * This program has become quite complex due to many options.\n");
    fprintf(stdout,"\t - It takes a single stokes & single channel file produced by\n");
    fprintf(stdout,"\t   the readfits and plots amplitude/phase  of the visibilties.\n");
    fprintf(stdout,"\t - One can limit the baselines to be plotted by using the options \n");
    fprintf(stdout,"\t   -ant & -base \n");
    fprintf(stdout,"\t - By pressing the key 't' one can switch between amplitude & phase \n");
    fprintf(stdout,"\t - By pressing the key '0','1' or '2' one can plot the data for all the sources\n");
    fprintf(stdout,"\t   or source 1, 2 respectively. \n");
    fprintf(stdout,"\t - By pressing 'h' one can get the help bar at the bottom \n");
    fprintf(stderr,"\t *** use the follwing options ***\n\n");
    fprintf(stderr,"\t ./ivplot <input file>  <imod> -src <isrc> -nplots <#nplots>  -ant <ant1> -base <ant2> \n");
    fprintf(stdout,"\t input file: file created by readfits \n");
    fprintf(stdout,"\t imod      : a=amp,p=phase <default=0>\n");
    fprintf(stdout,"\t ant       : first antenna \n");
    fprintf(stdout,"\t base      : second antennas  \n");
    fprintf(stderr,"\t Warning ! odering of the options is important !\n");
    return(-1);
  }else{
    for (i = 0; i < nargc; i++){
      if (!strcmp(nargv[i],"p"))
	imod = 1;
      if (!strcmp(nargv[i],"-nplots"))
	nplots = atoi(nargv[++i]);
      if (!strcmp(nargv[i],"-src"))
	isrc = atoi(nargv[++i]);
      
    }// for i
    return(0);
  }// for else 
}

int flip(int i){
  int j ;
  if (i==0) 
    j=1;
  else
    j=0;
  return(j);
}



