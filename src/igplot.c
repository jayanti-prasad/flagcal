#include <defs.h>
#include <flagcal.h>
#include<cpgplot.h>
/*---------------------------------------------------------------------------
  This program has been greatly improved and have been made user friendly 
  by adding the help bar and parameter labels. At present I do not think
  if any thing extra can be added. However, suggestions are welcome. 
                             --------- Jayanti Prasad
  program last modified on   -  Thu Jul  8 11:55:43 IST 2010       
 ---------------------------------------------------------------------------*/

// default values 
int ichan,istock=0,imod=0,nplots=4,ihelp=0;

cmplx *g; 
float  *tdata; 
int flip(int );
int ftime(float, float[],float *);
int read_input(int , char *[]);
int read_gain(char []);

int main(int argc, char*argv[]){
  int i,j,i1,l,l1,npage;
  int p = 5,justdoit = 1,nxsub = 10,nysub = 10,*isymb;
  float *x,*y,tmm[4],xtick=0.0,ytick=0.0,scl;
  float x1,x2,y1,y2,tx1,tx2,ty1,ty2,a1,b1; 
  char ch='A', label[max_len],xopt[]="BSTNZH",yopt[]="BSTN";
  char  help_text[max_len],para_text[max_len],extra_text[max_len];
 
  isymb =&p;
  sprintf(help_text,"<n:next><p:previous><=:zoom in><-:zoom out><t:amp/phase><b:reset><d:quit><h:help>");
  if (read_input(argc,argv) != 0)
    return(-1);
  npage  = nants/nplots;
  if (nants - npage * nplots > 0)
    npage++;
  sprintf(para_text,"channel=%d, stokes=%d",ichan,istock);
  sprintf(extra_text,"For help press h ");
  if(read_gain(argv[1]) !=0)
    return(-1);
  fprintf(stdout,"nchans=%ld nstokes=%ld ntimes=%d\n",nchans,nstokes,ntimes);
  x = (float *)malloc(ntimes*sizeof(float));
  y = (float *)malloc(ntimes*sizeof(float));
  scl = 1.0;
  i   = 0;   
  while (i < npage ){
  loop1:
    ftime(tdata[0],tmm,&x1); 
    ftime(tdata[ntimes-1],tmm,&x2); 
    //initlization for every page
    ch='A';
    while(ch=='A' && x1 !=x2){
      cpgbeg(0,"/xs",1,nplots);
      cpgpap(scl * 10.0,0.8);
      if (nplots < 4)
        cpgsch(1.0); 
      else
        cpgsch(2.0);
      cpgscr(0,0.0,0.0,0.0);//black 
      cpgscr(1,1.0,1.0,1.0);//white
      cpgscr(2,1.0,0.0,0.0);//red  
      cpgscr(3,0.0,1.0,0.0);//green
      cpgscr(4,1.0,1.0,0.0);//yellow 
      cpgsci(1);
      for(j = 0; j < nplots; j++){
	l = j + nplots*i; 
	if (l < nants){
	  y1 = 0.0;
	  y2 = 0.0;
	  sprintf(label,"Antennas %d",l+1);
	  for(i1 = 0; i1 < ntimes; i1++){
	    ftime(tdata[i1],tmm,&x[i1]);
	    l1 = istock + nstokes *(ichan+nchans*i1);
	    switch (imod) {
	    case 1:
	      y[i1] = cphase(g[l+nants*l1]);
	      break;
	    case 2:
	      y[i1] = g[l+nants*l1].re;
	      break;
	    case 3:
	    y[i1] = g[l+nants*l1].im;
	    break;
	    default:
	      y[i1] = cmod(g[l+nants*l1]);
	    break;
	    }//for switch 
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
	  if (j==0){
	    if (nplots >=3)
	      cpgsch(2);
	    else
	      cpgsch(1);
	    cpgsci(2);
	    cpgmtxt("T",1,0.9,0.5,para_text);
	    if(ihelp==0) 
	      cpgmtxt("T",1,0.1,0.5,extra_text);
	  }// for if
	  cpgsci(4); 
	  switch (imod) {
	  case 1:
	    cpglab("", "Deg",label);
	    break;
	  default:
	    cpglab("", "",label);
	    break;
	  }//for switch 
	  cpgsci(3); 
	  for(i1 = 0; i1 < ntimes; i1++)
	    cpgpnts(1,&x[i1],&y[i1],isymb,1);
	}// for if
	cpgsci(2);
	if(j==(nplots-1) && ihelp==1){
	   if (nplots >=3)
	     cpgsch(2);
	    else
	      cpgsch(1);
	   cpgmtxt("B",3,0.5,0.5,help_text);
	}// for of
      }// for j 
      cpgsci(1);
      cpgcurs(&tx1,&ty1,&ch);
      if (ch=='a' || ch=='A'){
	cpgband(2,1,tx1,ty1,&tx2,&ty2,&ch);
        x1=min(tx1,tx2);
	x2=max(tx1,tx2);
      }// for if
      switch (ch){
      case 'b':
	goto loop1;
        break;
      case 'p':
        i--;
        if ( i < 0) i=npage-1;
        goto loop1;
        break;
      case 'n':
        i++;
        if (i > npage-1) i=0;
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
      case 't':
        imod = flip(imod);
        goto loop1;
        break;
      case 'h':
        ihelp=flip(ihelp);
        goto loop1;
        break;
      }//for switch
      
    }// for while 
    cpgend();
  }// for i
  return(0);
}// main program ends here 

int read_gain(char gfile[]){
  FILE *inpf;
  int i1,j1,k1,l,l1;
  inpf = fopen(gfile,"rb");
  fread(&ntimes,sizeof(int),1,inpf);
  fread(&nchans,sizeof(long),1,inpf);
  fread(&nstokes,sizeof(long),1,inpf);
  g = (cmplx *)malloc(ntimes*nchans*nstokes*nants*sizeof(cmplx));
  tdata = (float *)malloc(ntimes*sizeof(float));
  for(i1 = 0; i1 < ntimes; i1++)
    tdata[i1] = 0.0; 
  for(i1 = 0; i1 < ntimes; i1++){
    fread(&tdata[i1],sizeof(float),1,inpf);
    for(j1 = 0; j1 < nchans; j1++){
      for(k1 = 0; k1 < nstokes; k1++){
        l1 = k1 + nstokes *(j1+nchans*i1);
        for(l = 0; l < nants; l++){
          fread(&g[l+nants*l1],sizeof(cmplx),1,inpf);
	}// for l
      }// for k1
    } // for j1
  }// for i1
  fclose(inpf);
  return(0);
}//read gain 

int read_input(int nargc, char *nargv[]){
  int i;
  if (nargc < 4){
    fprintf(stderr,"\t use the following options: \n");
    fprintf(stderr,"\t ./igplot  <inputfile> <mod> -c <chan> -s <stoke> -n <nplots> \n");
    fprintf(stderr,"\t input file: file created by flagcal (gain.dat)  \n");
    fprintf(stderr,"\t mod       : a=amp,p=phase,r=real,i=imag <default=a>\n");
    fprintf(stderr,"\t chan      : channel for which gains are plotted\n");
    fprintf(stderr,"\t stoke     : stoke parameter for which gains are plotted <default=0>\n");
    fprintf(stderr,"\t nplots    : number of plots per page <default>\n"); 
    return(-1);
  }// for if
  for(i = 0; i < nargc; i++){
    if (!strcmp(nargv[i],"a"))
      imod  = 0;
    if (!strcmp(nargv[i],"p"))
      imod  = 1;
    if (!strcmp(nargv[i],"r"))
      imod  = 2;
    if (!strcmp(nargv[i],"i"))
      imod  = 3;
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
