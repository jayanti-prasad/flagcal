/*-------------------------------------------------------------------------
- This is the main flagging program. It computes   (a) vsr,
  (b) bad_ant,  (c) bad_chan and (d) bad_base. Which are used
  later on for flagging. 

- This program also writes the above flagging information into files 
  which can be plotted using the plotting programs provided with 
  the package or with gnuplot. 

- Jayanti Prasad
Sat Jul 24 10:04:21 IST 2010
--------------------------------------------------------------------------*/

#include <defs.h>
#include <flagcal.h>

int longest_block(int , float , float  []);

int   flag_data(){
  int i1,i2,i3,j,j1,k1,l,l1,l2,m,n,n1,icond;
  float x[max_nsb],lblock,X,Y,S,V;
  cmplx z; 
      
  if((vsr = (float *)malloc(nsblocks*nchans*nstokes*nbaselines*sizeof(float)))==NULL)
    {perror("ERROR ! memory allocation for  vsr in flag_data()\n"); return(-1);}  
  if((dant = (int   *)malloc(nants*sizeof(int)))==NULL)
    {perror("ERROR ! memory allocation for  dant in flag_data()\n"); return(-1);}
  if((bad_ant  = (float *)malloc(nscans*nants*sizeof(float)))==NULL)
    {perror("ERROR ! memory allocation for  bad_ant in flag_data()\n"); return(-1);}
  if((bad_base  = (float *)malloc(nscans*nbaselines*sizeof(float)))==NULL)
    {perror("ERROR ! memory allocation for  bad_base in flag_data()\n"); return(-1);}
  if((bad_chans = (float *)malloc(nscans*nchans*sizeof(float)))==NULL)
    {perror("ERROR ! memory allocation for  bad_CHANS flag_data()\n"); return(-1);}
  
  for(i1 = 0; i1 < nscans; i1++){
    for(i2 = 0; i2 < nsb[i1]; i2++){
      j = sbid_in[i2+max_nsb*i1];  
      for(j1 = 0; j1 < nchans; j1++){
	for(k1 = 0; k1 < nstokes; k1++){
	  l1 = k1 + nstokes *(j1+nchans*j);
	  for(n = 0; n < nbaselines; n++){
            vsr[n+nbaselines*l1] = 0.0;
	  }// for n
	} // for k1
      } // for j1
    }// for i2
    for(l=0; l < nants; l++)
      bad_ant[l+nants*i1] = 0.0; 
    for(n=0; n < nbaselines; n++)
      bad_base[n+nbaselines*i1] = 0.0;
    for(j1=0; j1 < nchans; j1++)
      bad_chans[j1+nchans*i1] = 0.0;
  } // for i1 
  
  for(i1 = 0; i1 < nscans; i1++){
    if (src[idmap[sid[i1]]].type != 2){
      for(i2 = 0; i2 < nsb[i1]; i2++){
	j = sbid_in[i2+max_nsb*i1];  
	for(j1 = 1; j1 < nchans; j1++){
	  for(k1 = 0; k1 < nstokes; k1++){
	    l1  =  k1 + nstokes * (j1+nchans*j);
	    for(n = 0; n < nbaselines; n++){
              l = (id[n]+1)/nants;
              m = (id[n]+1) - nants *l;
              S = 0.0;V = 0.0;X = 0.0; Y = 0.0;
	      for(i3  = sb_start[j]; i3 < sb_end[j]; i3++){
		l2    =  k1 + nstokes * (j1+nchans*i3);
		z     = Cmplx(Visib[n+nbaselines*l2].re,Visib[n+nbaselines*l2].im);
		if (Visib[n+nbaselines*l2].w > 0.0){
		  X += z.re; 
		  Y += z.im; 
		  S +=  cmod(z);
		}// for if
   	      }// for i3
	      if (S > 0.0){
		vsr[n+nbaselines*l1] = sqrt(X*X+Y*Y)/S;
	      }// for if
            } // for n 
	  } // for k1
	} // for j1
      } // for i2
    } // for if 
  }// for i1
  
  for(i1 = 0; i1 < nscans; i1++){
    if (src[idmap[sid[i1]]].type != 2){
      for(j1 = 0; j1 < nchans; j1++){
	for(k1 = 0; k1 < nstokes_half; k1++){
	  for(n = 0; n < nbaselines; n++){
	    l = (id[n]+1)/nants;
	    m = (id[n]+1) - nants *l;
	    n1= 0;
	    icond = 0;
	    for(i2 = 0; i2 < nsb[i1]; i2++){
	      j = sbid_in[i2+max_nsb*i1];  
	      l1  =  k1 + nstokes * (j1+nchans*j);
	      x[i2] = vsr[n+nbaselines*l1];
	      if (x[i2] > thrs_vsr) n1++;
	    }// for i2;
	    if ((float) n1/nsb[i1] < thrs_ngood) 
	      icond = 1; 
	    if(!icond){	    
	      lblock = (float) longest_block(nsb[i1],thrs_vsr,x)/nsb[i1];
	      if(lblock < thrs_block)
		icond = 1;
	    }//for if
	    if(icond){
	      bad_ant[l+nants*i1]+=1.0;
	      bad_ant[m+nants*i1]+=1.0;
	      bad_chans[j1+nchans*i1]+=1.0;
	      bad_base[n+nbaselines*i1]+=1.0;
	    } // for if
	  }// for n
	} // for k1
      }// for j1
      for(l = 0; l  < nants; l++)
      	bad_ant[l+nants*i1]/=(nchans*nstokes_half*(nants-1));
      for(n = 0; n  < nbaselines; n++)
      	bad_base[n+nbaselines*i1]/=(nchans*nstokes_half); 
      for(j1= 0; j1 < nchans; j1++)
	bad_chans[j1+nchans*i1]/=nbaselines*nstokes_half;
    }// for if 
  }// for i1
  
  // this is for the source scan  

  for(i1 = 0; i1 < nscans; i1++){
    if (src[idmap[sid[i1]]].type == 2){
      for(l= 0; l < nants; l++){
	if (i1==0)
	  bad_ant[l+nants*i1] =  bad_ant[l+nants*(i1+1)];
	if (i1==nscans-1)
	  bad_ant[l+nants*i1] =  bad_ant[l+nants*(i1-1)];
	if (i1 > 0 && i1 < nscans-1)
	  bad_ant[l+nants*i1] = 
	    max(bad_ant[l+nants*(i1-1)],bad_ant[l+nants*(i1+1)]);
      }// for l
      for(j1= 0; j1 < nchans; j1++){
	if (i1==0)
	  bad_chans[j1+nchans*i1] =    bad_chans[j1+nchans*(i1+1)];
	if (i1==nscans-1)
	  bad_chans[j1+nchans*i1] =    bad_chans[j1+nchans*(i1-1)];
	if (i1 > 0 && i1 < nscans-1)
	  bad_chans[j1+nchans*i1] =
	    max(bad_chans[j1+nchans*(i1-1)],bad_chans[j1+nchans*(i1+1)]);
      }// for j1
      for(n=0; n < nbaselines; n++){
       	if (i1==0)
	  bad_base[n+nbaselines*i1] =    bad_base[n+nbaselines*(i1+1)];
	if (i1==nscans-1)
	  bad_base[n+nbaselines*i1] =    bad_base[n+nbaselines*(i1-1)];
	if (i1 > 0 && i1 < nscans-1)
	  bad_base[n+nbaselines*i1] = 
	    max(bad_base[n+nbaselines*(i1-1)],bad_base[n+nbaselines*(i1+1)]);
      }// for n
    } // for if
  } // for i1
  
  // completely dead antenna
  for(l = 0; l < nants; l++)
    dant[l] = 0;
  for(i1 = 0; i1 < nscans; i1++){
    for(l = 0; l < nants; l++){
      if (bad_ant[l+nants*i1] > thrs_ant)
	dant[l]++;
    }// for l 
  }// for i1



  return(0);
}

