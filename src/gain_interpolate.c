#include <defs.h>
#include <flagcal.h>

cmplx interpol(cmplx ,cmplx ,float);
float slope(int,int);
float median (int , float []);

int gain_intrapolate(){
  int i1, i2, j, j1,k1, l,np,l1, l2, l3,ii,jj;
  cmplx gtmp[nscans*nants*nchans*nstokes];
  float grad,x[ntimes],y[ntimes]; 
  
  fprintf(stdout,"\t  INTRAPOLATING gains for source\n");
  
  for(i1 = 0; i1 < nscans; i1++){
    for(j1 = 0; j1 < nchans; j1++){
      for(k1 = 0; k1 < nstokes; k1++){
	l1 = k1 + nstokes*(j1+nchans*i1);
	for(l = 0; l < nants; l++)
	  gtmp[l+nants*l1] = Cmplx(0.0,0.0);
      }// for k1
    }// for j1
  }// for i1 


  // gain for every scan is initilized to zero.
  
  for(i1 = 0; i1 < nscans; i1++){
    if (src[idmap[sid[i1]]].type != 2){
      for(j1 = 0; j1 < nchans; j1++){
	for(k1 = 0; k1 < nstokes; k1++){
	  l2 = k1 + nstokes*(j1+nchans*i1); 
	  for(l = 0; l < nants; l++){
	    np = 0;
	    for(i2 = 0; i2 < nab[i1]; i2++){
	      j = abid_in[i2+max_nab*i1]; 
	      l1 = k1 + nstokes*(j1+nchans*j);
	      if (cmod(g1[l+nants*l1]) < 0.0)  continue;
	      x[np]  = g1[l+nants*l1].re;
	      y[np]  = g1[l+nants*l1].im;
	      np++;
	    }// for i2 
	    if (np > 0)
	      gtmp[l+nants*l2] = Cmplx(median(np,x),median(np,y));
	    else
	      gtmp[l+nants*l2] = Cmplx(0.0,0.0);
	  }// for l
	}// for k1
      }// for j1
    } // for if
  }// for i1
  
  fprintf(stdout,"\t average computed\n");
  
  // now even the first scan can be that of a the science source 


  for(i1=0; i1 < nscans; i1++){
    if (src[idmap[sid[i1]]].type == 2 && bs[i1] > 0){ // if the scan is of the science source

      if (lc[i1]> 0)
	ii=lc[i1]; //  if the left scan is sepecified by the user
      else
	if (i1==0) // if the left scan is the first scan 
	  ii = i1+1;
	else
	  ii=i1-1;
      
      if (rc[i1] > 0)
	jj=rc[i1]; // if the right scan is specified by the user
      else
	if (i1==nscans-1) // if the right scan is the last scan 
	  jj=i1-1;
	else
	  jj=i1+1;
     
        if (ii > jj){
          fprintf(stdout,"ERROR ! left scan for interpolating is right of the right scan !\n");
	  return(-1);  
	}
      
      for(i2 = 0; i2 < nab[i1]; i2++){ // and is not flagged. 
	j  = abid_in[i2+max_nab*i1]; 
	grad = slope(i1,i2);
	for(j1 = 0; j1 < nchans; j1++){
	  for(k1 = 0; k1 < nstokes; k1++){
	    l1 = k1 + nstokes*(j1+nchans*j);
	    l2 = k1 + nstokes*(j1+nchans*ii);
	    l3 = k1 + nstokes*(j1+nchans*jj);
	    for(l=0; l < nants; l++)
	      g1[l+nants*l1] = interpol(gtmp[l+nants*l2],gtmp[l+nants*l3],grad);
	  }// for k1
	}// for j1
      }// for i2
    }// for if
  }// for i1

  
  return(0);
}


float slope(int i1, int i2){
  int j, j2, j3;
  float grad; 
  
  j  = abid_in[i2+max_nab*i1]; 
  // this is the current block in soource scan 
  
  j2 = abid_in[(nab[i1-1]-1)+ max_nab*(i1-1)]; 
  // this is the last  block of the left scan
  
  j3 = abid_in[0 + max_nab*(i1+1)]; 
  // this is the first block of the right scan 
  
  grad = (abmt[j]-abmt[j2])/(abmt[j3]-abmt[j2]);
  
  return(grad); 
}


 cmplx interpol(cmplx z1,cmplx z2,float m){
 cmplx z;

  /*     
  A   = cmod(z1)   + (cmod(z2)  -cmod(z1)  ) * m; 
  phi = cphase(z1) + (cphase(z2)-cphase(z1))* m;
  phi *= M_PI/180.0;
  z.re = A*cos(phi);
  z.im = A*sin(phi);
  */
 z.re = z1.re + (z2.re-z1.re) * m ;
 z.im = z1.im + (z2.im-z1.im) * m ;
 
  return(z);
}

