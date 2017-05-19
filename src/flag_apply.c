/* -------------------------------------------------------------------
   - This program applies flagging on the basis of vsr, bad_ant etc., which 
   are computed by the program flag_data.
   - Flagging for the original data and average data (used for computing gains)
   is done separately (please check the mapping in both the case) 
   - Note that flagging is done only for the first two stokes (if there are more than two)
   RR & LL and is exported to RL & LR on the basis of that.
   
   
   - Jayanti Prasad
   Sat Jul 24 10:09:59 IST 2010
   ---------------------------------------------------------------------*/
#include  <defs.h>
#include <flagcal.h>

int flag_apply_all(){
  int j,l,m,n,n1,i1,i2,i3,j1,k1,l1,l2,ll[nstokes],npflag[5]; 
  
  for(i1 = 0; i1 < 5; i1++)
    npflag[i1] = 0;
  
  for(n1=0,i1 = 0; i1 < nscans; i1++){
    if (gs[i1] > 0){
      for(i2 = 0; i2 < nsb[i1]; i2++){ 
	j  = sbid_in[i2+max_nsb*i1];
	for(i3 = sb_start[j]; i3 < sb_end[j]; i3++){
	  for(j1 = 0; j1 < nchans; j1++){
	    for(k1 = 0; k1 < nstokes_half; k1++){
	      l1   =  k1 + nstokes * (j1+nchans*j);
	      l2   =  k1 + nstokes * (j1+nchans*i3); 
	      for(n = 0; n < nbaselines; n++){
		l = (id[n]+1)/nants;
		m = (id[n]+1) - nants * l;
		if (src[idmap[sid[i1]]].type == 0 || src[idmap[sid[i1]]].type == 1){
		  // VSR flagging is only for calibrators 
		  if(vsr_mode ==1 && Visib[n+nbaselines*l2].w > 0.0){
		    if(vsr[n+nbaselines*l1] < thrs_vsr){ 
		      Visib[n+nbaselines*l2].w = -1.0;   
		      npflag[0]++; 
		    }// for if
		  }// for if
		} if (aflag_mode == 1  && Visib[n+nbaselines*l2].w > 0.0){
		  if ((bad_ant[l+nants*i1] > thrs_ant) || (bad_ant[m+nants*i1] > thrs_ant)){
		  Visib[n+nbaselines*l2].w  = -1.0;
		  npflag[1]++;
		  }// for if 
		}if(bflag_mode==1 && Visib[n+nbaselines*l2].w > 0.0) {
		  if (bad_base[n+nbaselines*i1] > thrs_base ){
		    Visib[n+nbaselines*l2].w  = -1.0;
		    npflag[2]++;
		  }//for if
		}if (cflag_mode ==1 && Visib[n+nbaselines*l2].w > 0.0){ 
		  if (bad_chans[j1+nchans*i1] > thrs_chan){ 
		    Visib[n+nbaselines*l2].w  = -1.0;
		  npflag[3]++;
		  }// for if
		}// for if
		n1++; 
	      } // for n
	    } // for k1
	  }// for j1
	} // for i2
      }// for i3
    }//for if
  } // for i1
  
  if (nstokes_half !=  nstokes){
    for(i1 = 0; i1 < nscans; i1++){
      if (gs[i1] > 0){
	for(i2 = 0; i2 < nsb[i1]; i2++){ 
	  j  = sbid_in[i2+max_nsb*i1];
	  for(i3 = sb_start[j]; i3 < sb_end[j]; i3++){
	    for(j1 = 0; j1 < nchans; j1++){
	      for(k1 = 0; k1 < nstokes; k1++)
		ll[k1] =  k1 + nstokes * (j1+nchans*i3);
	      for(n = 0; n < nbaselines; n++){
		if(Visib[n+nbaselines*ll[0]].w  < 0.0 ||  Visib[n+nbaselines*ll[1]].w  < 0.0){  
		  for(k1 = nstokes_half; k1 < nstokes; k1++)
		    Visib[n+nbaselines*ll[k1]].w =  -1.0;
		}// for if    
	      } // for n
	    } // for k1
	  }// for j1
	} // for i2
      }// for i3
    }// for if
  } // for i1
  
  fprintf(stdout,"\t Flagging [org] on the basis of VSR      = %6.6f\n",(float)npflag[0]/n1);
  fprintf(stdout,"\t Flagging [org] on the basis of bad ANT  = %6.3f\n",(float)npflag[1]/n1);
  fprintf(stdout,"\t Flagging [org] on the basis of bad BASE = %6.3f\n",(float)npflag[2]/n1);
  fprintf(stdout,"\t Flagging [org] on the basis of bad CHAN = %6.3f\n",(float)npflag[3]/n1);
  
  return(0);
  
}  

