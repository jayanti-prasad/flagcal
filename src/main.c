#include <defs.h>
#include <flagcal.h>

/*--------------------------------------------------------------
 
   This is the main program which does the following things:
  
  - reads a multi-source multi-channel FITS file (created by GMRT gvfits)
  - flags the visibility data
  - computes antenna based complex gains 
  - calibrates the visibility data
  - writes  output in a user defined fits file (single/multi channel, single/multi source)

   Check README for detail and defs.h for the definition of the variables.
  
  Jayanti Parasd 

  Tue Jun  1 10:45:53 IST 2010
 ------------------------------------------------------------- */

struct timeval  t1,t2;
struct timezone tzp;
char tfname[100];
FILE *tf; 

int main(int argc, char *argv[])
{
  
  gettimeofday (&t1, &tzp);      /* starting time                   */
  
  if (read_input(argc,argv))     /*  Read command line options  */
    return(0);
 
  
  sprintf(tfname,"cpu_%d.txt",nthreads); 
  
  tf = fopen(tfname,"w");  
  
  printf("\t time estimate will be written in %s\n",tfname);
  
  if (read_parameters(parafile)) /* Read parameters from a file */
    return(0); 
     
  check_badtimes_only=1;
  nbadtimes=0;
  if(access(flgfile, F_OK) !=-1)
    read_infile(2); // read flag file so that one can flag bad time ranges jnc 13/sep/11
  check_badtimes_only=0;

  if(read_fits(infits))          /* Read multi-source FITS file */
    return(0); 
   
  if(initilise())                /* indexing, memory allocation &  initaliaztion */
    return(0);	 
  
  if(init_flag())                 /* Initial flagging (Input) */ 
    return(0); 
  
   time_estimate("Reading and Initilization");               /* keep track of time */
   
  if(clip())                     /* clip very high values using MAX_VAL (see defs.h) */
    return(0); 
  
  if(mad_filter_global())        /*  This is the first mad filter  */
    return(0);
  
  time_estimate("Global Mad Filetering"); 
  
  if(setjy())                    /* Set the flux for flux calibrators */
    return(0);
  
  if(pre_mad_filter())           /* Remove high peaks from every scan  */
    return(0);
  
  time_estimate("Pre Mad Filtering"); 
  
  if(flag_data())                /* Flagging on the basis of VSR & A, B, C */
    return(0);

  write_flag_info(); 
  
  if(flag_apply_all())           /* Flagging on the basis of VSR & A, B, C */
  return(0);

  if (imode)                      /* If one  flagging information is needed */
    return(0);
  
  if(smooth_data())               /* Smooth the data for gain computation */
    return(0);  
 
  time_estimate("Smoothing data"); 
  
  if(gain_compute())              /* Call the main gain computing module */
    return(0); 
  
  if(getjy())                     /* compute the flux of phase cal */
    return(0); 
  
  if(gain_correct())              /* normalize the phase of gains wrt ref_ant */
    return(0);
  
  if(gain_intrapolate())          /* interpolate the gains for source */
    return(0);  
  
  if(gain_apply())                /* Apply the gains  */
    return(0);
  
  if(write_gains())               /* write gain file */
    return(0);

  time_estimate("Gain computation"); 
  
  if(remove_peaks_time())         /* Remove RFI peaks from time */
    return(0); 
  
  if(time_estimate("Removing peak in time"))
    return(0);
  
  if(remove_peaks_freq())         /* Remove RFI peaks from frequency */
    return(0);
  
  time_estimate("Removing peak in frequency"); 
  
  if(post_mad_filter())           /* Post mad filtering */
    return(0);
  
  time_estimate("Post mad filtering"); 
  

  write_fits(infits,outfits);     /* Write output FITS file */
  
  fprintf(stdout,"\t  everything DONE !\n");
  
  time_estimate("Writing output"); 
 
  fclose(tf); 
  
  return(0);
  
} // main programs ends 



int time_estimate(char text[]){
  float t; 
  gettimeofday (&t2, &tzp);
  t = (t2.tv_sec-t1.tv_sec) + (t2.tv_usec-t1.tv_usec)/1000000.0;
  fprintf(stdout,"\t %s time = %2.2f sec \n",text,t);
  fprintf(tf,"\t %s : time = %2.2f sec \n",text,t);
    
  t1 = t2; 
  return(0);
}

