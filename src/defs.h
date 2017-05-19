#ifndef _DEFS_H
#define _DEFS_H
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<omp.h>
#include <ctype.h>
#include <sys/time.h>
#include<sys/stat.h>
#include <sys/unistd.h>
#include<complex.h>

#define nants 30
#define nbaselines  435
#define nl 8
#define max_len  200 
#define MAX_TIMES 5000
#define MAX_SCANS 40
#define MAX_SRC   128 
#define MAX_BAD_TIMES 100

typedef struct bad_times{
  double start;
  double end;
}BadTimes;

char err_msg[1000]; 

/* parameters related to input fits file */

int  max_times; 
int  nsrc;    // number of sources in the multi-source file 
int  nhdu;    //   ''      HDUs    ''         ''            
int  ntimes;  // number of time samples  in the input FITS file  
int  nscans;  //   ''        scans               ''              


long nchans;  //   ''        chans               ''               
long nstokes; //   ''        stokes              ''              
long ncmplx ; //   ''         cmplex             ''                
long gcount,gcount2 ; //   ''        visbilities         ''              
long pcount ; //   ''        random group para   ''              
long nbit,naxis ; // standard FITS keyword                            



int  nthreads;    // number of openMp threads  
long nstokes_half;// half of the number of stokes 
long ggcount[MAX_SRC]; // number of visibilties for a scan 
float crval4,freq,ra,dec; // FITS keywords
float *pscal,*pzero;
/* Parameter for switching a module on or of */
int run_mode,clip_mode,gmad_mode,pmad_mode,vsr_mode,aflag_mode;
int bflag_mode,cflag_mode,qmad_mode,rfi_mode_t,rfi_mode_f; 
int range_t_mode;
float thrs_clip,thrs_gmad,thrs_pmad,thrs_vsr,thrs_ant;
float thrs_base,thrs_chan,thrs_comp,thrs1_comp;
float thrs_qmad,thrs_qmadlow,thrs_qmadmean,thrs_qmadb,thrs_qmadblow,thrs_qmadbmean;
float thrs_peak_t,thrs_peak_tlow,thrs_peak_f,thrs_peak_flow; 
float dthrs_peak_t,dthrs_peak_f;
float thrs_block,thrs_ngood;
int  nriter_t,nriter_f, win_min_t,win_min_f,dwin_t,dwin_f;
float thrs_range_t,thrs_low_range_t,thrs_mean_range_t;
int   win_min_rt,dwin_rt,nriter_rt;
float dthrs_rt;
int *tgrp,*bgrp;

/* parameters for output fits file  */
int echan;        // ending channel for output
int bchan;        // starting channel for output
int imode; // options for running the code
int src_out;      // output source id
/* parameters used for indexing the visibility data */
long   tid[MAX_SCANS]; // starting time for a scan 
long   gid[MAX_SCANS]; // starting group for a scan 
long   sid[MAX_SCANS]; // source id of scan 
long   idmap[MAX_SRC]; // map from id in scan to src no in src struct (sors)
                         // jnc 5/sep/11
int    id[nbaselines];  // id for a baseline 
int    dc[nants*nants];  // inver id for a baselin 
int    tpnt[MAX_TIMES];// number of time samples in a scan 
float  time_data[MAX_TIMES]; // time for a time sample 
int nablocks; // number of bolcks for computing gains (averaging blocks)
int nsblocks; // number of blocks for flagging        (flagging blocks)
int max_nab;  // number of time samples in the longest averaging block
int *nab;     // number of averaging blocks in a scan 
int *last_abl;// number of time samples in the last averaging block
int *abid_in; // inverse id of an averaging block  
int *ab_start;// starting time sample for an averaging block  
int *ab_end;  // ending time samples  for an averaging block   
float *abmt; // average time for an averaging block 
int max_nsb;  // number of time samples in the  longest flagging block
int *nsb;     // number of flagging blocks in a scan 
int *last_sbl;// number of time samples in the last flagging block in a scan 
int *sbid_in; // inverse id of a flagging  block  
int *sb_start;// starting time sample for a flagging block  
int *sb_end;  // ending time samples  for an flagging block   
/* parameters uded for flagging */ 
int    niter; // number of iterations for computing solutions
int    nflag;  // number of time samples in a flagging  block 
int    nsol; // number of time samples in an averaging block 
int    ref_ant;    // ference antenna for calibration 
float *bad_ant;
BadTimes badtimes[MAX_BAD_TIMES];
int    nbadtimes,check_badtimes_only;
float  dthrs_peak; // factor by which threshold is lowered 
int *ba,*bb,*bc,*bs,*gs,*rc,*lc;

float  eps;        // a small number 
float  alpha;      // parameter use for computing gain  
int    nriter;     // number of iterations for removing rfi peaks
int    w_smooth;   // shortest smoothing window for removing rfi peaks 
int    dw_smooth;  // factor used to increare smoothing window

cmplx3 *Visib,*V1; // input & average visibilties 
cmplx *g1;        //  average gain 

char infits[max_len],outfits[max_len],parafile[max_len],srcfile[max_len],flgfile[max_len],calfile[max_len];
char gainfile[max_len];

float  *vsr;       // an indicator for the badness of a flagging block
float  *bad_base;  // an indicator for the badness of a baseline
float  *bad_chans; // an indicator for the badness of a channel
int    *dant;      // an indicator for the overall (all scans) badness of an anteena

typedef struct source{char name[30]; int id; int type; float flux; float ra; float dec;} sors; 
typedef struct weight{char val[10];} weight; 
typedef struct sarray{char  val[32];}sar;


sar srtype[3],pname[30];
sors src[MAX_SRC];


#endif /*_defs_H */
