#ifndef FLAGCAL_H
#define FLAGCAL_H

// io functions
float gettime(float , float );
int  read_input(int ,char *[]);
int  read_parameters(char []);
int  read_fits(char[] );
int  write_fits(char [],char[]);
int  write_gains();
int  write_flag_info();
// flagging functions
int clip(); 
int init_flag();
int  pre_mad_filter();
int  flag_data();
int flag_apply_all();
int  remove_peaks_time();
int  remove_peaks_freq();
int post_mad_filter(); 
int mad_filter_global();
// indexing/monitoring  function
int  initilise();
int  smooth_data();
int  time_estimate(char[]);
// calibration function
int  gain_compute();
int  setjy();
int  getjy();
int  gain_correct();
int  gain_intrapolate();
int  gain_apply();


void  printerror( int status);
int   group_data(float *,float [],int *, int *,int*);
int   ftime(float, float[],float *);

int init_flag();
int access(const char *, int );
int print_message();
int read_infile(int );
int read_gains();
int write_index();
int string_decode(char s[], int *, int *);

float  max(float , float );
float  min(float , float );
#endif /* FLAGCAL_H */

