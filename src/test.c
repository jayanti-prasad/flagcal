#include <defs.h>
#include<flagcal.h>
#include <sys/unistd.h>
#include<string.h>

int main(){
  int i,n1=23,n2=28;
  FILE *inpfile;
  char string1[max_len],c[max_len], *p,*p1;
  int  iry[n1];
  float fry[n2];
  struct sarray   switch_val[2];

  char *ikwrd[23]={"nsol","niter","ref_ant","clip_mode","gmad_mode","pmad_mode","nflag","aflag_mode","bflag_mode","cflag_mode","qmad_mode","rfi_mode_t","nriter_t","win_min_t","dwin_t","rfi_mode_f","nriter_f","win_min_f","dwin_f","range_t_mode","win_min_rt","dwin_rt","nriter_rt"};
  char *fkwrd[28]={"alpha","eps","thrs_clip","thrs_gmad","thrs_pmad","thrs_vsr","thrs_ngood","thrs_ant","thrs_base","thrs_chan","thrs_qmad","thrs_low_qmad","thrs_mean_qmad","thrs_b_qmad","thrs_b_low_qmad","thrs_b_mean_qmad","thrs_comp","thrs1_comp","thrs_peak_t","thrs_low_peak_t","dpt","thrs_peak_f","thrs_low_peak_f","dpf","thrs_range_t","thrs_low_range_t","thrs_mean_range_t","dprt"};

  sprintf(switch_val[0].val,"OFF");
  sprintf(switch_val[1].val,"ON");

  for(i=0; i < n1; i++)
    printf("%s \n",ikwrd[i]);

  for(i=0; i < n2; i++)
    printf("%s \n",fkwrd[i]);


}
