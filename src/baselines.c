#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>


float min(float,float);
float max(float,float);
//int strcmp ( const char * , const char * );

int  bresolve(int nargc, char *nargv[],int *nbl, float base[]){
  int i,j,k,m,i1=0,i2=0,l,*ifound;
  int   *ant1,*ant2,nants1,nants2;
  int nants = 30;
  
  for (i = 0; i < nargc; i++){
    if (!strcmp(nargv[i],"-ant")){
      i1 = i;
    } if (!strcmp(nargv[i],"-base")){
      i2=i;
    } 
  }// for i 
  
  if (i1==0) nants1 = nants;
  if (i2==0) nants2 = nants;
  
  if (i1 !=0 && i2 !=0) nants1 = i2-i1-1;
  if (i1 !=0 && i2 ==0) nants1 = nargc-i1-1;
  
  if(i2 !=0) nants2 = nargc -i2-1;

  ant1=(int *)malloc(nants1*sizeof(int));
  ant2=(int *)malloc(nants2*sizeof(int));
  ifound=(int *)malloc(nants*nants*sizeof(int));
  
 
  if (i1==0){
    for(i=0; i < nants1; i++){
      ant1[i] = i+1;
    }
  }else{
    for(i = i1;i < i1+nants1; i++){
      ant1[i-i1] = atoi(nargv[i+1]);
    } 
  }
  
  if(i2==0){
    for(i = 0; i < nants2; i++){
      ant2[i] = i+1;
    }
  }else{
    for(i = i2; i < i2+nants2; i++){
      ant2[i-i2] = atoi(nargv[i+1]);
    }
  }

  for(i=1; i <= nants*nants ; i++)
    ifound[i] = 0;
  
  k = 0;
  for(i=0; i < nants1; i++){
    for(j=0; j < nants2; j++){
      if (ant1[i] != ant2[j]){
	l = (int) (min((float) ant1[i],(float)ant2[j]));
         m = (int) (max((float) ant1[i],(float)ant2[j]));
         if (ifound[m+nants*l] == 0){
	  base[k] = (float) (m+nants*l);
	  k++;
	  ifound[m+nants*l] = 1;
	 }// for if
      }// if ant1 !=ant2 
    }// for j
  }// for i 

  *nbl=k;


   return(0);
}
