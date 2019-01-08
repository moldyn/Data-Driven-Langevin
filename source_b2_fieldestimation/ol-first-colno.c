/*
 *   This file is part of OLANGEVIN
 *
 *   Copyright (c) 2013 Rainer Hegger
 *
 *   OLANGEVIN is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   OLANGEVIN is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with OLANGEVIN; if not, see <http://www.gnu.org/licenses/>.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <time.h>
#include <math.h>
#include <ctype.h>
#include "routines/olang.h"

unsigned int WHICHFUTURE=0;
unsigned int VERBOSITY=1;
unsigned long FLENGTH=1000,EXCLUDE=0;
double INNOISE=0.0;
char *COLUMN=NULL,*OUTFILE=NULL,*arfile=NULL;


struct param pars={.LENGTH=ULONG_MAX,.DIM=1,.EMB=1,.DELAY=1,
		   .MINN=50,.SIGMA=1.0,.AR_SIZE=0,.maxr=0.0};
char dimset=0,nofields=0;
unsigned long seed=0x4973592378L;

void show_options(char *progname)
{

  fprintf(stderr," Usage: %s [Options]\n",progname);
  fprintf(stderr," Options:\n");
  fprintf(stderr,"Everything not being a valid option will be interpreted"
          " as a possible"
          " datafile.\nIf no datafile is given stdin is read. Just - also"
          " means stdin\n");
  fprintf(stderr,"\t-l # of data to be used [default whole file]\n");
  fprintf(stderr,"\t-x # of lines to be ignored [default %lu]\n",EXCLUDE);
  fprintf(stderr,"\t-c column [default 1,...,# of components]\n");
  fprintf(stderr,"\t-m # of components [default %u]\n",pars.DIM);
  fprintf(stderr,"\t-M Embedding dimension [default 1 (no embedding)]\n");
  fprintf(stderr,"\t-d delay for the embedding [default %u]\n",pars.DELAY);
  fprintf(stderr,"\t-k # of neighbors  [default %u]\n",pars.MINN);
  fprintf(stderr,"\t-F # of component which defines future [not set]\n");
  fprintf(stderr,"\t-L # of iterations [default %lu]\n",FLENGTH);
#if defined(WEIGHTS)
  fprintf(stderr,"\t-r radius for the weights [default not set]\n");
#endif
  fprintf(stderr,"\t-N dont write field-file [not set]\n");
  fprintf(stderr,"\t-A name of ar-model file [default none]\n");
  fprintf(stderr,"\t-%% add initial noise (absolute amplitude) [none]\n");
  fprintf(stderr,"\t-I seed for the rnd-generator (If seed=0, the time\n"
          "\t\tcommand is used to set the seed) [Default: fixed]\n");
  fprintf(stderr,"\t-o output file [default 'datafile'.lang;"
	  " no -o means write to stdout]\n");
  fprintf(stderr,"\t-V verbosity level [default 1]\n\t\t"
          "0='only panic messages'\n\t\t"
          "1='+ input/output messages'\n");
  fprintf(stderr,"\t-h  show these options\n");
  exit(0);
}

void scan_options(int n,char **in)
{
  char *out;

  if ((out=check_option(in,n,'l','u')) != NULL)
    sscanf(out,"%lu",&(pars.LENGTH));
  if ((out=check_option(in,n,'x','u')) != NULL)
    sscanf(out,"%lu",&EXCLUDE);
  if ((out=check_option(in,n,'c','s')) != NULL) {
    COLUMN=out;
    dimset=1;
  }
  if ((out=check_option(in,n,'m','u')) != NULL)
    sscanf(out,"%u",&(pars.DIM));
  if ((out=check_option(in,n,'M','u')) != NULL)
    sscanf(out,"%u",&(pars.EMB));
  if ((out=check_option(in,n,'F','u')) != NULL)
    sscanf(out,"%u",&WHICHFUTURE);
  if ((out=check_option(in,n,'d','u')) != NULL)
    sscanf(out,"%u",&(pars.DELAY));
  if ((out=check_option(in,n,'L','u')) != NULL)
    sscanf(out,"%lu",&FLENGTH);
#if defined(WEIGHTS)
  if ((out=check_option(in,n,'r','f')) != NULL)
    sscanf(out,"%lf",&(pars.maxr));
#endif
  if ((out=check_option(in,n,'N','n')) != NULL)
    nofields=1;
  if ((out=check_option(in,n,'k','u')) != NULL)
    sscanf(out,"%u",&(pars.MINN));
  if ((out=check_option(in,n,'A','s')) != NULL)
    arfile=out;
  if ((out=check_option(in,n,'I','u')) != NULL) {
    sscanf(out,"%lu",&seed);
    if (seed == 0)
      seed=(unsigned long)time((time_t*)&seed);
  }
  if ((out=check_option(in,n,'V','u')) != NULL)
    sscanf(out,"%u",&VERBOSITY);
  if ((out=check_option(in,n,'%','f')) != NULL)
    sscanf(out,"%lf",&INNOISE);
  if ((out=check_option(in,n,'o','o')) != NULL)
    if (strlen(out) > 0)
      OUTFILE=out;
}


int main(int argc,char **argv)
{
  unsigned int outstrlen=0,offset,hdim,wdim=1;
  int i,j;
  unsigned int *future,has_future;
  unsigned long count;
  long li,lj;
  char *outstring,*infile,*fieldfile,*wcol;
  double **series,*sermin,*serinter,**dfuture,max,**cast,sweights;
  double *new,*drift,**diffusion;
  struct sfound sf;
  FILE *fout,*ffield;

  if (scan_help(argc,argv))
    show_options(argv[0]);

  for (i=0;i<argc;i++) {
    outstrlen += strlen(argv[i]);
    outstrlen++;
  }
  check_alloc(outstring=(char*)malloc(sizeof(char)*(outstrlen+8)));
  sprintf(outstring,"#Prog: ");
  offset=7;
  for (i=0;i<argc;i++) {
    sprintf(outstring+offset,"%s ",argv[i]);
    offset += (strlen(argv[i])+1);
  }

  scan_options(argc,argv);
  if (arfile == NULL) {
    exit(READ_AR_FILE_NO_NAME);
  }
  read_ar_file(&pars,arfile);
  init_noise(pars,seed);

  if (pars.MINN < (pars.DIM+1)) {
    fprintf(stderr,"Too few neighbors. System is underdetermined.\n");
    exit(LANGEVIN_MAIN_TOO_FEW_MINN);
  }

  infile=search_datafile(argc,argv,NULL,VERBOSITY);
  if (infile == NULL) {
    fprintf(stderr,"No input datafile found.\n");
    exit(LANGEVIN_MAIN_NO_INPUTFILE);
  }
  
  if (OUTFILE == NULL) {
    check_alloc(OUTFILE=(char*)calloc(strlen(infile)+6,(size_t)1));
    sprintf(OUTFILE,"%s.lang",infile);
    check_alloc(fieldfile=(char*)calloc(strlen(infile)+12,(size_t)1));
    sprintf(fieldfile,"%s.lang.field",infile);
  }
  else {
    check_alloc(fieldfile=(char*)calloc(strlen(OUTFILE)+7,(size_t)1));
    sprintf(fieldfile,"%s.field",OUTFILE);
  }
  OUTFILE=test_outfile(OUTFILE);
  if (!nofields)
    fieldfile=test_outfile(fieldfile);

  if (COLUMN == NULL)
    series=(double**)get_multi_series(infile,&(pars.LENGTH),EXCLUDE,
				      &(pars.DIM),"",dimset,VERBOSITY);
  else
    series=(double**)get_multi_series(infile,&(pars.LENGTH),EXCLUDE,
				      &(pars.DIM),COLUMN,dimset,VERBOSITY);

  if (INNOISE > 0.0)
    for (li=0;li<pars.LENGTH;li++)
      for (i=0;i<pars.DIM;i++)
	series[i][li] += 2.0*INNOISE*
	  ((double)rnd_long()/(double)ULONG_MAX-0.5);

  check_alloc(sermin=malloc(sizeof(double)*pars.DIM));
  check_alloc(serinter=malloc(sizeof(double)*pars.DIM));
  rescale_data(series,pars,sermin,serinter);

  max=serinter[0];
  for (i=1;i<pars.DIM;i++)
    if (serinter[i] > max)
      max=serinter[i];

  for (i=0;i<pars.DIM;i++) {
    for (j=0;j<pars.LENGTH;j++)
      series[i][j] /= max;
  }
  if (pars.maxr > 0.0)
    pars.maxr /= max;

  hdim=(pars.EMB-1)*pars.DELAY;
  pars.hdim=hdim;

  check_alloc(future=(unsigned int*)malloc(sizeof(int)*pars.LENGTH));
  if (WHICHFUTURE>0) {
    check_alloc(wcol=(char*)calloc((size_t)10,(size_t)1));
    sprintf(wcol,"%u",WHICHFUTURE);
    dfuture=(double**)get_multi_series(infile,&(pars.LENGTH),EXCLUDE,
				       &wdim,wcol,1,0L);
    for (li=hdim;li<pars.LENGTH-1;li++) {
      has_future= (dfuture[0][li]>0.0);
      for (lj=1;lj<=hdim;lj++)
	has_future &= (dfuture[0][li-lj]>0.0);
      future[li]=has_future;
    }
    future[pars.LENGTH-1]=0;
    free(dfuture);
  }
  else {
    for (li=0;li<pars.LENGTH-1;li++)
      future[li]=1;
    future[pars.LENGTH-1]=0;
  }

  check_alloc(sf.found=(unsigned long*)
	      malloc(sizeof(unsigned long)*pars.LENGTH));
  check_alloc(sf.distance=(double*)malloc(sizeof(double)*pars.LENGTH));
  check_alloc(sf.weight=(double*)malloc(sizeof(double)*pars.LENGTH));
  check_alloc(sf.count=(unsigned long*)malloc(sizeof(unsigned long)));
  check_alloc(sf.aveps=(double*)malloc(sizeof(double)));
  sf.count[0]=0;
  sf.aveps[0]=0.0;
  check_alloc(cast=(double**)malloc(sizeof(double*)*pars.DIM));
  for (i=0;i<pars.DIM;i++)
    check_alloc(cast[i]=(double*)malloc(sizeof(double)*(hdim+1)));
  check_alloc(new=(double*)malloc(sizeof(double)*pars.DIM));

  check_alloc(drift=(double*)malloc(sizeof(double)*pars.DIM));
  check_alloc(diffusion=(double**)malloc(sizeof(double*)*pars.DIM));
  for (i=0;i<pars.DIM;i++) {
    check_alloc(diffusion[i]=(double*)malloc(sizeof(double)*pars.DIM));
  }

  for (i=0;i<=hdim;i++)
    for (j=0;j<pars.DIM;j++)
      cast[j][hdim-i]=series[j][pars.LENGTH-1-i];

  pars.minminn=(unsigned long)((pars.DIM+3.0)/2.0+0.5);
  init_neighbor_search(series,pars,future);

  /* print program call and column labels */
  fout=fopen(OUTFILE,"w");
  fprintf(fout,"%s\n",outstring);
#if defined(MAXNORM)
  fprintf(fout,"#Norm: Maxnorm\n");
#else
  fprintf(fout,"#Norm: L2 Norm\n");
#endif
#if defined(WEIGHTS)
  if (pars.maxr > 0.0)
    fprintf(fout,"#Weights: on WFACT = %lf maxr = %lf\n",WFACT,pars.maxr*max);
  else
    fprintf(fout,"#Weights: on WFACT = %lf maxr not set\n",WFACT);
#else
  fprintf(fout,"#Weights: off\n");
#endif
  fprintf(fout,"#AR_SIZE: %u\n",pars.AR_SIZE);
  fprintf(fout,"#AR_SIGMAS: ");
  for (i=0;i<pars.DIM;i++)
    fprintf(fout,"%e ",pars.AR_SIGMA[i]);
  fprintf(fout,"\n");
  if (pars.AR_SIZE > 0) {
    for (i=0;i<pars.AR_SIZE;i++) {
      fprintf(fout,"#AR_MOD_%u: ",i+1);
      for (j=0;j<pars.DIM;j++)
	fprintf(fout,"%e ",pars.AR_MOD[j][i]);
      fprintf(fout,"\n");
    }
  }
  fprintf(fout,"#Content: ");
  for (i=0;i<pars.DIM;i++)
    fprintf(fout,"x%d ",i+1);
  fprintf(fout,"\n");
  fflush(fout);
  if (!nofields) {
    ffield=fopen(fieldfile,"w");
    fprintf(ffield,"%s\n",outstring);
#if defined(WEIGHTS)
  if (pars.maxr > 0.0)
    fprintf(ffield,"#Weights: on WFACT = %lf maxr = %lf\n",WFACT,pars.maxr*max);
  else
    fprintf(ffield,"#Weights: on WFACT = %lf maxr not set\n",WFACT);
#else
    fprintf(ffield,"#Norm: L2 Norm\n");
#endif
#if defined(WEIGHTS)
    fprintf(ffield,"#Weights: on WFACT = %lf\n",WFACT);
#else
    fprintf(ffield,"#Weights: off\n");
#endif
    fprintf(ffield,"#AR_SIZE: %u\n",pars.AR_SIZE);
    fprintf(ffield,"#AR_SIGMAS: ");
    for (i=0;i<pars.DIM;i++)
      fprintf(ffield,"%e ",pars.AR_SIGMA[i]);
    fprintf(ffield,"\n");
    if (pars.AR_SIZE > 0) {
      for (i=0;i<pars.AR_SIZE;i++) {
	fprintf(ffield,"#AR_MOD_%u: ",i+1);
	for (j=0;j<pars.DIM;j++)
	  fprintf(ffield,"%e ",pars.AR_MOD[j][i]);
	fprintf(ffield,"\n");
      }
    }
    fprintf(ffield,"#Content: ");
    for (i=0;i<pars.DIM;i++)
      for (j=0;j<pars.EMB;j++)
	if (pars.EMB>1) {
	  fprintf(ffield,"x%d_%d ",i+1,j+1);
	}
	else {
	  fprintf(ffield,"x%d ",i+1);
	}
    for (i=0;i<pars.DIM;i++)
      fprintf(ffield,"h%d ",i+1);
    for (i=0;i<pars.DIM;i++)
      for (j=0;j<=i;j++)
	fprintf(ffield,"D_%d_%d ",i+1,j+1);
#if defined(WEIGHTS)
    fprintf(ffield,"distance sweights\n");
#else
    fprintf(ffield,"distance\n");
#endif
    fflush(ffield);
  }
  count=0;
  while (count <FLENGTH) {
    search_neighbors(series,cast,pars,sf);
    get_fields(series,pars,sf,drift,diffusion);
    make_cast(cast,new,pars,drift,diffusion);

    for (i=0;i<pars.DIM;i++)
      fprintf(fout,"%lf ",new[i]*max+sermin[i]);
    fprintf(fout,"\n");
    fflush(fout);
    if (!nofields) {
      for (i=0;i<pars.DIM;i++)
	for (j=0;j<pars.EMB;j++)
	  fprintf(ffield,"%e ",cast[i][hdim-j*pars.DELAY]*max+sermin[i]);
      for (i=0;i<pars.DIM;i++)
	fprintf(ffield,"%e ",(drift[i]-cast[i][hdim])*max);
      for (i=0;i<pars.DIM;i++)
	for (j=0;j<=i;j++)
	  fprintf(ffield,"%e ",diffusion[i][j]*max);
#if defined(WEIGHTS)
      sweights=sf.weight[0];
      for (i=1;i<pars.MINN;i++)
	sweights += sf.weight[i];
      fprintf(ffield,"%e %lf\n",sf.distance[pars.MINN-1]*max,sweights);
#else
      fprintf(ffield,"%e\n",sf.distance[pars.MINN-1]*max);
#endif
      fflush(ffield);
    }
    
    for (i=1;i<=hdim;i++)
      for (j=0;j<pars.DIM;j++)
	cast[j][i-1]=cast[j][i];
    for (i=0;i<pars.DIM;i++)
      cast[i][hdim]=new[i];
    count++;
  }
  fclose(fout);
  if (!nofields) 
    fclose(ffield);

  return 0;
}
