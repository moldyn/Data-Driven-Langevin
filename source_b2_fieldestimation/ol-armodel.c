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

//#define RESCALE_DATA

unsigned int WHICHFUTURE=0;
unsigned int VERBOSITY=1;
unsigned long EXCLUDE=0;
double INNOISE=0.0;
char *COLUMN=NULL,*OUTFILE=NULL;


struct param pars={.LENGTH=ULONG_MAX,.DIM=1,.EMB=1,.DELAY=1,.MINN=50,
		   .SIGMA=1.0,.AR_SIZE=0};
char dimset=0;
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
  fprintf(stderr,"\t-p Order of the ar model [default %u]\n",pars.AR_SIZE);
  fprintf(stderr,"\t-F # of component which defines future [not set]\n");
  fprintf(stderr,"\t-%% add initial noise (absolute amplitude) [none]\n");
  fprintf(stderr,"\t-I seed for the rnd-generator (If seed=0, the time\n"
          "\t\tcommand is used to set the seed) [Default: fixed]\n");
  fprintf(stderr,"\t-o output file [default 'datafile-pxxx'.oar]");
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
  if ((out=check_option(in,n,'F','u')) != NULL)
    sscanf(out,"%u",&WHICHFUTURE);
  if ((out=check_option(in,n,'p','u')) != NULL)
    sscanf(out,"%u",&pars.AR_SIZE);
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
  unsigned int outstrlen=0,wdim=1,offset;
  unsigned int *future;
  long li,i,j;
  char *outstring,*infile,*wcol;
  double **series,**dfuture,*sigma,*aver,h;
  FILE *fout;

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

  infile=search_datafile(argc,argv,NULL,VERBOSITY);
  if (infile == NULL) {
    fprintf(stderr,"No input datafile found.\n");
    exit(LANGEVIN_MAIN_NO_INPUTFILE);
  }
  
  if (OUTFILE == NULL) {
    check_alloc(OUTFILE=(char*)calloc(strlen(infile)+11,(size_t)1));
    if (pars.AR_SIZE < 10)
      sprintf(OUTFILE,"%s-p00%u.oar",infile,pars.AR_SIZE);
    else {
      if (pars.AR_SIZE < 100)
	sprintf(OUTFILE,"%s-p0%u.oar",infile,pars.AR_SIZE);
      else
	sprintf(OUTFILE,"%s-p%u.oar",infile,pars.AR_SIZE);
    }
  }
  OUTFILE=test_outfile(OUTFILE);

  if (COLUMN == NULL)
    series=(double**)get_multi_series(infile,&(pars.LENGTH),EXCLUDE,
				      &(pars.DIM),"",dimset,VERBOSITY);
  else
    series=(double**)get_multi_series(infile,&(pars.LENGTH),EXCLUDE,
				      &(pars.DIM),COLUMN,dimset,VERBOSITY);


  if (INNOISE > 0.0) {
    rnd_init(seed);
    for (i=0;i<10000;i++)
      rnd_long();
    for (li=0;li<pars.LENGTH;li++)
      for (i=0;i<pars.DIM;i++)
	series[i][li] += 2.0*INNOISE*
	  ((double)rnd_long()/(double)ULONG_MAX-0.5);
  }

  check_alloc(future=(unsigned int*)malloc(sizeof(int)*pars.LENGTH));
  if (WHICHFUTURE>0) {
    check_alloc(wcol=(char*)calloc((size_t)10,(size_t)1));
    sprintf(wcol,"%u",WHICHFUTURE);
    dfuture=(double**)get_multi_series(infile,&(pars.LENGTH),EXCLUDE,
				       &wdim,wcol,1,0L);
    for (li=0;li<pars.LENGTH;li++) {
      future[li] = (unsigned int)(dfuture[0][li]>0.0);
    }
    free(dfuture);
    future[pars.LENGTH-1]=0;
  }
  else {
    for (li=0;li<pars.LENGTH-1;li++)
      future[li]=1;
    future[pars.LENGTH-1]=0;
  }

  check_alloc(pars.AR_MOD=(double**)malloc(sizeof(double*)*pars.DIM));
  for (i=0;i<pars.DIM;i++)
    check_alloc(pars.AR_MOD[i]=(double*)malloc(sizeof(double)*pars.AR_SIZE));
  check_alloc(sigma=(double*)malloc(sizeof(double)*pars.DIM));
  check_alloc(aver=(double*)malloc(sizeof(double)*pars.DIM));

  for (i=0;i<pars.DIM;i++) {
    aver[i]=sigma[i]=0.0;
    for (j=0;j<pars.LENGTH;j++) {
      h=series[i][j];
      aver[i] += h;
      sigma[i] += h*h;
    }
    aver[i] /= (double)pars.LENGTH;
    sigma[i]=sqrt(sigma[i]/((double)pars.LENGTH-1.0)-aver[i]*aver[i]);
  }
#if defined(RESCALE_DATA)
  for (i=0;i<pars.DIM;i++) {
    for (j=0;j<pars.LENGTH;j++)
      series[i][j]=(series[i][j]-aver[i])/sigma[i];
  }
#endif

  if (pars.AR_SIZE > 0)
    make_corr(series,future,pars,sigma);

  fout=fopen(OUTFILE,"w");
  fprintf(fout,"%u\n",pars.AR_SIZE);
  for (i=0;i<pars.DIM;i++)
    fprintf(fout,"%lf ",sigma[i]);
  fprintf(fout,"\n");
  for (j=0;j<pars.AR_SIZE;j++) {
    for (i=0;i<pars.DIM;i++)
      fprintf(fout,"%.10e ",pars.AR_MOD[i][j]);
    fprintf(fout,"\n");
  }
  fclose(fout);

  return 0;
}
