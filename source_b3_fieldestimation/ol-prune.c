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

unsigned int WHICHFUTURE=0,compdim=1;
unsigned int VERBOSITY=1,*PART,DEPTH=1,ODEPTH=1,MAXOUT=1;
unsigned long EXCLUDE=0;
char *COLUMN=NULL,*OUTFILE=NULL,*SPART,compdimset=0,scrambleset=0;


struct param pars={ULONG_MAX,1,1,1,50,1.0,.AR_SIZE=0};
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
  fprintf(stderr,"\t-C # of components to compare [default -m]\n");
  fprintf(stderr,"\t-F # of component which defines future [not set]\n");
  fprintf(stderr,"\t-P # of partitions [default 2,2,2,... (C-times)]\n");
  fprintf(stderr,"\t-D pattern depth [default %u]\n",ODEPTH);
  fprintf(stderr,"\t-M max pattern out per box [default %u]\n",MAXOUT);
  fprintf(stderr,"\t-s scramble output patterns [default no]\n");
  fprintf(stderr,"\t-o output file [default 'datafile'.pru;"
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
  if ((out=check_option(in,n,'C','u')) != NULL) {
    sscanf(out,"%u",&compdim);
    compdimset=1;
  }
  if ((out=check_option(in,n,'F','u')) != NULL)
    sscanf(out,"%u",&WHICHFUTURE);
  if ((out=check_option(in,n,'P','s')) != NULL)
    SPART=out;
  if ((out=check_option(in,n,'D','u')) != NULL)
    sscanf(out,"%u",&ODEPTH);
  if ((out=check_option(in,n,'M','u')) != NULL)
    sscanf(out,"%u",&MAXOUT);
  if ((out=check_option(in,n,'s','n')) != NULL)
    scrambleset=1;
  if ((out=check_option(in,n,'V','u')) != NULL)
    sscanf(out,"%u",&VERBOSITY);
  if ((out=check_option(in,n,'o','o')) != NULL)
    if (strlen(out) > 0)
      OUTFILE=out;
}


int main(int argc,char **argv)
{
  unsigned int outstrlen=0;
  int i,j,k,n;
  unsigned int offset,wdim=1,alldim;
  unsigned int has_future,*future,**useries,h;
  unsigned long count;
  long *segment;
  char *outstring,*infile,*wcol,*use_array;
  double **series,*sermin,*serinter,**dfuture;
  ptree *root;
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
 
  if (!compdimset) 
    compdim=pars.DIM;
  else
    if (compdim > pars.DIM) {
      fprintf(stderr,"-C value larger -m Value.\n");
      exit(LANGEVIN_MAIN_C_TOO_LARGE);
    }

  infile=search_datafile(argc,argv,NULL,VERBOSITY);
  if (infile == NULL) {
    fprintf(stderr,"No input datafile found.\n");
    exit(LANGEVIN_MAIN_NO_INPUTFILE);
  }
  
  if (OUTFILE == NULL) {
    check_alloc(OUTFILE=(char*)calloc(strlen(infile)+5,(size_t)1));
    sprintf(OUTFILE,"%s.pru",infile);
  }
  OUTFILE=test_outfile(OUTFILE);

  if (COLUMN == NULL)
    series=(double**)get_multi_series(infile,&(pars.LENGTH),EXCLUDE,
				      &(pars.DIM),"",dimset,VERBOSITY);
  else
    series=(double**)get_multi_series(infile,&(pars.LENGTH),EXCLUDE,
				      &(pars.DIM),COLUMN,dimset,VERBOSITY);

  check_alloc(sermin=malloc(sizeof(double)*pars.DIM));
  check_alloc(serinter=malloc(sizeof(double)*pars.DIM));
  rescale_data(series,pars,sermin,serinter);

  alldim=pars.DIM;
  pars.DIM=compdim;

  check_alloc(useries=(unsigned int**)malloc(sizeof(int*)*pars.DIM));
  check_alloc(PART=(unsigned int*)malloc(sizeof(int)*pars.DIM));
  if (SPART == NULL)
    for (i=0;i<pars.DIM;i++)
      PART[i]=2;
  else
    set_part(PART,SPART,pars);

  for (i=0;i<pars.DIM;i++) {
    check_alloc(useries[i]=(unsigned int*)malloc(sizeof(int)*pars.LENGTH));
    for (j=0;j<pars.LENGTH;j++) {
      h=(unsigned int)((series[i][j]*(double)PART[i])/serinter[i]);
      useries[i][j]=((h<PART[i])?h:(PART[i]-1));
    }
  }

  check_alloc(future=(unsigned int*)malloc(sizeof(int)*pars.LENGTH));
  if (WHICHFUTURE>0) {
    check_alloc(wcol=(char*)calloc((size_t)10,(size_t)1));
    sprintf(wcol,"%u",WHICHFUTURE);
    dfuture=(double**)get_multi_series(infile,&(pars.LENGTH),EXCLUDE,
				       &wdim,wcol,1,0L);
    for (i=0;i<pars.LENGTH-1;i++) {
      future[i] = (dfuture[0][i]>0.0);
    }
    future[pars.LENGTH-1]=1;
    free(dfuture);
  }
  else {
    for (i=0;i<pars.LENGTH;i++)
      future[i]=1;
  }
  root=make_ptree(PART[0]);

  count=0;
  check_alloc(segment=(long*)malloc(sizeof(long)*pars.LENGTH));
  for (i=0;i<pars.LENGTH;i += ODEPTH) {
    has_future=1;
    for (j=0;j<ODEPTH;j++)
      has_future &= future[i+j];
    if (has_future) {
      fill_tree(root,useries,pars,i,0,DEPTH,PART);
    }
    segment[count++]=i;
  }

  if (scrambleset)
    lscramble(segment,count);

  check_alloc(use_array=(char*)malloc(sizeof(char)*pars.LENGTH));
  for (i=0;i<count;i++) {
    n=segment[i];
    has_future=1;
    for (j=0;j<ODEPTH;j++)
      has_future &= future[n+j];
    if (has_future) {
      h=read_tree(root,useries,pars,n,0,DEPTH);
      if (h <= MAXOUT)
	use_array[i]=1;
      else 
	use_array[i]=0;
    }
    else 
      use_array[i]=2;
  }
  pars.DIM=alldim;

  fout=fopen(OUTFILE,"w");
  fprintf(fout,"%s\n",outstring);
  fprintf(fout,"#Content: ");
  for (i=0;i<pars.DIM;i++)
    fprintf(fout,"x%d ",i+1);
  fprintf(fout,"future time_index\n");
  fflush(fout);

  for (i=0;i<count;i++) {
    n=segment[i];
    if (use_array[i] == 1) {
      for (j=0;j<(ODEPTH-1);j++) {
	for (k=0;k<pars.DIM;k++)
	  fprintf(fout,"%f ",series[k][n+j]+sermin[k]);
	fprintf(fout,"1 %u\n",n+j);
      }
      for (k=0;k<pars.DIM;k++)
	fprintf(fout,"%f ",series[k][n+ODEPTH-1]+sermin[k]);
      if (i<(count-1)) {
	if (segment[i+1] == (n+ODEPTH) && use_array[i+1])
	  fprintf(fout,"1 %u\n",n+ODEPTH-1);
	else 
	  fprintf(fout,"0 %u\n",n+ODEPTH-1);
      }
      else 
	fprintf(fout,"0 %u\n",n+ODEPTH-1);
    }
  }

  fclose(fout);
}
