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
unsigned long FLENGTH=ULONG_MAX,EXCLUDE=0,ignore=0;
double INNOISE=0.0;
char *COLUMN=NULL,*OUTFILE=NULL,*refname=NULL;


struct param pars={ULONG_MAX,1,1,1,50,1.0,.AR_SIZE=0,.maxr=0.0};
char dimset=0,reference=0,noiseonly=0,neighborhood=0;
unsigned long seed=0x4973592378L;

void show_options(char *progname)
{

  fprintf(stderr," Usage: %s [Options]\n",progname);
  fprintf(stderr," Options:\n");
  fprintf(stderr,"Everything not being a valid option will be interpreted"
          " as a possible"
          " datafile.\nIf no datafile is given stdin is read. Just - also"
          " means stdin\n");
  fprintf(stderr,"\t-l # of lines to be used [default whole file]\n");
  fprintf(stderr,"\t-x # of lines to be ignored [default %lu]\n",EXCLUDE);
  fprintf(stderr,"\t-c column [default 1,...,# of components]\n");
  fprintf(stderr,"\t-m # of components [default %u]\n",pars.DIM);
  fprintf(stderr,"\t-M Embedding dimension [default 1 (no embedding)]\n");
  fprintf(stderr,"\t-d delay for the embedding [default %u]\n",pars.DELAY);
  fprintf(stderr,"\t-k # of neighbors  [default %u]\n",pars.MINN);
  fprintf(stderr,"\t-F # of component which defines future [not set]\n");
#if defined(WEIGHTS)
  fprintf(stderr,"\t-r radius for the weights [default not set]\n");
#endif
  fprintf(stderr,"\t-H output additional neighborhood info\n");
  fprintf(stderr,"\t-N output noise only [fields+noise]\n");
  fprintf(stderr,"\t-X ignore lines for test [default: 0]\n");
  fprintf(stderr,"\t-L # of lines for test [default all input data]\n");
  fprintf(stderr,"\t-t test points for field estimation [default same as data set]\n");
  fprintf(stderr,"\t-%% add initial noise (absolute amplitude) [none]\n");
  fprintf(stderr,"\t-I seed for the rnd-generator (If seed=0, the time\n"
          "\t\tcommand is used to set the seed) [Default: fixed]\n");
  fprintf(stderr,"\t-o output file [default 'datafile'.ltm;"
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
  if ((out=check_option(in,n,'H','n')) != NULL)
    neighborhood=1;
  if ((out=check_option(in,n,'N','n')) != NULL)
    noiseonly=1;
  if ((out=check_option(in,n,'d','u')) != NULL)
    sscanf(out,"%u",&(pars.DELAY));
  if ((out=check_option(in,n,'L','u')) != NULL)
    sscanf(out,"%lu",&FLENGTH);
  if ((out=check_option(in,n,'X','u')) != NULL)
    sscanf(out,"%lu",&ignore);
#if defined(WEIGHTS)
  if ((out=check_option(in,n,'r','f')) != NULL)
    sscanf(out,"%lf",&(pars.maxr));
#endif
  if ((out=check_option(in,n,'t','s')) != NULL) {
    reference=1;
    refname=out;
  }
  if ((out=check_option(in,n,'k','u')) != NULL)
    sscanf(out,"%u",&(pars.MINN));
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
  long li,lj,n,refn;
  char *outstring,*infile,*wcol;
  char ntoolarge;
  double **series,**rseries,*sermin,*serinter,**dfuture,max,**cast,sweights;
  double *new,*drift,**diffusion;
  double *eccentricity;
  double abs_ecc,var_ratio_fut,var_ratio_past;
  struct sfound sf;
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

  if (pars.MINN < (pars.DIM+1)) {
    fprintf(stderr,"Too few neighbors. System is underdetermined.\n");
    exit(LANGEVIN_MAIN_TOO_FEW_MINN);
  }

  infile=search_datafile(argc,argv,NULL,VERBOSITY);
  if (infile == NULL) {
    fprintf(stderr,"No input datafile found.\n");
    exit(LANGEVIN_MAIN_NO_INPUTFILE);
  }
  if (reference) {
    fout=fopen(refname,"r");
    if (fout == NULL) {
      fprintf(stderr,"No test point datafile found.\n");
      exit(LANGEVIN_MAIN_NO_INPUTFILE);
    }
    else
      fclose(fout);
  }

  if (OUTFILE == NULL) {
    check_alloc(OUTFILE=(char*)calloc(strlen(infile)+6,(size_t)1));
    sprintf(OUTFILE,"%s.ltm",infile);
  }
  OUTFILE=test_outfile(OUTFILE);

  if (!reference) {
    refname=infile;
    ignore += EXCLUDE;
  }

  if (COLUMN == NULL)
    series=(double**)get_multi_series(infile,&(pars.LENGTH),EXCLUDE,
				      &(pars.DIM),"",dimset,VERBOSITY);
  else
    series=(double**)get_multi_series(infile,&(pars.LENGTH),EXCLUDE,
				      &(pars.DIM),COLUMN,dimset,VERBOSITY);

  hdim=(pars.EMB-1)*pars.DELAY;
  pars.hdim=hdim;

  if (COLUMN == NULL)
    rseries=(double**)get_multi_series(refname,&FLENGTH,ignore,
				       &(pars.DIM),"",dimset,VERBOSITY);
  else
    rseries=(double**)get_multi_series(refname,&FLENGTH,ignore,
				       &(pars.DIM),COLUMN,dimset,VERBOSITY);

  if (VERBOSITY)
    fprintf(stderr,"%lf\n",rseries[0][0]);

  init_noise(pars,seed);

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
    for (j=0;j<FLENGTH;j++)
      rseries[i][j]=(rseries[i][j]-sermin[i])/max;
  }
  if (pars.maxr > 0.0)
    pars.maxr /= max;

  check_alloc(future=(unsigned int*)malloc(sizeof(int)*pars.LENGTH));
  if (WHICHFUTURE>0) {
    check_alloc(wcol=(char*)calloc((size_t)10,(size_t)1));
    sprintf(wcol,"%u",WHICHFUTURE);
    dfuture=(double**)get_multi_series(infile,&(pars.LENGTH),EXCLUDE,
				       &wdim,wcol,1,0L);
    for (li=hdim;li<pars.LENGTH;li++) {
      has_future= (dfuture[0][li]>0.0);
      for (lj=1;lj<=hdim;lj++)
	has_future &= (dfuture[0][li-lj]>0.0);
      future[li]=has_future;
    }
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
  check_alloc(eccentricity=(double*)malloc(sizeof(double)*pars.DIM));

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
  fprintf(fout,"#Content: ");
  for (i=0;i<pars.DIM;i++) {
    for (j=0;j<pars.EMB;j++) {
      if (pars.EMB>1) {
	fprintf(fout,"x%d_%d ",i+1,j+1);
      }
      else {
	fprintf(fout,"x%d ",i+1);
      }
    }
  }
  if (!noiseonly) {
    for (i=0;i<pars.DIM;i++)
      fprintf(fout,"h%d ",i+1);
    for (i=0;i<pars.DIM;i++)
      for (j=0;j<=i;j++)
	fprintf(fout,"D_%d_%d ",i+1,j+1);
    if (!reference) {
      for (i=0;i<pars.DIM;i++)
        fprintf(fout,"xi%d ",i+1);
    }
#if defined(WEIGHTS)
    fprintf(fout,"distance sweights ");
#else
    fprintf(fout,"distance ");
#endif
    if (neighborhood) {
      fprintf(fout,"abs_ecc ");
      for (i=0;i<pars.DIM;i++)
        fprintf(fout,"ecc%d ",i+1);
      fprintf(fout,"var_ratio_fut ");
    }
    fprintf(fout,"future\n");
  }
  else {
    if (!reference) {
      for (i=0;i<pars.DIM;i++)
        fprintf(fout,"xi%d ",i+1);
    }
#if defined(WEIGHTS)
    fprintf(fout,"distance sweights future\n");
#else
    fprintf(fout,"distance future\n");
#endif
  }
  fflush(fout);

  n=0;
  while (n < FLENGTH) {

    refn=n+ignore-EXCLUDE;
    if (!reference) {
      while (!future[refn]) {
        n++;
        refn++;
        if (n>=(pars.LENGTH-1))
          goto ntoolarge;
      }
    }
    for (i=0;i<=hdim;i++)
      for (j=0;j<pars.DIM;j++)
	cast[j][hdim-i]=rseries[j][n-i];

    /* search neighbors */
    if (!reference)
      pars.MINN += 2;

    search_neighbors(series,cast,pars,sf);

    if (!reference) {
      j=0;
      for (i=0;i<pars.MINN;i++) {
        if ((sf.found[i] == refn) || (sf.found[i] == (refn-1))) {
          j++;
        }
        sf.found[i]=sf.found[j];
        sf.distance[i]=sf.distance[j];
        j++;
      }
      pars.MINN -= 2;
    }

    /* calculate fields */
    get_fields(series,pars,sf,drift,diffusion);
    for (i=0;i<pars.DIM;i++)
      drift[i] -= cast[i][hdim];

    /* calculate noise */
    if (!reference)
      make_test(series,refn,new,pars,drift,diffusion);

    /* calculate neighborhood information */
    if (neighborhood)
      neighborhood_info(series,pars,sf,cast,eccentricity,&abs_ecc,
				      &var_ratio_fut,&var_ratio_past);

    /* print coordinate */
    for (i=0;i<pars.DIM;i++)
      for (j=0;j<pars.EMB;j++)
        fprintf(fout,"%e ",cast[i][hdim-j*pars.DELAY]*max+sermin[i]);

    /* print fields */
    if (!noiseonly) {
      for (i=0;i<pars.DIM;i++)
        fprintf(fout,"%e ",(drift[i])*max);
      for (i=0;i<pars.DIM;i++)
        for (j=0;j<=i;j++)
          fprintf(fout,"%e ",diffusion[i][j]*max);
    }

    /* print noise */
    if (!reference) {
      for (i=0;i<pars.DIM;i++)
        fprintf(fout,"%e ",new[i]);
    }

    /* print distance, weight information and future */
#if defined(WEIGHTS)
    sweights=sf.weight[0];
    for (i=1;i<pars.MINN;i++)
      sweights += sf.weight[i];
    fprintf(fout,"%e %lf ",sf.distance[pars.MINN-1]*max,sweights);
#else
    fprintf(fout,"%e ",sf.distance[pars.MINN-1]*max);
#endif
    if (neighborhood) {
      fprintf(fout,"%e ",abs_ecc*max);
      for (i=0;i<pars.DIM;i++)
        fprintf(fout,"%e ",(eccentricity[i])*max);
      fprintf(fout,"%e ",var_ratio_fut);
    }
    fprintf(fout,"%u\n",future[refn+1]);
    fflush(fout);
    n++;
  }
 ntoolarge: fclose(fout);

  return 0;
}
