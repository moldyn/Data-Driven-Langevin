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

#include <stdlib.h>
#include <string.h>
#include "olang.h"

void read_ar_file(struct param *p,char *fname)
{
  FILE *file=NULL;
  char **format,*input;
  unsigned int i,j;

  file=fopen(fname,"r");

  if (file==NULL) {
    fprintf(stderr,"ar-model file does not exist!\n");
    exit(READ_AR_FILE_NOT_EXIST);
  }

  check_alloc(input=(char*)malloc(sizeof(char)*1024));
  check_alloc(format=(char**)malloc(sizeof(char*)*p->DIM));
  for (i=0;i<p->DIM;i++) {
    check_alloc(format[i]=(char*)malloc(sizeof(char)*p->DIM*4));
    format[i][0]=0;
    for (j=0;j<i;j++)
      strcat(format[i],"%*lf");
    strcat(format[i],"%lf");
  }

  if (fgets(input,1024,file) != NULL) {
    sscanf(input,"%u\n",&(p->AR_SIZE));
  }
  else {
    fprintf(stderr,"Not enough lines in ar-model file!\n");
    exit(READ_AR_FILE_NOT_ENOUGH_LINES);
  }

  check_alloc((p->AR_SIGMA)=(double*)malloc(sizeof(double)*(p->DIM)));
  if (fgets(input,1024,file) != 0) {
    for (i=0;i<p->DIM;i++) {
      sscanf(input,format[i],&(p->AR_SIGMA[i]));
    }
  }
  else {
    fprintf(stderr,"Not enough lines in ar-model file!\n");
    exit(READ_AR_FILE_NOT_ENOUGH_LINES);
  }

  if (p->AR_SIZE > 0) {
    check_alloc((p->AR_MOD)=(double**)malloc(sizeof(double*)*(p->DIM)));
    for (j=0;j<p->DIM;j++)
      check_alloc((p->AR_MOD[j])=(double*)malloc(sizeof(double)*(p->AR_SIZE)));

    for (j=0;j<p->AR_SIZE;j++) {
      if (fgets(input,1024,file) != 0) {
	for (i=0;i<p->DIM;i++)
	  sscanf(input,format[i],&(p->AR_MOD[i][j]));
      }
      else {
	fprintf(stderr,"Not enough lines in ar-model file!\n");
	exit(READ_AR_FILE_NOT_ENOUGH_LINES);
      }
    }

    free(input);
    for (i=0;i<p->DIM;i++)
      free(format[i]);
    free(format);

    check_alloc(p->noise_hist=(double**)malloc(sizeof(double*)*p->DIM));
    for (i=0;i<p->DIM;i++)
      check_alloc(p->noise_hist[i]=(double*)malloc(sizeof(double)*p->AR_SIZE));
  }
}
