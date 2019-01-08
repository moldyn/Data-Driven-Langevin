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
#include <math.h>
#include "olang.h"

void make_corr(double **x,unsigned int *fut,
	       struct param p,double *s)
{
  unsigned int cl,cl1;
  unsigned long i,j,k,n,count,start,stop;
  double *h,**mat,**imat,*vec,*hvec,sigma;

  cl=p.AR_SIZE;
  cl1=cl+1;

  check_alloc(mat=(double**)malloc(sizeof(double*)*cl));
  check_alloc(vec=(double*)malloc(sizeof(double)*cl));
  check_alloc(hvec=(double*)malloc(sizeof(double)*cl1));
  for (i=0;i<cl;i++) {
    check_alloc(mat[i]=(double*)malloc(sizeof(double)*cl));
  }

  for (i=0;i<p.DIM;i++) {

    for (j=0;j<cl1;j++)
	hvec[j]=0.0;

    h=x[i];

    for (j=0;j<cl1;j++) {
      count=0;
      start=stop=0;
      while (start < (p.LENGTH-1)) {
	while (fut[stop] == 0) {
	  start=stop;
	  stop++;
	}
	while (fut[stop] != 0)
	    stop++;
	for (n=start+j;n<=stop;n++) {
	  hvec[j] += h[n]*h[n-j];
	  count++;
	}
	start=stop;
      }
      hvec[j] /= count;
    }

    for (j=0;j<cl;j++) {
      vec[j]=hvec[j+1];
      for (k=0;k<=j;k++) {
	mat[j][k]=mat[k][j]=hvec[j-k];
      }
    }
    imat=invert_matrix(mat,cl);
    for (j=0;j<cl;j++) {
      p.AR_MOD[i][j]=0.0;
      for (k=0;k<cl;k++)
	p.AR_MOD[i][j] += imat[j][k]*vec[k];
    }
    for (j=0;j<cl;j++)
      free(imat[j]);
    free(imat);
    sigma=hvec[0];
    for (j=0;j<cl;j++)
      sigma -= p.AR_MOD[i][j]*hvec[j+1];
    if (sigma < 0.0) {
      fprintf(stderr,"Standard deviation becomes negative. Exiting\n");
      exit(MAKE_CORR_NEGATIVE_SIGMA);
    }
    s[i]=sqrt(sigma);
  }

  free(hvec);
  for (i=0;i<cl;i++) 
    free(mat[i]);
  free(mat);
}
