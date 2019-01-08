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

char isinitgetfields=0;
double **mat,*av;

void mat_init(struct param p)
{
  unsigned int i;

  check_alloc(mat=(double**)malloc(sizeof(double*)*p.DIM));
  for (i=0;i<p.DIM;i++)
    check_alloc(mat[i]=(double*)malloc(sizeof(double)*p.DIM));
  check_alloc(av=(double*)malloc(sizeof(double)*p.DIM));
  isinitgetfields=1;
}

void get_fields(double **series,struct param p,struct sfound sf,
		double *drift,double **dif)
{
  unsigned int number;
  long i,j,k,d,d1,fi,fi1;
  double *sd,*sd1,sum,num,h,sweights,*weight;
  unsigned long *f;

  f=sf.found;

  if (!isinitgetfields)
    mat_init(p);

  number=p.MINN;  
  num=(double)number;
  weight=sf.weight;
#if defined(WEIGHTS)
  sweights=weight[0];
  for (i=1;i<p.MINN;i++)
    sweights += weight[i];
#else
  sweights=(double)p.MINN;
#endif

  /* average successor drift = <x_{n+1}> */
  for (d=0;d<p.DIM;d++) {
    h=0.0;
    sd=series[d];
    for (i=0;i<number;i++)
      h += sd[f[i]+1]*weight[i];
    av[d] = h/sweights;
    drift[d]=av[d];
  }

  /* correlation matrix mat = C(x_{n+1},x_{n+1}) = D D^T */
  for (d=0;d<p.DIM;d++) {
    sd=series[d];
    for (d1=0;d1<=d;d1++) {
      sd1=series[d1];
      sum=0.0;
      for (i=0;i<number;i++) {
	fi=f[i];
	fi1=fi+1;
#if defined(WEIGHTS)
	sum += (sd[fi1])*(sd1[fi1])*weight[i];
      }
      mat[d][d1]=(sum/sweights-av[d]*av[d1]);
      mat[d1][d]=mat[d][d1];
#else
	sum += (sd[fi1])*(sd1[fi1]);
      }
      mat[d][d1]=num/(num-1.0)*(sum/num-av[d]*av[d1]);
      mat[d1][d]=mat[d][d1];
#endif
    }
  }

#ifdef SYM_DIFFOP
  get_diffusion(mat,p);
  for (i=0;i<p.DIM;i++)
    for (j=0;j<p.DIM;j++)
      dif[i][j]=mat[i][j];
#else
  /* Cholesky decomposition to obtain dif = D */

  for (i=0;i<p.DIM;i++) {
    for (j=0;j<i;j++) {
      sum=mat[i][j];
      for (k=0;k<j;k++) {
	sum -= mat[i][k]*mat[j][k];
      }
      mat[i][j]=sum/mat[j][j];
    }

    sum=mat[i][i];
    for (j=0;j<i;j++) {
      sum -= mat[i][j]*mat[i][j];
    }
    if (sum >= 0.0) {
      mat[i][i]=sqrt(sum);
    }
    else {
      printf("%ld %e\n",i,sum);
      exit(GET_FIELDS_CHOLESKY);
    }
  }
  for (i=0;i<p.DIM;i++)
    for (j=0;j<=i;j++)
      dif[i][j]=mat[i][j];
#endif
}
