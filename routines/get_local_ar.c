/*
 *   This file is part of OLANGEVIN
 *
 *   Copyright (c) 2013-2014 Rainer Hegger
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

char isinitgetarcoeff=0;
double **mat,*av,*vec;

void mat_init_ar(struct param p)
{
  unsigned int i;

  check_alloc(mat=(double**)malloc(sizeof(double*)*(p.DIM*p.EMB)));
  for (i=0;i<p.DIM*p.EMB;i++)
    check_alloc(mat[i]=(double*)malloc(sizeof(double)*(p.DIM*p.EMB)));
  check_alloc(av=(double*)malloc(sizeof(double)*(p.DIM*(p.EMB+1))));
  check_alloc(vec=(double*)malloc(sizeof(double)*(p.DIM*p.EMB)));

  isinitgetarcoeff=1;
}

void get_ar_coeff(double **series,struct param p,struct sfound sf,
		  double *coef0,double **coef,double *dif)
{
  double **imat,*sd,*sd1,h,sum,sweight,*weight;
  unsigned long *f,which,count;
  unsigned int num,i,i1,i2,j,j1,j2,n;

  if (!isinitgetarcoeff) {
    mat_init_ar(p);
  }
  f=sf.found;
  num=p.MINN;
  weight=sf.weight;
#if defined(WEIGHT)
  sweight=weight[0];
  for (i=1;i<p.MINN;i++)
    sweight += weight[i];
#else
  sweight=(double)p.MINN;
#endif
  
  for (i=0;i<p.DIM*p.EMB;i++) {
    sum=0.0;
    i1=i/p.EMB;
    i2=i%p.EMB;
    for (n=0;n<num;n++) {
      which=f[n];
#if defined(WEIGHT)
      sum += series[i1][which-i2]*weight[n];
#else
      sum += series[i1][which-i2];
#endif
    }
    av[i]=sum/sweight;
  }
  for (i=0;i<p.DIM;i++) {
    sum=0.0;
    for (n=0;n<num;n++) {
      which=f[n];
#if defined(WEIGHT)
      sum += series[i][which+1]*weight[n];
#else
      sum += series[i][which+1];
#endif
    }
    av[p.DIM*p.EMB+i]=sum/sweight;
  }

  for (i=0;i<p.DIM*p.EMB;i++) {
    i1=i/p.EMB;
    i2=i%p.EMB;
    sd=series[i1];
    for (j=i;j<p.DIM*p.EMB;j++) {
      j1=j/p.EMB;
      sd1=series[j1];
      j2=j%p.EMB;
      sum=0.0;
      for (n=0;n<num;n++) {
	which=f[n];
#if defined(WEIGHT)
	sum += sd[which-i2]*sd1[which-j2]*weight[n];
#else
	sum += sd[which-i2]*sd1[which-j2];
#endif
      }
      mat[i][j]=sweight/(sweight-1.0)*(sum/sweight-av[i]*av[j]);
      mat[j][i]=mat[i][j];
    }
  }

  imat=invert_matrix(mat,p.DIM*p.EMB);

  for (i=0;i<p.DIM;i++) {
    sd=series[i];
    for (j=0;j<p.DIM*p.EMB;j++) {
      j1=j/p.EMB;
      sd1=series[j1];
      j2=j%p.EMB;
      sum=0.0;
      for (n=0;n<num;n++) {
	which=f[n];
#if defined(WEIGHT)
	sum += sd[which+1]*sd1[which-j2]*weight[n];
#else
	sum += sd[which+1]*sd1[which-j2];
#endif
      }
      vec[j]=sweight/(sweight-1.0)*(sum/sweight-av[p.DIM*p.EMB+i]*av[j]);
    }
    for (j=0;j<p.DIM*p.EMB;j++) {
      coef[i][j]=0.0;
      for (n=0;n<p.DIM*p.EMB;n++)
	coef[i][j] += imat[j][n]*vec[n];
    }

    coef0[i]=av[p.DIM*p.EMB+i];
    for (j=0;j<p.DIM*p.EMB;j++)
      coef0[i] -= coef[i][j]*av[j];

    dif[i]=0.0;
    for (n=0;n<num;n++) {
      which=f[n];
      h=0.0;
      for (j=0;j<p.DIM*p.EMB;j++) {
	j1=j/p.EMB;
	j2=j%p.EMB;
	h +=  coef[i][j]*series[j1][which-j2];
      }
      dif[i] += sqr(series[i][which+1]-coef0[i]-h);
    }
    dif[i]=sqrt(dif[i]/(double)num);
  }

  for (i=0;i<p.DIM*p.EMB;i++)
    free(imat[i]);
  free(imat);
}

void make_ar_cast(double **x,double *new,struct param p,double *coeff0,
		  double **coeff,double *dif)
{
  unsigned int i,i1,i2,d;

  get_noise(p,new);

  for (d=0;d<p.DIM;d++) {
    new[d] *= dif[d];
    new[d] += coeff0[d];
    for (i=0;i<p.DIM*p.EMB;i++) {
      i1=i/p.EMB;
      i2=i%p.EMB;
      new[d] += coeff[d][i]*x[i1][p.hdim-i2];
    }
  }
}
