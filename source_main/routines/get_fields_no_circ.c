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

char isinitgetfieldsnocirc=0;
double *av_fut,*av_past,**mat_pastpast,**mat_futpast,**mat_futfut,**y,**yf,**yp;

void mat_init_no_circ(struct param p)
{
  unsigned int i;

  check_alloc(mat_pastpast=(double**)malloc(sizeof(double*)*p.DIM));
  check_alloc(mat_futpast=(double**)malloc(sizeof(double*)*p.DIM));
  check_alloc(mat_futfut=(double**)malloc(sizeof(double*)*p.DIM));
  check_alloc(y=(double**)malloc(sizeof(double*)*p.DIM));
  check_alloc(yf=(double**)malloc(sizeof(double*)*p.DIM));
  check_alloc(yp=(double**)malloc(sizeof(double*)*p.DIM));
  for (i=0;i<p.DIM;i++) {
    check_alloc(mat_pastpast[i]=(double*)malloc(sizeof(double)*p.DIM));
    check_alloc(mat_futpast[i]=(double*)malloc(sizeof(double)*p.DIM));
    check_alloc(mat_futfut[i]=(double*)malloc(sizeof(double)*p.DIM));
    check_alloc(y[i]=(double*)malloc(sizeof(double)*p.MINN));
    check_alloc(yf[i]=(double*)malloc(sizeof(double)*p.MINN));
    check_alloc(yp[i]=(double*)malloc(sizeof(double)*p.MINN));
  }
  check_alloc(av_fut=(double*)malloc(sizeof(double)*p.DIM));
  check_alloc(av_past=(double*)malloc(sizeof(double)*p.DIM));

  isinitgetfieldsnocirc=1;
}

void get_fields_no_circ(double **series,struct param p,struct sfound sf,
		   double *drift,double **gamma,double **dif,double **x)
{
  unsigned int number;
  long i,j,k,l,d,d1,fi,fi1,fim1;
  double *sd,*sd1,sum,num,**imat,hf,hp,sweights,*weight,hw,shift;
  unsigned long *f,hdim;

  f=sf.found;

  if (!isinitgetfieldsnocirc)
    mat_init_no_circ(p);

  number=p.MINN;  
  num=(double)number;
  hdim=p.hdim;
  weight=sf.weight;
#if defined(WEIGHTS)
  sweights=weight[0];
  for (i=1;i<p.MINN;i++)
    sweights += weight[i];
#else
  sweights=(double)p.MINN;
#endif

  for (i=0;i<p.DIM;i++) {
    shift=x[i][hdim];
    for (j=0;j<p.MINN;j++) {
      y[i][j]=series[i][f[j]]-shift;
      yf[i][j]=series[i][f[j]+1]-shift;
      yp[i][j]=series[i][f[j]-1]-shift;
      if (y[i][j] < -M_PI) y[i][j]=pi2+y[i][j];
      if (y[i][j] > M_PI) y[i][j]=y[i][j]-pi2;
      if (yf[i][j] < -M_PI) yf[i][j]=pi2+yf[i][j];
      if (yf[i][j] > M_PI) yf[i][j]=yf[i][j]-pi2;
      if (yp[i][j] < -M_PI) yp[i][j]=pi2+yp[i][j];
      if (yp[i][j] > M_PI) yp[i][j]=yp[i][j]-pi2;
    }
  }
  
  /* average predecessor av_fut = <x_{n-1}> and
   *          successor av_past = <x_{n+1}> */
  for (d=0;d<p.DIM;d++) {
    hp=hf=0.0;
    for (i=0;i<number;i++) {
#if defined(WEIGHTS)
      hw=weight[i];
      hf += yf[d][i]*hw;
      hp += yp[d][i]*hw;
    }
    av_fut[d] = hf/sweights;
    av_past[d] = hp/sweights;
#else
      hf += yf[d][i];
      hp += yp[d][i];
    }
    av_fut[d] = hf/num;
    av_past[d] = hp/num;
#endif
  }

  /* correlation matrices mat_pastpast = C(x_{n-1},x_{n-1})
   *                    and mat_futfut = C(x_{n+1},x_{n+1}) */
  for (d=0;d<p.DIM;d++) {
    for (d1=0;d1<=d;d1++) {
      hp=hf=0.0;
      for (i=0;i<number;i++) {
#if defined(WEIGHTS)
	hw=weight[i];
	hp += yp[d][i]*yp[d1][i]*hw;
	hf += yf[d][i]*yf[d1][i]*hw;
      }
      mat_pastpast[d][d1]=(hp/sweights-av_past[d]*av_past[d1]);
      mat_pastpast[d1][d]=mat_pastpast[d][d1];
      mat_futfut[d][d1]=(hf/sweights-av_fut[d]*av_fut[d1]);
      mat_futfut[d1][d]=mat_futfut[d][d1];
#else
	hp += yp[d][i]*yp[d1][i];
	hf += yf[d][i]*yf[d1][i];
      }
      mat_pastpast[d][d1]=num/(num-1.0)*(hp/num-av_past[d]*av_past[d1]);
      mat_pastpast[d1][d]=mat_pastpast[d][d1];
      mat_futfut[d][d1]=num/(num-1.0)*(hf/num-av_fut[d]*av_fut[d1]);
      mat_futfut[d1][d]=mat_futfut[d][d1];
#endif
    }
}

  /* correlation matrix mat_futpast = C(x_{n+1},x_{n-1}) */
  for (d=0;d<p.DIM;d++) {
    for (d1=0;d1<p.DIM;d1++) {
      hp=0.0;
      for (i=0;i<number;i++) {
#if defined(WEIGHTS)
	hp += yf[d][i]*yp[d1][i]*weight[i];
      }
      mat_futpast[d][d1]=(hp/sweights-av_fut[d]*av_past[d1]);
#else
      hp += yf[d][i]*yp[d1][i];
      }
      mat_futpast[d][d1]=num/(num-1.0)*(hp/num-av_fut[d]*av_past[d1]);
#endif
    }
  }

  /* friction matrix gamma = C(x_{n+1},x_{n-1}) C^{-1}(x_{n-1},x_{n-1}) */
  imat=invert_matrix(mat_pastpast,p.DIM);

  for (d=0;d<p.DIM;d++) {
    for (d1=0;d1<p.DIM;d1++) {
      hf= mat_futpast[d][0]*imat[0][d1];
      for (k=1;k<p.DIM;k++) {
	hf += mat_futpast[d][k]*imat[k][d1];
      }
      gamma[d][d1]=hf;
    }
  }

  /* "squared" diffusion matrix now mat_futfut = K K^T
   * = C(x_{n+1},x_{n+1}) - gamma C(x_{n-1},x_{n-1}) gamma^T */
  for (i=0;i<p.DIM;i++) {
    for (j=0;j<=i;j++) {
      for (k=0;k<p.DIM;k++) {
	mat_futfut[i][j] -= mat_futpast[i][k]*gamma[j][k];
      }
      mat_futfut[j][i]=mat_futfut[i][j];
    }
  }

  for (d=0;d<p.DIM;d++)
    free(imat[d]);
  free(imat);
  
  /* Cholesky decomposition to obtain K */
  for (i=0;i<p.DIM;i++) {
    for (j=0;j<i;j++) {
      sum=mat_futfut[i][j];
      for (k=0;k<j;k++) {
	sum -= mat_futfut[i][k]*mat_futfut[j][k];
      }
      mat_futfut[i][j]=sum/mat_futfut[j][j];
    }

    sum=mat_futfut[i][i];
    for (j=0;j<i;j++) {
      sum -= mat_futfut[i][j]*mat_futfut[i][j];
    }
    if (sum >= 0.0) {
      mat_futfut[i][i]=sqrt(sum);
    }
    else {
      printf("%ld %e\n",i,sum);
      exit(GET_FIELDS_NO_CHOLESKY);
    }
  }

  /* diffusion matrix dif = K*/
  for (i=0;i<p.DIM;i++)
    for (j=0;j<=i;j++)
      dif[i][j]=mat_futfut[i][j];

  /* drift matrix drift = <x_{n+1}> - y + gamma*(y - <x_{n-1}>) */
  for (d=0;d<p.DIM;d++) {
    drift[d]=av_fut[d];
    for (k=0;k<p.DIM;k++)
      drift[d] -= gamma[d][k]*av_past[k];
  }
}
