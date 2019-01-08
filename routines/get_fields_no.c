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

char isinitgetfieldsno=0;
double *av_fut,*avv_fut,*av_past,*avv_past,*av_now,
  **mat_pastpast,**mat_futpast,**mat_futfut;

void mat_init_no(struct param p)
{
  unsigned int i;

  check_alloc(mat_pastpast=(double**)malloc(sizeof(double*)*p.DIM));
  check_alloc(mat_futpast=(double**)malloc(sizeof(double*)*p.DIM));
  check_alloc(mat_futfut=(double**)malloc(sizeof(double*)*p.DIM));
  for (i=0;i<p.DIM;i++) {
    check_alloc(mat_pastpast[i]=(double*)malloc(sizeof(double)*p.DIM));
    check_alloc(mat_futpast[i]=(double*)malloc(sizeof(double)*p.DIM));
    check_alloc(mat_futfut[i]=(double*)malloc(sizeof(double)*p.DIM));
  }
  check_alloc(av_fut=(double*)malloc(sizeof(double)*p.DIM));
  check_alloc(av_past=(double*)malloc(sizeof(double)*p.DIM));
  check_alloc(av_now=(double*)malloc(sizeof(double)*p.DIM));
  check_alloc(avv_fut=(double*)malloc(sizeof(double)*p.DIM));
  check_alloc(avv_past=(double*)malloc(sizeof(double)*p.DIM));

  isinitgetfieldsno=1;
}

void get_fields_no(double **series,struct param p,struct sfound sf,
		   double *drift,double **gamma,double **dif,double **x)
{
  unsigned int number;
  long i,j,k,l,d,d1,fi,fi1,fim1;
  double *sd,*sd1,sum,num,**imat,hf,hp,hn,sweights,*weight,hw;
  unsigned long *f,hdim;

  f=sf.found;

  if (!isinitgetfieldsno)
    mat_init_no(p);

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

  /* average predecessor av_fut = <x_{n-1}> and
   *          successor av_past = <x_{n+1}> and
   *            neighbor av_now =  <x_{n}> */
  for (d=0;d<p.DIM;d++) {
    hp=hf=hn=0.0;
    sd=series[d];
    for (i=0;i<number;i++) {
      fi=f[i];
#if defined(WEIGHTS)
      hw=weight[i];
      hf += sd[fi+1]*hw;
      hp += sd[fi-1]*hw;
      hn += sd[fi]*hw;
    }
    av_fut[d] = hf/sweights;
    av_past[d] = hp/sweights;
    av_now[d] = hn/sweights;
    avv_fut[d] = (hf-hn)/sweights;
    avv_past[d] = (hn-hp)/sweights;
#else
      hf += sd[fi+1];
      hp += sd[fi-1];
      hn += sd[fi];
    }
    av_fut[d] = hf/num;
    av_past[d] = hp/num;
    av_now[d] = hn/num;
    avv_fut[d] = (hf-hn)/num;
    avv_past[d] = (hn-hp)/num;
#endif
  }

  /* correlation matrices mat_pastpast = C(x_{n}-x_{n-1},x_{n}-x_{n-1})
   *                    and mat_futfut = C(x_{n+1}-x_{n},x_{n+1}-x_{n}) */
  for (d=0;d<p.DIM;d++) {
    sd=series[d];
    for (d1=0;d1<=d;d1++) {
      sd1=series[d1];
      hp=hf=0.0;
      for (i=0;i<number;i++) {
	fi=f[i];
	fi1=fi+1;
	fim1=fi-1;
#if defined(WEIGHTS)
	hw=weight[i];
	hp += (sd[fi]-sd[fim1])*(sd1[fi]-sd1[fim1])*hw;
	hf += (sd[fi1]-sd[fi])*(sd1[fi1]-sd1[fi])*hw;
      }
      mat_pastpast[d][d1]=(hp/sweights-avv_past[d]*avv_past[d1]);
      mat_pastpast[d1][d]=mat_pastpast[d][d1];
      mat_futfut[d][d1]=(hf/sweights-avv_fut[d]*avv_fut[d1]);
      mat_futfut[d1][d]=mat_futfut[d][d1];
#else
	hp += (sd[fi]-sd[fim1])*(sd1[fi]-sd1[fim1]);
	hf += (sd[fi1]-sd[fi])*(sd1[fi1]-sd1[fi]);
      }
      mat_pastpast[d][d1]=num/(num-1.0)*(hp/num-avv_past[d]*avv_past[d1]);
      mat_pastpast[d1][d]=mat_pastpast[d][d1];
      mat_futfut[d][d1]=num/(num-1.0)*(hf/num-avv_fut[d]*avv_fut[d1]);
      mat_futfut[d1][d]=mat_futfut[d][d1];
#endif
    }
}

  /* correlation matrix mat_futpast = C(x_{n+1}-x_{n},x_{n}-x_{n-1}) */
  for (d=0;d<p.DIM;d++) {
    sd=series[d];
    for (d1=0;d1<p.DIM;d1++) {
      sd1=series[d1];
      hp=0.0;
      for (i=0;i<number;i++) {
	fi=f[i];
#if defined(WEIGHTS)
	hp += (sd[fi+1]-sd[fi])*(sd1[fi]-sd1[fi-1])*weight[i];
      }
      mat_futpast[d][d1]=(hp/sweights-avv_fut[d]*avv_past[d1]);
#else
	hp += (sd[fi+1]-sd[fi])*(sd1[fi]-sd1[fi-1]);
      }
      mat_futpast[d][d1]=num/(num-1.0)*(hp/num-avv_fut[d]*avv_past[d1]);
#endif
    }
  }

  /* friction matrix gamma = -C(x_{n+1}-x_{n},x_{n}-x_{n-1}) C^{-1}(x_{n}-x_{n-1},x_{n}-x_{n-1}) */
  imat=invert_matrix(mat_pastpast,p.DIM);

  for (d=0;d<p.DIM;d++) {
    for (d1=0;d1<p.DIM;d1++) {
      hf= mat_futpast[d][0]*imat[0][d1];
      for (k=1;k<p.DIM;k++) {
	hf += mat_futpast[d][k]*imat[k][d1];
      }
      gamma[d][d1]=-hf;
    }
  }

  /* "squared" diffusion matrix now mat_futfut = K K^T
   * = C(x_{n+1}-x_{n},x_{n+1}-x_{n}) - gamma C(x_{n}-x_{n-1},x_{n}-x_{n-1}) gamma^T */ 
  for (i=0;i<p.DIM;i++) {
    for (j=0;j<=i;j++) {
      for (k=0;k<p.DIM;k++) {
	mat_futfut[i][j] += mat_futpast[i][k]*gamma[j][k]; /* gamma = -C(x_{n+1}-x_{n},x_{n}-x_{n-1}) C^{-1}(x_{n}-x_{n-1},x_{n}-x_{n-1}) inserted */
      }
      mat_futfut[j][i]=mat_futfut[i][j];
    }
  }

  for (d=0;d<p.DIM;d++)
    free(imat[d]);
  free(imat);

#ifdef SYM_DIFFOP
  get_diffusion(mat_futfut,p);

/* diffusion matrix dif = K*/
  for (i=0;i<p.DIM;i++)
    for (j=0;j<p.DIM;j++)
      dif[i][j]=mat_futfut[i][j];
#else
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
#endif

  /* drift matrix drift = <x_{n+1}-x_{n}> + gamma*<x_{n} - x_{n-1}> */
  for (d=0;d<p.DIM;d++) {
    drift[d]=av_fut[d]-av_now[d];
    for (k=0;k<p.DIM;k++)
      drift[d] += gamma[d][k]*(av_now[k]-av_past[k]);
  }
}
