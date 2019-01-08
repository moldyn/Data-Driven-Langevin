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

void neighborhood_info(double **series,struct param p,struct sfound sf,
                double **cast,double *eccentricity,double *abs_ecc,
                double *vratio_fut,double *vratio_past)
{
  unsigned int number;
  long i,d,fi,fi1,fim1;
  double *sd,num,h,sweights,*weight;
  double cur,ss,ss_fut,ss_past;
  unsigned long *f,hdim;

  f=sf.found;

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

  /* eccentricity = <x_{n}> - y_{m} */
  *abs_ecc = 0;
  for (d=0;d<p.DIM;d++) {
    h=0.0;
    sd=series[d];
    cur=cast[d][hdim];
    for (i=0;i<number;i++) {
      h += sd[f[i]]*weight[i];
    }
    h = h/sweights-cur;
    eccentricity[d]=h;
    *abs_ecc += pow(h,2);
  }
  *abs_ecc = sqrt(*abs_ecc);

  /* ratios vratio_fut = <(x_{n+1} - y_{m})^2> / <(x_{n} - y_{m})^2>
       and vratio_past = <(x_{n-1} - y_{m})^2> / <(x_{n} - y_{m})^2> */
  ss=ss_fut=ss_past=0;
  for (d=0;d<p.DIM;d++) {
    sd=series[d];
    cur=cast[d][hdim];
    for (i=0;i<number;i++) {
      fi=f[i];
      fi1=fi+1;
      fim1=fi-1;
      ss += pow(sd[fi]-cur,2)*weight[i];
      ss_fut += pow(sd[fi1]-cur,2)*weight[i];
      ss_past += pow(sd[fim1]-cur,2)*weight[i];
    }
  }
  *vratio_fut = ss_fut/ss;
  *vratio_past = ss_past/ss;
}
