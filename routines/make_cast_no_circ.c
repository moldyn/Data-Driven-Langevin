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
//#include <math.h>
#include "olang.h"
//#define TOO_BIG

char isinitmakecastnocirc=0;
double *vecno,*diffs;

void init_vec_no_circ(struct param p)
{
  check_alloc(vecno=(double*)malloc(sizeof(double)*p.DIM));
  check_alloc(diffs=(double*)malloc(sizeof(double)*p.DIM));
  isinitmakecastnocirc=1;
}

void make_cast_no_circ(double **x,double *new,struct param p,double *drift,
		       double **gamma,double **dif)
{
  unsigned int d,d1,hdim,toobig;
  double h;

  if (!isinitmakecastnocirc)
    init_vec_no_circ(p);

  hdim=p.hdim;
  for (d=0;d<p.DIM;d++) {
    h=x[d][hdim]-x[d][hdim-1];
    if (h > M_PI) h=h-pi2;
    if (h < -M_PI) h=h+pi2;
    diffs[d]=h;
  }
  toobig=1;
  while (toobig) {
    toobig=0;
    get_noise(p,vecno);
    for (d=0;d<p.DIM;d++) {
      if (vecno[d] >M_PI) toobig=1;
      if (vecno[d] < -M_PI) toobig=1;
    }
  }

  for (d=0;d<p.DIM;d++) {
    new[d]=x[d][hdim]+drift[d];
    for (d1=0;d1<=d;d1++) {
      new[d] += (dif[d][d1]*vecno[d1]/1.2-gamma[d][d1]*diffs[d1]);
    }
    for (d1=d+1;d1<p.DIM;d1++)
      new[d] -= gamma[d][d1]*diffs[d1];
    if (new[d] > M_PI) {
      new[d]=fmod(new[d],pi2);
      if (new[d] > M_PI) new[d]=new[d]-pi2;
    }
    if (new[d] < -M_PI) {
      new[d]=fmod(new[d],pi2);
      if (new[d] < -M_PI) new[d]=pi2+new[d];
    }
  }
}
