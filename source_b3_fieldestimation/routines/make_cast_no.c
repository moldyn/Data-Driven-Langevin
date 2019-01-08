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
#include "olang.h"
//#define TOO_BIG

char isinitmakecastno=0;
double *vecno;

void init_vec_no(struct param p)
{
  check_alloc(vecno=(double*)malloc(sizeof(double)*p.DIM));
  isinitmakecastno=1;
}

void make_cast_no(double **x,double *new,struct param p,double *drift,
	       double **gamma,double **dif)
{
  unsigned int d,d1,hdim;
  char too_big;

  if (!isinitmakecastno)
    init_vec_no(p);

  do {
    get_noise(p,vecno);

    hdim=p.hdim;

    for (d=0;d<p.DIM;d++) {
      new[d]=x[d][hdim]+drift[d];
      for (d1=0;d1<=d;d1++) {
	new[d] += (dif[d][d1]*vecno[d1]-gamma[d][d1]*
		   (x[d1][hdim]-x[d1][hdim-1]));
      }
      for (d1=d+1;d1<p.DIM;d1++)
	new[d] -= gamma[d][d1]*(x[d1][hdim]-x[d1][hdim-1]);
    }
#if defined(TOO_BIG)
    too_big=(new[0]<0.0)||(new[0]>1.0);
#else
    too_big=0;
#endif
  } while (too_big);
}
