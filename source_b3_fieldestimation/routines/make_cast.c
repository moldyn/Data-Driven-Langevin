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

char isinitmakecast=0;
double *vec;

void init_vec(struct param p)
{
  check_alloc(vec=(double*)malloc(sizeof(double)*p.DIM));
  isinitmakecast=1;
}

void make_cast(double **x,double *new,struct param p,double *drift,
	       double **dif)
{
  unsigned int d,d1;

  if (!isinitmakecast)
    init_vec(p);

  get_noise(p,vec);

  for (d=0;d<p.DIM;d++) {
    new[d]=drift[d];
    for (d1=0;d1<=d;d1++) {
      new[d] += dif[d][d1]*vec[d1];
    }
  }
}
