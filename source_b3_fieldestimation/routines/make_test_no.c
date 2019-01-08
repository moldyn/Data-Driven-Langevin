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

void make_test_no(double **x,unsigned long n,double *new,struct param p,double *drift,
	       double **gamma,double **dif)
{
  unsigned int d,d1;
  double help;

  for (d=0;d<p.DIM;d++)
    new[d]=0.0;

  for (d=0;d<p.DIM;d++) {
    help=0.0;
    for (d1=0;d1<d;d1++)
      help += dif[d][d1]*new[d1]-gamma[d][d1]*(x[d1][n]-x[d1][n-1]);
    for (d1=d;d1<p.DIM;d1++)
      help -= gamma[d][d1]*(x[d1][n]-x[d1][n-1]);
    new[d]=(x[d][n+1]-x[d][n]-drift[d]-help)/dif[d][d];
  }
}
