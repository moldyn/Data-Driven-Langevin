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

void get_ar_noise(struct param p,double *noi)
{
  unsigned int d,j;

  for (d=0;d<p.DIM;d++)
    noi[d]=gaussian(p.AR_SIGMA[d]);
  if (p.AR_SIZE > 0) {
    for (d=0;d<p.DIM;d++) {
      for (j=0;j<p.AR_SIZE;j++)
	noi[d] += p.AR_MOD[d][j]*p.noise_hist[d][j];
      for (j=p.AR_SIZE-1;j>=1;j--)
	p.noise_hist[d][j]=p.noise_hist[d][j-1];
      p.noise_hist[d][0]=noi[d];
      noi[d] *= p.SIGMA;
    }
  }
}

void init_noise(struct param p,unsigned long seed)
{
  unsigned long i,j;
  double *noi;

  rnd_init(seed);

  for (i=0;i<10000;i++) {
    rnd_long();
    gaussian(1.0);
  }

  if (p.AR_SIZE > 0) {
    check_alloc(noi=(double*)malloc(sizeof(double)*p.DIM));

    for (i=0;i<p.DIM;i++)
      for (j=0;j<p.AR_SIZE;j++)
	p.noise_hist[i][j]=gaussian(p.AR_SIGMA[i]);

    for (i=0;i<10000;i++)
      get_ar_noise(p,noi);

    free(noi);
  }
}

void get_noise(struct param p,double *noi)
{
  unsigned int d;

  if (p.AR_SIZE == 0) {
    for (d=0;d<p.DIM;d++)
      noi[d]=gaussian(p.SIGMA);
  }
  else {
    get_ar_noise(p,noi);
  }
}
