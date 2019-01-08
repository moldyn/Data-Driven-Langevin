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

#include "olang.h"

void rescale_data(double **data,struct param par,
		  double *min,double *inter)
{
  unsigned int i;
  unsigned long n;
  double *ser,hmin,h,hmax;

  for (i=0;i<par.DIM;i++) {
    ser=data[i];
    hmin=ser[0];
    hmax=ser[0];
    for (n=1;n<par.LENGTH;n++) {
      h=ser[n];
      if (hmin > h)
	hmin=h;
      else
	if (hmax < h)
	  hmax=h;
    }
    hmax -= hmin;
    for (n=0;n<par.LENGTH;n++)
      ser[n] -= hmin;
    min[i]=hmin;
    inter[i]=hmax;
  }
}
