/*
 *   This file is part of OLANGEVIN
 *
 *   Copyright (c) 2013- Rainer Hegger
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include "olang.h"

#define SIZE 2048

void lscramble(long *arr,unsigned long size)
{
  long i,j,k,m;
  unsigned long rnd,rndf,hlength,allscr=0;
  long *scfound,*scnhelp,*ret,scnfound;
  long scbox[SIZE],lswap,element,scbox1=SIZE-1;
  double *rz,*schelp,swap;
  
  rnd_init(0x45adf);

  check_alloc(rz=(double*)malloc(sizeof(double)*size));
  check_alloc(ret=(long*)malloc(sizeof(long)*size));
  check_alloc(scfound=(long*)malloc(sizeof(long)*size));
  check_alloc(scnhelp=(long*)malloc(sizeof(long)*size));
  check_alloc(schelp=(double*)malloc(sizeof(double)*size));

  for (i=0;i<size;i++)
    rz[i]=(double)(rnd_long())/ULONG_MAX;
  
  for (i=0;i<SIZE;i++)
    scbox[i]= -1;
  for (i=0;i<size;i++) {
    m=(int)(rz[i]*(double)SIZE)&scbox1;
    scfound[i]=scbox[m];
    scbox[m]=i;
  }
  for (i=0;i<SIZE;i++) {
    scnfound=0;
    element=scbox[i];
    while(element != -1) {
      scnhelp[scnfound]=element;
      schelp[scnfound++]=rz[element];
      element=scfound[element];
    }
    
    for (j=0;j<scnfound-1;j++)
      for (k=j+1;k<scnfound;k++)
	if (schelp[k] < schelp[j]) {
	  swap=schelp[k];
	  schelp[k]=schelp[j];
	  schelp[j]=swap;
	  lswap=scnhelp[k];
	  scnhelp[k]=scnhelp[j];
	  scnhelp[j]=lswap;
	}
    for (j=0;j<scnfound;j++)
      ret[allscr+j]=arr[scnhelp[j]];
    allscr += scnfound;
  }
  for (i=0;i<size;i++)
    arr[i]=ret[i];

  free(rz);
  free(ret);
  free(scfound);
  free(schelp);
  free(scnhelp);
}
