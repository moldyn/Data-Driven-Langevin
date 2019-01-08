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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "olang.h"
#include "olang_err.h"

#define NEIGH 7
#define EPSILONS 100
#define EPSFAC 2.0

struct s_neighbor {
  unsigned int n;
  long ***box;
  long *list;
} neigh[NEIGH];

double starteps[EPSILONS][EPSILONS];
unsigned long countstarteps[EPSILONS][EPSILONS],minminn;
unsigned int nsseconddim,nsthirddim;

void put_in_boxes_circ(int which,double **x,struct param p,unsigned int *fut)
{
  long n,i,j,k,N;
  unsigned int hdim;
  long ***box,*list;
  double scale;

  hdim=p.hdim;
  box=neigh[which].box;
  N=neigh[which].n;
  scale=(double)N;
  list=neigh[which].list;

  for (i=0;i<N;i++)
    for (j=0;j<N;j++)
      for (k=0;k<N;k++)
	box[i][j][k]=-1;;

  for (n=hdim;n<p.LENGTH;n++) {
    if (fut[n]) {
      i=(long)((x[0][n]+M_PI)/pi2*scale)%N;
      j=(long)((x[nsseconddim][n]+M_PI)/pi2*scale)%N;
      k=(long)((x[nsthirddim][n]+M_PI)/pi2*scale)%N;
      list[n]=box[i][j][k];
      box[i][j][k]=n;
    }
  }
}

void init_neighbor_search_circ(double **x,struct param p,unsigned int *fut)
{
  long n=3;
  int i,j,k;

  if (p.DIM == 1) {
    nsseconddim=0;
    nsthirddim=0;
  }
  else if (p.DIM == 2) {
    nsseconddim=1;
    nsthirddim=1;
  }
  else {
    nsseconddim=1;
    nsthirddim=2;
  }

  for (i=0;i<NEIGH;i++) {
    neigh[i].n=n;
    check_alloc(neigh[i].list=(long*)malloc(sizeof(long)*p.LENGTH));
    check_alloc(neigh[i].box=(long***)malloc(sizeof(long**)*n));
    for (j=0;j<n;j++) {
      check_alloc(neigh[i].box[j]=(long**)malloc(sizeof(long*)*n));
      for (k=0;k<n;k++)
	check_alloc(neigh[i].box[j][k]=(long*)malloc(sizeof(long)*n));
    }
    put_in_boxes_circ(i,x,p,fut);
    n=n*2;
  }
  for (i=0;i<EPSILONS;i++) {
    for (j=0;j<EPSILONS;j++) {
      starteps[i][j]=0.0;
      countstarteps[i][j]=0;
    }
  }
  minminn=p.minminn;
}

unsigned int find_neighbors_circ(double **x,double **y,unsigned int box,
				 struct param p,struct sfound sf,double eps)
{
  long i,j,k,i1,i2,j1,j2,k1,k2,ld,le,x2,y2;
  long element;
  unsigned long nfound=0,*f;
  double hmax,dx,h,eps2,*distance;
  long ***hfbox,*hflist;
  unsigned int hdel,del,dim,emb,N,hdim;

  N=neigh[box].n;
  del=p.DELAY;
  dim=p.DIM;
  emb=p.EMB;
  hfbox=neigh[box].box;
  hflist=neigh[box].list;

  hdim=p.hdim;
  eps*=pi2;
  eps2=eps*eps;

  distance=sf.distance;
  f=sf.found;

  h=(y[0][hdim]+M_PI)/pi2;
  if (h>=1.0) 
    i=N-1;
  else {
    if (h<0.0) i=0;
    else i=((long)(h*N));
  }
  
  h=(y[nsseconddim][hdim]+M_PI)/pi2;
  if (h>=1.0) 
    j=N-1;
  else {
    if (h<0.0) j=0;
    else j=((long)(h*N));
  }
  
  h=(y[nsthirddim][hdim]+M_PI)/pi2;
  if (h>=1.0) 
    k=N-1;
  else {
    if (h<0.0) k=0;
    else k=((long)(h*N));
  }

  for (i1=i-1;i1<=i+1;i1++) {
    i2=i1;
    if (i1 < 0) i2=N-1;
    if (i1 >= N) i2=0;
    for (j1=j-1;j1<=j+1;j1++) {
      j2=j1;
      if (j1 < 0) j2=N-1;
      if (j1 >= N) j2=0;
      for (k1=k-1;k1<=k+1;k1++) {
	k2=k1;
	if (k1 < 0) k2=N-1;
	if (k1 >= N) k2=0;

	element=hfbox[i2][j2][k2];
	while (element != -1) {
	  hmax=0.0;
	  for (le=0;le<emb;le++) {
	    hdel=le*del;
	    x2=element-hdel;
	    y2=hdim-hdel;
#if defined(MAXNORM)
	    for (ld=0;ld<dim;ld++) {
	      dx=fabs(x[ld][x2]-y[ld][y2]);
	      if (dx > M_PI) dx=pi2-dx;
	      if (dx > hmax) {
		hmax=dx;
		if (dx > eps)
		  goto toolarge;
	      }
	    }
#else
	    for (ld=0;ld<dim;ld++) {
	      dx=fabs(x[ld][x2]-y[ld][y2]);
	      if (dx > M_PI) dx=pi2-dx;
	      hmax += dx*dx;
	      if (hmax > eps2)
		goto toolarge;
	    }
#endif
	  }
#if !defined(MAXNORM)
	  hmax=sqrt(hmax);
#endif
	  f[nfound]=element;
	  distance[nfound++]=hmax;
	toolarge:	element=hflist[element];
	}
      }
    }
  }

  return nfound;
}

void search_neighbors_circ(double **x,double **cy,struct param p,struct sfound sf)
{
  int ei,ej,i,hdim,whichsize,whichbox;
  long found;
  double epsilon;

  hdim=p.hdim;

  ei=(int)(cy[0][hdim]*EPSILONS);
  if (ei < 0) ei=0; else if (ei>(EPSILONS-1)) ei=EPSILONS-1;
  ej=(int)(cy[nsseconddim][hdim]*EPSILONS);
  if (ej < 0) ej=0; else if (ej>(EPSILONS-1)) ej=EPSILONS-1;

  if (countstarteps[ei][ej] > 0)
    epsilon=starteps[ei][ej]/countstarteps[ei][ej];
  else
    epsilon=0.001;

  found=0;
  epsilon /= EPSFAC;

  while (found < p.MINN) {
    epsilon *= EPSFAC;
    whichsize=(int)(1.0/epsilon);
    if (whichsize > neigh[NEIGH-1].n)
      whichbox=NEIGH-1;
    else {
      whichbox=0;
      for (i=NEIGH-1;i>0;i--) {
	if ((whichsize > neigh[i-1].n) && (whichsize <= neigh[i].n))
	  whichbox=i-1;
      }
    }
    found=find_neighbors_circ(x,cy,whichbox,p,sf,epsilon);
  }
  sort(found,p,sf);
  if (sf.distance[p.MINN-1] == 0.0) {
    fprintf(stderr,"all neighbors collapse to one point. Maybe add\n"
	    "initial noise to the data\n");
    exit(SEARCH_NEIGHBORS_ZERO_DISTANCE);
  }
  starteps[ei][ej] += sf.distance[p.MINN-1];
  countstarteps[ei][ej]++;
}
