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
#include "olang.h"

void set_part(unsigned int *part,char *spart,struct param p)
{
  unsigned int i,j,last;

  sscanf(spart,"%u",&last);
  part[0]=last;
  for (i=1;i<p.DIM;i++) {
    while ((spart[0] != ',') && (spart[0] != '\0')) {spart+=1;}
    if (spart[0] != '\0') spart++;
    if ((unsigned int)strlen(spart) > 0) {
      sscanf(spart,"%u",&last);
      part[i]=last;
    }
    else
      break;
  }
  for (j=i;j<p.DIM;j++)
    part[j]=last;

}

ptree *make_ptree(unsigned int size)
{
  ptree *node;
  unsigned int i;
  check_alloc(node=(ptree*)malloc(sizeof(ptree)));
  node->num=0;
  node->used=0;
  node->which=NULL;
  check_alloc(node->next=(ptree**)malloc(sizeof(ptree*)*size));
  for (i=0;i<size;i++)
    node->next[i]=NULL;

  return node;
}

void fill_tree(ptree *node,unsigned int **x,struct param p,unsigned int act,
	       unsigned int lcount,unsigned int depth,unsigned int *part)
{
  unsigned int comp,n;
  unsigned int next;

  comp=lcount%p.DIM;
  n=lcount/p.DIM;
  next=x[comp][act+n];

  if ((n<(depth-1)) || (comp < (p.DIM-1))) {
    if (node->next[next] == NULL)
      node->next[next]=make_ptree(part[comp]);
    fill_tree(node->next[next],x,p,act,lcount+1,depth,part);
  }
  else {
    if (node->next[next] == NULL)
      node->next[next]=make_ptree(part[comp]);
    node=node->next[next];
    node->which=(unsigned long*)realloc(node->which,
					sizeof(long)*((node->num)+1));
    node->which[node->num]=act;
    node->num += 1;
  }
}

unsigned int read_tree(ptree *node,unsigned int **x,struct param p,
		       unsigned int act,unsigned int lcount,
		       unsigned int depth)
{
  unsigned int comp,n,next;

  comp=lcount%p.DIM;
  n=lcount/p.DIM;
  next=x[comp][act+n];

  if ((n<(depth-1)) || (comp < (p.DIM-1))) {
    return read_tree(node->next[next],x,p,act,lcount+1,depth);
  }
  else {
    node=node->next[next];
    node->used++;
    return node->used;
  }
}

