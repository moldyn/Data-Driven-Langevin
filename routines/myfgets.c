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
#include <string.h>
#include "olang.h"

char* myfgets(char *str,int *size,FILE *fin,unsigned int verbosity)
{
  char *ret;
  char *hstr=NULL;
  char last;

  ret=fgets(str,*size,fin);
  if (ret == NULL)
    return NULL;

  last=str[strlen(str)-1];

  while (last != '\n') {
    *size += INPUT_SIZE;
    check_alloc(hstr=(char*)calloc((size_t)INPUT_SIZE,(size_t)1));
    check_alloc(str=realloc(str,(size_t)*size));
    ret=fgets(hstr,INPUT_SIZE,fin);
    strcat(str,hstr);
    if (verbosity&VER_INPUT)
      fprintf(stderr,"Line in file too long. Increasing input size\n");
    last=str[strlen(str)-1];
    free(hstr);
  }

  if (ret == NULL)
    return NULL;
  else
    return str;
}
