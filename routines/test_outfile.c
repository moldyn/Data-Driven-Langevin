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
#include <unistd.h>
#include "olang.h"
#include "olang_err.h"

char *test_outfile(char *name)
{
  char *nname;
  unsigned int len,i=1;
  FILE *file;

  if (access(name,F_OK) != -1) {
    len=strlen(name);
    check_alloc(nname=(char*)malloc(sizeof(char)*(len+4)));
    sprintf(nname,"%s",name);
    while (access(nname,F_OK) != -1) {
      sprintf(nname,"%s.%u",name,i);
      i++;
    }
    if (i > 1) {
      free(name);
      name=nname;
    }
  }
  file=fopen(name,"a");
  if (file == NULL) {
    fprintf(stderr,"Couldn't open %s for writing. Exiting\n",name);
    exit(TEST_OUTFILE_NO_WRITE_ACCESS);
  }
  fclose(file);

  return name;
}
