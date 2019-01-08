#
#   This file is part of OLANGEVIN
#
#   Copyright (c) 2013- Rainer Hegger
#
#   OLANGEVIN is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; either version 3 of the License, or
#   (at your option) any later version.
#
#   OLANGEVIN is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with OLANGEVIN; if not, see <http://www.gnu.org/licenses/>.
#
SHELL = /bin/sh

AR = ar
ARFLAGS = r
CC = gcc
CFLAGS = -O3
RANLIB = ranlib

ALL = check_alloc.o check_option.o scan_help.o myfgets.o\
      get_multi_series.o search_datafile.o test_outfile.o\
      rescale_data.o rand.o neighbor_search.o get_fields.o\
      make_cast.o make_test.o invert_matrix.o get_fields_no.o\
      make_cast_no.o get_noise.o read_ar_file.o make_test_no.o make_corr.o\
      prune.o scramble.o eigen.o get_local_ar.o get_fields_no_circ.o\
      make_cast_no_circ.o neighbor_search_circ.o neighborhood_info.o

libolang.a: $(ALL)
	$(AR) $(ARFLAGS) libolang.a $?
	$(RANLIB) libolang.a

clean:
	@rm -f *.a *.o *~ #*#
