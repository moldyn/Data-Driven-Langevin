#
#   This file is part of OLANGEVIN
#
#   Copyright (c) 2013 Rainer Hegger
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
prefix = /home/bl20/Arbeitsfläche/dLE2backup/olangevin_b2
BINDIR = ${prefix}/bin
SUFFIX = "_b2"

CC = gcc
CFLAGS = -O3
AR = ar
ARFLAGS = r
INSTALL = /usr/bin/install -c

LOADLIBS = routines/libolang.a -lm

# list of executables we want to produce
   ALL =   ol-first ol-first-tm ol-second ol-second-tm ol-first-temp \
	ol-second-temp ol-first-colno ol-second-colno ol-check-xdot2 \
	ol-armodel ol-prune ol-second-colno-temp ol-first-colno-temp \
        ol-search-neighbors ol-local-ar ol-second-circ

all: $(ALL) 

routines/libolang.a:
	(cd routines && $(MAKE))

$(ALL): routines/libolang.a *.c
	-$(CC) $(CFLAGS) $(COPTS) -o $@ $@.c $(LOADLIBS)

install: all
	-for bin in $(ALL); do $(INSTALL) $$bin $(BINDIR)/$$bin$(SUFFIX); done

clean:
	@rm -f *.o *~ #*#
	@rm -f $(ALL)
	-(cd routines && $(MAKE) clean)

uninstall:
	-@for bin in $(ALL); do rm -f $(BINDIR)/$$bin; done
