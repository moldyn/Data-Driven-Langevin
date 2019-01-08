# Quick Start

0) Download the files in source_main to some folder datadrivenlangevin
1) Create inside datadrivenlangevin the subfolder bin
2) Adjust the Makefile inside datadrivenlangevin and inside datadrivenlangevin/routines
(was downloaded) as you need it
3) Run make clean in the shell
4) Run make install in the shell

After this you get variouse programs inside of datadrivenlangevin/bin. Important are:

- ol-second_b1: generates a Langevin model trajectory based on the Markovian Langevin
equation by estimating the Langevin field from some data set, see for example
Schaudinnus et al.,  Phys. Rev. Lett. 115, 050602 (2015)

- ol-second-tm_b1: generates a so-called "testmodel"-run which means that fields and 
noise which would have produced a given trajectory if it would have been produced based
on the Markovian langevin equation are calculated. 

- ol-first_b1: does the same as ol-second_b1 but uses the OVERDAMPED Langevin equation
as model, see Hegger and Stock, J. Chem. Phys. 130, 034106 (2009).

- ol-second-tm_b1: the same as ol-second-tm_b1 but uses the OVERDAMPED Langevin equation
as model.

# Description of Source Folders

source_main contains the normally used c-code. Stick to it for the start.

source_b2_fieldestimation contains the code for an alternative estimation of the fields

source_b3_fieldestimation contains the code for another alternative estimation of the fields

FOR FURTHER INFOMATION PLEASE HAVE A LOOK AT THE READMES OF THE DIFFERENT source_... FOLDERS.

# Licensing
This file is part of OLANGEVIN

Copyright (c) 2013 Rainer Hegger

OLANGEVIN is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

OLANGEVIN is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with OLANGEVIN; if not, see <http://www.gnu.org/licenses/>.
