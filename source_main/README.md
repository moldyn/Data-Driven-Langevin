# Installation
To install the programs adjust the makefile in the main directory and the one 
in the routines subdirectory to your needs. If you are using a linux system, you should
be fine with the default settings. Maybe you have to change the BINDIR variable.


# Programs:
  1. General remarks
  2. General command line switches
  3. The programs
    3.1 ol-second_b1
    3.2 ol-second-tm_b1
    3.3 ol-first_b1
    3.4 ol-first-tm_b1
 
1. General remarks
The programs are written to be used within a shell. There is nothing similar to a
GUI and there will (most probably) never be such a thing. This especially means that,
if you are using Windows or similar, you have to install an environment that allows you
to compile and run the programs. Since I'm not using Windows at all, I can not
(and I will not) assist you at all.

2. General command line switches
The programs use command line options to change parameters. Some of the switches are 
general, in the sense that they are equal for all programs. Some are only used for 
some of the programs. The main general switches are:

	-h	shows a help page with all possible switches available for this program.

	-l num	read 'num' lines from the input file. If 'num' is bigger than the number
		of lines in the file, the file is read to the end. If not set, the whole
		file is read. Lines starting with a '#' are treated as comments are
		ignored.

	-x num	The first 'num' lines of the input file are ignored, regardless of 
		whether these lines are comments or data lines.

	-m num	read 'num' columns from the input files.

	-c str	specifies which columns are read from the input file. If not set,
		the first 'num' (as specified with -m ) columns are read. If you
		would like to change give set 'str' to a comma separated list of
		unsigned integers. Example '-m 2 -c 1,3' means: read columns 1 
		and 3. The order is kept. That means '-c 2,1' reads column 2 as the
		first component and column 1 as the second. If the number given with
		-m is larger than the number of columns given with -c, the columns
		succeeding the largest column given with -c are read in addition.
		This means: '-m 3 -c 3,1' will read columns 3, 1 and 4.

	-M num	Form delay vectors of dimension 'num' in addition to the number
		of components given with -m. Checkout the literature for the
		meaning of delay embedding.

	-d num	The time delay with which the delay vectors are build. Checkout
		the literature for the meaning.

	-k num	This specifies the neighborhood size for the calculation of the
		fields. 'num' means that the neighborhood is chosen that big that
		it exactly contains the 'num' closest spatial neighbors of the
		point to be forecasted.

	-L num	Iterate a trajectory of length 'num'.

	-F num	'num' is the column of the file that specifies whether the actual 
		time step (the row of the input file) has a temporal successor. This
		column of the file should only contain the values '0' and '1'.
		'1' means there is a successor
		'0' means there is none.
		This is used if the input files contains several pieces of trajectories
		glued together.

	-o str	'str' is the name of the output datafile.

	-V num	Verbosity level. The programs may produce some informational messages
		which can be turned off by this switch.

		
3. The programs

The programs whose names contain "second" use the Markovian Langevin equation, the programs
whose names contain "first" use the Overdamped Langevin equation. WE RECOMMENT TO USE THE
PROGRAMS WHICH USE THE MARKOVIAN LANGEVIN EQUATION SINCE ITS BROADER APPLICABLE THAN THE
OVERDAMPED LANGEVIN EQUATION.

3.1 ol-second_b1

This program simply takes the last point of the input datafile and iterates it into
the future, thus creating an artificial trajectory. The program creates two datafiles.
The first two lines of both output files are comments containing the call of the
program with all switches given on the first and column labels on the second line.

The first one simply contains the iterated trajectory consisting of the 'num' columns
given by '-m'. The second file contains the fields. The filename is 
'name of outputfile'.field

The first m columns are the components of the vector where the field was estimated.

The next m columns are the components of the drift field.

The next m*m columns are the components of the friction field

The next m*(m+1)/2 are the relevant components of the diffusion field (the lower
triangular matrix achieved from the Cholesky decomposition).

The last component is the distance of the 'k'-th neighbor from the actual point,
which means it is the neighborhood size used for the local averages.

Additional command line switches:
	-% num	Add noise with absolute amplitude 'num' to the data before	
		starting the analysis. This is meant for cases where the
		data has such low resolution, that many points collapse to
		a single position given rise to ill defined distances.

	-I num	Set the seed for the internal random number generator to
		'num'. This might be useful if a ensemble of trajectories
		shall be created.

3.2 ol-second-tm_b1
This program does not iterate a trajectory but simply creates the fields
at the positions of the existing data. This is done by taking points from
the given input data and fitting the field as if the point was not part of
the input data file.

One output file is created that contains
m+m+m*m+m*(m+1)/2+m+1 columns.

The first m columns are the components of the vector where the field was 
estimated.

The next m columns are the components of the drift field.

The next m*m columns are the components of the friction field

The next m*(m+1)/2 are the relevant components of the diffusion field (the
lower triangular matrix achieved from the Cholesky decomposition).

The next m columns contain the noise which would be needed to reach the
successor of the point, which is known.

The last component is the distance of the 'k'-th neighor from the actual point, 
which means it is the neighborhood size used for the local averages.

Additional command line switches:
	-% num	Add noise with absolute amplitude 'num' to the data before	
		starting the analysis. This is meant for cases where the
		data has such low resolution, that many points collapse to
		a single position given rise to ill defined distances.

	-I num	Set the seed for the internal random number generator to
		'num'. This might be useful if a ensemble of trajectories
		shall be created.


3.3 ol-first_b1

This program simply takes the last point of the input datafile and iterates it into
the future, thus creating an artificial trajectory. The program creates two datafiles.
The first two lines of both output files are comments containing the call of the
program with all switches given on the first and column labels on the second line.

The first one simply contains the iterated trajectory consisting of the 'num' columns
given by '-m'. The second file contains the fields. The filename is 
'name of outputfile'.field

Say the following switches were used '-m m -M M -k k'
Then the file contains m*M+m+m*(m+1)/2+1 columns.

The first m*M columns are the components of the vector where the field was estimated.
The order is: First inputfile column + all its delays, second inputfile column + all
its delays, ...

The next m columns are the components of the drift field.

The next m*(m+1)/2 are the relevant components of the diffusion field (the lower
triangular matrix achieved from the Cholesky decomposition).

The last component is the distance of the 'k'-th neighbor from the actual point,
which means it is the neighborhood size used for the local averages.

Additional command line switches:
	-% num	Add noise with absolute amplitude 'num' to the data before	
		starting the analysis. This is meant for cases where the
		data has such low resolution, that many points collapse to
		a single position given rise to ill defined distances.

	-I num	Set the seed for the internal random number generator to
		'num'. This might be useful if a ensemble of trajectories
		shall be created.
		
3.4 ol-first-tm_b1

This program does not iterate a trajectory but simply creates the fields
at the positions of the existing data. This is done by taking points from
the given input data and fitting the field as if the point was not part of
the input data file.

One output file is created that contains
m*M+m+m*(m+1)/2+m+1 columns.

The first m*M columns are the components of the vector where the field was 
estimated.
The order is: First inputfile column + all its delays, second inputfile
column + all its delays, ...

The next m columns are the components of the drift field.

The next m*(m+1)/2 are the relevant components of the diffusion field (the
lower triangular matrix achieved from the Cholesky decomposition).

The next m columns contain the noise which would be needed to reach the
successor of the point, which is known.

The last component is the distance of the 'k'-th neighor from the actual point, 
which means it is the neighborhood size used for the local averages.

Additional command line switches:
	-% num	Add noise with absolute amplitude 'num' to the data before	
		starting the analysis. This is meant for cases where the
		data has such low resolution, that many points collapse to
		a single position given rise to ill defined distances.

	-I num	Set the seed for the internal random number generator to
		'num'. This might be useful if a ensemble of trajectories
		shall be created.
		
# References

To have a look at the basic ideas and some applications of the data driven
Langevin equation please have a look at 

- Hegger and Stock, J. Chem. Phys. 130, 034106 (2009): here, the approach
is developed based on the overdamped Langevin equation

- Schaudinnus et al.,  Phys. Rev. Lett. 115, 050602 (2015): here, the 
Markovian Langevin equation is used instead of the overdamped equation

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
