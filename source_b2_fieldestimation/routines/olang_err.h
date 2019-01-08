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

/* These definitions give the exit codes for the C part of the olangevin 
   package.
   Typically the name is build up of, first, the name of the routine creating
   the exception, secondly, sort of an description of the exception.
   */

#ifndef _OLANG_ERR_H
#define _OLANG_ERR_H

#define CHECK_OPTION_NOT_UNSIGNED 1
#define CHECK_OPTION_NOT_INTEGER 2
#define CHECK_OPTION_NOT_FLOAT 3
#define CHECK_OPTION_NOT_TWO 4
#define CHECK_OPTION_NOT_THREE 5
#define CHECK_OPTION_C_NO_VALUE 6

#define CHECK_ALLOC_NOT_ENOUGH_MEMORY 7

#define GET_MULTI_SERIES_INPUT_SIZE 8
#define GET_MULTI_SERIES_WRONG_TYPE_OF_C 9
#define GET_MULTI_SERIES_NO_LINES 10

#define TEST_OUTFILE_NO_WRITE_ACCESS 11

#define LANGEVIN_MAIN_TOO_FEW_MINN 12
#define LANGEVIN_MAIN_NO_INPUTFILE 13

#define GET_FIELDS_CHOLESKY 14

#define SEARCH_NEIGHBORS_ZERO_DISTANCE 15

#define SOLVELE_SINGULAR_MATRIX 16

#define GET_FIELDS_NO_CHOLESKY 17

#define READ_AR_FILE_NO_NAME 18
#define READ_AR_FILE_NOT_EXIST 19
#define READ_AR_FILE_NOT_ENOUGH_LINES 20

#define MAKE_CORR_NEGATIVE_SIGMA 21

#define LANGEVIN_MAIN_C_TOO_LARGE 22

#define EIG2_TOO_MANY_ITERATIONS 23

#define TICA_DELAY_TOO_LARGE 24
#endif
