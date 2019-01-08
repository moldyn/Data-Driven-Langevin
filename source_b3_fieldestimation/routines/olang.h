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

#ifndef _OLANG_ROUTINES_H
#define _OLANG_ROUTINES_H

#include <stdio.h>
#include <stdlib.h>

#ifndef _OLANG_ERRORS_H
#include "olang_err.h"
#endif

/* Instead of using the Euclidean Norm, use the Max-Norm */
//#define MAXNORM

//#define WEIGHTS
#define WFACT 0.25

struct param {
  unsigned long LENGTH;
  unsigned int DIM;
  unsigned int EMB;
  unsigned int DELAY;
  unsigned int MINN;
  double SIGMA;
  unsigned int AR_SIZE;
  double *AR_SIGMA;
  double **AR_MOD;
  double **noise_hist;
  unsigned int hdim;
  double maxr;
  unsigned long minminn;
};

struct sfound {
  unsigned long *found;
  double *distance;
  double *weight;
  double *aveps;
  unsigned long *count;
};

typedef struct ptree {
  struct ptree **next;
  unsigned long num;
  unsigned long used;
  unsigned long *which;
} ptree;

#define sqr(x) ((x)*(x))

#ifdef __cplusplus
extern "C" {
#endif

/* The possible names of the verbosity levels */
#define VER_INPUT 0x1
#define VER_USR1 0x2
#define VER_USR2 0x4
#define VER_USR3 0x8
#define VER_USR4 0x10
#define VER_USR5 0x20
#define VER_USR6 0x40
#define VER_FIRST_LINE 0x80

#define INPUT_SIZE 1024

  extern double **get_multi_series(char *,unsigned long *,unsigned long,
				   unsigned int *,char *,char,unsigned int);
  extern char* check_option(char **in,int n,int which,int type);
  extern int scan_help(int n,char **in);
  extern char* myfgets(char *str,int *size,FILE *fin,unsigned int verbosity);
  extern char* search_datafile(int n,char **names,unsigned int *col,
			       unsigned int verbosity);
  extern char* test_outfile(char *name);
  extern void rescale_data(double **ser,struct param par,
			   double *min,double *inter);
  extern void check_alloc(void *pnt);

  extern void get_fields(double **series,struct param p,struct sfound sf,
			 double *drift,double **dif);
  extern void make_cast(double **x,double *new,struct param p,double *drift,
			double **dif);
  extern void make_test(double **x,unsigned long n,double *new,struct param p,
			double *drift,double **dif);
  extern void read_ar_file(struct param *p,char *fname);
  extern void get_fields_no(double **series,struct param p,struct sfound sf,
			    double *drift,double **gamma,double **dif,
			    double **x);
  extern void make_cast_no(double **x,double *new,struct param p,double *drift,
			   double **gamma,double **dif);
  extern void make_test_no(double **x,unsigned long n,double *new,struct param p,
			   double *drift,double **gamma,double **dif);
  extern void make_corr(double **x,unsigned int *fut,
			struct param p,double *sigma);
  extern void neighborhood_info(double **series,struct param p,struct sfound sf,
			double **cast,double *eccentricity,double *abs_ecc,
                        double *vratio_fut,double *vratio_past);


  /* routines from rand.c */
  extern void rnd_init(unsigned long);
  extern unsigned long rnd_long();
  extern unsigned long rnd_1279();
  extern unsigned long rnd69069();
  extern double gaussian(double);

  /* routines from neighbor_search.c */
  extern void put_in_boxes(int which,double **x,struct param p,
			   unsigned int *fut);
  extern void init_neighbor_search(double **x,struct param p,
				   unsigned int *fut);
  extern void sort(unsigned long nfound,struct param p,struct sfound sf);
  extern unsigned int find_neighbors(double **x,double **y,unsigned int box,
				     struct param p,struct sfound sf,
				     double eps);
  extern void search_neighbors(double **x,double **y,struct param p,
				 struct sfound sf);

  /*routines from invert_matrix.c */
  extern void solvele(double **mat,double *vec,unsigned int n);
  extern double **invert_matrix(double **mat,unsigned int size);

  /* routines from get_noise.c */
  extern void init_noise(struct param p,unsigned long seed);
  extern void get_noise(struct param p,double *noi);
  extern void get_ar_noise(struct param p,double *noi);

  /* routines from prune.c */
  extern void set_part(unsigned int* part,char *spart,struct param p);
  extern ptree *make_ptree(unsigned int size);
  extern void fill_tree(ptree *node,unsigned int **x,struct param p,
			unsigned int act,unsigned int count,
			unsigned int depth,unsigned int *part);
  extern unsigned int read_tree(ptree *node,unsigned int **x,struct param p,
				unsigned int act,unsigned int count,
				unsigned int depth);

  /* routines from scramble.c */
  extern void lscramble(long *arr,unsigned long size);

  /* routines from eigen.c */
  extern void my_eigen(double **mat,unsigned long dim,double *eig);

  /* routines from get_local_ar.c */
  extern void get_ar_coeff(double **series,struct param p,struct sfound sf,
			   double *av,double **coef,double *dif);
  extern void make_ar_cast(double **x,double *new,struct param p,double *av,
			   double **coeff,double *dif);
#ifdef __cplusplus
}
#endif

#endif
