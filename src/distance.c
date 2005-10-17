/*! \file distance.c
 * \brief all functions requiered for R dist function and C hcluster function.
 *
 *  \date Created: probably in 1995
 *  \date Last modified: Time-stamp: <2005-10-09 13:12:06 antoine>
 *
 *  \author R core members, and lately: Antoine Lucas 
 *
 *  
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1995, 1996  Robert Gentleman and Ross Ihaka
 *  Copyright (C) 1998, 2001  Robert Gentleman, Ross Ihaka and the
 *                            R Development Core Team
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <float.h>
#include <R_ext/Arith.h>
#include <R_ext/Error.h>
#include "mva.h"

/** \brief Distance euclidean (i.e. sqrt(sum of square) )
 */
double R_euclidean(double *x, int nr, int nc, int i1, int i2,int * flag, void ** opt)
{
    double dev, dist;
    int count, j;

    count= 0;
    dist = 0;
    for(j = 0 ; j < nc ; j++) {
	if(R_FINITE(x[i1]) && R_FINITE(x[i2])) {
	    dev = (x[i1] - x[i2]);
	    dist += dev * dev;
	    count++;
	}
	i1 += nr;
	i2 += nr;
    }
    if(count == 0)
      { 
	*flag = 0;
	return NA_REAL;
      }

    if(count != nc) dist /= ((double)count/nc);
    return sqrt(dist);
}

/** \brief Distance maximum (supremum norm)
 */
double R_maximum(double *x, int nr, int nc, int i1, int i2, int * flag, void ** opt)
{
    double dev, dist;
    int count, j;

    count = 0;
    dist = -DBL_MAX;
    for(j = 0 ; j < nc ; j++) {
	if(R_FINITE(x[i1]) && R_FINITE(x[i2])) {
	    dev = fabs(x[i1] - x[i2]);
	    if(dev > dist)
		dist = dev;
	    count++;
	}
	i1 += nr;
	i2 += nr;
    }
    if(count == 0)
      {
	*flag = 0;
	return NA_REAL;
      }
    return dist;
}

/** \brief Distance manhattan
 */
double R_manhattan(double *x, int nr, int nc, int i1, int i2, int * flag, void ** opt)
{
    double dist;
    int count, j;

    count = 0;
    dist = 0;
    for(j = 0 ; j < nc ; j++) {
	if(R_FINITE(x[i1]) && R_FINITE(x[i2])) {
	    dist += fabs(x[i1] - x[i2]);
	    count++;
	}
	i1 += nr;
	i2 += nr;
    }
    if(count == 0)
      {
	*flag = 0;
	return NA_REAL;
      }
    if(count != nc) dist /= ((double)count/nc);
    return dist;
}

/** \brief Distance canberra
 */
double R_canberra(double *x, int nr, int nc, int i1, int i2, int * flag, void ** opt)
{
    double dist, sum, diff;
    int count, j;

    count = 0;
    dist = 0;
    for(j = 0 ; j < nc ; j++) {
	if(R_FINITE(x[i1]) && R_FINITE(x[i2])) {
	    sum = fabs(x[i1] + x[i2]);
	    diff = fabs(x[i1] - x[i2]);
	    if (sum > DBL_MIN || diff > DBL_MIN) {
		dist += diff/sum;
		count++;
	    }
	}
	i1 += nr;
	i2 += nr;
    }
    if(count == 0)
      {
	*flag = 0;
	return NA_REAL;
      }
    if(count != nc) dist /= ((double)count/nc);
    return dist;
}

/** \brief Distance binary
 */
double R_dist_binary(double *x, int nr, int nc, int i1, int i2,int * flag, void ** opt)
{
    int total, count, dist;
    int j;

    total = 0;
    count = 0;
    dist = 0;

    for(j = 0 ; j < nc ; j++) {
	if(R_FINITE(x[i1]) && R_FINITE(x[i2])) {
	    if(x[i1] || x[i2]){
		count++;
		if( ! (x[i1] && x[i2]) ) dist++;
	    }
	    total++;
	}
	i1 += nr;
	i2 += nr;
    }

    if(total == 0)
      {
	*flag = 0;
	return NA_REAL;
      }
    if(count == 0) return 0;
    return (double) dist / count;
}

/** \brief Pearson / Pearson centered (correlation)
 *  \note Added by A. Lucas
 */
double R_pearson(double *x, int nr, int nc, int i1, int i2,int * flag, void ** opt)
{
    double num,sum1,sum2, dist;
    int count,j;

    count= 0;
    num = 0;
    sum1 = 0;
    sum2 = 0;

    for(j = 0 ; j < nc ; j++) {
	if(R_FINITE(x[i1]) && R_FINITE(x[i2])) {
	    num += (x[i1] * x[i2]);
	    sum1 += (x[i1] * x[i1]);
	    sum2 += (x[i2] * x[i2]);
	    count++;
	}
	i1 += nr;
	i2 += nr;
    }
    if(count == 0) 
      {
	*flag = 0;
	return NA_REAL;
      }
    dist = 1 - ( num / sqrt(sum1 * sum2) );
    return dist;
}


/** \brief Distance correlation (Uncentered Pearson)
 *  \note Added by A. Lucas
 */
double R_correlation(double *x, int nr, int nc, int i1, int i2,int * flag, void ** opt)
{
    double num,denum,sumx,sumy,sumxx,sumyy,sumxy;
    int count,j;

    count= 0;
    sumx=0;
    sumy=0;
    sumxx=0;
    sumyy=0;
    sumxy=0;


    for(j = 0 ; j < nc ; j++) {
	if(R_FINITE(x[i1]) && R_FINITE(x[i2])) {
	    sumxy += (x[i1] * x[i2]);
	    sumx += x[i1];
	    sumy += x[i2];
	    sumxx += x[i1] * x[i1];
	    sumyy += x[i2] * x[i2];
	    count++;
	}
	i1 += nr;
	i2 += nr;
    }
    if(count == 0)
      {
	*flag = 0;
	return NA_REAL;
      }
    num = sumxy - ( sumx*sumy /count );
    denum = sqrt( (sumxx - (sumx*sumx /count ) )* (sumyy - (sumy*sumy /count ) ) );
    return 1 - (num / denum);
}

/** \brief Spearman distance (rank base metric)
 *  \note Added by A. Lucas
 */
double R_spearman(double *x, int nr, int nc, int i1, int i2,int * flag, void ** opt)
{
  int j;
  double * data_tri;
  int * order_tri;
  int * rank_tri;
  int n;
  double diffrang=0;

  /* initialisation of variables */
  data_tri = (double * ) opt[0];
  order_tri = (int * )   opt[1];
  rank_tri  = (int * )   opt[2];
  for(j = 0 ; j < nc ; j++) {
    if(!(R_FINITE(x[i1]) && R_FINITE(x[i2])))
      {
	*flag = 0;
	return NA_REAL;	
      }
    order_tri[j] = rank_tri[j] = j;
    order_tri[j+nc] = rank_tri[j+nc] = j;
    data_tri[j] = x[i1];
    data_tri[j+nc] = x[i2];
    i1 += nr;
    i2 += nr;
  }

  n  = nc;
  /* sort and compute rank */
  /* First list */
  rsort_rank_order(data_tri, order_tri,rank_tri, &n);
  /* Second list */
  rsort_rank_order(&(data_tri[nc]),&( order_tri[nc]),&(rank_tri[nc]), &n);

  for(j=0;j<nc;j++)
    {
      diffrang += pow((double) ( rank_tri[j] - rank_tri[j+nc]),2);
    }
  return(  diffrang );

  /*
   * verifcation in R:
   * Dist(x,method='spearman') ; n =dim(x)[2]
   * l=c(x[3,],x[4,]); sum((rank(l[1:n])-rank(l[(n+1):(2*n)]))^2)
   * cor.test(x[3,],x[4,],method="spearm")$statistic
   */

}

enum { EUCLIDEAN=1, MAXIMUM, MANHATTAN, CANBERRA, BINARY ,PEARSON, CORRELATION, SPEARMAN};
/* == 1,2,..., defined by order in the r function dist */


/**
 * R_distance: compute distance. Function called direclty by R
 * \brief compute distance and call one of function R_euclidean or R_...
 * \param x input matrix
 * \param nr,nc number of row and columns
 *        nr individuals with nc values.
 * \param d distance half matrix.
 * \param diag if we compute diagonal of dist matrix (usualy: no).
 * \param method 1, 2,... method used
 * \param ierr error return; 1 good; 0 missing values
 */
void R_distance(double *x, int *nr, int *nc, double *d, int *diag, 
		int *method, int * ierr)
{
    int dc, i, j, ij;
    /* for spearman dist */
    void ** opt;
    double * data_tri;
    int * order_tri;
    int * rank_tri;
    double (*distfun)(double*, int, int, int, int, int *, void **) = NULL;

    switch(*method) {
    case EUCLIDEAN:
	distfun = R_euclidean;
	break;
    case MAXIMUM:
	distfun = R_maximum;
	break;
    case MANHATTAN:
	distfun = R_manhattan;
	break;
    case CANBERRA:
	distfun = R_canberra;
	break;
    case BINARY:
	distfun = R_dist_binary;
	break;
    case PEARSON:
	distfun = R_pearson;
	break;
    case CORRELATION:
	distfun = R_correlation;
	break;
    case SPEARMAN:
        distfun = R_spearman;
	opt = (void * ) malloc (  3 * sizeof(void));
	data_tri  = (double * ) malloc (2*  (*nc) * sizeof(double));
	order_tri  = (int * ) malloc (2 * (*nc) * sizeof(int));
	rank_tri  = (int * ) malloc (2 * (*nc) * sizeof(int));
	if( (data_tri == NULL) || (order_tri == NULL) || (rank_tri == NULL)) 
	  error("distance(): unable to alloc memory");
	opt[0] = (void *) data_tri;
	opt[1] = (void *) order_tri;
	opt[2] = (void *) rank_tri;
	break;

    default:
	error("distance(): invalid distance");
    }


    *ierr = 1; /* res = 1 => no missing values
	          res = 0 => missings values */

    dc = (*diag) ? 0 : 1; /* diag=1:  we do the diagonal */
    ij = 0;
    for(j = 0 ; j <= *nr ; j++)
	for(i = j+dc ; i < *nr ; i++)
	    d[ij++] = distfun(x, *nr, *nc, i, j,ierr,opt);


    if((*method) == SPEARMAN)
      {
	free(opt);
	free(data_tri);
	free(order_tri);
	free(rank_tri);
      }
    /* If spearman: free(vecteur pour tri) */
}



/**
 * Sort, and return order and rank
 * \brief This function is used by R_spearman.
 * order and rank must be initialised with 0:(n-1)
 * make sort, return x sorted
 * order = order(x) - 1 
 * rank = rank(x) -1
 * 
 */

void rsort_rank_order(double *x, int *order, int *rank, int * n)
{
    double v;
    int i, j, h, iv;


    for (h = 1; h <= *n / 9; h = 3 * h + 1);
    for (; h > 0; h /= 3)
	for (i = h; i < *n; i++) {
	    v = x[i]; iv = order[i];
	    j = i;
	    while (j >= h && (x[j - h] > v))
	      {
		x[j] = x[j - h]; order[j] = order[j-h];
		rank[order[j]] = j; 
		j -= h; 
	      }
	    x[j] = v; order[j] = iv; 
	    rank[order[j]] = j; 
	}

}
