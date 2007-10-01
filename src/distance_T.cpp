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


#define _AMAP_DISTANCE_TEMPLATE_CPP 1

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "distance_T.h"
#include "distance.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <R_ext/Arith.h>
#include <R_ext/Error.h>
#include <limits>

#ifndef __MINGW_H
#include <pthread.h>
#endif

#define MAX( A , B )  ( ( A ) > ( B ) ? ( A ) : ( B ) )
#define MIN( A , B )  ( ( A ) < ( B ) ? ( A ) : ( B ) )


/** \brief Distance euclidean (i.e. sqrt(sum of square) )
 */
template<class T> T  distance_T<T>::R_euclidean(double * x, double * y , int nr_x, int nr_y, int nc, 
						int i1, int i2,
						int * flag, void ** opt)
{
  T dev, dist;
  int count, j;

  count= 0;
  dist = 0;
  for(j = 0 ; j < nc ; j++) {
    if(R_FINITE(x[i1]) && R_FINITE(y[i2])) {
      dev = (x[i1] - y[i2]);
      dist += dev * dev;
      count++;
    }
    i1 += nr_x;
    i2 += nr_y;
  }
  if(count == 0)
    { 
      *flag = 0;
      return NA_REAL;
    }

  if(count != nc) dist /= ((T)count/nc);
  return sqrt(dist);
}

/** \brief Distance maximum (supremum norm)
 */
template<class T> T  distance_T<T>::R_maximum(double * x, double * y , int nr_x, int nr_y, int nc, 
					      int i1, int i2,
					      int * flag, void ** opt)
{
  T dev, dist;
  int count, j;

  count = 0;
  dist = std::numeric_limits<T>::min();
  for(j = 0 ; j < nc ; j++) {
    if(R_FINITE(x[i1]) && R_FINITE(y[i2])) {
      dev = fabs(x[i1] - y[i2]);
      if(dev > dist)
	dist = dev;
      count++;
    }
    i1 += nr_x;
    i2 += nr_y;
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
template<class T> T  distance_T<T>::R_manhattan(double * x, double * y , int nr_x, int nr_y, int nc, 
						int i1, int i2,
						int * flag, void ** opt)
{
  T dist;
  int count, j;

  count = 0;
  dist = 0;
  for(j = 0 ; j < nc ; j++) {
    if(R_FINITE(x[i1]) && R_FINITE(y[i2])) {
      dist += fabs(x[i1] - y[i2]);
      count++;
    }
    i1 += nr_x;
    i2 += nr_y;
  }
  if(count == 0)
    {
      *flag = 0;
      return NA_REAL;
    }
  if(count != nc) dist /= ((T)count/nc);
  return dist;
}

/** \brief Distance canberra
 */
template<class T> T  distance_T<T>::R_canberra(double * x, double * y , int nr_x, int nr_y, int nc, 
					       int i1, int i2,
					       int * flag, void ** opt)
{
  T dist, sum, diff;
  int count, j;

  count = 0;
  dist = 0;
  for(j = 0 ; j < nc ; j++) {
    if(R_FINITE(x[i1]) && R_FINITE(y[i2])) {
      sum = fabs(x[i1] + y[i2]);
      diff = fabs(x[i1] - y[i2]);
      if (sum > DBL_MIN || diff > DBL_MIN) {
	dist += diff/sum;
	count++;
      }
    }
    i1 += nr_x;
    i2 += nr_y;
  }
  if(count == 0)
    {
      *flag = 0;
      return NA_REAL;
    }
  if(count != nc) dist /= ((T)count/nc);
  return dist;
}

/** \brief Distance binary
 */
template<class T> T  distance_T<T>::R_dist_binary(double * x, double * y , int nr_x, int nr_y, int nc, 
						  int i1, int i2,
						  int * flag, void ** opt)
{
  int total, count, dist;
  int j;

  total = 0;
  count = 0;
  dist = 0;

  for(j = 0 ; j < nc ; j++) {
    if(R_FINITE(x[i1]) && R_FINITE(y[i2])) {
      if(x[i1] || y[i2]){
	count++;
	if( ! (x[i1] && y[i2]) ) dist++;
      }
      total++;
    }
    i1 += nr_x;
    i2 += nr_y;
  }

  if(total == 0)
    {
      *flag = 0;
      return NA_REAL;
    }
  if(count == 0) return 0;
  return (T) dist / count;
}

/** \brief Pearson / Pearson centered (correlation)
 *  \note Added by A. Lucas
 */
template<class T> T  distance_T<T>::R_pearson(double * x, double * y , int nr_x, int nr_y, int nc, 
					      int i1, int i2,
					      int * flag, void ** opt)
{
  T num,sum1,sum2, dist;
  int count,j;

  count= 0;
  num = 0;
  sum1 = 0;
  sum2 = 0;

  for(j = 0 ; j < nc ; j++) {
    if(R_FINITE(x[i1]) && R_FINITE(y[i2])) {
      num += (x[i1] * y[i2]);
      sum1 += (x[i1] * x[i1]);
      sum2 += (y[i2] * y[i2]);
      count++;
    }
    i1 += nr_x;
    i2 += nr_y;
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
template<class T> T  distance_T<T>::R_correlation(double * x, double * y , int nr_x, int nr_y, int nc, 
						  int i1, int i2,
						  int * flag, void ** opt)
{
  T num,denum,sumx,sumy,sumxx,sumyy,sumxy;
  int count,j;

  count= 0;
  sumx=0;
  sumy=0;
  sumxx=0;
  sumyy=0;
  sumxy=0;


  for(j = 0 ; j < nc ; j++) {
    if(R_FINITE(x[i1]) && R_FINITE(y[i2])) {
      sumxy += (x[i1] * y[i2]);
      sumx += x[i1];
      sumy += y[i2];
      sumxx += x[i1] * x[i1];
      sumyy += y[i2] * y[i2];
      count++;
    }
    i1 += nr_x;
    i2 += nr_y;
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
template<class T> T  distance_T<T>::R_spearman(double * x, double * y , int nr_x, int nr_y, int nc, 
					       int i1, int i2,
					       int * flag, void ** opt)
{
  int j;
  double * data_tri;
  int * order_tri;
  int * rank_tri;
  int n;
  T diffrang=0;

  // initialisation of variables they are 3 vectors of size 2*nc 
  data_tri = (double * ) opt[0];
  order_tri = (int * )   opt[1];
  rank_tri  = (int * )   opt[2];
  for(j = 0 ; j < nc ; j++) {
    if(!(R_FINITE(x[i1]) && R_FINITE(y[i2])))
      {
	*flag = 0;
	return NA_REAL;	
      }
    order_tri[j] = rank_tri[j] = j;
    order_tri[j+nc] = rank_tri[j+nc] = j;
    data_tri[j] = x[i1];
    data_tri[j+nc] = y[i2];
    i1 += nr_x;
    i2 += nr_y;
  }

  n  = nc;
  /* sort and compute rank */
  /* First list */
  rsort_rank_order(data_tri, order_tri,rank_tri, &n);
  /* Second list */
  rsort_rank_order(&(data_tri[nc]),&( order_tri[nc]),&(rank_tri[nc]), &n);

  for(j=0;j<nc;j++)
    {
      diffrang += pow((T) ( rank_tri[j] - rank_tri[j+nc]),2);
    }
  return(  diffrang );

  /*
   * verification in R:
   * Dist(x,method='spearman') ; n =dim(x)[2]
   * l=c(x[3,],x[4,]); sum((rank(l[1:n])-rank(l[(n+1):(2*n)]))^2)
   * cor.test(x[3,],x[4,],method="spearm")$statistic
   */

}

/**
 * R_distance: compute parallelized distance. Function called direclty by R
 * \brief compute distance and call function thread_dist
 * that call one of function R_euclidean or R_...
 * \param x input matrix
 * \param nr,nc number of row and columns
 *        nr individuals with nc values.
 * \param d distance half matrix.
 * \param diag if we compute diagonal of dist matrix (usualy: no).
 * \param method 1, 2,... method used
 * \param nbprocess: number of threads to create
 * \param ierr error return; 1 good; 0 missing values
 */
template<class T> void  distance_T<T>::distance(double *x, int *nr,
						int *nc, T *d, int *diag, 
						int *method,int *nbprocess, 
						int * ierr)
{



  int  i;
  T_argument * arguments;

  bool dc = (*diag) ? 0 : 1; /* diag=1:  we do the diagonal */ 
  
  /*
   * Arguments sent to thread (adress):
   * number of thread
   * nr
   * nc 
   * dc
   * *x
   * *d
   * *method
   * *ierr
   */ 

  arguments = (T_argument * ) malloc ((*nbprocess) * sizeof( T_argument ));


  //printf("nb processs %d\n",*nbprocess);

  for(i=0; i< *nbprocess; ++i)
    {
      arguments[i].id =i;
      arguments[i].x=x;
      arguments[i].nr = nr;
      arguments[i].nc = nc;
      arguments[i].dc = dc;
      arguments[i].d = d;
      arguments[i].method = method;
      arguments[i].nbprocess= *nbprocess;
      arguments[i].ierr=ierr;
    }
  *ierr = 1; /* res = 1 => no missing values
		res = 0 => missings values */


#ifndef __MINGW_H
  pthread_t * th = (pthread_t *) malloc ( *nbprocess * sizeof(pthread_t));

  for (i=0; i < *nbprocess ; i++)
    {
      pthread_create(th+i,0,distance_T<T>::thread_dist,(void *)(arguments+i));
    }

  /* Attends la fin    */
  for (i=0; i < *nbprocess ; i++)
    {      
      pthread_join(*(th+i),NULL);
    }      
  free( th);

#else

  // p_thread not yet used on windows.

  arguments[0].nbprocess = 1;
  thread_dist((void *)arguments);
#endif

  free( arguments );





}


/** thread_dist function that compute distance.
 *
 */
template <class T> void* distance_T<T>::thread_dist(void* arguments_void)
{

  int nbprocess,nr,nc,i,j,dc,ij;
  T_argument * arguments = static_cast<T_argument*>(arguments_void); 
  T * d;
  double * x;
  int * method;
  int * ierr;
  /* for spearman dist */
  void *opt[3] ;
  double * data_tri = 0;
  int * order_tri = 0;
  int * rank_tri = 0;
  T (*distfun)(double*,double*,int, int, int, int, int, int *, void **) = NULL;


  short int no = arguments[0].id;
  nr = *arguments[0].nr;
  nc = *arguments[0].nc;
  dc = arguments[0].dc;
  x  = arguments[0].x;
  d  = arguments[0].d;
  method =  arguments[0].method;
  nbprocess = arguments[0].nbprocess;
  ierr =  arguments[0].ierr;

    
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
    data_tri  = (double * ) malloc (2*  (nc) * sizeof(double));
    order_tri  = (int * ) malloc (2 * (nc) * sizeof(int));
    rank_tri  = (int * ) malloc (2 * (nc) * sizeof(int));
    if( (data_tri == NULL) || (order_tri == NULL) || (rank_tri == NULL)) 
      error("distance(): unable to alloc memory");
    opt[0] = (void *) data_tri;
    opt[1] = (void *) order_tri;
    opt[2] = (void *) rank_tri;
    break;

  default:
    error("distance(): invalid distance");
  }
    

  /*
    debut = ((nr+1) / nbprocess + 1 ) * no ;
    fin =  min ( ((nr+1) / nbprocess + 1) * ( no + 1 ) , (nr+1));
  */

  /* debut des boucles 0
     fin: nr+1 */

    
  int debut = (int) floor( ((nr+1.)*nbprocess - sqrt( (nr+1.)*(nr+1.) * nbprocess * nbprocess - (nr+1.)*(nr+1.) * nbprocess * no  ) )/nbprocess);
  int fin = (int) floor(((nr+1.)*nbprocess - sqrt( (nr+1.)*(nr+1.) * nbprocess * nbprocess - (nr+1.)*(nr+1.) * nbprocess * (no+1.)  ) )/nbprocess);

    
  //printf("Thread %d debut %d fin %d\n",no,debut,fin);    


  // here: the computation !
  //    for(j = 0 ; j <= nr ; j++)
  for(j = debut ; j < fin ; j++)
    {
      ij = (2 * (nr-dc) - j +1) * (j) /2 ;
      for(i = j+dc ; i < nr ; i++)
	{
	  d[ij++] = distfun(x,x,nr, nr, nc, i, j,ierr,opt);
	}
    }

  if((*method) == SPEARMAN)
    {
      free(data_tri);
      free(order_tri);
      free(rank_tri);
    }

  return (void*)0;
    
}









/**
 * R_distance_kms: compute distance between individual i1 and
 * centroid i2
 * \brief compute distance and call one of function R_euclidean or R_...
 * \brief This function is called by kmeans_Lloyd2 
 * \param x input matrix (individuals)
 * \param y input matrix (centroids)
 * \param nr1,nr2,nc number of row (nr1:x, nr2:y) and columns
 *        nr individuals with nc values.
 * \param i1, i2: indice of individuals (individual i1, centroid i2)
 * \param method 1, 2,... method used
  * \param ierr for NA 0 if no value can be comuted due to NA
  * \param opt optional parameter send to spearman dist.
 */
template <class T> T distance_T<T>::distance_kms(double *x,double *y, int nr1,int nr2, int nc,int i1,int i2, int *method, 
						 int * ierr, void ** opt)
{
  /*
   * compute distance x[i1,*] - y[i2,*]
   * x matrix n x p
   * y matrix m x p
   * nr1 = n; nr2 = m; nc =p
   */
  
  T res;

  T (*distfun)(double*,double*,int, int, int, int, int, int *, void **) = NULL;

  
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
    break;

  default:
    error("distance(): invalid distance");
  }
    
  res = distfun(x,y, nr1,nr2, nc, i1, i2,ierr, opt);
  return( res);
}
