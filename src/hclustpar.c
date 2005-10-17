/** \file hclustpar.c
 * \brief Parallelized version of hcluster (use distpar routine)
 *   
 *                                                          
 *   \date Created       : 14/11/02 
 *   \date Last Modified : Time-stamp: <2005-11-12 15:54:58 antoine>
 *
 * \author Antoine Lucas.
 */


#include <stdlib.h>
#include <math.h>
#include "mva.h"


#define max( A , B )  ( ( A ) > ( B ) ? ( A ) : ( B ) )
#define min( A , B )  ( ( A ) < ( B ) ? ( A ) : ( B ) )

/** Paralelized hierarchical clustering
 * \brief allocate distance matrix execute function R_distancepar, launch 
 * hclust on this distance matrix, and free memory of distance matrix.
 * \param x: data nr x nc
 * \param nr,nc number of row and columns		  
 * \param membr: member, vector  1:nr		  
 * \param method  integer -> distance method	  
 * \param diag integer: if we compute diagonal in distance matrix (usually no)
 * \param iopt    integer -> link used 
 * \param nbprocess nb of process for parallelization
 * \param ia, ib   result (merge)			  
 * \param iorder  result (order)			  
 * \param crit    result (height)			  
 * \param result  flag  0 => correct		  
 *               1 => Pb			  
 *               2 => Cannot allocate memory 
 *               3 => Pb with distance matrix
 */
void hclusterpar(double *x, int *nr, int *nc, int *diag, int *method, int *iopt ,int *ia , int *ib,int *iorder,double *crit,double *membr,int *nbprocess,int * result)
{
  int  len;
  int flag;
  double *d;

  *result = 1;

#ifndef __MINGW_H

  len = (*nr * (*nr-1) )/2; 
  d = (double*) malloc (len * sizeof(double));
  if(d == NULL )
    {
      printf("AMAP: Not enought system memory to allocate matrix distance\n"); 
      *result = 2;
      return;
    }
  /*
   * Calculate d: distance matrix
   */
  R_distancepar(x,nr,nc,d,diag,method,nbprocess,&flag);
  if(flag == 0)
    {
      printf("AMAP: Unable to compute Hierarchical Clustering: missing values in distance matrix\n"); 
      *result = 3;
      return;
    }

  /*
   *  Hierarchical clustering
   */
  hclust(nr,&len,iopt ,ia ,ib,iorder,crit,membr,d,result);
  free(d); 

  *result = 0;

#else
  /* No threads on Windows yet */
  hcluster(x, nr, nc, diag, method, iopt ,ia ,ib,iorder, crit, membr, result);
  
#endif


}

