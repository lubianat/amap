/*! \file : kmeans.c
 * 
 *
 * \brief  K-means clustering 
 *
 * \date Created       : before 2005
 * \date Last Modified : Time-stamp: <2005-10-09 14:06:38 antoine>
 *
 * \author R core team. Modified by A. Lucas for distance choice.
 *
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 2004   The R Development Core Team.
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
 */

#include <R.h>
/*#include "modreg.h" */ /* for declarations for registration */ 
#include "kmeans.h" 

#include "distance_T.h" 



/** K-means clustering using Lloyd algorithm.
 * \brief compute k-nearest centroid of our dataset.
 * \param x matrix of size nxp: input data
 * \param pn nb of individual (pn=n)
 * \param pp number of observation by individual (pp=p)
 * \param cen matrix of size k*p: centroids
 * \param pk number of centroids (k)
 * \param cl vector of flag of size n
 * \param pmaxiter integer: maximum iteration
 * \param nc vector of size k: number of individuals in cluster k.
 * \param wss vector of size k: sum of square in each cluster.
 * \param method: which method to use.
 */
void kmeans_Lloyd2(double *x, int *pn, int *pp, double *cen, int *pk, int *cl, 
		  int *pmaxiter, int *nc, double *wss, int * method)
{
  /* n: nb of individuals 
   * k: nb of clusters
   * p: number of abservations by individuals
   * x: matrix of size nxp
   * cen: matrix of size kxp
   */
    int n = *pn, k = *pk, p = *pp, maxiter = *pmaxiter;
    int iter, i, j, c, it, inew = 0;
    double best, dd;
    Rboolean updated;
    distance_T<double>::T_tri opt;
    int  ierr[0];
    double * data_tri;
    int * order_tri;
    int * rank_tri;
    
    if( (*method == distance_T<double>::SPEARMAN) ||  (*method == distance_T<double>::KENDALL))
      {
	opt.data_tri_x  = (double * ) malloc ( p * sizeof(double));
	opt.order_tri_x  = (int * ) malloc ( p * sizeof(int));
	opt.rank_tri_x  = (int * ) malloc ( p * sizeof(int));
	opt.data_tri_y  = (double * ) malloc ( p * sizeof(double));
	opt.order_tri_y  = (int * ) malloc ( p * sizeof(int));
	opt.rank_tri_y  = (int * ) malloc ( p * sizeof(int));
	if( (opt.data_tri_x == NULL) || (opt.order_tri_x == NULL) || (opt.rank_tri_x == NULL) ||
	    (opt.data_tri_y == NULL) || (opt.order_tri_y == NULL) || (opt.rank_tri_y == NULL)) 
	  error("distance(): unable to alloc memory");

      }

    for(i = 0; i < n; i++) cl[i] = -1;
    for(iter = 0; iter < maxiter; iter++) {
	updated = FALSE;
	for(i = 0; i < n; i++) {
	    /* find nearest centre for each point */
	    best = R_PosInf;
	    for(j = 0; j < k; j++) {

  	        dd = distance_T<double>::distance_kms(x,cen,n,k,p,i,j,method,ierr,opt);
		/*printf("| %f",dd);
		 */
		if(dd < best) {
		    best = dd;
		    inew = j+1;
		}
	    }
	    if(cl[i] != inew) {
		updated = TRUE;
		cl[i] = inew;
	    }
	}
	if(!updated) break;
	/* update each centre */
	for(j = 0; j < k*p; j++) cen[j] = 0.0;
	for(j = 0; j < k; j++) nc[j] = 0;
	for(i = 0; i < n; i++) {
	    it = cl[i] - 1; nc[it]++;
	    for(c = 0; c < p; c++) cen[it+c*k] += x[i+c*n];
	}
	for(j = 0; j < k*p; j++) cen[j] /= nc[j % k];
    }

    *pmaxiter = iter + 1;
    /*    for(j = 0; j < k; j++) wss[j] = 0.0; */
    for(i = 0; i < n; i++) {
	it = cl[i] - 1;
	wss[it] = distance_T<double>::distance_kms(x,cen,n,k,p,i,it,method,ierr,opt);
	wss[it] = wss[it] * (wss[it]) ;
	/*
	for(c = 0; c < p; c++) {
	    tmp = x[i+n*c] - cen[it+k*c];
	    wss[it] += tmp * tmp;
	    }*/
    }


    if( (*method == distance_T<double>::SPEARMAN) ||  (*method == distance_T<double>::KENDALL))
      {
	free(opt.data_tri_x);
	free(opt.rank_tri_x);
	free(opt.order_tri_x);	
	free(opt.data_tri_y);
	free(opt.rank_tri_y);
	free(opt.order_tri_y);	
      }
}
