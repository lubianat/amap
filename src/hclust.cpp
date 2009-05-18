/** \file hclust.c
 *                                                            
 *   \brief Hierarchical Clustering.
 *                                                          
 *   \date Created       : 14/11/02 
 *   \date Last Modified : Time-stamp: <2007-10-03 20:33:11 antoine>
 *

 *  \author F. Murtagh, ESA/ESO/STECF, Garching, February 1986. 
 *  \author Modifications for R: Ross Ihaka, Dec 1996                 
 *                        Fritz Leisch, Jun 2000               
 *   all vars declared:   Martin Maechler, Apr 2001            
 *  \author C adaptation:        A. Lucas, Dec 2002
 */

/*
 *   Parameters:                                               
 *                                                          
 *   N                 the number of points being clustered     
 *   DISS(LEN)         dissimilarities in lower half diagonal   
 *                     storage; LEN = N.N-1/2,                  
 *   IOPT              clustering criterion to be used,         
 *   IA, IB, CRIT      history of agglomerations; dimensions    
 *                     N, first N-1 locations only used,        
 *   MEMBR, NN, DISNN  vectors of length N, used to store       
 *                     cluster cardinalities, current nearest   
 *                     neighbour, and the dissimilarity assoc.  
 *                     with the latter.                         
 *                     MEMBR must be initialized by R to the    
 *                     default of  rep(1, N)                    
 *   FLAG              boolean indicator of agglomerable obj.  
 *                     clusters.                                
 *                                                          
 *  
 *   Fortran algorithme:
 *   F. Murtagh, ESA/ESO/STECF, Garching, February 1986.        
 *   Modifications for R: Ross Ihaka, Dec 1996                 
 *                        Fritz Leisch, Jun 2000               
 *   all vars declared:   Martin Maechler, Apr 2001            
 *   C adaptation:        A. Lucas, Dec 2002
 * ------------------------------------------------------------ */



#include <stdlib.h>
#include <math.h>
#include "hclust.h"
#include "distance.h"
#include "hclust_T.h"

#define max( A , B )  ( ( A ) > ( B ) ? ( A ) : ( B ) )
#define min( A , B )  ( ( A ) < ( B ) ? ( A ) : ( B ) )

/** Hierachical clustering
 * \brief compute hierachical clustering from a distance matrix 
 * hclust can be called directly by R, or with hcluster C function.
 * \param n number of individuals
 * \param len = n (n-1) / 2 (size of distance matrix)
 * \param *iopt integer -> link used
 * \param ia, ib result (merge)
 * \param iorder result (order)
 * \param crit result (height)
 * \param membr number of individuals in each cluster.
 * \param diss distance matrix (size len)
 * \param result  flag  0 => correct              
 *               1 => Pb                          
 *               2 => Cannot allocate memory 
 *               3 => Pb with distance matrix
 *
 * \note this is an adaptation of the fortran function designed from the
 * R core team.
 */
void hclust(int *n,int *len, int *iopt ,int *ia , int *ib,int *iorder,double *crit,double *membr,double *diss,int *result)
{
  hclust_T::hclust<double>(n,len,iopt ,ia ,ib,iorder,
			   crit,membr,diss, result);
} /* end function hclust */




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
 * \param precision 1: float; 2: double
 * \param ia, ib   result (merge)			  
 * \param iorder  result (order)			  
 * \param crit    result (height)			  
 * \param result  flag  0 => correct		  
 *               1 => Pb			  
 *               2 => Cannot allocate memory 
 *               3 => Pb with distance matrix
 *
 * \note this is an adaptation of the fortran function designed from the
 * R core team.
 */
void hcluster(double *x, int *nr, int *nc, int *diag, int *method, int *iopt ,int *ia , int *ib,int *iorder,double *crit,double *membr,int *nbprocess,int * precision, int * result)
{
  

  if(*precision == 1)
    
    hclust_T::hcluster<float>(x,nr,nc,diag,method, iopt ,
				  ia ,ib,iorder,crit,membr,
				  nbprocess,  result);
  else
    hclust_T::hcluster<double>(x,nr,nc,diag,method, iopt ,
				   ia ,ib,iorder,crit,membr,
				   nbprocess,  result);
  return;
    
}






/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*                                                             */
/* Given a HIERARCHIC CLUSTERING, described as a sequence of   */
/* agglomerations, prepare the seq. of aggloms. and "horiz."   */
/* order of objects for plotting the dendrogram using S routine*/
/* 'plclust'.                                                  */
/*                                                             */
/* Parameters:                                                 */
/*                                                             */
/* IA, IB:       vectors of dimension N defining the agglomer- */
/*                ations.                                      */
/* IIA, IIB:     used to store IA and IB values differently    */
/*               (in form needed for S command `plclust`       */
/* IORDER:       "horiz." order of objects for dendrogram      */
/*                                                             */
/* F. Murtagh, ESA/ESO/STECF, Garching, June 1991              */
/* C adaptation:              A. Lucas, Nov 2002               */
/*                                                             */
/* HISTORY                                                     */
/*                                                             */
/* Adapted from routine HCASS, which additionally determines   */
/*  cluster assignments at all levels, at extra comput. expense*/
/*                                                             */
/*-------------------------------------------------------------*/

void hierclust::hcass2( int *n, int *ia,  int *ib,int *iorder, int *iia, int *iib)
{
  int i,j,k,k1,k2,loc;

  /*  Following bit is to get seq. of merges into format acceptable to plclust
   *  I coded clusters as lowest seq. no. of constituents; S's `hclust' codes
   *  singletons as -ve numbers, and non-singletons with their seq. nos.
   */

  for (i=0; i< *n; i++ ) 
    {
      iia[i]= - ia[i];
      iib[i]= - ib[i];
    }

  for (i=0; i< (*n-2); i++ ) 
    /* In the following, smallest (+ve or -ve) seq. no. wanted */
    {
      k = min ( ia[i] , ib[i] );
      for (j=i+1; j< (*n-1); j++ ) 
	{
	  if( ia[j] == k ) iia[j]= i+1;
	  if( ib[j] == k ) iib[j]= i+1;
	}
    }  

   for (i=0; i< (*n-1); i++ ) 
     {
       if( (iia[i] > 0) && (iib[i] < 0) )
	 {
	   k= iia[i];
	   iia[i] = iib[i];
	   iib[i] = k;
	 }
       if( (iia[i] > 0) && (iib[i] > 0))
	 {
	   k1= min ( iia[i], iib[i] );
	   k2= max ( iia[i], iib[i] );
	   iia[i] = k1;
	   iib[i] = k2;
	 }
     }



   /* New part for 'order' */

   iorder[0] = - iia[*n-2];
   iorder[1] = - iib[*n-2];
   loc = 2;
   for (i=(*n-3); i>= 0; i-- ) 
     {
       for (j=0; j< loc; j++ ) 
	 {
	   if ( -iorder[j] == i +1)
	     {
	       /* REPLACE IORDER(J) WITH IIA(I) AND IIB(I) */ 
	       iorder[j] = - iia[i];
	       if (j == (loc-1)) 
		 {
		   loc++;
		   iorder[loc-1]= - iib[i];

		   break; /* for j */

		 }
	       loc++;
	       for (k=loc-1; k >= (j+1); k--)
		 {
		   iorder[k]=iorder[k-1];
		 }
	       iorder[j+1]= - iib[i];
	       break; /* for j */

	     }
	 }     
     }
}
