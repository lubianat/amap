/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
 *                                                            
 *   HIERARCHICAL CLUSTERING using (user-specified) criterion. 
 *                                                          
 *   Created       : 14/11/02 
 *   Last Modified : Time-stamp: <2005-02-26 18:07:29 antoine>
 *
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
#include "mva.h"


#define max( A , B )  ( ( A ) > ( B ) ? ( A ) : ( B ) )
#define min( A , B )  ( ( A ) < ( B ) ? ( A ) : ( B ) )


void hclusterpar(double *x, int *nr, int *nc, int *diag, int *method, int *iopt ,int *ia , int *ib,int *iorder,double *crit,double *membr,int *nbprocess,int * result)
     /*
      * x: data nr x nc
      * membr: member, vector  1:nr
      * method  integer -> distance method
      * diag
      * iopt    integer -> link used
      * ia ib   result (merge)
      * iorder  result (order)
      * crit    result (height)
      * nbprocess nb of process for parallelization
      * result  flag  0 => correct
      *               1 => Pb
      *               2 => Cannot allocate memory
      *               3 => Pb with distance matrix
      */ 

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

