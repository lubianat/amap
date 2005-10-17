/** \file hclust.c
 *                                                            
 *   \brief Hierarchical Clustering.
 *                                                          
 *   \date Created       : 14/11/02 
 *   \date Last Modified : Time-stamp: <2005-10-09 14:43:14 antoine>
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
#include "mva.h"


#define max( A , B )  ( ( A ) > ( B ) ? ( A ) : ( B ) )
#define min( A , B )  ( ( A ) < ( B ) ? ( A ) : ( B ) )

/** hierarchical clustering
 * \brief allocate distance matrix execute function R_distance, launch 
 * hclust on this distance matrix, and free memory of distance matrix.
 * \param x: data nr x nc
 * \param nr,nc number of row and columns		  
 * \param membr: member, vector  1:nr		  
 * \param method  integer -> distance method	  
 * \param diag integer: if we compute diagonal in distance matrix (usually no)
 * \param iopt    integer -> link used 
 * \param ia, ib   result (merge)			  
 * \param iorder  result (order)			  
 * \param crit    result (height)			  
 * \param result  flag  0 => correct		  
 *               1 => Pb			  
 *               2 => Cannot allocate memory 
 *               3 => Pb with distance matrix
 */

void hcluster(double *x, int *nr, int *nc, int *diag, int *method, int *iopt ,int *ia , int *ib,int *iorder,double *crit,double *membr,int *result)
{
  int  len;
  int flag;
  double *d;

  *result = 1;

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

  R_distance(x,nr,nc,d,diag,method,&flag);
  
  if(flag == 0)
    {
      printf("AMAP: Unable to compute Hierarchical Clustering: missing values in distance matrix\n"); 
      *result = 3;
      return;
    }

  /*
   *  Hierarchical clustering
   */
  hclust (nr,&len,iopt ,ia ,ib,iorder,crit,membr,d,result);
  free(d); 
  *result = 0;
}


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

  int im,jm,jj,i,j,ncl,ind,i2,j2,k,ind1,ind2,ind3;
  double inf,dmin,x,xx;
  int *nn;
  double *disnn;
  short int *flag;
  int *iia;
  int *iib;


  *result = 1;
  nn    = (int*) malloc (*n * sizeof(int));
  disnn = (double*) malloc (*n * sizeof(double));
  flag  = (short int*) malloc (*n * sizeof(short int));
  if(nn == NULL || disnn == NULL || flag == NULL )
    {
      printf("AMAP: Not enought system memory \n"); 
      *result = 2;
      return;
    }


  /* Initialisation */
  for ( i=0; i<*n ; i++)
    flag[i]=1;
 
  ncl=*n;
  inf= 1e20;

  /*
   * Carry out an agglomeration - first create list of NNs
   */
  for ( i=0; i<(*n-1) ; i++)
    {
      dmin=inf;
      for (j=i+1; j<*n; j++)
	{
	  ind = ioffst(*n,i,j);
	  if (diss[ind]< dmin)
	    {
	      dmin = diss[ind];
	      jm=j;
	    }
	}
      nn[i]=jm;
      disnn[i]=dmin;
    }
    
  /*
   *  Repeat previous steps until N-1 agglomerations carried out.
   */
  while (ncl > 1)
    {
      /*
       * Next, determine least diss. using list of NNs
       */
      dmin = inf;
      for ( i=0; i<(*n-1) ; i++)
	{
	  if( flag[i] )
	    {
	      if (disnn[i] < dmin )
		{
		  dmin = disnn[i];
		  im=i;
		  jm=nn[i];
		}
	    }
	}
      ncl = ncl-1;

      /*
       * This allows an agglomeration to be carried out.
       * At step n-ncl, we found dmin=dist[i2,j2]
       */
      i2=min (im,jm);
      j2=max (im,jm);
      ia[*n-ncl-1]=i2+1;
      ib[*n-ncl-1]=j2+1;
      crit[*n-ncl-1]=dmin;
 
     
      /*
       * Update dissimilarities from new cluster.
       */
      flag[j2]=0;
      dmin=inf;
      for ( k=0; k< *n ; k++)
	{ 
	  if(flag[k] && (k != i2) )
	    {
	      x =  membr[i2] + membr[j2] + membr[k];
	      if (i2 < k)
		{
		  ind1 = ioffst(*n,i2,k);
		}
	      else
		{
		  ind1 = ioffst(*n,k,i2);
		}
	      if (j2 < k)
		{
		  ind2 = ioffst(*n,j2,k);
		}
	      else
		{
		  ind2 = ioffst(*n,k,j2);
		}
	      ind3=ioffst(*n,i2,j2);
	      xx=diss[ind3];
	      /*
	       * Gi & Gj are agglomerated => Gii
	       * We are calculating D(Gii,Gk) (for all k)
	       *
	       * diss[ind1] = D(Gi,Gk) (will be replaced by  D(Gii,Gk))
	       * diss[ind2] = D(Gj,Gk) 
	       * xx = diss[ind3] = D(Gi,Gj)
	       *
	       * membr[i2] = #Gi
	       * membr[j2] = #Gj
	       * membr[k]  = #Gk
	       * 
	       * x = #Gi + #Gj + #Gk
	       */
	      switch(*iopt)
		{
		  /*
		   * WARD'S MINIMUM VARIANCE METHOD - IOPT=1.
		   */	      
		case 1: 
		  {
		    diss[ind1] = (membr[i2]+membr[k])* diss[ind1] + 
		      (membr[j2]+membr[k])* diss[ind2] - 
		      membr[k] * xx;
		    diss[ind1] = diss[ind1] / x;
		    break; 
		  }
		  /*
		   * SINGLE LINK METHOD - IOPT=2.
		   */
		case 2: diss[ind1] = min (diss[ind1],diss[ind2]); break; 
		  /*
		   * COMPLETE LINK METHOD - IOPT=3.
		   */
		case 3: diss[ind1] = max (diss[ind1],diss[ind2]); break; 
		  /*
		   * AVERAGE LINK (OR GROUP AVERAGE) METHOD - IOPT=4.
		   */
		case 4:  diss[ind1] = ( membr[i2] * diss[ind1] +
					membr[j2] * diss[ind2] ) /
			   (membr[i2] + membr[j2]); 
		  break; 
		  /*
		   *  MCQUITTY'S METHOD - IOPT=5.
		   */
		case 5:  diss[ind1] = 0.5 * diss[ind1]+0.5*diss[ind2]; 
		  break;
		  /*
		   * MEDIAN (GOWER'S) METHOD - IOPT=6.
		   */
		case 6:  diss[ind1] = 0.5* diss[ind1]+0.5*diss[ind2] -0.25*xx;
		  break;
		  /*
		   * CENTROID METHOD - IOPT=7.
		   */
		case 7:
		  diss[ind1] = (membr[i2]*diss[ind1] + membr[j2]*diss[ind2] - 
				membr[i2] * membr[j2]*xx /
				(membr[i2] + membr[j2]) ) /
		    (membr[i2] + membr[j2]);
		  break;
		} /* end switch */
	      if( (i2 <= k) && ( diss[ind1] < dmin ) )
		{
		  dmin = diss[ind1];
		  jj=k;
		}
	    } /* 800 */
	}
      membr[i2] = membr[i2] +   membr[j2];
      disnn[i2] = dmin;

      nn[i2] = jj;


      /*
       *  Update list of NNs insofar as this is required.
       */
      for ( i=0; i< (*n-1) ; i++)
	{
	  if( flag[i] && ((nn[i] == i2 ) || (nn[i] == j2) ) )
	    {
	      /* (Redetermine NN of I:)   */
	      dmin = inf;
	      for  ( j=i+1; j< *n ; j++)
		{
		  ind= ioffst(*n,i,j);
		  if(flag[j] && (i != j) && (diss[ind] < dmin) )
		    {
		      dmin =diss[ind];
		      jj=j;
		    }
		  nn[i] = jj;
		  disnn[i] = dmin;
		}
	    }
	}
    } /* end of while */

  
  free(nn);
  free(disnn);
  free(flag);

  
  iia = (int*) malloc (*n * sizeof(int));
  iib = (int*) malloc (*n * sizeof(int));
  if(iia == NULL || iib == NULL )
    {
      printf("AMAP: Not enought system memory \n"); 
      *result = 2;
      return;
    }


  
  hcass2(n,ia,ib,iorder,iia,iib);
  

  /*
   * Cette boucle devrait etre evitee...
   */
  
  for (i=0; i< *n; i++ ) 
    {
      ia[i]= iia[i];
      ib[i]= iib[i];
    }

  free(iia);
  free(iib);
  
  *result = 0;
} /* end function hclust */



/** Return indice 
 * \brief The upper half diagonal distance matrix is stored as a vector...
 * so distance between individual i and j is stored at postion ioffst(i,j)
 * \param n number of individuals (distance matrix is nxn)
 * \param i,j: indices in matrix
 */
int ioffst(int n,int i,int j)
     /* Map row I and column J of upper half diagonal symmetric matrix 
      * onto vector.  i < j 
      */
{
  return j+i*n-(i+1)*(i+2)/2 ;
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

/** Hierachical clustering subroutine
 * \brief compute hierachical clustering from a distance matrix 
 * This routine is called by hclust
 * \param n number of individuals
 * \param ia, ib result (merge)
 * \param iia, iib result (merge)
 * \param iorder result (order)
 *
 * \note this is an adaptation of the fortran function designed from the
 * R core team.
 */
void hcass2( int *n, int *ia,  int *ib,int *iorder, int *iia, int *iib)
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
