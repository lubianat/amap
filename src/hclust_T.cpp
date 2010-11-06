
#define _AMAP_HCLUST_TEMPLATE_CPP 1


#include "hclust.h"
#include "distance_T.h"
#include "distance.h"
#include "hclust_T.h"
#include <stdio.h>



namespace hclust_T
{

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
   * \param nbprocess number of processes (thread) used
   * \param result  flag  0 => correct		  
   *               1 => Pb			  
   *               2 => Cannot allocate memory 
   *               3 => Pb with distance matrix
   *
   */
  template <class T> void hcluster(double *x, int *nr, int *nc, int *diag, int *method, int *iopt ,int *ia , int *ib,int *iorder,double *crit,double *membr,int *nbprocess, int * result)
  {
    int  len;
    int flag;
    T *d;
    
    *result = 1;
    

    len = (*nr * (*nr-1) )/2; 
    d = (T*) malloc (len * sizeof(T));
    if(d == NULL )
      {
	printf("AMAP: Not enought system memory to allocate matrix distance\n"); 
	*result = 2;
	return;
      }
    /*
     * Calculate d: distance matrix
     */
    distance_T<T>::distance(x,nr,nc,d,diag,method,nbprocess,&flag,-1);
    if(flag == 0)
      {
	printf("AMAP: Unable to compute Hierarchical Clustering: missing values in distance matrix\n"); 
	*result = 3;
	return;
      }


    /*
     *  Hierarchical clustering
     */
      hclust_T::hclust<T>(nbprocess,x,*nr,*nc,method,nr,&len,iopt ,ia ,ib,iorder,crit,membr,d,result);
		free(d); 
    *result = 0;
    
  }


	
  template <class T> void hclust(int * nbprocess,double *mx,int nr, int nc,
				 int *method,int *n,int *len, int *iopt ,int *ia , 
				 int *ib,int *iorder,double *crit,double *membr,T *diss,int *result)
  {
    int im,jm,jj,i,j,ncl,ind,i2,j2,k,ind1,ind2,ind3;
    double inf,dmin,x,xx;
    int *nn;
    int *items = NULL;
    double *disnn;
    short int *flag;
    int *iia;
    int *iib;
    int count,h,idx1,idxitem1,idx2;

    
    *result = 1;
    nn    = (int*) malloc (*n * sizeof(int));
    if(*iopt==CENTROID2) 
      {
	items   = (int*) malloc (*n*nc * sizeof(int));
      }
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
      {
	flag[i]=1;
	idxitem1=i;
	if(items != NULL) 
	  {
	    for ( j=0; j<nc ; j++) 
	      {
		items[idxitem1]=1;  /* only one item per cluster and per coordinate for iopt==8*/
		idxitem1+=nr;
	      }
	  }
      }
    
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
	i2=MIN (im,jm);
	j2=MAX (im,jm);
	ia[*n-ncl-1]=i2+1;
	ib[*n-ncl-1]=j2+1;
	crit[*n-ncl-1]=dmin;
	
	
	/*
	 * Update dissimilarities from new cluster.
	 */
	flag[j2]=0;
	dmin=inf;
	jj=0;

	/*
	 * Cluster3 CENTROID METHOD - IOPT=8.
	 */
	if (*iopt==CENTROID2) 
	  {
	    /* compute centroind coordinates*/
	    idx1=i2;
	    idx2=j2;
	    ind3=ioffst(*n,i2,j2);
	    //printf("Aggregate %d-%d %d-%d (method=%d, dmin=%f (%f))\n",i2,j2,im,jm,*method,dmin,diss[ind3]);
	    for(h = 0 ; h< nc ; h++) 
	      {
		if(R_FINITE(mx[idx1]) && R_FINITE(mx[idx2])) 
		  {
		    mx[idx1] = (items[idx1]*mx[idx1] + items[idx2]*mx[idx2])/(items[idx1]+items[idx2]);
		    items[idx1]+=items[idx2];
		    
		  }
		else
		  {
		    if(!R_FINITE(mx[idx1]) && R_FINITE(mx[idx2])) 
		      {
			mx[idx1] = mx[idx2];
			items[idx1]=items[idx2];
		      }
		  }
		idx1 += nr;
		idx2 += nr;
	      }
	    int flg;
	    int dg=0;
	    /* recompute all distances in parallel */
	    distance_T<T>::distance(mx,&nr,&nc,diss,&dg,method,nbprocess,&flg, i2);
	    /*update disnn and nn*/
	    for  ( k=0; k< *n ; k++) 
	      {		
		if(flag[k] && (k != i2) )
		  {
		    if (i2 < k)
		      {
			ind1 = ioffst(*n,i2,k);
		      }
		    else
		      {
			ind1 = ioffst(*n,k,i2);
		      }

		    if( (i2 < k) && ( diss[ind1] < dmin ) )
		      {
			dmin = diss[ind1];
			jj=k;
		      }
		    if( (i2 > k) && (nn[k]!=j2) && ( diss[ind1] < disnn[k] ) )
		      {
			disnn[k]=diss[ind1];
			nn[k]=i2;
		      }
		  }
	      }
	  } 
	else 
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
		  case WARD: 
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
		  case SINGLE: diss[ind1] = MIN (diss[ind1],diss[ind2]); break; 
		    /*
		     * COMPLETE LINK METHOD - IOPT=3.
		     */
		  case COMPLETE: diss[ind1] = MAX (diss[ind1],diss[ind2]); break; 
		    /*
		     * AVERAGE LINK (OR GROUP AVERAGE) METHOD - IOPT=4.
		     */
		  case AVERAGE:  diss[ind1] = ( membr[i2] * diss[ind1] +
					  membr[j2] * diss[ind2] ) /
			     (membr[i2] + membr[j2]); 
		    break; 
		    /*
		     *  MCQUITTY'S METHOD - IOPT=5.
		     */
		  case MCQUITTY:  diss[ind1] = 0.5 * diss[ind1]+0.5*diss[ind2]; 
		    break;
		    /*
		     * MEDIAN (GOWER'S) METHOD - IOPT=6.
		     */
		  case MEDIAN:  diss[ind1] = 0.5* diss[ind1]+0.5*diss[ind2] -0.25*xx;
		    break;
		    /*
		     * CENTROID METHOD - IOPT=7.
		     */
		  case CENTROID:
		    diss[ind1] = (membr[i2]*diss[ind1] + membr[j2]*diss[ind2] - 
				  membr[i2] * membr[j2]*xx /
				  (membr[i2] + membr[j2]) ) /
		      (membr[i2] + membr[j2]);
		    break;
		  } /* end switch */
					
	
		if( (i2 < k) && ( diss[ind1] < dmin ) )
		  {
		    dmin = diss[ind1];
		    jj=k;
		  }
		if( (i2 > k) && (nn[k]!=j2) && ( diss[ind1] < disnn[k] ) )
		  {
		    disnn[k]=diss[ind1];
		    nn[k]=i2;
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
	/*printf("%d/%d\n",ncl,*n);*/
      } /* end of while */

  
    free(nn);
    free(disnn);
    free(flag);

    if(items != NULL )
      free(items);
  
    iia = (int*) malloc (*n * sizeof(int));
    iib = (int*) malloc (*n * sizeof(int));
    if(iia == NULL || iib == NULL )
      {
	printf("AMAP: Not enought system memory \n"); 
	*result = 2;
	return;
      }


  
    hierclust::hcass2(n,ia,ib,iorder,iia,iib);
  

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
  }


}
