#include <stdlib.h>
#include <math.h>

#include "mva.h" 
#ifndef __MINGW_H
#include <pthread.h>
#endif

#define NBPROCESS  2
 
/*  int NBPROCESS = 2; */



#define max( A , B )  ( ( A ) > ( B ) ? ( A ) : ( B ) )
#define min( A , B )  ( ( A ) < ( B ) ? ( A ) : ( B ) )

enum { EUCLIDEAN=1, MAXIMUM, MANHATTAN, CANBERRA, BINARY ,PEARSON, CORRELATION, SPEARMAN};
/* == 1,2,..., defined by order in the R function dist */



void R_distancepar(double *x, int *nr, int *nc, double *d, int *diag, int *method,int *nbprocess, int * ierr)
{
#ifndef __MINGW_H
    int dc,numero_thread;
    short int * jobs;
    int  i;
    double (*distfun)(double*, int, int, int, int, int) = NULL;
    pthread_t * th;
    void ** arguments;

    dc = (*diag) ? 0 : 1; /* diag=1:  we do the diagonal */ 

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
    /*    numero_thread=0;*/
    arguments = malloc ( 9 * sizeof( void *));
    jobs      = malloc ( *nr * sizeof(short int));
    th = (pthread_t *) malloc ( *nbprocess * sizeof(pthread_t));

    for(i=0; i< *nr; i++){jobs[i]=0;}  

    /*    arguments[0] =  &numero_thread ;*/
    arguments[0] =  jobs;
    arguments[1] =  nr ;
    arguments[2] =  nc ;
    arguments[3] =  &dc ;
    arguments[4] =  x ;
    arguments[5] =  d ;
    arguments[6] =  method ;
    arguments[7] =  nbprocess;
    arguments[8] =  ierr;

    *ierr = 1; /* res = 1 => no missing values
	          res = 0 => missings values */


    for (i=0; i < *nbprocess ; i++)
    {
      pthread_create(th+i,0,thread_dist,(void *)arguments);
    }

   /* Attends la fin    */
  for (i=0; i < *nbprocess ; i++)
    {      
      pthread_join(*(th+i),NULL);
    }      
#else
  R_distance(x, nr, nc, d,diag,method,ierr);

#endif

}

#ifndef __MINGW_H
void* thread_dist(void* arguments)
{

  long int nbprocess,no,nr,nc,i,j,dc,debut,fin,ij;
  int * tmp;
  short int *jobs;
  void ** arguments2; 
  double * d;
  double * x;
  int * method;
  int * ierr;
  /* for spearman dist */
  void ** opt;
  double * data_tri;
  int * order_tri;
  int * rank_tri;
  double (*distfun)(double*, int, int, int, int, int *, void **) = NULL;
  
  

  arguments2 = (void **) arguments;

  jobs = (short int *)  arguments2[0];
  nr = * (int *) arguments2[1];
  nc = * (int *) arguments2[2];
  dc = * (int *) arguments2[3];
  x  = (double *)  arguments2[4];
  d  = (double *)  arguments2[5];
  method = (int *) arguments2[6];
  nbprocess = * (int *) arguments2[7];
  ierr =  (int *) arguments2[8];

  /*
    no = * (int *) arguments2[0];*/
  /* Increment du thread */ 
  /*  tmp = (int *) arguments2[0];
  *tmp = no+1;
  */

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

  /*
    debut = (long int) floor( ((nr+1.)*nbprocess - sqrt( (nr+1.)*(nr+1.) * nbprocess * nbprocess - (nr+1.)*(nr+1.) * nbprocess * no  ) )/nbprocess);
    fin = (long int) floor(((nr+1.)*nbprocess - sqrt( (nr+1.)*(nr+1.) * nbprocess * nbprocess - (nr+1.)*(nr+1.) * nbprocess * (no+1.)  ) )/nbprocess);
    printf("Thread %d debut %d fin %d\n",no,debut,fin);

  */


    for(j = 0 ; j <= nr ; j++)
      if(jobs[j]==0)
	{
	  jobs[j]=1;
	  ij = (2 * (nr-dc) - j +1) * (j) /2 ;
	  for(i = j+dc ; i < nr ; i++)
	    {
	      d[ij++] = distfun(x, nr, nc, i, j,ierr,opt);
	    }
	}

    if((*method) == SPEARMAN)
      {
	free(opt);
	free(data_tri);
	free(order_tri);
	free(rank_tri);
      }

    return (void*)0;

}
#endif

