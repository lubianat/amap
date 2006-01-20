

#include <stdio.h>


/** \fn diss compute a dissimilarity matrix
 * \brief diss(i,j) = number of values in individual i that
 *                    are equal to values in individual j 
 *                    - number of values that are different
 */


void diss (int * data, double * res, int * n, int * p,double * w)
{
  int k,i,j=0;
 


  for(i=0; i< (*n); i++)
    {
      /* follow lines: case when i=j */
      j = i;
      res[i + (j * (*n))] = *p;

      for(j=i+1; j < (*n); j++)
	{
	  
	  res[j + (i * (*n))] = 0;
	  for(k=0; k< (*p);k++)
	    {
	      if(data[i + (k * (*n))] == data[j + (k * (*n))])
		(res[j + (i * (*n))]) +=w[k] ;
	      else
		(res[j + (i * (*n))]) -=w[k] ;
	    }

	  res[i + (j * (*n))] = res[j + (i * (*n))] ;

	}

    }

}
