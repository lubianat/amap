
#include "distance.h"
#include "distance_T.h"



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
void R_distance(double *x, int *nr, int *nc, double *d, int *diag, int *method,int *nbprocess, int * ierr)
{


  distance_T<double>::distance(x,nr,nc, d,diag,method,
			       nbprocess, ierr);
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
