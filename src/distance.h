
#ifdef __cplusplus
extern "C" {
#endif


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
  void R_distance(double *x, int *nr, int *nc, double *d, int *diag, 
		  int *method,int * nbprocess, int * ierr);



  void rsort_rank_order(double *x, int *order, int *rank, int * n);

#ifdef __cplusplus
};
#endif


