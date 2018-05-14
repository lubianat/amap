#ifdef __cplusplus
extern "C" {
#endif


  
  void hclust(int *n,int *len, int *iopt ,int *ia , int *ib,int *iorder,double *crit,double *membr,double *diss, int *result);

  void hcluster(double *x, int *nr, int *nc, int *diag, int *method, int *iopt ,int *ia , int *ib,int *iorder,double *crit,double *membr,int *nbprocess,int * precision, int * result);

#ifdef __cplusplus
};




namespace hierclust
{
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
  void hcass2( int *n, int *ia,  int *ib,int *iorder, int *iia, int *iib);

};


#endif
