extern "C"
{

  
  void hclust(int *n,int *len, int *iopt ,int *ia , int *ib,int *iorder,double *crit,double *membr,double *diss, int *result);

  void hcluster(double *x, int *nr, int *nc, int *diag, int *method, int *iopt ,int *ia , int *ib,int *iorder,double *crit,double *membr,int *nbprocess,int * precision, int * result);

  int ioffst(int n,int i,int j);


  void hcass2( int *n, int *ia,  int *ib,int *iorder, int *iia, int *iib);
}
