

#ifndef _AMAP_HCLUST_TEMPLATE
#define _AMAP_HCLUST_TEMPLATE 1

namespace hclust_T
{

  template <class T> void hcluster(double *x, int *nr, int *nc, 
				   int *diag, int *method, int *iopt ,
				   int *ia , int *ib,int *iorder,
				   double *crit,double *membr,
				   int *nbprocess, int * result);
  
  template <class T> void hclust(int *n,int *len, int *iopt ,int *ia ,
				 int *ib,int *iorder,double *crit,
				 double *membr,T *diss,int *result);
  

}

#endif


#ifndef _AMAP_HCLUST_TEMPLATE_CPP
#define _AMAP_HCLUST_TEMPLATE_CPP 1
#include "hclust_T.cpp"
#endif
