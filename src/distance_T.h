
#ifndef _AMAP_DISTANCE_TEMPLATE
#define _AMAP_DISTANCE_TEMPLATE 1

/* == 1,2,..., defined by order in the r function dist */
enum { EUCLIDEAN=1, MAXIMUM, MANHATTAN, CANBERRA, BINARY ,PEARSON, CORRELATION, SPEARMAN};

template<class T> class distance_T
{

 private:

  struct T_argument
  {
    short int id;
    double * x;
    int * nr;
    int * nc;
    bool dc;
    T * d;
    int * method;
    int  nbprocess;
    int * ierr;
  };


  
  // only static functions; no attributes
  distance_T();
    

  ~distance_T();

 public:

  static void distance(double *x, int *nr, int *nc, T *d, int *diag,
		       int *method,int *nbprocess, int * ierr);

 private:
  
  static void* thread_dist(void* arguments);

  /** \brief Distance euclidean (i.e. sqrt(sum of square) )
   */
  static T  R_euclidean(double * x, int nr, int nc, int i1, int i2,int * flag, void ** opt);

  /** \brief Distance maximum (supremum norm)
   */
  static T R_maximum(double * x, int nr, int nc, int i1, int i2, int * flag, void ** opt);

  /** \brief Distance manhattan
   */
  static T R_manhattan(double * x, int nr, int nc, int i1, int i2, int * flag, void ** opt);

  /** \brief Distance canberra
   */
  static T R_canberra(double * x, int nr, int nc, int i1, int i2, int * flag, void ** opt);

  /** \brief Distance binary
   */
  static T R_dist_binary(double * x, int nr, int nc, int i1, int i2,int * flag, void ** opt);

  /** \brief Pearson / Pearson centered (correlation)
   *  \note Added by A. Lucas
   */
  static T R_pearson(double * x, int nr, int nc, int i1, int i2,int * flag, void ** opt);

  /** \brief Distance correlation (Uncentered Pearson)
   *  \note Added by A. Lucas
   */
  static T R_correlation(double * x, int nr, int nc, int i1, int i2,int * flag, void ** opt);

  /** \brief Spearman distance (rank base metric)
   *  \note Added by A. Lucas
   */
  static T R_spearman(double * x, int nr, int nc, int i1, int i2,int * flag, void ** opt);


};

#endif

#ifndef _AMAP_DISTANCE_TEMPLATE_CPP 
#define _AMAP_DISTANCE_TEMPLATE_CPP 1
#include "distance_T.cpp"
#endif


