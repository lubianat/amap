
#ifndef _AMAP_DISTANCE_TEMPLATE
#define _AMAP_DISTANCE_TEMPLATE 1



template<class T> class distance_T
{

 public:
  /* == 1,2,..., defined by order in the r function dist */
  enum { EUCLIDEAN=1, MAXIMUM, MANHATTAN, CANBERRA, BINARY ,PEARSON, CORRELATION, SPEARMAN};

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


  /** \brief R_distance_kms: compute distance between individual i1 and
   * centroid i2
   *
   * compute distance and call one of function R_euclidean or R_...
   * This function is called by kmeans_Lloyd2 
   *
   * \param x input matrix (individuals)
   * \param y input matrix (centroids)
   * \param nr1,nr2,nc number of row (nr1:x, nr2:y) and columns
   *        nr individuals with nc values.
   * \param i1, i2: indice of individuals (individual i1, centroid i2)
   * \param method 1, 2,... method used
   * \param ierr for NA 0 if no value can be comuted due to NA
   * \param opt optional parameter required for spearman
   */
  static T distance_kms(double *x,double *y, int nr1,int nr2, 
			int nc,int i1,int i2, int *method,
			int * ierr, void ** opt);


 private:
  
  static void* thread_dist(void* arguments);

  /** \brief Distance euclidean (i.e. sqrt(sum of square) )
   */
  static T  R_euclidean(double * x, double * y , int nr_x, int nr_y, int nc, 
			int i1, int i2,
			int * flag, void ** opt);

  /** \brief Distance maximum (supremum norm)
   */
  static T R_maximum(double * x, double * y , int nr_x, int nr_y, int nc, 
		     int i1, int i2,
		     int * flag, void ** opt);

  /** \brief Distance manhattan
   */
  static T R_manhattan(double * x, double * y , int nr_x, int nr_y, int nc, 
		       int i1, int i2,
		       int * flag, void ** opt);

  /** \brief Distance canberra
   */
  static T R_canberra(double * x, double * y , int nr_x, int nr_y, int nc, 
		      int i1, int i2,
		      int * flag, void ** opt);

  /** \brief Distance binary
   */
  static T R_dist_binary(double * x, double * y , int nr_x, int nr_y, int nc, 
			 int i1, int i2,
			 int * flag, void ** opt);

  /** \brief Pearson / Pearson centered (correlation)
   *  \note Added by A. Lucas
   */
  static T R_pearson(double * x, double * y , int nr_x, int nr_y, int nc, 
		     int i1, int i2,
		     int * flag, void ** opt);

  /** \brief Distance correlation (Uncentered Pearson)
   *  \note Added by A. Lucas
   */
  static T R_correlation(double * x, double * y , int nr_x, int nr_y, int nc, 
			 int i1, int i2,
			 int * flag, void ** opt);

  /** \brief Spearman distance (rank base metric)
   *  \note Added by A. Lucas
   */
  static T R_spearman(double * x, double * y , int nr_x, int nr_y, int nc, 
		      int i1, int i2,
		      int * flag, void ** opt);


};

#endif

#ifndef _AMAP_DISTANCE_TEMPLATE_CPP 
#define _AMAP_DISTANCE_TEMPLATE_CPP 1
#include "distance_T.cpp"
#endif


