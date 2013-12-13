#ifndef ACPROB_H
#define ACPROB_H 1


extern "C" {


  /*! noyau: base function for kernel computation
   * \brief compute a kernel. called by W or VarRob
   * \param u: input data (scalar)
   * \param k: char, g, q, t, e, c, u for:
   *        gaussien = (2*pi)^(-1/2) * exp(-u^2/2),
   *        quartic   = 15/16 * (1-u^2)^2 * (abs(u)<1),
   *        triweight = 35/32 * (1-u^2)^3 * (abs(u)<1),
   *        epanechikov = 3/4 * (1-u^2) *   (abs(u)<1),
   *        cosinus = pi/4 * cos (u*pi/2) * (abs(u)<1),
   *        uniform = 1/2 * (abs(u)<1),
   * \param res: output data: one value of type double
   */
  void noyau(double *u, char **k,double *res) ;

  
  /** W compute a "local" variance This function is called directly by R
   * \brief this functions compute a "local" variance, using a specific kernel.
   * This function depends on function mult, norm and noyau
   * \param x: matrix; input data n x p
   * \param h: window size (scalar)
   * \param d: scalar product matrix p x p
   * \param n: length x
   * \param p: number of columns of x
   * \param kernel: kernel utilised
   * \param res: matrix returned (allocated in R)
   *
   * \note
   *   matrix x[i,j] (n x p) is substitute
   *   by vector vecteur x[i + j x n]  (i=0..n-1 ; j=0..p-1) 
   *
   * \param result  flag  0 => correct
   *               1 => Pb
   *               2 => Cannot allocate memory
   */
  void W(double *x,double *h,double *d,int *n,int *p,char **kernel,double *res, int * result);


  /** VarRob: compute robust variance. Function called birecly by R
   * \brief Robust variance... gives a low ponderation to isolated 
   * values. This ponderation is determined with kernel and window size.
   * This function depends on function mult, norm and noyau
   * \param x: data matrix n x p
   * \param h: kernel window size (scalar)
   * \param d: matrix of scalar product p x p
   * \param n: lenght of x
   * \param p: length of x
   * \param kernel: kernel used
   * \param res result (matrix)
   * \param result  flag  0 => correct
   *               1 => Pb
   *               2 => Cannot allocate memory
   */
  void VarRob(double *x,double *h,double *d,int *n,int *p,char **kernel,double *res, int * result);

};

#endif
