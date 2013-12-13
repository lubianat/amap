/*! \file : acprob.c
 * 
 *
 * \brief  Robust principal component analysis
 *
 * \date Created       : 06/11/02 
 * \date Last Modified : Time-stamp: <2013-12-03 22:05:39 antoine>
 *
 * This Message must remain attached to this source code at all times.
 * This source code may not be redistributed, nor may any executable
 * based upon this source code.
 *
 * \author  Antoine Lucas (antoinelucas@libertysurf.fr)
 *
 * Please notify me of any changes made to the code
 * 
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "R.h"
#ifndef M_PIl
#define M_PIl          3.1415926535897932384626433832795029L  /* pi */
#endif
#include "smartPtr.h"
#include "acprob.h"

/* Compiliation: */
/* R CMD SHLIB acprob.c */


/* Fonction noyau dans R: */
/* dyn.load(paste("prog/acp/src/acprob", .Platform$dynlib.ext, sep = "")) */
/* .C("noyau",as.double(0.6),as.character('g'),res=double(1)) */


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
void noyau(double *u, char **k,double *res)  
{
  double pi= M_PIl;

  switch (**k)
	{
	case 'g' : *res = pow(2 * pi,-0.5) * exp(-pow(*u ,2)/2)   ; break;
	case 'q' : *res = 15.0/16 * pow(1- pow(*u,2),2) * (abs(*u)<1); break;
	case 't' : *res = 35.0/32 * pow(1- pow(*u,2),3) * (abs(*u)<1); break;
	case 'e' : *res =  3.0/4  * (1- pow(*u,2))   * (abs(*u)<1); break;
	case 'c' : *res = pi/4  *cos(*u * pi/2) * (abs(*u)<1); break;
	case 'u' : *res = 1.0/2 * (abs(*u)<1)              ; break;
	}
  /*  return *res; */
}

  /*
  dans R:
  switch(kernel,
        gaussien = (2*pi)^(-1/2) * exp(-u^2/2),
        quartic   = 15/16 * (1-u^2)^2 * (abs(u)<1),
        triweight = 35/32 * (1-u^2)^3 * (abs(u)<1),
        epanechikov = 3/4 * (1-u^2) *   (abs(u)<1),
        cosinus = pi/4 * cos (u*pi/2) * (abs(u)<1),
        uniform = 1/2 * (abs(u)<1),
  */

/** norm
 * \brief compute norm: sqrt(x'.d.x). called by W or VarRob
 * \param x: vector of size p
 * \param p: size of vector x and matrix d
 * \param d: matrix of size pxp
 */
double norm(double *x,int *p,double *d)
     /*
      * x:     vecteur p:1
      * d:     matrice p x p
      * On calcule sqrt( x'.d.x ) 
      */
{
  int i,j;
  double res=0;
  for (i=0; i < *p ; i++)
    for (j=0; j < *p ; j++)
      res += d[i+ j * *p]* x[i]*x[j];
  return sqrt ( res );
}

/** mult
 * \brief multiplication x.x' (return a matrix). called by W or VarRob
 * \param x: vector of size p
 * \param p: size of vector x and matrix d
 * \param res: matrix of result
 */
void mult(double *x,int *p,double *res)
     /*
      * x: vecteur p:1
      * On calcule la matrice x.x' 
      */
{
  int i,j;
  for (i=0; i < *p ; i++)
    for (j=0; j < *p ; j++)
      res[i+ j * *p] = x [i] * x [j]  ;
}



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
void W(double *x,double *h,double *d,int *n,int *p,char **kernel,double *res, int * result)
{
  SmartPtr<double> delta (*p);
  SmartPtr<double> temp (*p * *p);
  double N=0,K=0,som=0;
  int i,j,k,l;
  
  *result = 1; 

  for (l=0; l < (*p * *p) ; l++)
    res[l]=0;
  for (i=0; i < (*n-1) ; i++)
    {
      for (j=i+1; j < *n ; j++)
	{
	  /* delta = Xi-Xj  (ligne i et j de la matrice X) */
	  for (k=0; k < *p ; k++)
	    delta[k]=x[i+(k * *n)]- x[j+(k * *n)];
	  N = norm(delta.get(),p,d)/ *h;

	  /* tmp2 = K ( |delta/h|^2 )  */
	  noyau(&N,kernel,&K);
	  som += K;

	 /*   temp = delta * delta'  (matrice) */
	  mult(delta.get(),p,temp.get());
	  for (l=0; l < (*p * *p) ; l++)
	    res[l] += temp[l] * K ;
	}
    }

  for (l=0; l < (*p * *p) ; l++)
    res[l] = res[l] / som ;


  *result = 0; 
}


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
void VarRob(double *x,double *h,double *d,int *n,int *p,char **kernel,double *res, int * result)
{
  int i,j;
  SmartPtr<double> temp (*p * *p);
  SmartPtr<double> Xi (*p);
  double N=0,K=0,som=0;

  *result = 1;


  som = 0;
  for (i=0; i < *n ; i++)
    {
      for (j=0; j < *p ; j++)
	Xi[j]=x[i+(j * *n)];
     
      N = norm(Xi.get(),p,d)/ *h;
      noyau(&N,kernel,&K);

      mult(Xi.get(),p,temp.get());
      for (j=0; j < (*p * *p) ; j++)
	res[j] += temp[j] * K ;
      som += K;
    }

  for (j=0; j < (*p * *p) ; j++)
    res[j] = res[j] / som ;


  *result = 0;
}


