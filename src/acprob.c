/*
 * File : acprob.c
 *
 * Description   : Robust principal component analysis
 *
 * Created       : 06/11/02 
 * Last Modified : Time-stamp: <2003-02-12 15:20:33 lucas>
 *
 * This Message must remain attached to this source code at all times.
 * This source code may not be redistributed, nor may any executable
 * based upon this source code.
 *
 * Authored by Antoine Lucas (lucas@toulouse.inra.fr)
 *
 * Please notify me of any changes made to the code
 * 
 * 
 */

#include <stdlib.h>
#include <math.h>
#ifndef M_PIl
#define M_PIl          3.1415926535897932384626433832795029L  /* pi */
#endif


/* Compiliation: */
/* R CMD SHLIB acprob.c */


/* Fonction noyau dans R: */
/* dyn.load(paste("prog/acp/src/acprob", .Platform$dynlib.ext, sep = "")) */
/* .C("noyau",as.double(0.6),as.character('g'),res=double(1)) */


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



int W(double *x,double *h,double *d,int *n,int *p,char **kernel,double *res, int * result)
     /*
      * x: matrice des données n x p
      * h: largeure de la fenetre (scalaire)
      * d: matrice de produit scalaire p x p
      * n: longueur de x
      * p: largeure de x
      * k: noyau utilise
      *
      * Convention:
      *   la matrice x[i,j] (n x p) est remplacée par 
      *   le vecteur x[i + j x n]  (i=0..n-1 ; j=0..p-1) 
      *
      * result  flag  0 => correct
      *               1 => Pb
      *               2 => Cannot allocate memory
      */
{
  double *delta;
  double *temp;
  double N=0,K=0,som=0;
  int i,j,k,l;
  
  *result = 1; 
  delta = (double*) malloc (*p * sizeof(double));
  temp  = (double*) malloc (*p * *p * sizeof(double));

  if(delta == NULL || temp == NULL )
    {
      printf("AMAP: Not enought system memory \n"); 
      *result = 2; 
      return(0);
    }


  for (l=0; l < (*p * *p) ; l++)
    res[l]=0;
  for (i=0; i < (*n-1) ; i++)
    {
      for (j=i+1; j < *n ; j++)
	{
	  /* delta = Xi-Xj  (ligne i et j de la matrice X) */
	  for (k=0; k < *p ; k++)
	    delta[k]=x[i+(k * *n)]- x[j+(k * *n)];
	  N = norm(delta,p,d)/ *h;

	  /* tmp2 = K ( |delta/h|^2 )  */
	  noyau(&N,kernel,&K);
	  som += K;

	 /*   temp = delta * delta'  (matrice) */
	  mult(delta,p,temp);
	  for (l=0; l < (*p * *p) ; l++)
	    res[l] += temp[l] * K ;
	}
    }

  for (l=0; l < (*p * *p) ; l++)
    res[l] = res[l] / som ;

  free(delta);
  free(temp);
  *result = 0; 
}

int VarRob(double *x,double *h,double *d,int *n,int *p,char **kernel,double *res, int * result)
     /*
      * x: matrice des données n x p
      * h: largeure de la fenetre (scalaire)
      * d: matrice de produit scalaire p x p
      * n: longueur de x
      * p: largeure de x
      * k: noyau utilise
      * result  flag  0 => correct
      *               1 => Pb
      *               2 => Cannot allocate memory
      */
{
  int i,j;
  double *temp, *Xi;
  double N=0,K=0,som=0;

  *result = 1;
  temp  = (double*) malloc (*p * *p * sizeof(double));
  Xi = (double*) malloc (*p * sizeof(double));
  if(temp == NULL || Xi == NULL )
    {
      printf("AMAP: Not enought system memory \n"); 
      *result = 2;
      return(0);
    }


  som = 0;
  for (i=0; i < *n ; i++)
    {
      for (j=0; j < *p ; j++)
	Xi[j]=x[i+(j * *n)];
     
      N = norm(Xi,p,d)/ *h;
      noyau(&N,kernel,&K);

      mult(Xi,p,temp);
      for (j=0; j < (*p * *p) ; j++)
	res[j] += temp[j] * K ;
      som += K;
    }

  for (j=0; j < (*p * *p) ; j++)
    res[j] = res[j] / som ;

  free(Xi);
  free(temp);
  *result = 0;
}


/*
main()
{
  double m[8]={1,2,3,4,5,6,7,8};
  double d[4]={1,0,0,1};
  int p=2,n=4;
  double h=1;
  char k1='k';
  char *k2;
  char **k3;
  double res[4];

  k2 = (char*) malloc ( sizeof(char));
  k3 = (char**) malloc ( sizeof(char* ));
  if(k2 == NULL || k3 == NULL )
    {
      printf("AMAP: Not enought system memory \n"); 
      return(0);
    }


  k2 = &k1;
  k3 = &k2;
  W(m,&h,d,&n,&p,k3,res);

}
*/
