#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>
#include <R.h>
#include <Rinternals.h>

#include "distance.h"
#include "acprob.h"
#include "hclust.h"
#include "kmeans.h"

void diss (int * data, double * res, int * n, int * p,double * w);
void matind(int * nblev,int * data,int * res, int * n, int * m,int * k);
void checkMatrix(double * values, int * nrow, int * ncol);
void checkMatrixTriangle(double * values, int * nrow, int * isDiagonal) ;

/*
static const  R_CMethodDef cMethods[]  = {
   {"R_distance", (DL_FUNC) &R_distance, 8,{REALSXP,INTSXP, INTSXP,REALSXP,INTSXP, INTSXP,INTSXP, INTSXP}},
   {"rsort_rank_order", (DL_FUNC) &rsort_rank_order, 4,{REALSXP,INTSXP, INTSXP,INTSXP}},
   {"noyau",  (DL_FUNC) &noyau,  3, {REALSXP, STRSXP, REALSXP}},
   {"W",      (DL_FUNC) &W,      8, {REALSXP,REALSXP,REALSXP,INTSXP, INTSXP,STRSXP,REALSXP,INTSXP}},
   {"VarRob", (DL_FUNC) &VarRob, 8, {REALSXP,REALSXP,REALSXP,INTSXP, INTSXP,STRSXP,REALSXP,INTSXP}},
   {"hclust", (DL_FUNC)&hclust, 10, {INTSXP, INTSXP, INTSXP,INTSXP, INTSXP, INTSXP ,REALSXP ,REALSXP ,REALSXP,INTSXP}},
   {"hcluster",(DL_FUNC) &hcluster, 14,{REALSXP,INTSXP, INTSXP, INTSXP,INTSXP, INTSXP, INTSXP,INTSXP, INTSXP,REALSXP ,REALSXP, INTSXP,INTSXP, INTSXP}},
   {"kmeans_Lloyd2",(DL_FUNC) &kmeans_Lloyd2,10, {REALSXP, INTSXP, INTSXP,REALSXP, INTSXP, INTSXP, INTSXP, INTSXP,REALSXP, INTSXP}},
   {NULL, NULL, 0 }
   };*/

static const  R_CMethodDef cMethods[]  = {
  {"matind",           (DL_FUNC) &matind,6},
  {"diss",             (DL_FUNC) &diss, 5},
  {"checkMatrixTriangle", (DL_FUNC)&checkMatrixTriangle, 3},
  {"checkMatrix", (DL_FUNC)&checkMatrix, 3},
  {"R_distance",       (DL_FUNC) &R_distance, 8},
  {"rsort_rank_order", (DL_FUNC) &rsort_rank_order, 4},
  {"noyau",            (DL_FUNC) &noyau,  3},
  {"W",                (DL_FUNC) &W,      8},
  {"VarRob",           (DL_FUNC) &VarRob, 8},
  {"hclust",           (DL_FUNC) &hclust, 10},
  {"hcluster",         (DL_FUNC) &hcluster, 14},
  {"kmeans_Lloyd2",    (DL_FUNC) &kmeans_Lloyd2,10},
  {NULL, NULL, 0 }
};


void
F77_NAME(pnkfmb)(int *fmbvr, int *triabs, int *allsol, int *n, double *couts,
                 int * yasve, int * y, int * renum, double* bornth, int * nbcl0,
		 double* z0,int * nbcl,  double* z, int * nbemp,
		 int * nbdep, int * nbsol, int * nazp);



static const R_FortranMethodDef fMethods[] = {
  {"pnkfmb", (DL_FUNC) &F77_NAME(pnkfmb),  17},
  {NULL, NULL, 0} };




void attribute_visible R_init_amap(DllInfo *info)
{
  R_registerRoutines(info, cMethods,NULL, fMethods, NULL);
  R_useDynamicSymbols(info, FALSE);
  R_forceSymbols(info, TRUE);

}
