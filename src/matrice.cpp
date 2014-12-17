#include "matrice.h"
#include "matriceTriangle.h"
#include <R.h>

extern "C" {
  /**
   * purpose of this function is to display a matrix, for unitary testing.
   */
  void checkMatrix(double * values, int * nrow, int * ncol)  {
    amap::matrice<double> myMatrix(values, *nrow, *ncol);

    Rprintf("\n");
    for (int i = 0 ; i < *nrow; i++) {
      amap::vecteur<double> row = myMatrix.getRow(i);

      for (int j  = 0 ; j < row.size(); j++) {
	Rprintf("%f\t", row[j]);
      }

      Rprintf("\n");      
    }

  }


  void checkMatrixTriangle(double * values, int * nrow, int * isDiagonal)  {
    amap::matriceTriangle<double> myMatrix(values, *nrow, *isDiagonal);

    Rprintf("\n");
    for (int i = 0 ; i < *nrow; i++) {
      amap::vecteur<double> row = myMatrix.getRow(i);

      for (int j  = 0 ; j < row.size(); j++) {
	Rprintf("%f\t", row[j]);
      }

      Rprintf("\n");      
    }

  }

};
