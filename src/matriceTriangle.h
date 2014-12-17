#ifndef AMAP_TRIANGULAR_MATRIX
#define AMAP_TRIANGULAR_MATRIX

#include "matrice.h"

namespace amap {
  /**
   * Matrix data.
   *
   * a triangular matrix with 4 row and 4 cols
   * +---+---+---+---+
   * |   |   |   |   |
   * +---+---+---+---+
   * | 1 |   |   |   |
   * +---+---+---+---+
   * | 2 | 4 |   |   |
   * +---+---+---+---+
   * | 3 | 5 | 6 |   |
   * +---+---+---+---+
   * 
   * matrix is  4x4 but dataSet is of a size 6 (when not diagonal).
   * dataSet is of size 10 if matrix is diagonal.
   * 
   */

  template<class T> class matriceTriangle : public matrice<T> {

  private:
    T blankValue;
    bool isDiagonal;
  public:
    /**
     * Contructor.
     * \param values_p the data matrix
     * \param nrows_p the number of rows
     */
    matriceTriangle(T * values_p, int nrow_p, bool isDiagonal_p) 
      : matrice<T>(values_p, nrow_p, nrow_p),
	isDiagonal(isDiagonal_p)
      {
      };
  
    
    /**
     * Accessor on data.
     */
    virtual T & operator[] (int index) {
      blankValue = 0;


      int i = index % this->getNrow();
      int j = index / this->getNrow();

      if (i == j && !isDiagonal) {
	return blankValue;
      }

      if (i <= j)
	{
	  int temp = i;
	  i = j;
	  j = temp;
	}
      int maxRowSize = this->getNrow();
      if (!isDiagonal) {
	i--;
	maxRowSize--;
      }

      //
      // j * *(j + 1) / 2 : nombres de points qui sont dans "le triangle du haut".
      //

      int ij = i + j * maxRowSize - j * (j + 1) / 2;


      return matrice<T>::operator[](ij);
    };


  };


};


#endif
