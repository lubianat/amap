

/** matind  
 * \brief create indicatrices instead of factors.
 * 
 * \param nblev number of levels for each variables (size 1xm)
 * \param data: input data n individuals with m variables
 * \param res: return matrix, indicatrics of size nxsum(nblev)
 * \param n,m size of matrix data.
 * 
 *
  color   size           color.blue color.red size.large size.medium size.small
   blue  large                    1         0          1           0          0
    red  large		          0         1          1           0          0
    red  small     => 	          0         1          0           0          1
   blue medium		          1         0          0           1          0
    red  large		          0         1          1           0          0

    nblev = 2,3

    x  = 1 1 
         2 1
         2 2  
         1 3
         2 1
 */

void matind(int * nblev,int * data,int * res, int * n, int * m,int * k)
{
  int i,j,curr_col_in_res=0;

  /* for all variables ... */
  for(j = 0; j< *m; j++)
    {
      /* for all individuals ... */
      for(i = 0; i< *n; i++)
	{
	  if((i + (*n)*(curr_col_in_res+ data[i+ (*n)*j] -1)) >= 0)
	    {
	      /*	  res[i,(curr_col_in_res+ data[i,j] -1)] ++;*/
	      res[i + (*n)*(curr_col_in_res+ data[i+ (*n)*j] -1)] ++;	  
	    }
	}
      curr_col_in_res  = curr_col_in_res + nblev[j];
    }

}
