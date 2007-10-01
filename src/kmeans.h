
#ifndef KMEANS_HEADER_AMAP
#define KMEANS_HEADER_AMAP 1


extern "C" 
{
  
  void kmeans_Lloyd2(double *x, int *pn, int *pp, double *cen, int *pk, int *cl, 
		     int *pmaxiter, int *nc, double *wss, int * method);
    
}





#endif
