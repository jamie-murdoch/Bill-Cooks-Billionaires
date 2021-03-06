/****************************************************************************/
/*                                                                          */
/*                      Utility Functions for CO759                         */
/*                                                                          */
/****************************************************************************/

#ifndef __CO759_UTIL_H
#define __CO759_UTIL_H

#define LP_EPSILON 0.000001

#define CO759_SWAP(x,y,temp) {temp = x; x = y; y = temp;}

#include <vector>

double CO759_zeit (void);
double CO759_real_zeit (void);

int CO759_build_xy (int ncount, std::vector<double>& xlist, std::vector<double>& ylist, int gridsize);
int CO759_build_xy (int ncount, double *xlist, double *ylist, int gridsize);

bool is_almost_integral(double x);

#endif  /* __CO759_UTIL_H */
