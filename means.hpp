/*We compute and return the arithemetic mean, the geometric mean, and the harmonic mean*/

#ifndef _StatsLibExtended_means_HPP
#define _StatsLibExtended_means_HPP

template <typename T>
double mean(T &x);

template <typename T>
double mean(T &x, const T &b);


#include "means.ipp"
#endif
