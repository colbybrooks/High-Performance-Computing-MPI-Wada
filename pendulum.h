#ifndef PENDULUM_H_INCLUDED
#define PENDULUM_H_INCLUDED
#include <cmath>
#include "rk4.h"

const real_t PI = 3.1415926535897932;
const real_t TWOPI = 2*PI;

// The Poincare' map for the forced damped pendulum.
extern void poincare(array_t& y, size_t niter);
#endif // PENDULUM_H_INCLUDED
