#ifndef RK4_H_INCLUDED
#define RK4_H_INCLUDED
#include <valarray>

using real_t = double;  // class definition for creating doubles i.a float, single
using array_t = std::valarray<real_t>;
// class definition for creating a array of doubles

// Calling interface for the vector field that defines
// the ordinary differential equation.  Y.SIZE() defines the number of equations
// in the system.  Ensure that DY.SIZE() equals Y.SIZE() for the derivatives.
using vecfield_t = void (*)(real_t t, const array_t& y, array_t& dy);

// Calling interfaces for the fourth-order RK method
extern void rkstep(vecfield_t f, real_t t, real_t h, array_t& y);
extern void rk4(vecfield_t f, array_t& y, real_t tstart, real_t tend,
  size_t nsteps);


#endif // RK4_H_INCLUDED
