#include <iostream>
#include <cmath>
#include "rk4.h"
#include <valarray>
using namespace std;

// MYODE is the name of your actual vector field.  It must have the same
// arguments in the same order as the prototype in rk4.h. Void returns no outputs

// RK4 - integrate the vector field F from TSTART to TEND in NSTEPS equal time steps.

void
rk4(vecfield_t f, array_t& y, real_t tstart, real_t tend, size_t nsteps)
{
  real_t h = (tend-tstart)/nsteps;  // Declaring and evaluating h as a double
  real_t t = tstart;
  size_t n = y.size();  // number of equations
  if(n == 0) return;
  array_t dy1(n), ytemp(n), ytemp3(n), ytemp4(n), k1(n), k2(n),k3(n),k4(n);

  for (size_t i = 0; i < nsteps; i++){

    // Declaring the following variables a double array of size(n)
    f(t, y, k1);    // Evaluating the function for k1

    for(size_t j = 0; j < n; j++){  // first partial RK step
        ytemp[j] = y[j] + (h*0.5)*k1[j];    // Evaluating ytemp for the following equation for k2
   }
    f(t + 0.5 * h ,ytemp ,k2);

    for(size_t j = 0; j < n; j++){
        ytemp3[j] = y[j] + (h*0.5)*k2[j];
   }
    f(t + 0.5 * h ,ytemp3 ,k3);

    for(size_t j = 0; j < n; j++){
        ytemp4[j] = y[j] + h*k3[j];
   }
    f(t + h ,ytemp4 ,k4);

    for(size_t j = 0; j < n; j++){
        y[j] += (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j]) * h / 6.0;
    }
        t = t+h;
  } //iterates rkstep for nsteps, changing t from tstart to tend by +h
  return;
}

