#include <iostream>
#include <cmath>
#include "rk4.h"
#include <valarray>
using namespace std;
#include "pendulum.h"

// This code defines the Poincare' map for the forced damped nonlinear
// pendulum, x'' + damping*x' + sin(x) = force*cos(t).
// Of course, you may use the same sampled implementation as in the C
// code, because C++ is largely upwardly compatible with C.
// This version illustrates the use of a "static class" to structure the
// map and its parameters.

// There is only ever one instance of a static class.
// Static classes do not have constructors; instead, their
// data elements are initialized outside of the class as shown below.

// Any class (not just static classes) can declare a STATIC member function.
// The difference between static member functions and ordinary ones is that
// the static member functions do not have a THIS pointer passed implicitly
// among their arguments.  Consequently, static member functions cannot
// access private data elements unless they are also STATIC.

// The restrictions nevertheless work for this application because
// (1) we integrate exactly one forced, damped pendulum and
// (2) the calling interface to the RK4 procedure does not admit a function
//     containing a pointer to an arbitrary data structure.

// One potential advantage of this arrangement over the C version is that
// the parameters DAMPING and FORCE are in their own namespace.  Outside
// of the member functions, they can be accessed only with the qualified
// names FORCED_DAMPED::DAMPING and FORCED_DAMPED::FORCE and so cannot be
// confused with other variables named DAMPING and FORCE that may be defined
// elsewhere in the program.

class forced_damped {
   public:
   static real_t damping;  // map parameters, specified below.
   static real_t force;
};

real_t forced_damped::damping = 0.2;
real_t forced_damped::force = 1.66;

//  The first-order system that defines the forced, damped pendulum
   static void
   pendulum(real_t t, const array_t& y, array_t& dy) {
      size_t neq = y.size();
      for(size_t i = 0; i < neq; i +=2 ){
         dy[i] = y[i+1];   // derivative of position is velocity
         dy[i+1] = forced_damped::force*std::cos(t) - std::sin(y[i]) - forced_damped::damping*y[i+1]; // equation defining pendulum velocity
      }
      return;

   }

   // Member function that defines the Poincare' map
   void poincare(array_t& y, size_t niter)  {
      const size_t NSTEPS = 256;
      size_t neq = y.size();
      for(size_t j = 0; j < niter; j++ )  {
         rk4(pendulum, y, 0.0, TWOPI, NSTEPS); // iterate rk4 solver niter times
            for (size_t i = 0; i < neq; i+=2){
            y[i] = remainder(y[i],TWOPI);
            }
      }

// calculate remainder for all even elements
    //for (size_t i = 0; i < neq; i+=2){
         //y[i] = remainder(y[i],TWOPI); // wrap position back into [-PI,PI] if outside this range
      //}
      return;
   }
