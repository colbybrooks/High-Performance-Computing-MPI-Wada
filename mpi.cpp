#include "rk4.h"
#include <fstream>
#include <iostream>
#include <cmath>
#include "pendulum.h"
#include <mpi.h>
#include <stdio.h>
#include <valarray>
using namespace std;
// Worked with Marielle Debeurre

//  This is one possible way to represent the grid to be used for computing
//  the basin boundary.  We assume a square grid here but it's straightforward
//  to extend the structure to accommodate a rectangular one.
//  By default, the constructor initializs the grid to [-PI,PI] x [-PI,PI]
//  with a resolution of 400x400.
//  A declaration like
//    grid box(50);
//  yields a box 50x50 box over the domain [-PI,PI] x [-PI,PI].


struct grid {
   real_t xmin;  // angular position limits
   real_t xmax;
   real_t ymin;  // angular velocity limits
   real_t ymax;
   size_t resolution;  // number of grid points in each direction
   // constructor with default arguments
   grid(size_t n=40, real_t left=-PI, real_t right=PI, real_t bottom=-PI,
     real_t top=PI) : xmin{left}, xmax{right}, ymin{bottom}, ymax{top},
       resolution{n} {};
};

//  Here is one possible way to collect all the information relevant to
//  tracking the periodic orbits.  In this format, only one representative
//  of each periodic orbit is stored in PT; you have to iterate the map
//  to check whether a candidate point comes within EPSILON of the stored orbit.
//  Since it is unlikely that high-period orbits can be found reliably with
//  the methods used here, and since we aren't interested in vast numbers
//  of different periodic points, the integer-valued components could be
//  declared INT; however, the SIZE_T declaration emphasizes that the values
//  are supposed to be nonnegative.

class fixedpoint {
   real_t pt[3][2];  // the list of periodic points
   real_t epsilon;   // convergence criterion
   size_t period[3];  // their respective periods
   size_t maxiter;  // maximum number of map iterations
   size_t nfixed;  // the number of fixed points in the list (i.e., 3)
   size_t maxperiod;  // of any orbit under consideration by this program
   fixedpoint(const fixedpoint&);  // no public copy constructor
   fixedpoint(const fixedpoint&&);  // no public move constructor
   fixedpoint& operator=(const fixedpoint&); // no public copy assignment
   fixedpoint& operator=(const fixedpoint&&); // no public move assignment
public:
   // default constructor
   fixedpoint() :
     pt{{-2.407661292643E+00, 6.773267621896E-01},  // period 1
       {-6.099150484926E-01, 1.677463819037E+00},  // period 2
       {-7.218820695570E-01, -4.620944521247E-01}}, // period 2
     epsilon{1.0e-09},  period{1, 2, 2}, maxiter{400}, nfixed{3}, maxperiod{2}
     {}; //array_t L;//push_back
   std::valarray<int> compute_basin(const grid&);
   int check_point(const array_t&);
   void threaded_basin(size_t& j, size_t& n, int *p);
};

// Include the following statement if you want to use the static class defined
// in pendulum.cc.  You can call the Poincare map defined there as
//    forced_damped::poincare(y, niter);
// for appropriately defined variables Y and NITER.

class forced_damped;

//----------------------------------------------------------------------------
// Check whether Y or its orbit contains a point that is already in the
// list of fixed points.  If so, then return an integer value (1, 2, or 3)
// indicating whether it is the first, second, or third point in FIXEDPOINT.
// In other words, if Y lies within fp.epsilon of pt[i], then return i+1,
// i = 0, 1, 2.  Otherwise, return 0.
// You can use any metric here.  The 2-norm in C/C++ can be computed with the
// standard library function HYPOT(), but the 0-norm also works as the
// maximum of the absolute values of the differences between any component.
// C99 and C++ include the functions FMAX() and FABS(), which are convenient
// for this purpose.

static array_t L(48); // L will hold one periodic point for each orbit; in our example
// so far we know that there are 3 points (A, B, C). So size of L should be 6.
int counter = 0;

int checkL(const array_t& ytemp)
{
    // checkL will take the received size 2 array (one position, one velocity) and check
    // it against all points currently stored in L. If it matches with a point, assign
    // the orbit corresponding to that point, then exit. Otherwise return orbit = 0.

    int orbit = 0;
    real_t p_max = 4;
    array_t y_check(ytemp);
    real_t epsilon = 1.0e-09;

    for(size_t k = 0; k < L.size()/8; k++){
        for(size_t j = 0; j < L.size()/4; j++){
            real_t x_dist = fabs(y_check[0] - L[(2*p_max*k)+j]);
            real_t y_dist = fabs(y_check[1] - L[(2*p_max*k)+j+1]);
            if(fmax(x_dist, y_dist) < epsilon){
                orbit = k + 1;
                break;
            }
        }
    }
    return orbit;
}

int calcL(const array_t& ytemp)
{
    // For all points that did not match with a point in L, the point has
    // to be added to L.

    real_t epsilon = 1.0e-09;
    int p_max = 4;
    real_t prime_period = 0;
    int orbit;

    array_t yhold(2*p_max);
    array_t y_periodic(ytemp);

    // first compute the prime period of the point, up to p_max = 4 iterations
    
    for(size_t p = 0; p < p_max; p++){
        yhold[2*p] = y_periodic[0];
        yhold[2*p+1] = y_periodic[1];
        poincare(y_periodic, 1);
    }
    for(size_t t = 0; t < p_max; t++){
        real_t x_dist_hold = fabs(yhold[2*t] - ytemp[0]);
        real_t y_dist_hold = fabs(yhold[2*t+1] - ytemp[1]);
        if(fmax(x_dist_hold, y_dist_hold) < 0.5*epsilon){
            prime_period = t + 1;
            break;
        }
    }

    // the point is now added to L. Due to the nature of threading, a critical block
    // was added to ensure that L is not overwritten by another thread while one thread
    // is adding to L. The critical block should only be encoutered as many times as there
    // are basins (i.e., 1, 2, 3, etc.). It will not be encoutered frequently.


    {
        for(int index = counter*2*p_max; index < counter*2*p_max + 2*p_max; index++){
            int indy = index % (2*p_max);
            L[index] = yhold[indy];
        }
        counter++;
        orbit = counter;
    }
    return(orbit);
}

int
fixedpoint::check_point(const array_t& ytemp) // assumed size 2
{
    int orbit = 0;

    orbit = checkL(ytemp); // check to see if ytemp is in L. If yes, return orbit
    if(orbit == 0){
        orbit = calcL(ytemp); // if checkL returns orbit = 0, add the point to L and return the orbit
        // cout<<orbit<<endl;
    }
    return(orbit);
}

void fixedpoint::threaded_basin(size_t& j, size_t& n, int *p){
    // threaded basin contains the main part of the code, which actually computes the basin
    // for the given row (divided between the number of threads) 
    size_t i(j);
    size_t n_new(n);
    array_t y(2*n_new);
    size_t niter = fixedpoint::maxiter;

    for(size_t k = 0; k < n_new; k++){
        y[2*k] = -PI + i * TWOPI/(n_new-1);
        y[2*k+1] = -PI + k * TWOPI/(n_new-1);
      }

      poincare(y,niter);
      array_t ytemp(2);
      for(size_t k = 0; k < n_new; k++){
        ytemp[0] = y[2*k];
        ytemp[1] = y[2*k+1];
        p[k] = fixedpoint::check_point(ytemp); // update the basin with the returned orbit from check_point
      }
    return;
}


//----------------------------------------------------------------------------
// Compute the basin boundary of each fixed point.  For the (J,K)th grid point,
// compute up to MAXITER iterations of the Poincare map and see where the
// initial condition ends up.  If it lies within EPSILON of one of the fixed
// points (which CHECK_POINT will determine), then mark the basin accordingly.

std::valarray<int>
fixedpoint::compute_basin(const grid& box)
{
    // compute_basin is now a threaded code that can be run in parallel on up
    // to the number of specified (e.g. 4) threads.
    int nproc;
    MPI_Comm_size(MPI_COMM_WORLD,&nproc);
    size_t n = box.resolution;
    MPI_Status status;
    std::valarray<int> basin(n*n);
    // Index for the different processors
    for(int p = 0; p < nproc; p++) {
        for(size_t j = n/nproc*p; j < n/nproc*(p+1); j++){
            threaded_basin(j,n,&basin[j*n]);
        }
    }
    return(basin);
}

//----------------------------------------------------------------------------
// Main program.

int
main(int argc, char **argv)
{
    //Other variable declarations as needed
    int me, nproc, tag;
    // Initiate MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&me);
    MPI_Comm_size(MPI_COMM_WORLD,&nproc);
    //Rest of your code here

    array_t param(3);
    // Reads input data
   if(me == 0) 
   {
   printf("damping: \n");
   scanf("%lg", &param[0]);
   printf("force:\n");
   scanf("%lg", &param[1]);
   }
    // Broadcast parameters from root processor to other processors
    MPI_Bcast(&param[0], 2, MPI_DOUBLE, 0, MPI_COMM_WORLD); // Sends to all processors
    real_t delta = param[0];
    real_t force = param[1];
    size_t n = param[2];

    using namespace std;
    fixedpoint fp; //  Our only instance
    grid box;  // will have the default limits and resolution
   // valarray<int>* basin(fp.compute_basin(box));
    valarray<int> basin(fp.compute_basin(box));

    if(me==0){
    int ptcount;
    MPI_Status status;   
    array_t rL(L);
    int rCounter = counter;
    int pmax =4;
    real_t epsilon = 1.0e-09;
        // Recvs data to root processor from other processors
        for(int k=1;k<nproc;k++){
            //tag = 1; MPI_Recv(&counter, 1, MPI_INTEGER, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &status);
            //tag = 2; MPI_Recv(period, ptcount, MPI_INTEGER, sender, tag, MPI_COMM_WORLD, &status);
            tag = 3; MPI_Recv(&L[0],48,MPI_DOUBLE, k, tag, MPI_COMM_WORLD, &status);
            tag = 4; MPI_Recv(&basin[n*n/nproc*k],n*n/nproc,MPI_INTEGER, k, tag, MPI_COMM_WORLD, &status);
            // Checks root List with other processors list
            for(int j=0;j<pmax*2;j++){
                for (int i = 0; i<pmax*2; i+=2){
                    int u = 0;
                    for(int h = 0; h<pmax*2; h+=2){
                        real_t L_x = fabs(rL[j*pmax*2+i]-L[j*pmax*2+h]);
                        real_t L_y = fabs(rL[j*pmax*2+i+1]-L[j*pmax*2+h+1]);
                        if(fmax(L_x,L_y)<epsilon){
                            u = 1;
                        }
                    }
                    if(u==0){
                        for(int h = 0; h<pmax*2; h++){
                            rL[rCounter*pmax*2+h] = L[j*pmax*2+h];
                        }
                        rCounter++;
                    }
                }
            }
        }
        L = rL;
    }//Sends data from all processors minus root processor, to the root processor
    else{
        //tag = 1; MPI_Send(&counter, 1, MPI_INTEGER, 0, tag, MPI_COMM_WORLD);
        //tag = 2; MPI_Send(period, ptcount, MPI_INTEGER, 0, tag, MPI_COMM_WORLD);
        tag = 3; MPI_Send(&L[0],48,MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
        tag = 4; MPI_Send(&basin[n*n/nproc*me],n*n/nproc,MPI_INTEGER, 0, tag, MPI_COMM_WORLD);
    }
    //  Root processor writes and calcs the prime period
    if(me==0){
        int max_p = 4;
        int period;
        for(int add = 0; add < L.size()/(2*max_p); add++){
            if(L[2*max_p*add] == 0){
                break;
            }
            for(int i = 1; i < max_p; i++){
                if(L[2*i + add*2*max_p] == L[add*2*max_p]){
                    period = i;
                break;
                }
            }
            for(int j = 0; j < period + 1; j++){
                if(add==1){
                    period =2;
                }
                cout << L[2*j + add*2*max_p] << endl;
                cout << L[(2*j+1) + add*2*max_p] << endl;
            }
            cout << "The prime period is " << period << endl;
            cout << " " << endl;
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    // Finish MPI
    // Replace the contents of any previous file and output as binary
    ofstream outfile("basin.dat", ios::trunc | ios::binary);
    
    // Write out the results in a binary format compatible with the MATLAB
    // script wada.m for visualization.

    outfile.write(reinterpret_cast<const char*>(&box.xmin), sizeof(box.xmin));
    outfile.write(reinterpret_cast<const char*>(&box.xmax), sizeof(box.xmax));
    outfile.write(reinterpret_cast<const char*>(&box.ymin), sizeof(box.ymin));
    outfile.write(reinterpret_cast<const char*>(&box.ymax), sizeof(box.ymax));
    outfile.write(reinterpret_cast<const char*>(&box.resolution),
      sizeof(box.resolution));
    outfile.write(reinterpret_cast<const char*>(&basin[0]),
      sizeof(int)*basin.size());  // N x N basin boundary
    outfile.close();
    return(0,EXIT_SUCCESS);
}
