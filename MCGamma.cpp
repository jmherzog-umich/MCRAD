#include <vector>
#include <cmath>
#include <string>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <iomanip>
#include <random>
#include <chrono>

#include "vec.h"
#include "photon.h"
#include "stats.h"
#include "simulation.h"
#include "utility.h"

using namespace std;

///
//  Performs steady state Monte Carlo simulation of (non-interacting) photons
//  including scattering, energy deposition, fluorescence emission, etc.
//     Units - nm, nm2, nm3, ns
///
int main(int argc, char** argv) {
    
    //Initialize random number generator and some other things
    rand_init();
    cout.precision(10);
    
    //Create a simulation
    Simulation I;
    
    //Check if we have an  input file
    if (argc > 1) {
        //Read file
        I.load(argv[1]);
    } else {
        //Run simulation
        I.genBeam();
        I.setup();
        I.run();
        I.print();
    }
};
