#if _MPI
#include <mpi.h>
#endif

#include <vector>
#include <cmath>
#include <string>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <iomanip>
#include <random>
#include <chrono>
#include <sstream>

#include "vec.h"
#include "photon.h"
#include "stats.h"
#include "simulation.h"
#include "utility.h"

#if _DEBUG
#include <fenv.h> 
#endif

using namespace std;

///
//  Performs steady state Monte Carlo simulation of (non-interacting) photons
//  including scattering, energy deposition, fluorescence emission, etc.
//     Units - nm, nm2, nm3, ns
///
int main(int argc, char** argv) {
    
    //Debug setup
    #if _DEBUG
    feenableexcept(FE_INVALID | FE_OVERFLOW);
    #endif
    
    //Initialize random number generator and some other things
    rand_init();
    cout.precision(10);
    
    //Create a simulation
    Simulation I;
    
    //Initialize MPI, load data in process 0, pass to all others, then run and collect
    #if _MPI
    int pSize, pRank, nameLen;
    char name[MPI_MAX_PROCESSOR_NAME];
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &pSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &pRank);
    MPI_Get_processor_name(name, &nameLen);
    #else
    int pRank = 0;
    #endif
    
    //Allocate space to hold the settings file
    string args, inputs;
    vector<string> settings;
    
    #if _MPI
    int len1, len2;
    #endif
    
    //Check if we have an input file and process any other arguments in main process
    if (pRank == 0) {
        args = readfile(argv[1]);
        for (int i = 2; i < argc; i ++) {
            inputs.append(argv[i]);
            inputs.append("\n");
        }
        #if _MPI
        len1 = args.size();
        len2 = inputs.size();
        #endif
    }
    
    //Broadcast args file
    #if _MPI
    MPI_Bcast(&len1, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (pRank > 0) args.resize(len1);
    MPI_Bcast(const_cast<char*>(args.data()), len1, MPI_CHAR, 0, MPI_COMM_WORLD);
    
    //Broadcast inputs and create settings
    MPI_Bcast(&len2, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (pRank > 0) inputs.resize(len2);
    MPI_Bcast(const_cast<char*>(inputs.data()), len2, MPI_CHAR, 0, MPI_COMM_WORLD);
    #endif
    
    //Create stringstream and separate key value args into vector
    stringstream ss(args), isstream(inputs);
    string buf;
    while (getline(isstream, buf)) {
        settings.push_back(buf);
    }
    
    //Now load settings from file/string
    I.setid(pRank);
    I.load(ss);
    I.readArgs(settings);
    I.exec();
    
};
