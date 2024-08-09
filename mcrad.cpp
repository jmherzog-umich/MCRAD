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

#if _MPI
#include <mpi.h>
#endif

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
    int pSize = 1;
    #endif
    
    //Allocate space to hold the settings file
    string args, inputs;
    vector<string> settings;
    
    #if _MPI
    int len1, len2;
    #endif
    
    //Check if we have an input file and process any other arguments in main process
    if (pRank == 0) {
        if (argc > 0)
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
    
    //Store list of strings to write to db at end
    string dbentries;
    
    //Store lists of settings that are changed on subsequent iterations
    vector<vector<string> > ranges;
    vector<int> rangeids;
    
    //Create stringstream and separate key value args into vector, with settings set for the first iteration
    unsigned long int simmax = 1;
    stringstream ss(args), isstream(inputs);
    string buf;
    while (getline(isstream, buf)) {
        if (buf.find(":") != string::npos) {
            //Store all settings in an array
            ranges.push_back(evalrange(buf, pRank, pSize));
            rangeids.push_back(settings.size());
            
            //Set the buffer to the zero index setting
            buf = ranges.back().at(0);
            
            //And update the max number of simulations for this node
            if (ranges.back().size() > simmax)
                simmax = ranges.back().size();
        }
        settings.push_back(buf);
    }
    
    //Run the simulation, and repeat 
    for (unsigned long int simnum = 0; simnum < simmax; simnum ++) {
        //Set this iteration's values
        for (unsigned long int j = 0; j < rangeids.size(); j ++)
            settings.at(rangeids.at(j)) = ranges.at(j).at(simnum);
    
        //Run the simulation
        I.setid(pRank+simnum*pSize);
        I.load(ss);
        I.readArgs(settings);
        I.exec();
        
        //Collect data from results
        dbentries.append(I.getdbline());
    }
    
    #if _MPI
    //First get the size of all the individual strings from each node
    int len3;
    vector<int> dbentryLens(pSize);
    vector<int> offsets(pSize,0);
    
    //Compute and send string lengths
    len3 = dbentries.size();
    MPI_Gather(&len3, 1, MPI_INT, &dbentryLens[0], 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    //Now make a receive buffer string that's long enough based on this
    string recvb;
    if (pRank == 0) {
        int len4 = 1;
        for (unsigned long int i = 0; i < dbentryLens.size(); i ++) {
            if (i > 0)
                offsets.at(i) = offsets.at(i-1) + dbentryLens.at(i-1);
            len4 *= dbentryLens.at(i);
        }
        recvb.resize(len4);
    }
    MPI_Gatherv(const_cast<char*>(dbentries.data()), len3, MPI_CHAR,
                const_cast<char*>(recvb.data()), &dbentryLens[0], &offsets[0],
                MPI_CHAR, 0, MPI_COMM_WORLD);
    
    I.write("", recvb);
    #else
    I.write("", dbentries);
    #endif
    
};
