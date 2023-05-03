#include <sstream>
#include <string>
#include <fstream>
#include <iostream>

#ifndef _SIMINFO_H_
#define _SIMINFO_H_

using namespace std;

struct SimInfo {
    //Incident beam
    double Ss=5e5;                  //Scattering cross-section
    double Sa=0e0;                  //absorption cross-section
    double g=0.99;                  //Scattering anisotropy
    double n0=1.600;                //Refractive index of front surface (SiO2)
    double n=1.330;                 //Refractive index of medium (water)
    double nx=1.600;                //Refractive index of back surface (SiO2)
	double sin0=0.0;                //Incident sin(theta) = 0 for normal incidence (-1 to 1)
    
    //Sample info
    int N0 = 5e6;           //Number of initial photon groups
    double dens=1e-13;      //Number density in bulk material
    double L=1e7;           //Length of sample [nm]
    double R=1e7;           //Radius of sample [nm]
    double T=300;           //Length of time to integrate [ns]
    
    //Beam parameters
    double Rb=5e5;          //Beam radius [nm]
    
    //Convergence and other settings
    double Wmin = 1e-10;    //Start terminating packets if they get this small
    double Wm = 0.1;        //Weighting for Roulette termination procedure
    int maxstep = 5e6;      //Maximum number of steps before terminating program
    
    //Calculated values
    double WSCALE;
    
    void load(const string& fname);
};

void SimInfo::load(const string& fname) {

        //Allocate some memory and open the file
        string cmd, key;
        stringstream cmdstr;
        ifstream ifile(fname);
        
        //Process line by line
        while (!ifile.eof()) {
            //Get the command
            getline(ifile, cmd);
            cmdstr.str(cmd);
            cmdstr >> key;
            
            cerr << cmd << "   KEY = " << key << endl;
            
            //Store value
            if (!key.compare("Ss"))
                cmdstr >> Ss;
            else if (!key.compare("Sa"))
                cmdstr >> Sa;
            else if (!key.compare("g"))
                cmdstr >> g;
            else if (!key.compare("n0"))
                cmdstr >> n0;
            else if (!key.compare("n"))
                cmdstr >> n;
            else if (!key.compare("nx"))
                cmdstr >> nx;
            else if (!key.compare("sin0"))
                cmdstr >> sin0;
            else if (!key.compare("N0"))
                cmdstr >> N0;
            else if (!key.compare("dens"))
                cmdstr >> dens;
            else if (!key.compare("L"))
                cmdstr >> L;
            else if (!key.compare("R"))
                cmdstr >> R;
            else if (!key.compare("T"))
                cmdstr >> T;
            else if (!key.compare("Rb"))
                cmdstr >> Rb;
            else if (!key.compare("Wm"))
                cmdstr >> Wm;
            else if (!key.compare("Wmin"))
                cmdstr >> Wmin;
            else if (!key.compare("maxstep"))
                cmdstr >> maxstep;
                
            //Clear values
            key.clear();
            cmd.clear();
            cmdstr.clear();
        }
        
        //Close file
        ifile.close();
        return;
}

#endif
