#include <sstream>
#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <random>
#include <chrono>

#include "stats.h"
#include "vec.h"
#include "photon.h"
#include "utility.h"

#ifndef _SIMINFO_H_
#define _SIMINFO_H_

#define CONST_C  2.998e8
#define CONST_PI 3.1415926535

using namespace std;

class Simulation {
    public:
        Simulation();
        void load(const string& fname);
        void setup();
        void run();
        void genBeam();
        void genFluor();
        void print();
        
        enum SimFlags {
            BackWall      = 0b00000001,     //Whether there is a backwall in the simulation or semi-infinite
            FrontWall     = 0b00000010,     //Whether there is a front wall or semi-infinite
            RadialWall    = 0b00000100,     //Whether there is a cylindrical wall or infinite
            GaussBeam     = 0b00001000     //Whether beam is Gaussian or a top-hat
        };
    
    private:
        //Incident beam
        double Ss=5e5;                  //Scattering cross-section
        double Sa=5e5;                  //absorption cross-section
        double g=0.98;                  //Scattering anisotropy
        double n0=1.600;                //Refractive index of front surface (glass)
        double n=1.330;                 //Refractive index of medium (water)
        double nx=1.600;                //Refractive index of back surface (glass)
        double sin0=0.0;                //Incident sin(theta) = 0 for normal incidence (-1 to 1)
        double FQY=1;                   //Fluorescence quantum yield
        
        //Sample info
        int N0 = 5e6;           //Number of initial photon groups
        double dens=1e-13;      //Number density in bulk material
        double L=1e7;           //Length of sample [nm]
        double R=1e7;           //Radius of sample [nm]
        double T=300;           //Length of time to integrate [ns]
        
        //Sim Flags
        SimFlags flags = (Simulation::SimFlags)3;
        
        //Beam parameters
        double Rb=5e5;          //Beam radius [nm]
        
        //Convergence and other settings
        double Wmin = 1e-10;    //Start terminating packets if they get this small
        double Wm = 0.1;        //Weighting for Roulette termination procedure
        int maxstep = 5e6;      //Maximum number of steps before terminating program
        
        //Stats settings
        int Zres = 100;
        int Rres = 100;
        int Tres = 1;
        int THres = 50;
        int momentlvl = 4;
        
        //Calculated values
        double WSCALE;
        Stats OUT;
        vector<Photon> PHOTONS;
        
        //Random values
        default_random_engine _GEN;
        uniform_real_distribution<double> _DIST;
        
        //Methods
        double roll();
};

Simulation::Simulation() {
    Ss=5e5; Sa=5e5; g=0.98; n0=1.600; n=1.330; nx=1.600; sin0=0.0; FQY=1; 
    N0 = 5e6; dens=1e-13; L=1e7; R=1e7; T=300; 
    Rb=5e5;
    Wmin = 1e-10; Wm = 0.1; maxstep = 5e6; 
    Zres = 100; Rres = 100; Tres = 1; THres = 50; momentlvl = 4;
    flags = (Simulation::SimFlags)3;
}

template<class T>
void remove(vector<T>& v, unsigned int i) {
    if (i != v.size())
        v.at(i) = v.back();
    v.pop_back();
}

double Simulation::roll() {
    double eps = 0;
    while (eps == 0 || eps == 1) {
        eps = _DIST(_GEN);
    }
    return eps;
};

void Simulation::print() {
    OUT.print();
}

void Simulation::genBeam() {
    //Alert the user
    cout << "Generating incident beam photons" << endl;
    
    //Initialize some values
    PHOTONS = vector<Photon>(N0, Photon(vec(),vec(sin0, 0, sqrt(1-sin0*sin0)), g));
    if (flags & SimFlags::FrontWall)
        WSCALE = (1.0 - (n0-n)*(n0-n)/(n0+n)/(n0+n));
    else
        WSCALE = 1.0;
    double eps, eps2;
    
    //Loop through and generate beam if needed
    if (Rb > 1) {
        for (int i = 0; i < N0; i ++) {
            if (flags & SimFlags::GaussBeam)
                eps = Rb * erfinvf((float)roll());
            else
                eps = Rb * sqrt(roll());
            eps2 = roll() * 2.0 * CONST_PI;
            PHOTONS.at(i).x.X += eps * cos(eps2);
            PHOTONS.at(i).x.Y += eps * sin(eps2);
        }
    }
}

void Simulation::genFluor() {
    //Alert the user
    cout << "Generating fluorescence photons" << endl;
    
    //Initialize some values
    PHOTONS = vector<Photon>(N0, Photon(vec(),vec(), g));
    WSCALE = FQY;
    double eps, eps2, eps3, r1, r2;
    int Nrz;
    int i = 0;
    
    //Loop through each bin
    for (int iR = 0; iR < OUT.Rres; iR ++) {
        for (int iZ = 0; iZ < OUT.Zres; iZ ++) {
            
            //Figure out how many photons inside this element
            Nrz = (int) (N0*OUT.DEP.at(iZ*OUT.Rres + iR)/OUT._ADEP)+0.5;
            
            //Loop through each photon in the bin
            for (int ii = 0; ii < Nrz; ii ++) {
            
                //Randomize direction
                eps = roll();
                eps2 = roll() * 2.0 * CONST_PI;
                PHOTONS.at(i).mu.X = sqrt(1-eps*eps) * cos(eps2);
                PHOTONS.at(i).mu.Y = sqrt(1-eps*eps) * sin(eps2);
                PHOTONS.at(i).mu.Z = eps * cos(eps2);
                
                //Get coordinates
                r1 = iR * OUT.dR * R;
                r2 = (iR+1) * OUT.dR * R;
                
                //Randomize position inside volume element
                eps = sqrt(roll()*(r2*r2-r1*r1) + r2*r2);
                eps2 = roll() * 2.0 * CONST_PI;
                eps3 = iZ * OUT.dZ * L + (OUT.dZ * L) * roll();
                PHOTONS.at(i).x.X = eps * cos(eps2);
                PHOTONS.at(i).x.Y = eps * sin(eps2);
                PHOTONS.at(i).x.Z = eps3;
                
                //Increment photon count
                i++;
            }
        }
    }
}

void Simulation::setup() {
    //Create output stats block
    OUT = Stats(N0, Zres, Rres, Tres, THres, momentlvl);
    
    //Setup random number generator
    _DIST = uniform_real_distribution<double>(0.0,1.0);
    _GEN.seed(time(0));
    WSCALE = 1;
    
    //Fresnel coefficients
    double Rspec = (n0-n)*(n0-n)/(n0+n)/(n0+n);
    double Rx = (nx-n)*(nx-n)/(nx+n)/(nx+n); // I.n/I.nx = m -> Rx = (1-m)^2/(1+m)^2 
    double Tspec = 1.0 - Rspec;
    
    //Print some settings
    cout << endl << endl;
    cout << "==================================================================" << endl;
    cout << "Settings" << endl;
    cout << "==================================================================" << endl;
    cout << "Phase function: ";
    switch (Photon::phase) {
        case Photon::PhaseFunction::HenyeyGreenstein:
            cout << "Henyey-Greenstein" << endl;
            break;
        case Photon::PhaseFunction::Rayleigh:
            cout << "Rayleigh" << endl;
            break;
        default:
            cout << "Other" << endl;
            break;
    }
    cout << "g: " << g << endl;
    cout << "Back Wall: " << ((flags & SimFlags::BackWall) ? "True" : "False" ) << endl;
    cout << "Front Wall: " << ((flags & SimFlags::FrontWall) ? "True" : "False" ) << endl;
    cout << "Radial Wall: " << ((flags & SimFlags::RadialWall) ? "True" : "False" ) << endl;
    cout << "Gaussian Beam Profile: " << ((flags & SimFlags::GaussBeam) ? "True" : "False" ) << endl;
    cout << "Photon packets: " << N0 << endl;
    cout << "Scattering cross-section: " << Ss << endl;
    cout << "Absorption cross-section: " << Sa << endl;
    cout << "Front Refractive Index: " << n0 << endl;
    cout << "Medium Refractive Index: " << n << endl;
    cout << "Back Refractive Index: " << nx << endl;
    cout << "Incident sin(theta): " << sin0 << endl;
    cout << endl << endl;
    
    //Warn the user
    cout << "==================================================================" << endl;
    cout << "Starting simulation" << endl;
    cout << "==================================================================" << endl;
    
    //Define output arrays
    cout << "Simulation size: " << R << " x " << L << endl;
    cout << endl << "Z: ";
    for (int i = 0; i < OUT.Zres; i ++)
        cout << scientific << setw(14) << OUT.dZ * i * L;
    cout << endl << "R: ";
    for (int i = 0; i < OUT.Rres; i ++)
        cout << scientific << setw(14) << OUT.dR * i * R;
    cout << endl << endl;
        
    //Fresnel coefficients
    cout <<endl <<endl;
    cout << "Fresnel coefficients:" << endl;
    cout << "  R0  = " << Rspec << endl;
    cout << "  Rx  = " << Rx << endl;
    cout << "  T0  = " << Tspec << endl;
        
    //Print MFPs
    cout << "Scattering mean free path: " << 1.0/Ss/dens << "nm" << endl;
    cout << "Absorption mean free path: " << 1.0/Sa/dens << "nm" << endl;
}

void Simulation::run() {
    
    //Setup clocks
    auto t0 = chrono::system_clock::now();
    auto t1 = chrono::system_clock::now();
    auto t2 = chrono::system_clock::now();
    std::chrono::duration<double> ELAPSED;
    time_t STARTTIME = chrono::system_clock::to_time_t(t0);
    string STARTTIMESTR(ctime(&STARTTIME));
    STARTTIMESTR = STARTTIMESTR.erase(STARTTIMESTR.back());
    
    //Initialize a whole lot of temporary variables
    double eps, s, ds, mu, newZ, dt, newW, tempR, oldW;
    bool done = false, donestep = false;
    const vec NORM(R, R, L);
    
    //loop through steps
    int STEP = 0; 
    while (STEP <= maxstep) {
        
        //Prepare for loop
        done = false;
        mu = (Sa+Ss)*dens;
        
        //Loop through each photon and update (go forward to delete if needed)
        for (unsigned int i = 0; true; i ++) {   //Don't check vector size cause it will change
        
            //Check that it's a good vector
            if (i >= PHOTONS.size())                                   //Exit if we've covered the whole vector
                break;
            while (PHOTONS.at(i).W <= Wmin) {                          //Loop until W!=0
                if (i < PHOTONS.size()-1) {                            //If not end...
                    OUT.terminate(PHOTONS.at(i).x/NORM, PHOTONS.at(i).t/T);
                    remove(PHOTONS, i);                                // delete and continue
                } else {                                                //Otherwise
                    if (i == 0)
                        PHOTONS.clear();                               // If i==0 and we can't remove any more, then we're at the end
                    done = true;                                        // propagate "we're done" downstream
                    break;                                              // stop looking for photons
                }
            }
            if (done)
                break;                              // and stop simulation
        
            //Roll the dice
            eps = roll();
            
            //Calculate new displacement then time
            s = -log(eps)/mu; dt = s/CONST_C*n;
            
            //Check for boundary collision - only Z-direction for simplicity
            newZ = s*PHOTONS.at(i).mu.Z + PHOTONS.at(i).x.Z;
            donestep = false;
            while (!donestep) {
                //Look for transmission events (hitting back surface)
                if ( (flags & SimFlags::BackWall) and (newZ>L)) {
                    //Find incremental shift, move to boundary and reflect
                    ds = (L-PHOTONS.at(i).x.Z)/PHOTONS.at(i).mu.Z;
                    s -= ds;
                    PHOTONS.at(i).x = PHOTONS.at(i).mu * ds;
                    PHOTONS.at(i).mu.Z = -PHOTONS.at(i).mu.Z;
                    
                    //Calculate reflection coefficient and count reflection in output
                    tempR = OUT.reflect(PHOTONS.at(i).mu, PHOTONS.at(i).W*WSCALE, n/nx, true);
                        
                    //update the newZ value
                    PHOTONS.at(i).W = PHOTONS.at(i).W * tempR;
                    newZ = s*PHOTONS.at(i).mu.Z + PHOTONS.at(i).x.Z;
                }
                
                //Look for reflection events (front surface interaction)
                else if ( (flags & SimFlags::FrontWall) and (newZ < 0)) {
                    //Find incremental shift, move to boundary and reflect
                    ds = -PHOTONS.at(i).x.Z/PHOTONS.at(i).mu.Z;
                    s -= ds;
                    PHOTONS.at(i).x = PHOTONS.at(i).x + PHOTONS.at(i).mu * ds;
                    PHOTONS.at(i).mu.Z = -PHOTONS.at(i).mu.Z;
                    
                    //Calculate reflection coefficient and count reflection in output
                    tempR = OUT.reflect(PHOTONS.at(i).mu, PHOTONS.at(i).W*WSCALE, n/n0, false);
                        
                    //update the newZ value
                    PHOTONS.at(i).W = PHOTONS.at(i).W * tempR;
                    newZ = s*PHOTONS.at(i).mu.Z + PHOTONS.at(i).x.Z;
                }
                
                //No walls or we haven't hit a wall, so continue
                else
                    donestep = true;
            }
            
            //Normal position update, absorption and scattering
            oldW = PHOTONS.at(i).W;
            newW = oldW * Sa/(Sa+Ss);
            PHOTONS.at(i).x = PHOTONS.at(i).x + PHOTONS.at(i).mu * s;
            PHOTONS.at(i).W -= newW;
            PHOTONS.at(i).t += dt;
            
            //Account for absorbed photons
            OUT.scatter(PHOTONS.at(i).x/NORM, PHOTONS.at(i).t/T, oldW*WSCALE, newW*WSCALE);
                        
            //Set new direction
            PHOTONS.at(i).Scatter(roll(), roll());
            
            //Terminate old packets
            if (PHOTONS.at(i).W <= Wmin) {
                eps = roll();
                if (eps <= Wm)
                    PHOTONS.at(i).W /= Wm;
                else
                    PHOTONS.at(i).W = 0;
            };
        }
        
        //PRINT DEBUG INFO:
        if (STEP%1 == 0) {
            cerr<<"Step "<<STEP<<"/"<<maxstep;
            if (STEP == 0) {
                cerr << "   Start time: " << STARTTIMESTR <<"     "<<N0<<" particles";
            } else {
                t2 = chrono::system_clock::now();
                ELAPSED = t2 - t0;
                cerr << "   Time: " << ELAPSED.count() << "/";
                ELAPSED = t2 - t1;
                cerr << ELAPSED.count() << " s";
                cerr << "     "<<PHOTONS.size() <<" particles remain";
                t1 = t2;
            } 
            cerr<<endl;
        }
        STEP += 1;
        
        //Check if program done
        if (PHOTONS.size() == 0) {
            cerr<<"Finished calculating energy deposition (in "<<STEP<<" steps)"<<endl;
            break;
        }
    }
    
    //Report status
    cout<<"Calculation complete in "<<STEP-1<<" steps ("<<maxstep<<" allowed)"<<endl;
}

void Simulation::load(const string& fname) {

        //Allocate some memory and open the file
        string cmd, key;
        stringstream cmdstr;
        ifstream ifile(fname);
        int tmpint, tmpflags = (int)flags;
        bool tmpbool;
        
        //Process line by line
        while (!ifile.eof()) {
            //Get the command
            getline(ifile, cmd);
            cmdstr.str(cmd);
            cmdstr >> key;
            
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
            else if (!key.compare("Zres"))
                cmdstr >> Zres;
            else if (!key.compare("Rres"))
                cmdstr >> Rres;
            else if (!key.compare("Tres"))
                cmdstr >> Tres;
            else if (!key.compare("THres"))
                cmdstr >> THres;
            else if (!key.compare("moments"))
                cmdstr >> momentlvl;
            else if (!key.compare("FQY"))
                cmdstr >> FQY;
            //More complicated settings
            else if (!key.compare("phase")) {
                cmdstr >> tmpint;
                Photon::phase = (Photon::PhaseFunction)tmpint;
            } else if (!key.compare("backwall")) {
                cmdstr >> tmpbool;
                if (tmpbool)
                    tmpflags |= SimFlags::BackWall;
                else
                    tmpflags &= ~(SimFlags::BackWall);
                flags = (SimFlags) tmpflags;
            } else if (!key.compare("frontwall")) {
                cmdstr >> tmpbool;
                if (tmpbool)
                    tmpflags |= SimFlags::FrontWall;
                else
                    tmpflags &= ~(SimFlags::FrontWall);
                flags = (SimFlags) tmpflags;
            } else if (!key.compare("cylwall")) {
                cmdstr >> tmpbool;
                if (tmpbool)
                    tmpflags |= SimFlags::RadialWall;
                else
                    tmpflags &= ~(SimFlags::RadialWall);
                flags = (SimFlags) tmpflags;
            }
            else if (!key.compare("gaussbeam")) {
                cmdstr >> tmpbool;
                if (tmpbool)
                    tmpflags |= SimFlags::GaussBeam;
                else
                    tmpflags &= ~(SimFlags::GaussBeam);
                flags = (SimFlags) tmpflags;
            }
            //Commands to run the program
            else if (!key.compare("run"))
                run();
            else if (!key.compare("setup"))
                setup();
            else if (!key.compare("genBeam"))
                genBeam();
            else if (!key.compare("genFluor"))
                genFluor();
            else if (!key.compare("print"))
                print();
                
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
