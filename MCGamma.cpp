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

#define CONST_C  2.998e8
#define CONST_PI 3.1415926535

using namespace std;

//Some global defs for random number generation
default_random_engine _GEN;
uniform_real_distribution<double> _DIST(0.0,1.0);
double roll() {
    double eps = 0;
    while (eps == 0 || eps == 1) {
        eps = _DIST(_GEN);
    }
    return eps;
};

//Function to 'remove' element from vector efficiently
template<class T>
void remove(vector<T>& v, unsigned int i) {
    if (i != v.size())
        v.at(i) = v.back();
    v.pop_back();
}

///
//  Performs steady state Monte Carlo simulation of (non-interacting) photons
//  including scattering, energy deposition, fluorescence emission, etc.
//     Units - nm, nm2, nm3, ns
//     Use ion density and convert particle scattering cross-section to per ion basis
///
int main(int argc, char** argv) {
    
    //Define parameters - Incident beam
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
    
    //Handle input - for now just a few shortcuts
    if (argc > 1)
        dens = atof(argv[1]);
    if (argc > 2)
        N0 = atof(argv[2]);
    if (argc > 3)
        Ss = atof(argv[3]);
    if (argc > 4)
        Sa = atof(argv[4]);
    if (argc > 5)
        g = atof(argv[5]);
    if (argc > 6)
        sin0 = atof(argv[6]);
     
    //Fresnel coefficients
    double Rspec = (n0-n)*(n0-n)/(n0+n)/(n0+n);
    double Rx = (nx-n)*(nx-n)/(nx+n)/(nx+n); // n/nx = m -> Rx = (1-m)^2/(1+m)^2 
    double Tspec = 1.0 - Rspec;
    const vec NORM(R, R, L);
            
    //Define output arrays
    cout << "Simulation size: " << R << " x " << L << endl;
    Stats OUT(N0);
    cout << endl << "Z: ";
    for (int i = 0; i < OUT.Zres; i ++)
        cout << scientific << setw(14) << OUT.dZ * i * L;
    cout << endl << "R: ";
    for (int i = 0; i < OUT.Rres; i ++)
        cout << scientific << setw(14) << OUT.dR * i * R;
    cout << endl << endl;
        
    //Print MFPs
    cout << "Scattering mean free path (INCIDENT): " << 1.0/Ss/dens << "nm" << endl;
    cout << "Absorption mean free path (INCIDENT): " << 1.0/Sa/dens << "nm" << endl;
    
    //Initialize a whole lot of temporary variables
    double eps, eps2, s, ds, mu, newZ, dt, newW, tempR, oldW;
    bool done = false;
    
    //Initialize photons
    _GEN.seed(time(0));
    vector<Photon> INCIDENT(N0, Photon(vec(),vec(sin0, 0, sqrt(1-sin0*sin0)), g));
    if (Rb > 1) {
        for (int i = 0; i < N0; i ++) {
            eps = Rb * sqrt(roll());
            eps2 = roll() * 2.0 * CONST_PI;
            INCIDENT.at(i).x.X += eps * cos(eps2);
            INCIDENT.at(i).x.Y += eps * sin(eps2);
        }
    }
    
    //Setup clocks
    auto t0 = chrono::system_clock::now();
    auto t1 = chrono::system_clock::now();
    auto t2 = chrono::system_clock::now();
    std::chrono::duration<double> ELAPSED;
    time_t STARTTIME = chrono::system_clock::to_time_t(t0);
    string STARTTIMESTR(ctime(&STARTTIME));
    STARTTIMESTR = STARTTIMESTR.erase(STARTTIMESTR.back());
    
    //loop through steps
    int STEP = 0; 
    while (STEP <= maxstep) {
        
        //Prepare for loop
        done = false;
        mu = (Sa+Ss)*dens;
        
        //Loop through each photon and update (go forward to delete if needed)
        for (unsigned int i = 0; true; i ++) {   //Don't check vector size cause it will change
        
            //Check that it's a good vector
            if (i >= INCIDENT.size())                                   //Exit if we've covered the whole vector
                break;
            while (INCIDENT.at(i).W <= Wmin) {                          //Loop until W!=0
                if (i < INCIDENT.size()-1) {                            //If not end...
                    OUT.terminate(INCIDENT.at(i).x/NORM, INCIDENT.at(i).t/T);
                    remove(INCIDENT, i);                                // delete and continue
                } else {                                                //Otherwise
                    if (i == 0)
                        INCIDENT.clear();                               // If i==0 and we can't remove any more, then we're at the end
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
            newZ = s*INCIDENT.at(i).mu.Z + INCIDENT.at(i).x.Z;
            while ((newZ>L) || (newZ < 0)) {
                //Look for transmission events (hitting back surface)
                if (newZ>L) {
                    //Find incremental shift, move to boundary and reflect
                    ds = (L-INCIDENT.at(i).x.Z)/INCIDENT.at(i).mu.Z;
                    s -= ds;
                    INCIDENT.at(i).x = INCIDENT.at(i).mu * ds;
                    INCIDENT.at(i).mu.Z = -INCIDENT.at(i).mu.Z;
                    
                    //Calculate reflection coefficient and count reflection in output
                    tempR = OUT.reflect(INCIDENT.at(i).mu, INCIDENT.at(i).W*Tspec, n/nx, true);
                        
                    //update the newZ value
                    INCIDENT.at(i).W = INCIDENT.at(i).W * tempR;
                    newZ = s*INCIDENT.at(i).mu.Z + INCIDENT.at(i).x.Z;
                }
                
                //Look for reflection events (front surface interaction)
                else {
                    //Find incremental shift, move to boundary and reflect
                    ds = -INCIDENT.at(i).x.Z/INCIDENT.at(i).mu.Z;
                    s -= ds;
                    INCIDENT.at(i).x = INCIDENT.at(i).x + INCIDENT.at(i).mu * ds;
                    INCIDENT.at(i).mu.Z = -INCIDENT.at(i).mu.Z;
                    
                    //Calculate reflection coefficient and count reflection in output
                    tempR = OUT.reflect(INCIDENT.at(i).mu, INCIDENT.at(i).W*Tspec, n/n0, false);
                        
                    //update the newZ value
                    INCIDENT.at(i).W = INCIDENT.at(i).W * tempR;
                    newZ = s*INCIDENT.at(i).mu.Z + INCIDENT.at(i).x.Z;
                }
            }
            
            //Normal position update, absorption and scattering
            oldW = INCIDENT.at(i).W;
            newW = oldW * Sa/(Sa+Ss);
            INCIDENT.at(i).x = INCIDENT.at(i).x + INCIDENT.at(i).mu * s;
            INCIDENT.at(i).W -= newW;
            INCIDENT.at(i).t += dt;
            
            //Account for absorbed photons
            OUT.scatter(INCIDENT.at(i).x/NORM, INCIDENT.at(i).t/T, oldW*Tspec, newW*Tspec);
                        
            //Set new direction
            INCIDENT.at(i).Scatter(roll(), roll());
            
            //Terminate old packets
            if (INCIDENT.at(i).W <= Wmin) {
                eps = roll();
                if (eps <= Wm)
                    INCIDENT.at(i).W /= Wm;
                else
                    INCIDENT.at(i).W = 0;
            };
        }
        
        //PRINT DEBUG INFO:
        if (STEP%1 == 0) {
            cerr<<"INC: Step "<<STEP<<"/"<<maxstep;
            if (STEP == 0) {
                cerr << "   Start time: " << STARTTIMESTR <<"     "<<N0<<" particles";
            } else {
                t2 = chrono::system_clock::now();
                ELAPSED = t2 - t0;
                cerr << "   Time: " << ELAPSED.count() << "/";
                ELAPSED = t2 - t1;
                cerr << ELAPSED.count() << " s";
                cerr << "     "<<INCIDENT.size() <<" particles remain";
                t1 = t2;
            } 
            cerr<<endl;
        }
        STEP += 1;
        
        //Check if program done
        if (INCIDENT.size() == 0) {
            cerr<<"Finished calculating incident energy deposition (in "<<STEP<<" steps)"<<endl;
            break;
        }
    }
    
    //Report status
    cout<<"INCIDENT calculation complete in "<<STEP-1<<" steps ("<<maxstep<<" allowed)"<<endl;
    //Fresnel coefficients
    cout <<endl <<endl;
    cout << "Fresnel coefficients:" << endl;
    cout << "  R0  = " << Rspec << endl;
    cout << "  Rx  = " << Rx << endl;
    OUT.print();

};
