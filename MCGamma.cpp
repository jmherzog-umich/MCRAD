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
#include "siminfo.h"

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

vector<Photon> genBeam(SimInfo& I) {
    vector<Photon> INCIDENT (I.N0, Photon(vec(),vec(I.sin0, 0, sqrt(1-I.sin0*I.sin0)), I.g));
    double eps, eps2;
    if (I.Rb > 1) {
        for (int i = 0; i < I.N0; i ++) {
            eps = I.Rb * sqrt(roll());
            eps2 = roll() * 2.0 * CONST_PI;
            INCIDENT.at(i).x.X += eps * cos(eps2);
            INCIDENT.at(i).x.Y += eps * sin(eps2);
        }
    }
    return INCIDENT;
}

vector<Photon> genFluor(SimInfo& I, Stats& DEP) {

}

void RunSim(SimInfo& I, Stats& OUT, vector<Photon>& PHOTONS) {
    //Setup clocks
    auto t0 = chrono::system_clock::now();
    auto t1 = chrono::system_clock::now();
    auto t2 = chrono::system_clock::now();
    std::chrono::duration<double> ELAPSED;
    time_t STARTTIME = chrono::system_clock::to_time_t(t0);
    string STARTTIMESTR(ctime(&STARTTIME));
    STARTTIMESTR = STARTTIMESTR.erase(STARTTIMESTR.back());
    
    //Initialize a whole lot of temporary variables
    double eps, eps2, s, ds, mu, newZ, dt, newW, tempR, oldW;
    bool done = false;
    const vec NORM(I.R, I.R, I.L);
    
    //loop through steps
    int STEP = 0; 
    while (STEP <= I.maxstep) {
        
        //Prepare for loop
        done = false;
        mu = (I.Sa+I.Ss)*I.dens;
        
        //Loop through each photon and update (go forward to delete if needed)
        for (unsigned int i = 0; true; i ++) {   //Don't check vector size cause it will change
        
            //Check that it's a good vector
            if (i >= PHOTONS.size())                                   //Exit if we've covered the whole vector
                break;
            while (PHOTONS.at(i).W <= I.Wmin) {                          //Loop until W!=0
                if (i < PHOTONS.size()-1) {                            //If not end...
                    OUT.terminate(PHOTONS.at(i).x/NORM, PHOTONS.at(i).t/I.T);
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
            s = -log(eps)/mu; dt = s/CONST_C*I.n;
            
            //Check for boundary collision - only Z-direction for simplicity
            newZ = s*PHOTONS.at(i).mu.Z + PHOTONS.at(i).x.Z;
            while ((newZ>I.L) || (newZ < 0)) {
                //Look for transmission events (hitting back surface)
                if (newZ>I.L) {
                    //Find incremental shift, move to boundary and reflect
                    ds = (I.L-PHOTONS.at(i).x.Z)/PHOTONS.at(i).mu.Z;
                    s -= ds;
                    PHOTONS.at(i).x = PHOTONS.at(i).mu * ds;
                    PHOTONS.at(i).mu.Z = -PHOTONS.at(i).mu.Z;
                    
                    //Calculate reflection coefficient and count reflection in output
                    tempR = OUT.reflect(PHOTONS.at(i).mu, PHOTONS.at(i).W*I.WSCALE, I.n/I.nx, true);
                        
                    //update the newZ value
                    PHOTONS.at(i).W = PHOTONS.at(i).W * tempR;
                    newZ = s*PHOTONS.at(i).mu.Z + PHOTONS.at(i).x.Z;
                }
                
                //Look for reflection events (front surface interaction)
                else {
                    //Find incremental shift, move to boundary and reflect
                    ds = -PHOTONS.at(i).x.Z/PHOTONS.at(i).mu.Z;
                    s -= ds;
                    PHOTONS.at(i).x = PHOTONS.at(i).x + PHOTONS.at(i).mu * ds;
                    PHOTONS.at(i).mu.Z = -PHOTONS.at(i).mu.Z;
                    
                    //Calculate reflection coefficient and count reflection in output
                    tempR = OUT.reflect(PHOTONS.at(i).mu, PHOTONS.at(i).W*I.WSCALE, I.n/I.n0, false);
                        
                    //update the newZ value
                    PHOTONS.at(i).W = PHOTONS.at(i).W * tempR;
                    newZ = s*PHOTONS.at(i).mu.Z + PHOTONS.at(i).x.Z;
                }
            }
            
            //Normal position update, absorption and scattering
            oldW = PHOTONS.at(i).W;
            newW = oldW * I.Sa/(I.Sa+I.Ss);
            PHOTONS.at(i).x = PHOTONS.at(i).x + PHOTONS.at(i).mu * s;
            PHOTONS.at(i).W -= newW;
            PHOTONS.at(i).t += dt;
            
            //Account for absorbed photons
            OUT.scatter(PHOTONS.at(i).x/NORM, PHOTONS.at(i).t/I.T, oldW*I.WSCALE, newW*I.WSCALE);
                        
            //Set new direction
            PHOTONS.at(i).Scatter(roll(), roll());
            
            //Terminate old packets
            if (PHOTONS.at(i).W <= I.Wmin) {
                eps = roll();
                if (eps <= I.Wm)
                    PHOTONS.at(i).W /= I.Wm;
                else
                    PHOTONS.at(i).W = 0;
            };
        }
        
        //PRINT DEBUG INFO:
        if (STEP%1 == 0) {
            cerr<<"Step "<<STEP<<"/"<<I.maxstep;
            if (STEP == 0) {
                cerr << "   Start time: " << STARTTIMESTR <<"     "<<I.N0<<" particles";
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
    cout<<"Calculation complete in "<<STEP-1<<" steps ("<<I.maxstep<<" allowed)"<<endl;
}

void PrintSimInfo(SimInfo& I, Stats& OUT, bool Fluor) {
    //Fresnel coefficients
    double Rspec = (I.n0-I.n)*(I.n0-I.n)/(I.n0+I.n)/(I.n0+I.n);
    double Rx = (I.nx-I.n)*(I.nx-I.n)/(I.nx+I.n)/(I.nx+I.n); // I.n/I.nx = m -> Rx = (1-m)^2/(1+m)^2 
    double Tspec = 1.0 - Rspec;
    if (Fluor)
        I.WSCALE = 1;
    else
        I.WSCALE = Tspec;
            
    //Warn the user
    cout << "==================================================================" << endl;
    if (Fluor)
        cout << "FLUORESCENCE EMISSION" << endl;
    else
        cout << "INCIDENT BEAM" << endl;
    cout << "==================================================================" << endl;
    
    //Define output arrays
    cout << "Simulation size: " << I.R << " x " << I.L << endl;
    cout << endl << "Z: ";
    for (int i = 0; i < OUT.Zres; i ++)
        cout << scientific << setw(14) << OUT.dZ * i * I.L;
    cout << endl << "R: ";
    for (int i = 0; i < OUT.Rres; i ++)
        cout << scientific << setw(14) << OUT.dR * i * I.R;
    cout << endl << endl;
        
    //Fresnel coefficients
    cout <<endl <<endl;
    cout << "Fresnel coefficients:" << endl;
    cout << "  R0  = " << Rspec << endl;
    cout << "  Rx  = " << Rx << endl;
        
    //Print MFPs
    string TAG;
    if (Fluor)
        TAG = "EMITTED";
    else
        TAG = "INCIDENT";
    cout << "Scattering mean free path (" << TAG << "): " << 1.0/I.Ss/I.dens << "nm" << endl;
    cout << "Absorption mean free path (" << TAG << "): " << 1.0/I.Sa/I.dens << "nm" << endl;
}

void RunIncidentSim(SimInfo& I) {
    
    //Print info
    Stats OUT(I.N0);
    PrintSimInfo(I, OUT, false);
    
    //Initialize photons
    _GEN.seed(time(0));
    vector<Photon> INCIDENT = genBeam(I);
        
    //Run simulation
    RunSim(I, OUT, INCIDENT);
    OUT.print();        
}

void RunFluorSim(SimInfo& I, Stats& DEP) {
    
    //Print info
    Stats OUT(I.N0);
    PrintSimInfo(I, OUT, true);
    
    //Initialize photons
    _GEN.seed(time(0));
    vector<Photon> EMITTED = genFluor(I, DEP);
        
    //Run simulation
    RunSim(I, OUT, EMITTED);
    OUT.print();        
}


///
//  Performs steady state Monte Carlo simulation of (non-interacting) photons
//  including scattering, energy deposition, fluorescence emission, etc.
//     Units - nm, nm2, nm3, ns
//     Use ion density and convert particle scattering cross-section to per ion basis
///
int main(int argc, char** argv) {
    
    //Read file
    SimInfo I;
    if (argc > 1)
        I.load(argv[1]);
     
    //Run simulation
    RunIncidentSim(I);

};
