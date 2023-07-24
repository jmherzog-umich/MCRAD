#include <sstream>
#include <algorithm>
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
#define CONST_EPS 1e-6

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
        
        enum struct SimFlags {
            BackWall      = 0b00000001,     //Whether there is a backwall in the simulation or semi-infinite
            FrontWall     = 0b00000010,     //Whether there is a front wall or semi-infinite
            RadialWall    = 0b00000100,     //Whether there is a cylindrical wall or infinite
            Interference  = 0b00001000,     //Whether we should include interference in the image calculations
            Saturation    = 0b00010000,     //Whether we should allow saturation of absorption in the simulation
        };
        
        //Definition for beam geometries and random spreading functions
        enum struct BeamSpread { Collimated = 0, Uniform = 1, Gaussian = 2, Lambertian = 3, Isotropic = 4 };
        enum struct BeamType { TopHat = 0, Gaussian = 1, TopHatEllipse = 2, GaussianEllipse = 3, TopHatAnnulus = 4, GaussianAnnulus = 5 };
        enum struct BeamDuration { TopHat = 0, Gaussian = 1, Cauchy = 2 };
    
    private:
        //Incident beam
        double Ss=5e5;                  //Scattering cross-section
        double Sa=5e5;                  //absorption cross-section
        double g=0.98;                  //Scattering anisotropy
        double w=5.5e5;                 //Angular frequency in GHz (5.5e5 ~ 545 nm light)
        double E = 1;                   //Total number of photons (not packets, but quanta of light)
        double n0=1.600;                //Refractive index of front surface (glass)
        double n=1.330;                 //Refractive index of medium (water)
        double nx=1.600;                //Refractive index of back surface (glass)
        double nr=1.600;                //Refractive index of outer radial surface (glass)
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
        BeamType beamprofile = Simulation::BeamType::TopHat;
        BeamSpread spreadfxn = Simulation::BeamSpread::Collimated;
        BeamDuration beamdur = Simulation::BeamDuration::TopHat;
        
        //Beam parameters
        double Rb=5e5;          //Beam radius [nm]
        double Pb=0;            //Beam spread parameter: standard deviation (gauss)/width (uniform)
        double Sb=2.5e5;        //Beam profile parameter: linewidth (TopHatEllipse, Guass1D), Inner radius (TopHatAnnulus, GaussianAnnulus)
        double Zb = 0;          //Focus location of beam (0 for infinity/unfocused)
        double Tb = 0;          //Laser pulse duration
        
        //Convergence and other settings
        double Wmin = 1e-10;    //Start terminating packets if they get this small
        double Wm = 0.1;        //Weighting for Roulette termination procedure
        unsigned int maxstep = 5e6;      //Maximum number of steps before terminating program
        
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
    N0 = 5e6; dens=1e-13; L=1e7; R=1e7; T=300; nr=1.600;
    Rb=5e5; E = 1e10; w=5.5e5; Zb = 0; Tb = 0;
    Wmin = 1e-10; Wm = 0.1; maxstep = 5e6; 
    Zres = 100; Rres = 100; Tres = 1; THres = 50; momentlvl = 4;
    flags = (Simulation::SimFlags)7;
    spreadfxn = Simulation::BeamSpread::Collimated;
    beamprofile = Simulation::BeamType::TopHat;
    beamdur = Simulation::BeamDuration::TopHat;
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
    cerr << "Generating incident beam photons" << endl;
    
    //Initialize some values
    PHOTONS = vector<Photon>(N0, Photon(vec(),vec(sin0, 0, sqrt(1-sin0*sin0)), g));
    if ((int)flags & (int)SimFlags::FrontWall)
        WSCALE = (1.0 - (n0-n)*(n0-n)/(n0+n)/(n0+n));
    else
        WSCALE = 1.0;
    double eps, eps2;
    const vec NORM(R, R, L);
    
    //Loop through and generate beam if needed
    double x,y;
    double Xf = Zb*tan(asin(sin0));
    for (int i = 0; i < N0; i ++) {
        
        ///
        //  Calculate beam profile
        ///
        if (Rb > 1) {        
            //Random numbers
            eps2 = roll();
            eps = roll();
        
            //Use the randos to generate x,y values based on model
            switch (beamprofile) {
                case Simulation::BeamType::TopHat :
                    x = Rb * sqrt(eps) * cos(2.0 * CONST_PI * eps2);
                    y = Rb * sqrt(eps) * sin(2.0 * CONST_PI * eps2);
                    break;
                    
                case Simulation::BeamType::Gaussian :
                    x = Rb * erfinvf((float)eps);
                    y = x * sin(2.0 * CONST_PI * eps2);
                    x = x * cos(2.0 * CONST_PI * eps2);
                    break;
                    
                case Simulation::BeamType::TopHatEllipse :
                    eps2 = atan(Sb/Rb*tan(2*CONST_PI*eps2));
                    eps = Sb*Rb/sqrt(pow(Sb*cos(eps2),2)+pow(Rb*sin(eps2),2)) * sqrt(eps);
                    x = eps*cos(eps2);
                    y = eps*sin(eps2);
                    break;
                    
                case Simulation::BeamType::GaussianEllipse :
                    x = Rb * erfinvf((float)eps);
                    y = Sb * erfinvf((float)eps2);
                    break;
                    
                case Simulation::BeamType::TopHatAnnulus :
                    x = sqrt(Sb*Sb + (Rb*Rb-Sb*Sb)*eps) * cos(2.0 * CONST_PI * eps2);
                    y = sqrt(Sb*Sb + (Rb*Rb-Sb*Sb)*eps) * sin(2.0 * CONST_PI * eps2);
                    break;
                    
                case Simulation::BeamType::GaussianAnnulus :
                    x = Rb + (Sb - Rb) * erfinvf((float)eps);
                    y = x * sin(2.0 * CONST_PI * eps2);
                    x = x * cos(2.0 * CONST_PI * eps2);
                    break;
                    
                default:
                    x = 0;
                    y = 0;
                    break;
            }
            
            //Set the positions
            PHOTONS.at(i).x.X += x;
            PHOTONS.at(i).x.Y += y;
        }
        
        ///
        //  Calculate photon packet time
        ///
        if (Tb > 0) {
            eps = roll();
            switch (beamdur) {
                case Simulation::BeamDuration::TopHat :
                    PHOTONS.at(i).t = eps * Tb;
                    break;
                case Simulation::BeamDuration::Gaussian :
                    PHOTONS.at(i).t = erfinvf((float)eps);
                    break;
                    
                case Simulation::BeamDuration::Cauchy :
                    PHOTONS.at(i).t = -Tb / tan(CONST_PI * eps);
                    break;
            }
        }
        
        ///
        //  Update directions
        ///
        //Focus the beam
        if (Zb > 0) {
            PHOTONS.at(i).mu = vec(Xf, 0, Zb) - PHOTONS.at(i).x;
            PHOTONS.at(i).mu /= PHOTONS.at(i).mu.norm();
        }
    
        //Implement diffraction here?
        ///
        // TODO: DIFFRACTION CALCULATION
        ///
    
        //Exit now if we aren't further spreading the beam
        if (spreadfxn == Simulation::BeamSpread::Collimated)
            continue;
            
        //Roll new random numbers
        eps2 = roll()*2*CONST_PI; //Azimuthal angle
        eps = roll();
        
        //Sample distributions for dispersion on top of focus
        switch (spreadfxn) {
            case Simulation::BeamSpread::Collimated :
                break;
            case Simulation::BeamSpread::Uniform :
                break;
            case Simulation::BeamSpread::Gaussian :
                break;
            case Simulation::BeamSpread::Lambertian :
                break;
            case Simulation::BeamSpread::Isotropic :
                break;
        }
    }   
}

void Simulation::genFluor() {
    //Alert the user
    cerr << "Generating fluorescence photons" << endl;
    
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
    OUT.saturate = (int)flags&(int)SimFlags::Saturation;
    OUT.interfere = (int)flags&(int)SimFlags::Interference;
    
    //Add to stats IC
    const vec NORM(R, R, L);
    for (unsigned int i = 0; i < PHOTONS.size(); i ++)
        OUT.initialize(PHOTONS.at(i), NORM);
    
    //Setup random number generator
    _DIST = uniform_real_distribution<double>(0.0,1.0);
    _GEN.seed(time(0));
    WSCALE = 1;
    
    //Fresnel coefficients
    double Rspec = (n0-n)*(n0-n)/(n0+n)/(n0+n);
    double Rx = (nx-n)*(nx-n)/(nx+n)/(nx+n); // I.n/I.nx = m -> Rx = (1-m)^2/(1+m)^2 
    double Tspec = 1.0 - Rspec;
    
    //Print some settings
    cout << "==================================================================" << endl;
    cout << "Settings" << endl;
    cout << "==================================================================" << endl;
    cout << "Phase function: ";
    switch (Photon::phase) {
        case Photon::PhaseFunction::HenyeyGreenstein: cout << "Henyey-Greenstein" << endl; break;
        case Photon::PhaseFunction::Rayleigh: cout << "Rayleigh" << endl; break;
        default: cout << "Other" << endl; break;
    }
    cout << "Beam profile: ";
    switch (beamprofile) {
        case Simulation::BeamType::TopHat : cout << "Top Hat (r = " << Rb << ")" << endl;  break;
        case Simulation::BeamType::Gaussian : cout << "Gaussian (sigma = " << Rb << ")" << endl; break;
        case Simulation::BeamType::TopHatEllipse : cout << "Elliptical Top Hat (a = " << Rb << ", b = " << Sb << ")" << endl; break;
        case Simulation::BeamType::GaussianEllipse : cout << "Elliptical Gaussian (sigmax = " << Rb << ", sigmay = " << Sb << ")" << endl; break;
        case Simulation::BeamType::TopHatAnnulus : cout << "Annular Top Hat (ri = " << Sb << ", ro = " << Rb << ")" << endl; break;
        case Simulation::BeamType::GaussianAnnulus : cout << "Annular Gaussian (mu = " << Rb << ", sigma = " << Sb << ")" << endl; break;
        default : cout << "Other" << endl; break;
    }
    
    cout << "Beam spread function: ";
    switch (spreadfxn) {
        case Simulation::BeamSpread::Collimated : cout << "Collimated" << endl; break;
        case Simulation::BeamSpread::Uniform : cout << "Uniform [0," << Pb << ")" << endl; break;
        case Simulation::BeamSpread::Gaussian : cout << "Gaussian (sigma = " << Pb << ")" << endl; break;
        case Simulation::BeamSpread::Lambertian : cout << "Lambertian [0, " << Pb << ")" << endl; break;
        case Simulation::BeamSpread::Isotropic : cout << "Isotropic [0, " << Pb << ")" << endl; break;
        default: cout << "Other" << endl; break;
    }
    
    cout << "Beam temporal profile: ";
    switch (beamdur) {
        case Simulation::BeamDuration::TopHat : cout << "TopHat (" << Tb << " ns)" << endl; break;
        case Simulation::BeamDuration::Gaussian : cout << "Gaussian (sigma = " << Tb << ")" << endl; break;
        case Simulation::BeamDuration::Cauchy : cout << "Gaussian (gamma = " << Tb << ")" << endl; break;
        default: cout << "Other" << endl; break;
    }
    
    cout << "g: " << g << endl;
    cout << "Focus at depth: " << ((Zb > 0) ? Zb : INFINITY) << endl;
    cout << "Back Wall: " << (((int)flags & (int)SimFlags::BackWall) ? "True" : "False" ) << endl;
    cout << "Front Wall: " << (((int)flags & (int)SimFlags::FrontWall) ? "True" : "False" ) << endl;
    cout << "Radial Wall: " << (((int)flags & (int)SimFlags::RadialWall) ? "True" : "False" ) << endl;
    cout << "Packet interference: " << (((int)flags & (int)SimFlags::Interference) ? "True" : "False" ) << endl;
    cout << "Time-dependent: " << ((T > 0) ? true : false) << endl;
    cout << "Saturation: " << (((int)flags & (int)SimFlags::Saturation) ? "True" : "False" ) << endl;
    cout << "Photon packets: " << N0 << endl;
    cout << "Scattering cross-section: " << Ss << endl;
    cout << "Absorption cross-section: " << Sa << endl;
    cout << "Front Refractive Index: " << n0 << endl;
    cout << "Medium Refractive Index: " << n << endl;
    cout << "Back Refractive Index: " << nx << endl;
    cout << "Side Refractive Index: " << nr << endl;
    cout << "Incident sin(theta): " << sin0 << endl;
    cout << "Albedo: " << Ss/(Ss+Sa) << endl;
    cout << "Theoretical Max Step Count: " << ceil(log(Wmin)/log(Ss/(Ss+Sa)) + 1.0/Wm) << endl;
    cout << endl << endl;
    
    //Define output arrays
    cout << "Simulation Domain (R,Z,T): " << R << ", " << L << ", " << T << endl;
    cout << "Grid points (nR,nZ,nT): " << Rres << ", " << Zres << ", " << Tres << endl;
    cout << endl << "Z: ";
    for (int i = 0; i < OUT.Zres; i ++)
        cout << scientific << setw(14) << OUT.dZ * i * L;
    cout << endl << "R: ";
    for (int i = 0; i < OUT.Rres; i ++)
        cout << scientific << setw(14) << OUT.dR * i * R;
    cout << endl << endl;
        
    //Fresnel coefficients
    cout << "Fresnel coefficients:" << endl;
    cout << "  R0  = " << Rspec << endl;
    cout << "  Rx  = " << Rx << endl;
    cout << "  T0  = " << Tspec << endl;
    
    //Print MFPs
    cout << "Scattering mean free path: " << 1.0/Ss/dens << "nm" << endl;
    cout << "Absorption mean free path: " << 1.0/Sa/dens << "nm" << endl;
    
    //Warn the user
    cout << "==================================================================" << endl;
    cout << "Starting simulation" << endl;
    cout << "==================================================================" << endl;
    
    //Sort the photons if we're time dependent
    if (T > 0)
        sort(begin(PHOTONS), end(PHOTONS));
}

void Simulation::run() {
    
    //Setup clocks
    auto t0 = chrono::system_clock::now();
    auto t1 = chrono::system_clock::now();
    auto t2 = chrono::system_clock::now();
    double tsim = 0;
    std::chrono::duration<double> ELAPSED;
    time_t STARTTIME = chrono::system_clock::to_time_t(t0);
    string STARTTIMESTR(ctime(&STARTTIME));
    STARTTIMESTR = STARTTIMESTR.erase(STARTTIMESTR.back());
    
    //Initialize a whole lot of temporary variables
    double eps, s, ds, ds2, ka, ka0, ks, k, tempR, Fabs, xf, m;
    bool done = false;
    int reflect;
    const vec NORM(R, R, L);
    vec x, u;
    
    //Scattering/absorption coefficients
    ks = Ss*dens;
    ka0 = Sa*dens;
    ka = ka0;
    k = ka + ks;
    
    //Some default values just in case
    m = 1;
    
    //loop through steps
    unsigned int STEP = 0; 
    while (STEP <= maxstep) {
        
        //Prepare for loop
        done = false;
        
        //Loop through each photon and update (go forward to delete if needed)
        for (unsigned int i = 0; true; i ++) {   //Don't check vector size cause it will change
        
            //Check that it's a good vector
            if (i >= PHOTONS.size())             //Exit if we've covered the whole vector
                break;
            while (PHOTONS.at(i).W <= Wmin) {    //Loop until W!=0
                if (i < PHOTONS.size()-1) {      //If not end...
                    remove(PHOTONS, i);          // delete and continue
                } else {                         //Otherwise
                    if (i == 0)
                        PHOTONS.clear();         // If i==0 and we can't remove any more, then we're at the end
                    done = true;                 // propagate "we're done" downstream
                    break;                       // stop looking for photons
                }
            }
            if (done)
                break;                           // and stop simulation
        
            //Roll the dice
            eps = roll();
            s = -log(eps);
            
            //Check for boundary collision and cell collisions
            while (s > CONST_EPS) {
                //Get the local k and calculate the expected displacement
                if ((int)flags & (int)SimFlags::Saturation) {
                    Fabs = OUT.E(PHOTONS.at(i).x / NORM) / dens *E/N0;
                    if (Fabs > 0)
                        ka = ka0 / Fabs * (1.0 - exp(-Fabs));
                    else
                        ka = ka0;
                }
                k = ka + ks;
                    
                //Calculate the distance to the next interface (step through
                // possible one and keep the shortest distance)
                //-If no obstruction, we just pass freely the whole distance (reflect = 0)
                ds = s; ds2 = s; reflect = 0;
                    
                //-Calculate distances to adjacent points (but increment by 
                // EPS to ensure that the new position will resolve to the
                // correct cell.
                if ((int)flags & (int)SimFlags::Saturation) {
                    
                    //Check next z-shell if abs(mu.Z > 0)
                    xf = PHOTONS.at(i).x.Z / (L/Zres);
                    if ( PHOTONS.at(i).mu.Z > 0 )
                        ds2 = ((ceil(xf)-xf) * (L/Zres) / PHOTONS.at(i).mu.Z) * k;
                    else if (PHOTONS.at(i).mu.Z < 0)
                        ds2 = (-(xf-floor(xf)) * (L/Zres) / PHOTONS.at(i).mu.Z) * k;
                    if (ds2 < ds and ds2 > 0)
                        ds = ds2;
 
                    //Check radial shell
                    xf = ceil(PHOTONS.at(i).x.r() / (R/Rres));
                    ds2 = PHOTONS.at(i).intersectR(xf) * k;
                    if (ds2 < ds and ds2 > 0)
                        ds = ds2;
                    if (xf > 0) {
                        ds2 = PHOTONS.at(i).intersectR(xf-R/Rres) * k;
                        if (ds2 < ds and ds2 > 0)
                            ds = ds2;
                    }
                }
                
                //-Calculate distance to back wall
                if ( ((int)flags & (int)SimFlags::BackWall) and (PHOTONS.at(i).mu.Z > 0)) {
                    //If the back wall is closer than the current ds value
                    ds2 = (L - PHOTONS.at(i).x.Z)/PHOTONS.at(i).mu.Z * k;
                    if (ds2 <= ds and ds2 >= 0) {
                        ds = ds2;
                        reflect = 1;
                    }
                }
                
                //-Calculate distance to front wall
                if ( ((int)flags & (int)SimFlags::FrontWall) and (PHOTONS.at(i).mu.Z < 0)) {
                    //If the front wall is closer than the current ds value
                    ds2 = -PHOTONS.at(i).x.Z/PHOTONS.at(i).mu.Z * k;
                    if (ds2 <= ds and ds2 >= 0) {
                        ds = ds2;
                        reflect = 2;
                    }
                }
                
                //-Calculate distance to radial wall
                if ((int)flags & (int)SimFlags::RadialWall) {
                    ds2 = PHOTONS.at(i).intersectR(R) * k;
                    if (ds2 <= ds and ds2 >= 0) {
                        ds = ds2;
                        reflect = 3;
                    }
                }
                
                //Move photon as necessary and update the incremental shift
                s -= ds;
                PHOTONS.at(i).x += PHOTONS.at(i).mu * ds / k;
                PHOTONS.at(i).t += ds / k / CONST_C * n;
                
                //Look for transmission events (hitting back surface) from the current cell
                if ( reflect == 1 ) {
                    //Set new direction and m value
                    PHOTONS.at(i).mu.Z = -PHOTONS.at(i).mu.Z;
                    PHOTONS.at(i).x.Z = L-CONST_EPS;
                    m = n / nx;
                }
                
                //Look for reflection events (front surface interaction) from the current cell
                else if ( reflect == 2 ) {
                    //Set new direction and m value
                    PHOTONS.at(i).mu.Z = -PHOTONS.at(i).mu.Z;
                    PHOTONS.at(i).x.Z = CONST_EPS;
                    m = n / n0;
                }
                
                //Check for reflections from radial wall
                else if ( reflect == 3) {
                    //Set new direction
                    u = PHOTONS.at(i).x; u.Z = 0; u = -u / u.norm();
                    PHOTONS.at(i).mu = PHOTONS.at(i).mu - u * u.dot(PHOTONS.at(i).mu) * 2;
                    
                    //Get new m value
                    m = n / nr;
                    
                    //Make sure we go just below r = R
                    tempR = R / PHOTONS.at(i).x.r();
                    PHOTONS.at(i).x.X *= tempR * (1 - CONST_EPS);
                    PHOTONS.at(i).x.Y *= tempR * (1 - CONST_EPS);
                }
                
                //And implement the reflections
                if (reflect) {
                    tempR = OUT.reflect(PHOTONS.at(i), NORM, T, WSCALE, m, reflect, STEP);
                    PHOTONS.at(i).W = PHOTONS.at(i).W * tempR;
                    if (m < 1)
                        PHOTONS.at(i).flipPhase = !PHOTONS.at(i).flipPhase;
                }
            }
            
            //Calculate new photon info and store scattering event
            OUT.scatter(PHOTONS.at(i), NORM, T, WSCALE, ka/k);
            PHOTONS.at(i).W -= PHOTONS.at(i).W * ka/k;
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
    cout<<"Calculation complete in "<<STEP-1<<" steps ("<<maxstep<<" allowed)"<<endl<<endl;
    OUT.finalize();
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
            else if (!key.compare("nr"))
                cmdstr >> nr;
            else if (!key.compare("E"))
                cmdstr >> E;
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
            else if (!key.compare("Pb"))
                cmdstr >> Pb;
            else if (!key.compare("Sb"))
                cmdstr >> Sb;
            else if (!key.compare("Zb"))
                cmdstr >> Zb;
            else if (!key.compare("Tb"))
                cmdstr >> Tb;
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
                    tmpflags |= (int)SimFlags::BackWall;
                else
                    tmpflags &= ~((int)SimFlags::BackWall);
                flags = (SimFlags) tmpflags;
            } else if (!key.compare("frontwall")) {
                cmdstr >> tmpbool;
                if (tmpbool)
                    tmpflags |= (int)SimFlags::FrontWall;
                else
                    tmpflags &= ~((int)SimFlags::FrontWall);
                flags = (SimFlags) tmpflags;
            } else if (!key.compare("cylwall")) {
                cmdstr >> tmpbool;
                if (tmpbool)
                    tmpflags |= (int)SimFlags::RadialWall;
                else
                    tmpflags &= ~((int)SimFlags::RadialWall);
                flags = (SimFlags) tmpflags;
            } else if (!key.compare("beamdiffraction")) {
                cmdstr >> tmpbool;
                if (tmpbool)
                    tmpflags |= (int)SimFlags::Interference;
                else
                    tmpflags &= ~((int)SimFlags::Interference);
                flags = (SimFlags) tmpflags;
            } else if (!key.compare("saturation")) {
                cmdstr >> tmpbool;
                if (tmpbool)
                    tmpflags |= (int)SimFlags::Saturation;
                else
                    tmpflags &= ~((int)SimFlags::Saturation);
                flags = (SimFlags) tmpflags;
            } else if (!key.compare("beamprofile")) {
                cmdstr >> tmpint;
                beamprofile = (Simulation::BeamType)tmpint;
            } else if (!key.compare("beamspread")) {
                cmdstr >> tmpint;
                spreadfxn = (Simulation::BeamSpread)tmpint;
            } else if (!key.compare("beamwidth")) {
                cmdstr >> tmpint;
                beamdur = (Simulation::BeamDuration)tmpint;
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
