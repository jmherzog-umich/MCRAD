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
#include <chrono>

#include "stats.h"
#include "vec.h"
#include "photon.h"
#include "utility.h"
#include "medium.h"
#include "beam.h"
#include "grid.h"
#include "image.h"
#include "camera.h"

#ifndef _SIMINFO_H_
#define _SIMINFO_H_

#define CONST_C  299.8
#define CONST_PI 3.1415926535
#define CONST_EPS 1e-10

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
            SinglePhoton  = 0b00100000,     //Whether to operate in single photon mode
            Fluorescence  = 0b01000000,     //Whether to generate fluorescence photons during the run
            Cartesian     = 0b10000000,     //Whether the radial wall is square or not
        };
    
    private:
        //Incident beam
        double n0=1.600;                //Refractive index of front surface (glass)
        double nx=1.600;                //Refractive index of back surface (glass)
        double nr=1.600;                //Refractive index of outer radial surface (glass)
        
        //Sample info
        int N0 = 5e6;                   //Number of initial photon groups
        double L=1e4;                   //Length of sample [um]
        double R=1e4;                   //Radius of sample [um]
        double dT=1000;                 //Time step in simulation (0 for static)
        
        //Sim Flags
        SimFlags flags = (Simulation::SimFlags)3;
        
        //Define medium, beam, and grids
        Medium medium;
        Beam beam;
        Grid grid;
        Stats stats;
        Image imgIC, imgR, imgT;
        Camera cam;
        
        //Convergence and other settings
        double Wmin = 1e-10;            //Start terminating packets if they get this small
        double Wm = 0.1;                //Weighting for Roulette termination procedure
        unsigned int maxstep = 5e6;     //Maximum number of steps before terminating program
        
        //Stats settings
        int Zres = 20;
        int Rres = 20;
        int Tres = 20;
        int THres = 10;
        int momentlvl = 4;
        int printSteps = 0;
        
        //Calculated values
        vector<Photon> PHOTONS;
        double tsim;
};

Simulation::Simulation() {
    //Set default parameters
    n0=1.600; nx=1.600; nr=1.600;
    L=1e4; R=1e4; dT=10; 
    N0 = 1e6; maxstep = 5e6;
    Wmin = 1e-10; Wm = 0.1;
    Zres = 20; Rres = 20; Tres = 20; THres = 10; momentlvl = 4;
    flags = (Simulation::SimFlags)135;
    beam = Beam();
    medium = Medium();
    imgIC = Image(5,5,2*beam.Rb/5, false);
    imgR = Image(10,10,2*beam.Rb/5, ((int)flags & (int)SimFlags::Interference));
    imgT = Image(10,10,2*beam.Rb/5, ((int)flags & (int)SimFlags::Interference));
    cam = Camera(5, 1000, -1, 1.2, 50000, vec(100000,0,5000), vec(-1,0,0));
    printSteps = 0;
    tsim = 0;
}

template<class T>
void remove(vector<T>& v, unsigned int i) {
    if (i != v.size())
        v.at(i) = v.back();
    v.pop_back();
}

void Simulation::print() {
    cout << "==================================================================" << endl;
    cout << "Transmitted and reflected beams" << endl;
    cout << "==================================================================" << endl;
    cout << "Reflected beam intensity profile [photons]" << endl;
    imgR.printGrid();
    imgR.print();
    
    cout << "Transmitted beam intensity profile [photons]" << endl;
    imgT.printGrid();
    imgT.print();

    cout << "==================================================================" << endl;
    cout << "Simulated image " << endl;
    cout << "==================================================================" << endl;
    cam.printGrid();
    cam.print();

    cout << "==================================================================" << endl;
    cout << "Statistics" << endl;
    cout << "==================================================================" << endl;
    stats.print();
}

void Simulation::genBeam() {
    //Alert the user
    cerr << "Generating incident beam photons..." << endl;
    PHOTONS = beam.sampleBeam(N0);
}

void Simulation::genFluor() {
    //Alert the user
    cerr << "FLUORESCENCE NOT IMPLEMENTED!" << endl;
    
}

void Simulation::setup() {
    //Create output stats block
    stats = Stats(Tres, THres, dT, momentlvl);
    stats.setup();
    
    //Setup grid
    if ((int)flags & (int)SimFlags::Cartesian)
        grid = Grid(R, R, L, Rres, 1, Zres);
    else
        grid = Grid(R, L, Rres, 1, Zres);
    grid.newval("Absorption");
    grid.newval("Incident");
    grid.clear();
    
    //Setup images and cameras
    imgIC.clear();
    imgT.clear();
    imgR.clear();
    cam.setup();
    
    //Add to stats IC
    double R;
    for (unsigned int i = 0; i < PHOTONS.size(); i ++) {
        //Store IC
        imgIC.image(PHOTONS.at(i).x.X, PHOTONS.at(i).x.Y, PHOTONS.at(i).W);
        
        //Reflect photons at interface
        R = stats.initialize(PHOTONS.at(i), n0/medium.n(PHOTONS.at(i).v));
        PHOTONS.at(i).W *= (1-R);
        PHOTONS.at(i).S = -log(roll());
    }
    
    //Fresnel coefficients
    double n = medium.n();
    double Rspec = (n0-n)*(n0-n)/(n0+n)/(n0+n);
    double Rx = (nx-n)*(nx-n)/(nx+n)/(nx+n); // I.n/I.nx = m -> Rx = (1-m)^2/(1+m)^2 
    double Tspec = 1.0 - Rspec;
    
    //Print some settings
    cout << "==================================================================" << endl;
    cout << "Settings [UNITS: um, ps, THz]" << endl;
    cout << "==================================================================" << endl;
        
    cout << "Back Wall: " << (((int)flags & (int)SimFlags::BackWall) ? "True" : "False" ) << endl;
    cout << "Front Wall: " << (((int)flags & (int)SimFlags::FrontWall) ? "True" : "False" ) << endl;
    cout << "Radial Wall: " << (((int)flags & (int)SimFlags::RadialWall) ? "True" : "False" ) << endl;
    cout << "Packet interference: " << (((int)flags & (int)SimFlags::Interference) ? "True" : "False" ) << endl;
    cout << "Time-dependent: " << ((dT > 0) ? "True" : "False") << endl;
    cout << "Saturation: " << (((int)flags & (int)SimFlags::Saturation) ? "True" : "False" ) << endl;
    cout << "Single-photon mode: " << (((int)flags & (int)SimFlags::SinglePhoton) ? "True" : "False" ) << endl;
    cout << "Fluorescence generation: " << (((int)flags & (int)SimFlags::Fluorescence) ? "True" : "False" ) << endl;
    cout << "Side wall geometry: " << (((int)flags & (int)SimFlags::Cartesian) ? "Cartesian" : "Cylindrical" ) << endl;
    cout << "Photon packets: " << N0 << endl;
    cout << "Total photons: " << stats.PHI << endl;
    cout << "Front Refractive Index: " << n0 << endl;
    cout << "Back Refractive Index: " << nx << endl;
    cout << "Side Refractive Index: " << nr << endl;
    cout << "Theoretical Max (Static) Step Count: " << ceil(log(Wmin)/log(medium.albedo())) + 1.0/Wm << endl;
    if (dT > 0)
        cout << "Theoretical Max (Dynamic) Step Count: " << ceil((ceil(log(Wmin)/log(medium.albedo())) + 1.0/Wm) / medium.we() / stats.dT) << endl;
    cout << endl;
    
    //Print medium settings
    medium.print();
    
    //Print beam info
    beam.print();
        
    //Print camera info
    cam.printSetup();
        
    //Define output arrays
    grid.printGrid();
        
    //Fresnel coefficients
    cout << "Fresnel coefficients:" << endl;
    cout << "  R0  = " << Rspec << endl;
    cout << "  Rx  = " << Rx << endl;
    cout << "  T0  = " << Tspec << endl << endl;
    
    //Print the initial condition
    cout << "==================================================================" << endl;
    cout << "Initial condition" << endl;
    cout << "==================================================================" << endl;
    cout << "Incident beam profile before specular reflection [photons]:" << endl;
    imgIC.printGrid();
    imgIC.print();
    
    //Warn the user
    cout << "==================================================================" << endl;
    cout << "Starting simulation" << endl;
    cout << "==================================================================" << endl;
    
    //Sort the photons if we're time dependent
    if (dT > 0) {
        cerr << "Sorting photon vector..." << endl;
        sort(begin(PHOTONS), end(PHOTONS));
        tsim = PHOTONS.at(0).t + dT;
    }
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
    double eps, ds, ds2, ka, ka0, ks, k, tempR, Fabs, m, c, sint;
    unsigned long int ii,jj,kk;
    bool xy = false, done = false;
    int reflect;
    vec norm, u;
    Photon pt;
    
    //Scattering/absorption coefficients
    ks = medium.ks();
    ka0 = medium.ka();
    ka = ka0;
    k = ka + ks;
    ii = jj = kk = 0;
    
    //Some default values just in case
    double n = medium.n();
    double Emin = Wmin * beam.E;
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
            if ((dT>0) and (PHOTONS.at(i).t >= tsim))       //If time-dependent, and tsim is earlier than photon time, skip to next iteration (photons are sorted)
                break;                           //
            while (PHOTONS.at(i).W <= Emin) {    //Loop until W!=0
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
            
            //Check for boundary collision and cell collisions
            while (true) {
            
                //Check our current s-value
                if (PHOTONS.at(i).S <= CONST_EPS) {
                    
                    //Calculate grid values
                    ii = grid.ind1(PHOTONS.at(i).x);
                    jj = grid.ind2(PHOTONS.at(i).x);
                    kk = grid.ind3(PHOTONS.at(i).x);
                    
                    //Decide if we're treating this as a packet or single-photon
                    if ((int)flags & (int)SimFlags::SinglePhoton) {
                        //Roll to determine scatter or absorb
                        eps = roll();
                        if (eps < ka/k) { //Absorb
                            grid.at(0, ii,jj,kk) += PHOTONS.at(i).W;
                            stats.absorb(PHOTONS.at(i), PHOTONS.at(i).W);
                            PHOTONS.at(i).W = -1;
                            PHOTONS.at(i).S = 0;
                        } else { //Scatter
                            grid.at(1, ii,jj,kk) += PHOTONS.at(i).W;
                            stats.scatter(PHOTONS.at(i), PHOTONS.at(i).W);
                            PHOTONS.at(i).Scatter(roll(), roll(), medium);
                            PHOTONS.at(i).S = -log(roll());
                        }
                    } else {
                        //Update grids
                        grid.at(0, ii,jj,kk) += PHOTONS.at(i).W * ka/k;
                        grid.at(1, ii,jj,kk) += PHOTONS.at(i).W;
                    
                        //Update stats
                        stats.scatter(PHOTONS.at(i), PHOTONS.at(i).W);
                        stats.absorb(PHOTONS.at(i), PHOTONS.at(i).W * ka / k);
                    
                        //-update photon properties
                        PHOTONS.at(i).W *= ks/k;
                        PHOTONS.at(i).Scatter(roll(), roll(), medium);
                        PHOTONS.at(i).S = -log(roll());
                    }
                }
            
                //Get the local k and calculate the expected displacement
                if ((int)flags & (int)SimFlags::Saturation) {
                    Fabs = grid.norm(0, PHOTONS.at(i).x) / medium.dens();  //E = photons absorbed/volume, dens = molecules/volume, E/dens = photons absorbed/molecule
                    if (Fabs > 0)
                        ka = ka0 / Fabs * (1.0 - exp(-Fabs));
                    else
                        ka = ka0;
                }
                k = ka + ks;
                    
                //Calculate the distance to the next interface (step through
                // possible one and keep the shortest distance)
                //-If no obstruction, we just pass freely the whole distance (reflect = 0)
                ds = PHOTONS.at(i).S; ds2 = PHOTONS.at(i).S; reflect = 0;
                    
                //-Calculate distances to adjacent points (but increment by 
                // EPS to ensure that the new position will resolve to the
                // correct cell.
                if ((int)flags & (int)SimFlags::Saturation) {
                    ds2 = grid.intersect(PHOTONS.at(i).x, PHOTONS.at(i).mu);
                    if (ds2 <= ds and ds2 > 0)
                        ds = ds2;
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
                    //Check cylindrical geometry
                    if ((int)flags & (int)SimFlags::Cartesian)
                        ds2 = PHOTONS.at(i).intersectXY(R,xy) * k;
                    else
                        ds2 = PHOTONS.at(i).intersectR(R) * k;
                    if (ds2 < ds and ds2 >= 0) {
                        ds = ds2;
                        reflect = 3;
                    }
                }
                
                //-Check time boundary
                if (dT > 0) {
                    ds2 = (tsim - PHOTONS.at(i).t) * CONST_C * k;
                    if (ds2 < ds and ds2 >= 0) {
                        ds = ds2;
                        reflect = 0;
                    }
                }
                
                //Move photon as necessary and update the incremental shift
                PHOTONS.at(i).S -= ds;
                PHOTONS.at(i).x += PHOTONS.at(i).mu * ds / k;
                PHOTONS.at(i).t += ds / k / CONST_C * n;
                
                //Look for transmission events (hitting back surface) from the current cell
                if ( reflect == 1 ) {
                    //Set new direction and m value
                    norm = vec(0,0,-1);
                    u = PHOTONS.at(i).mu; u.Z = -u.Z;
                    PHOTONS.at(i).x.Z = L-CONST_EPS;
                    m = n / nx;
                }
                
                //Look for reflection events (front surface interaction) from the current cell
                else if ( reflect == 2 ) {
                    //Set new direction and m value
                    norm = vec(0,0,1);
                    u = PHOTONS.at(i).mu; u.Z = -u.Z;
                    PHOTONS.at(i).x.Z = CONST_EPS;
                    m = n / n0;
                }
                
                //Check for reflections from radial wall
                else if ( reflect == 3) {
                    //If we have a cartesian system
                    if ((int)flags & (int)SimFlags::Cartesian) {
                        if (xy) {
                            norm = vec(copysign(1,-PHOTONS.at(i).mu.X),0,0);
                            u = PHOTONS.at(i).mu; u.X = -u.X;
                            PHOTONS.at(i).x.X = copysign(R-CONST_EPS, PHOTONS.at(i).mu.X);
                        } else {
                            norm = vec(0,copysign(1,-PHOTONS.at(i).mu.Y),0);
                            u = PHOTONS.at(i).mu; u.Y = -u.Y;
                            PHOTONS.at(i).x.Y = copysign(R-CONST_EPS, PHOTONS.at(i).mu.Y);
                        }
                    } else {
                        //Set new direction
                        norm = PHOTONS.at(i).x; norm.Z = 0; norm = -norm / norm.norm();
                        c = norm.dot(PHOTONS.at(i).mu);
                        u = PHOTONS.at(i).mu - norm * c * 2.0;
                        
                        //Make sure we go just below r = R
                        tempR = (R - CONST_EPS) / PHOTONS.at(i).x.r();
                        PHOTONS.at(i).x.X = PHOTONS.at(i).x.X * tempR;
                        PHOTONS.at(i).x.Y = PHOTONS.at(i).x.Y * tempR;
                    }
                    m = n / nr;
                }
                
                if ((PHOTONS.at(i).x.Z < 0) or (PHOTONS.at(i).x.Z > L) or (abs(PHOTONS.at(i).x.X) > R) or (abs(PHOTONS.at(i).x.Y) > R))
                    cerr << "WTF!?!" << endl;
                
                //And implement the reflections
                if (reflect) {
                    //Calculate reflection coefficient and sint
                    tempR = stats.getR(PHOTONS.at(i).mu, norm, m, sint);
                    c = norm.dot(PHOTONS.at(i).mu);
                    
                    //Single-photon mode or photon packet mode
                    if ((int)flags & (int)SimFlags::SinglePhoton) {
                        eps = roll();
                        if (eps < tempR) { //Reflect
                            //If reflecting off higher index, flip phase
                            if (m < 1)
                                PHOTONS.at(i).flipPhase = !PHOTONS.at(i).flipPhase;
                            
                            //Reflect the packet
                            PHOTONS.at(i).mu = u;
                            
                        } else { //Transmit
                            //Increment images
                            if (reflect == 1)
                                imgT.image(PHOTONS.at(i).x.X, PHOTONS.at(i).x.Y, PHOTONS.at(i).W, PHOTONS.at(i).phi());
                            else if (reflect == 2)
                                imgR.image(PHOTONS.at(i).x.X, PHOTONS.at(i).x.Y, PHOTONS.at(i).W, PHOTONS.at(i).phi());
                                
                            //Increment stats - only uses sint for angle, not mu so we don't need to update yet
                            stats.reflect(PHOTONS.at(i), reflect, 0, sint);
                                
                            //Now image transmitted photon - u is the incident normal, P.mu is the reflected normal
                            pt = Photon(PHOTONS.at(i));
                            pt.mu = pt.mu * m - norm*(m*c + sqrt(1.0-m*m*(1-c*c)));
                            cam.image(pt);
                            
                            //Kill photon
                            PHOTONS.at(i).W = -1;
                            PHOTONS.at(i).S = 0;
                        }
                    } else {
                        //Increment images
                        if (reflect == 1)
                            imgT.image(PHOTONS.at(i).x.X, PHOTONS.at(i).x.Y, PHOTONS.at(i).W * (1-tempR), PHOTONS.at(i).phi());
                        else if (reflect == 2)
                            imgR.image(PHOTONS.at(i).x.X, PHOTONS.at(i).x.Y, PHOTONS.at(i).W * (1-tempR), PHOTONS.at(i).phi());
                            
                        //Now image transmitted photon
                        pt = Photon(PHOTONS.at(i));
                        pt.mu = pt.mu * m + norm*(m*c - sqrt(1.0-m*m*(1-c*c)));
                        pt.W = pt.W * (1-tempR);
                        cam.image(pt);
                        
                        //Update stats and photon
                        stats.reflect(PHOTONS.at(i), reflect, tempR, sint);
                        PHOTONS.at(i).W *= tempR;
                        if (m < 1)
                            PHOTONS.at(i).flipPhase = !PHOTONS.at(i).flipPhase;
                            
                        //And reflect the remainder of the photon packet
                        PHOTONS.at(i).mu = u;
                    }
                }
                
                //Break here if we're dynamic and the photon time reaches tsim
                if ((dT > 0) and (PHOTONS.at(i).t >= tsim))
                    break;
                //... or if we're static and we reach a small enough S value  
                if ((dT <= 0) and (PHOTONS.at(i).S <= CONST_EPS)) {
                    PHOTONS.at(i).S = 0;
                    break;
                }
            }
            
            //Terminate old packets
            if (PHOTONS.at(i).W <= Emin) {
                eps = roll();
                if (eps <= Wm)
                    PHOTONS.at(i).W /= Wm;
                else
                    PHOTONS.at(i).W = 0;
            }
        }
        
        //Increment T-step
        tsim += dT;
        
        //PRINT DEBUG INFO:
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
        STEP += 1;
        
        //Check if program done
        if (PHOTONS.size() == 0) {
            cerr<<"Finished calculating energy deposition (in "<<STEP<<" steps)"<<endl;
            break;
        }
        
        //Print the grids if requested
        if ((printSteps) and ((STEP % printSteps) == 0)) {
            cout << "Step: " << ((dT>0)?tsim : (double)STEP) << " " << ((dT>0)?" ps":"") << endl;
            grid.print();
            cout << "--------------------" << endl;
        }
    }
    
    //Print final grid only if we don't ask for steps
    if (!printSteps) {
        cout << "Final state: " << ((dT>0)?tsim : (double)STEP) << " " << ((dT>0)?" ps":"") << endl;
        grid.print();
    }
    
    //Report status
    cout<<"Calculation complete in "<<STEP-1<<" steps ("<<maxstep<<" allowed)"<<endl<<endl;
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
                cmdstr >> medium._Ss;
            else if (!key.compare("Sa"))
                cmdstr >> medium._Sa;
            else if (!key.compare("g"))
                cmdstr >> medium._g;
            else if (!key.compare("g2"))
                cmdstr >> medium._g2;
            else if (!key.compare("n0"))
                cmdstr >> n0;
            else if (!key.compare("n"))
                cmdstr >> medium._n;
            else if (!key.compare("nx"))
                cmdstr >> nx;
            else if (!key.compare("nr"))
                cmdstr >> nr;
            else if (!key.compare("E"))
                cmdstr >> beam.E;
            else if (!key.compare("sin0"))
                cmdstr >> beam.sin0;
            else if (!key.compare("N0"))
                cmdstr >> N0;
            else if (!key.compare("dens"))
                cmdstr >> medium._dens;
            else if (!key.compare("L"))
                cmdstr >> L;
            else if (!key.compare("R"))
                cmdstr >> R;
            else if (!key.compare("dT"))
                cmdstr >> dT;
            else if (!key.compare("Rb"))
                cmdstr >> beam.Rb;
            else if (!key.compare("Pb"))
                cmdstr >> beam.Pb;
            else if (!key.compare("Sb"))
                cmdstr >> beam.Sb;
            else if (!key.compare("Zb"))
                cmdstr >> beam.Zb;
            else if (!key.compare("Tb"))
                cmdstr >> beam.Tb;
            else if (!key.compare("Wm"))
                cmdstr >> Wm;
            else if (!key.compare("Wmin"))
                cmdstr >> Wmin;
            else if (!key.compare("maxstep"))
                cmdstr >> maxstep;
            else if (!key.compare("printsteps"))
                cmdstr >> printSteps;
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
                cmdstr >> medium._FQY;
            //More complicated settings
            else if (!key.compare("phase")) {
                cmdstr >> tmpint;
                medium.phase = (Medium::PhaseFunction)tmpint;
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
            } else if (!key.compare("interference")) {
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
            } else if (!key.compare("singlephoton")) {
                cmdstr >> tmpbool;
                if (tmpbool)
                    tmpflags |= (int)SimFlags::SinglePhoton;
                else
                    tmpflags &= ~((int)SimFlags::SinglePhoton);
                flags = (SimFlags) tmpflags;
            } else if (!key.compare("fluorescence")) {
                cmdstr >> tmpbool;
                if (tmpbool)
                    tmpflags |= (int)SimFlags::Fluorescence;
                else
                    tmpflags &= ~((int)SimFlags::Fluorescence);
                flags = (SimFlags) tmpflags;
            } else if (!key.compare("cartesian")) {
                cmdstr >> tmpbool;
                if (tmpbool)
                    tmpflags |= (int)SimFlags::Cartesian;
                else
                    tmpflags &= ~((int)SimFlags::Cartesian);
                flags = (SimFlags) tmpflags;
            } else if (!key.compare("beamprofile")) {
                cmdstr >> tmpint;
                beam.beamprofile = (Beam::BeamType)tmpint;
            } else if (!key.compare("beamspread")) {
                cmdstr >> tmpint;
                beam.spreadfxn = (Beam::BeamSpread)tmpint;
            } else if (!key.compare("beamwidth")) {
                cmdstr >> tmpint;
                beam.beamdur = (Beam::BeamDuration)tmpint;
            } else if (!key.compare("beamspec")) {
                cmdstr >> tmpint;
                beam.beamspec = (Beam::BeamSpectrum)tmpint;
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
                
            else
                cerr << "Unknown input option: " << cmd << endl;
                
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
