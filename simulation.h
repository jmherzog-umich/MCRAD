#include <sstream>
#include <algorithm>
#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <chrono>
#include <memory>

#include "stats.h"
#include "vec.h"
#include "spectrum.h"
#include "photon.h"
#include "utility.h"
#include "medium.h"
#include "beam.h"
#include "grid.h"
#include "image.h"
#include "camera.h"
#include "raypath.h"

#include "modules/cylgrid.h"
#include "modules/rectgrid.h"
#include "modules/sphgrid.h"

#ifndef _SIMINFO_H_
#define _SIMINFO_H_

#define CONST_PI            3.1415926535897932
#define CONST_C             299.792458
#define CONST_HBAR          1.054571817679489e-22
#define CONST_HK            47.992430733662212
#define CONST_APERY         1.2020569031595942
#define CONST_WIEN          1.5936242600400401
#define CONST_EPS           1e-10
#define CONST_PLANCKMAX     0.6476102378919149

using namespace std;

class Simulation {
    public:
        Simulation();
        void load(stringstream& settings);
        void setup();
        void run();
        void genBeam();
        void genFluor();
        void readArgs(const vector<string>& args);
        void exec();
        void setid(int n);
        bool set(const string &key, const vector<string>& val);
        
        void emitPhoton(const vec& x, double W, double t0);
        
        void print();
        void printstatus(int STEP);
        void printstatusheader();
        void printsettings();
        
        string getdbline() const;
        void write(string s, string line) const;
        void writedb(ostream& OF) const;
        void writedbheader(ostream& OF) const;
        
        enum struct SimFlags {
            BackWall                = 0b000000000001,     //Whether there is a backwall in the simulation or semi-infinite
            FrontWall               = 0b000000000010,     //Whether there is a front wall or semi-infinite
            RadialWall              = 0b000000000100,     //Whether there is a cylindrical wall or infinite
            Interference            = 0b000000001000,     //Whether we should include interference in the image calculations
            Saturation              = 0b000000010000,     //Whether we should allow saturation of absorption in the simulation
            SinglePhoton            = 0b000000100000,     //Whether to operate in single photon mode
            Fluorescence            = 0b000001000000,     //Whether to generate fluorescence photons during the run
            Cartesian               = 0b000010000000,     //Whether the radial wall is square or not
            TimeResolved            = 0b000100000000,     //Whether the radial wall is square or not
            PeriodicXY              = 0b001000000000,     //Whether the radial boundary condition is periodic
            PeriodicZ               = 0b010000000000,     //Whether the Z-axis boundaries enforce periodic behavior
            FluorescenceTrapping    = 0b100000000000,     //Whether to let fluorescence emission be re-absorbed (and re-emitted)
        };
    
    private:
        //Incident beam
        double n0=1.600;                //Refractive index of front surface (glass)
        double nx=1.600;                //Refractive index of back surface (glass)
        double nr=1.600;                //Refractive index of outer radial surface (glass)
        
        //Sample info
        unsigned long int N0 = 5e6;     //Number of initial photon groups
        double L=1e4;                   //Length of sample [um]
        double R=1e4;                   //Radius of sample [um]
        double dT=1;                    //Time step in simulation
        double dTf = 100;               //Time step for fluorescence
        double Fmin = 400;              //Minimum frequency for the simulation
        double Fmax = 1000;             //Maximum frequency for the simulation
        
        int id;                         //Integer ID for the specific process to keep track of simultaneous runs
        
        //Sim Flags
        SimFlags flags = (Simulation::SimFlags)3;
        
        //Define medium, beam, and grids
        Medium medium;
        Beam beam;
        unique_ptr<Grid> grid;
        Stats stats;
        Image imgIC, imgR, imgT;
        Camera cam;
        
        //Convergence and other settings
        double Wmin = 1e-10;            //Start terminating packets if they get this small
        unsigned int maxstep = 5e6;     //Maximum number of steps before terminating program
        
        //Stats settings
        unsigned int Zres = 20;
        unsigned int Rres = 20;
        unsigned int Tres = 20;
        unsigned int THres = 10;
        unsigned int Fres = 20;
        unsigned int momentlvl = 4;
        unsigned int printSteps = 0;
        unsigned long int trackPhoton = 0;
        unsigned long int storepaths = 50;
        
        //Calculated values
        vector<Photon> PHOTONS;
        multimap<double, Photon> NEWPHOTONS;
        vector<RayPath> PATHS;
        unsigned long int storedFpaths = 0;
        double tsim, tfluor, Emin;
        
        //Output settings
        string dbfile;
        string basename;
        bool outputfile;
        
        ostream oout;
        
        ofstream ofile;
        streambuf* obuffer;
        
        //Other things to keep track of
        chrono::system_clock::time_point t0;
        chrono::system_clock::time_point t1;
        chrono::system_clock::time_point t2;
        time_t STARTTIME;
        string STARTTIMESTR;
};

Simulation::Simulation() : oout(NULL) {
    //Set default parameters
    n0=1.600; nx=1.600; nr=1.600;
    L=1e4; R=1e4; dT=1; dTf = 100;
    N0 = 1e6; maxstep = 5e6;
    Wmin = 1e-10;
    Zres = 20; Rres = 20; Tres = 20; THres = 10; momentlvl = 4;
    flags = (Simulation::SimFlags)135;
    printSteps = 0; storedFpaths = 0;
    tsim = 0; tfluor = 0; storepaths = 50; trackPhoton = 0;
    id = 0; outputfile = false;
}

template<class T>
void remove(vector<T>& v, unsigned int i) {
    if (i != v.size())
        v.at(i) = v.back();
    v.pop_back();
}

void Simulation::setid(int n) {
    id = n;
}

void Simulation::print() {
    //Write output file
    util::writeheader(oout, "Reflected beam intensity profile [photons]");
    imgR.printGrid(oout);
    imgR.print(oout);
    
    util::writeheader(oout, "Transmitted beam intensity profile [photons]");
    imgT.printGrid(oout);
    imgT.print(oout);

    util::writeheader(oout, "Simulated image");
    cam.printGrid(oout);
    cam.print(oout);

    util::writeheader(oout, "Statistics");
    stats.print(oout);
    
    //Store path data
    if (storepaths) {
        string fname = util::makefilename(RayPath::ofbasename, "xyz", id);
        cerr << "Writing " << PATHS.size() << " paths to xyz file: " << fname << endl;
        ofstream FILE(fname);
        for (unsigned long int j = 0; j < PATHS.size(); j ++)
            PATHS.at(j).print(FILE);
        FILE.close();
    }
    
    //Close file if needed
    if (outputfile)
        ofile.close();
    
} 

void Simulation::genBeam() {
    //Alert the user
    cerr << "Generating incident beam photons...";
    PHOTONS = beam.sampleBeam(N0);
    cerr << "done!" << endl;
}

void Simulation::emitPhoton(const vec& x, double W0, double t0) {
    //Calculate direction
    double eps = 2*util::roll() - 1;
    double sint = sin(2*CONST_PI*util::roll());
    Photon pt;
    
    //Create fluorescence photon
    pt.x = x;
    pt.mu = vec(eps * sint, sqrt(1-eps*eps)*sint, sqrt(1-sint*sint));
    pt.t = t0;
    pt.v = medium.emit_v(util::roll());
    pt.W = W0 * medium.FQY();
    pt.S = util::logroll();
    pt.flags = (Photon::PhotonFlags)5;
    
    //Store ray paths of fluorescence data 
    if (storedFpaths < storepaths) {
        pt.storeRayPath(PATHS.at(storepaths + storedFpaths));
        storedFpaths ++;
    }
    
    //Emit the photon in stats
    stats.emit(pt);
    NEWPHOTONS.insert({pt.t, pt});
};

void Simulation::genFluor() {
    //If time-resolved, emit FQY * ABS * (1-exp(-dT/tau)) photons per cell divided into Nf packets (obviously weight > WMIN)
    //-Otherwise, emit all light
    
    //Initialize some values
    int nf = 0;                         //Number of packets we'll generate in this cell
    double np = 0;
    double wp = 0;
    double ff;
    unsigned long int ii;
    vec x;
    
    //If time-resolved, use dT to determine emission, otherwise emit all
    if ((int)flags & (int)SimFlags::TimeResolved)
        ff = exp(-(tsim-tfluor)/medium.tau());
    else
        ff = 0;
    
    //Iterate through grid and generate photons
    for (auto it = grid->begin(); !it.end(); it ++) {
        //Cache index
        ii = it.index();
        
        //Figure out how many particles decay at this time step
        np = grid->at(4, ii) * (1.0 - ff);
        
        //See if this cell is OK
        if (np < Emin)
            continue;
        
        //Subtract weight from grid
        grid->at(4, ii) *= ff;
            
        //Figure out how many packets we CAN generate, and weight per packet
        nf = 1 + (int) (np / beam.E);
        wp = np / (double) nf;

        //Loop through packet count, generate, and add to photon vector
        for (int i = 0; i < nf; i ++) {
            
            //Calculate new position, direction, etc.
            x = grid->rand(ii);
            emitPhoton(x, wp, tsim);
        }
    }
    tfluor = tsim;
    
    //Add to vector 
    std::transform( NEWPHOTONS.begin(), NEWPHOTONS.end(), back_inserter( PHOTONS ), [](auto &kv){ return kv.second;});
    NEWPHOTONS.clear();
    
}

void Simulation::printsettings() {
    //Fresnel coefficients
    double n = medium.n();
    double Rspec = (n0-n)*(n0-n)/(n0+n)/(n0+n);
    double Rx = (nx-n)*(nx-n)/(nx+n)/(nx+n); // I.n/I.nx = m -> Rx = (1-m)^2/(1+m)^2 
    double Tspec = 1.0 - Rspec;
    
    //Print some settings
    util::writeheader(oout, "Settings [UNITS: um, ps, THz]");
    oout << "Back Wall: " << (((int)flags & (int)SimFlags::BackWall) ? "True" : "False" ) << endl;
    oout << "Front Wall: " << (((int)flags & (int)SimFlags::FrontWall) ? "True" : "False" ) << endl;
    oout << "Radial Wall: " << (((int)flags & (int)SimFlags::RadialWall) ? "True" : "False" ) << endl;
    oout << "Packet interference: " << (((int)flags & (int)SimFlags::Interference) ? "True" : "False" ) << endl;
    oout << "Time-dependent: " << (((int)flags & (int)SimFlags::TimeResolved) ? "True" : "False") << endl;
    oout << "Saturation: " << (((int)flags & (int)SimFlags::Saturation) ? "True" : "False" ) << endl;
    oout << "Single-photon mode: " << (((int)flags & (int)SimFlags::SinglePhoton) ? "True" : "False" ) << endl;
    oout << "Fluorescence generation: " << (((int)flags & (int)SimFlags::Fluorescence) ? "True" : "False" ) << endl;
    if ((int)flags & (int)SimFlags::Fluorescence)
        oout << "Fluorescence-trapping: " << (((int)flags & (int)SimFlags::FluorescenceTrapping) ? "True" : "False") << endl;
    oout << "Side wall geometry: " << (((int)flags & (int)SimFlags::Cartesian) ? "Cartesian" : "Cylindrical" ) << endl;
    oout << "Periodic XY: " << (((int)flags & (int)SimFlags::PeriodicXY) ? "True" : "False" ) << endl;
    if (((int)flags & (int)SimFlags::PeriodicXY) && !((int)flags & (int)SimFlags::Cartesian))
        oout << "   WARNING: PERIODIC RADIAL BC WITH CYLINDRICAL COORDINATES. TAKE CARE INTERPRETING RESULTS." << endl;
    oout << "Periodic Z: " << (((int)flags & (int)SimFlags::PeriodicZ) ? "True" : "False" ) << endl;
    oout << "Photon packets: " << N0 << endl;
    oout << "Total photons: " << stats.PHI << endl;
    oout << "Frequency range: " << Fmin << " - " << Fmax << " THz" << endl;
    oout << "Simulation time-step: " << dT << " ps" << endl;
    if (((int)flags & (int)SimFlags::TimeResolved) && ((int)flags & (int)SimFlags::Fluorescence))
        oout << "Fluorescence time-step: " << dTf << " ps" << endl;
    oout << "Front Refractive Index: " << n0 << endl;
    oout << "Back Refractive Index: " << nx << endl;
    oout << "Side Refractive Index: " << nr << endl;
    if (storepaths) {
        oout << "Storing " << storepaths << " ray paths to file: " << RayPath::ofbasename << ".xyz" << endl;
        if ((int)flags & (int)SimFlags::Fluorescence)
            oout << "Storing " << storepaths << " fluorescence ray paths to file: " << RayPath::ofbasename << ".xyz" << endl;
    }
    oout << "Theoretical max (static) step count for absorption: " << ceil(log(Wmin)/log(medium.albedo(beam.beamspec.peak())) + 1.0/Photon::Wm) << endl;
    if ((int)flags & (int)SimFlags::TimeResolved)
        oout << "Theoretical max (dynamic) step count for absorption: " << ceil((ceil(log(Wmin)/log(medium.albedo())) + 1.0/Photon::Wm) / medium.we() / stats.dT) << endl;
    oout << endl;
    
    //Fresnel coefficients
    oout << "Fresnel coefficients:" << endl;
    oout << "  R0  = " << Rspec << endl;
    oout << "  Rx  = " << Rx << endl;
    oout << "  T0  = " << Tspec << endl << endl;
    
    //Print medium settings
    util::writeheader(oout, "Medium properties");
    medium.print(oout);
    oout << "Peak incident frequency: " << beam.beamspec.peak() << " THz" << endl;
    medium.print_at_f(oout, beam.beamspec.peak());
    oout << "Peak fluorescence frequency: " << medium.peak_v() << " THz" << endl;
    medium.print_at_f(oout, medium.peak_v());
    
    //Print beam info
    util::writeheader(oout, "Source properties");
    beam.print(oout);
        
    //Print camera info
    util::writeheader(oout, "Camera properties");
    cam.printSetup(oout);
        
    //Define output arrays
    util::writeheader(oout, "Output grid settings");
    grid->printGrid(oout);
    
    //Print the initial condition
    util::writeheader(oout, "Initial condition");
    oout << "Incident beam profile before specular reflection [photons]:" << endl;
    imgIC.printGrid(oout);
    imgIC.print(oout);
    
    //Warn the user
    util::writeheader(oout, "Starting simulation");
}

void Simulation::setup() {
    //Create output stats block
    stats = Stats(Tres, THres, Fres, dT, dTf, Fmin, Fmax, momentlvl);
    stats.setup();
    medium.Fmin = Fmin;
    medium.Fmax = Fmax;
    Spectrum::setFreqLimits(Fmin, Fmax);
    
    //Create output stream
    if (outputfile) {
        ofile.open(util::makefilename(basename, "out", id));
        obuffer = ofile.rdbuf();
    } else
        obuffer = cout.rdbuf();
    oout.rdbuf(obuffer);
    
    //Create images and cameras
    imgIC = Image(10,10,2*beam.Rb/5, false);
    imgR = Image(10,10,2*beam.Rb/5, ((int)flags & (int)SimFlags::Interference));
    imgT = Image(10,10,2*beam.Rb/5, ((int)flags & (int)SimFlags::Interference));
    cam = Camera(10, L/50, -0.2, 1.2, 50000, vec(300000,0,L/2), vec(-1,0,0));   //M=-0.2
    
    //Setup grid
    if ((int)flags & (int)SimFlags::Cartesian)
        grid = make_unique<RectGrid>(2*R, 2*R, L, Rres, 1, Zres);
    else
        grid = make_unique<CylGrid>(R, L, Rres, 1, Zres);
    grid->newval("Absorption");
    grid->newval("Incident");
    if ((int)flags & (int)SimFlags::Fluorescence) {
        grid->newval("Reabsorption");
        grid->newval("IncidentEmission");
        grid->newval("ExcitedState");
    }
    grid->clear();
    
    //Setup images and cameras
    imgIC.clear();
    imgT.clear();
    imgR.clear();
    cam.setup();
    
    //Initialize PATHS
    PATHS = vector<RayPath>(storepaths * (((int)flags&(int)SimFlags::Fluorescence) ? 2 : 1));
    
    //Add to stats IC
    double R;
    for (unsigned int i = 0; i < PHOTONS.size(); i ++) {
        //Store IC
        imgIC.image(PHOTONS.at(i).x.X, PHOTONS.at(i).x.Y, PHOTONS.at(i).W);
        
        //Reflect photons at interface
        R = stats.initialize(PHOTONS.at(i), n0/medium.n(PHOTONS.at(i).v));
        PHOTONS.at(i).W *= (1-R);
        PHOTONS.at(i).S = util::logroll();
        
        //Generate raypath objects if needed
        if (i < storepaths)
            PHOTONS.at(i).storeRayPath(PATHS.at(i));
    }
    
    //Print the settings
    printsettings();
    
    //Sort the photons if we're time dependent
    if ((int)flags & (int)SimFlags::TimeResolved) {
        cerr << "Sorting photon vector...";
        sort(begin(PHOTONS), end(PHOTONS));
        tsim = PHOTONS.at(0).t + dT;
        tfluor = tsim;
        cerr << " Done!" << endl;
    }
    
}

void Simulation::printstatusheader() {
    cerr << "   Start time: " << STARTTIMESTR << "     " << N0 << " particles";
    cerr << "      Max steps: " << maxstep << endl;
    cerr << "-------------------------------------------------------------------------------" << endl;
    cerr << " Step [-]  sim t [ps]  total/step time [s]      particles       grid [photons] " << endl;
    cerr << "-------------------------------------------------------------------------------" << endl;
}

void Simulation::printstatus(int STEP) {
    //Step and sim time
    chrono::duration<double> ELAPSED;
    cerr << setw(10) << STEP;
    cerr << " ";
    if ((int)flags & (int)SimFlags::TimeResolved) {
        cerr << setw(10) << fixed << setprecision(3) << tsim;
        cerr << " ";
    } else
        cerr << "           ";
    
    //Real time
    ELAPSED = t2 - t0;
    cerr << setw(10) << fixed << setprecision(1) << ELAPSED.count();
    ELAPSED = t2 - t1;
    cerr << "/" << setw(10) << fixed << setprecision(3) << ELAPSED.count();
    cerr << " ";
    
    //Photons and energy
    cerr << setw(9) << PHOTONS.size();
    if ((int)flags & (int)SimFlags::Fluorescence)
        cerr << "+" << setw(9) << NEWPHOTONS.size();
    else
        cerr << "          ";
    if ((int)flags & (int)SimFlags::Fluorescence)
        cerr << " " << setw(10) << scientific << setprecision(3) << grid->sum(4);
    else
        cerr << " " << setw(10) << scientific << setprecision(3) << grid->sum(0);
    
    //And finally new line
    cerr << endl;
}

void Simulation::run() {
    //Setup clocks
    t0 = t1 = t2 = chrono::system_clock::now();
    STARTTIME = chrono::system_clock::to_time_t(t0);
    STARTTIMESTR = ctime(&STARTTIME);
    STARTTIMESTR = STARTTIMESTR.erase(STARTTIMESTR.back());
    
    //Initialize a whole lot of temporary variables
    double eps, ds, ds2, ka, ka0, ks, k, tempR, Fabs, m, c, sint;
    unsigned long int ii,jj,kk;
    bool done = false, END=false;
    int reflect;
    vec norm, u;
    Photon pt;
    ii = jj = kk = 0;
    
    //Some default values just in case
    double n = medium.n();
    Emin = Wmin * beam.E;
    m = 1;
    
    //loop through steps
    unsigned int STEP = 0; 
    while (STEP <= maxstep) {
        
        //Track the current photon
        if (trackPhoton)
            PHOTONS.at(trackPhoton).printstatus();
        
        //Prepare for loop
        done = false;
        
        //Loop through each photon and update (go forward to delete if needed)
        for (unsigned long int i = 0; true; i ++) {   //Don't check vector size cause it will change
        
            //Check that it's a good vector
            if (i >= PHOTONS.size())             //Exit if we've covered the whole vector
                break;
            if (((int)flags & (int)SimFlags::TimeResolved) and (PHOTONS.at(i).t >= tsim))       //If time-dependent, and tsim is earlier than photon time, skip to next iteration (photons are sorted)
                break;                           //
            while (PHOTONS.at(i).W <= Emin) {    //Loop until W!=0
                if (i < PHOTONS.size()-1)        //If not end, delete and continue
                    remove(PHOTONS, i);          //Swap with back and remove if we have enough elements to not screw up the path
                else {                           //Otherwise
                    if (i == 0)
                        PHOTONS.clear();         // If i==0 and we can't remove any more, then we're at the end
                    done = true;                 // propagate "we're done" downstream
                    break;                       // stop looking for photons
                }
            }
            if (done)
                break;                           // and stop simulation
            
            //Get scattering/abs coeffs for this photon
            ks = medium.ks(PHOTONS.at(i).v);
            ka0 = medium.ka(PHOTONS.at(i).v);
            ka = ka0;
            k = ka + ks;
            
            //Check for boundary collision and cell collisions
            while (true) {
            
                //Check our current s-value
                if (PHOTONS.at(i).S <= CONST_EPS) {
                    
                    //Calculate grid values
                    ii = grid->ind1(PHOTONS.at(i).x);
                    jj = grid->ind2(PHOTONS.at(i).x);
                    kk = grid->ind3(PHOTONS.at(i).x);
                    
                    // TODO: After changes below moving stats into grid, we can combine the grid->at and absorb calls into one.
                    // 
                    
                    //Decide if we're treating this as a packet or single-photon
                    if ((int)flags & (int)SimFlags::SinglePhoton) {
                        //Roll to determine scatter or absorb
                        eps = util::roll();
                        if (eps < ka/k) { //Absorb
                            //Calculate stats
                            grid->at((PHOTONS.at(i).isFluorescence()?2:0), ii,jj,kk) += PHOTONS.at(i).W;
                            stats.absorb(PHOTONS.at(i), PHOTONS.at(i).W);
                            
                            //Generate fluorescence photons
                            if ((int)flags & (int)SimFlags::Fluorescence)
                                if (((int)flags & (int)SimFlags::FluorescenceTrapping) || (!PHOTONS.at(i).isFluorescence())){
                                    grid->at(4, ii,jj,kk) += PHOTONS.at(i).W;
                                    emitPhoton(PHOTONS.at(i).x, PHOTONS.at(i).W, tsim + medium.emit_tau(util::logroll()));
                                }
                            
                            //And kill the photon
                            PHOTONS.at(i).Kill();
                            break;
                            
                        } else { //Scatter
                            grid->at((PHOTONS.at(i).isFluorescence()?3:1), ii,jj,kk) += PHOTONS.at(i).W;
                            stats.scatter(PHOTONS.at(i), PHOTONS.at(i).W);
                            PHOTONS.at(i).Scatter(util::roll(), util::roll(), medium, util::logroll());
                        }
                    } else {
                        //Update grids
                        grid->at((PHOTONS.at(i).isFluorescence()?2:0), ii,jj,kk) += PHOTONS.at(i).W * ka/k;
                        if ((int)flags & (int)SimFlags::Fluorescence)
                            if (((int)flags & (int)SimFlags::FluorescenceTrapping) || (!PHOTONS.at(i).isFluorescence()))
                                grid->at(4, ii,jj,kk) += PHOTONS.at(i).W * ka/k;
                        grid->at((PHOTONS.at(i).isFluorescence()?3:1), ii,jj,kk) += PHOTONS.at(i).W;
                    
                        //Update stats
                        stats.scatter(PHOTONS.at(i), PHOTONS.at(i).W);
                        stats.absorb(PHOTONS.at(i), PHOTONS.at(i).W * ka / k);
                    
                        //-update photon properties
                        PHOTONS.at(i).W *= ks/k;
                        PHOTONS.at(i).Scatter(util::roll(), util::roll(), medium, util::logroll());
                    }
                }
            
                //Get the local k and calculate the expected displacement
                if ((int)flags & (int)SimFlags::Saturation) {
                    //E = photons absorbed/volume, dens = molecules/volume, E/dens = photons absorbed/molecule
                    if (PHOTONS.at(i).isFluorescence())
                        Fabs = grid->dens(2, PHOTONS.at(i).x) / medium.dens();
                    else if ((int)flags & (int)SimFlags::Fluorescence)
                        Fabs = grid->dens(4, PHOTONS.at(i).x) / medium.dens();
                    else
                        Fabs = grid->dens(0, PHOTONS.at(i).x) / medium.dens();
                    if (Fabs > 0)
                        ka = ka0 / Fabs * (1.0 - exp(-Fabs));
                    else
                        ka = ka0;
                }
                k = ka + ks;
                    
                //Calculate the distance to the next interface (step through
                // possible one and keep the shortest distance)
                //-If no obstruction, we just pass freely the whole distance (reflect = 0)
                ds = PHOTONS.at(i).S/k;
                ds2 = ds;
                reflect = 0;
                
                //-Check adjacent cells if we need ot
                if ((int)flags & (int)SimFlags::Saturation)
                  grid->collideCell(PHOTONS.at(i).x, PHOTONS.at(i).mu, ds, reflect);
                  
                //-Check front/back walls
                if (!((int)flags & (int)SimFlags::PeriodicZ)) {
                    if ((int)flags & (int)SimFlags::FrontWall)
                        grid->collideFront(PHOTONS.at(i).x, PHOTONS.at(i).mu, ds, reflect);
                    if ((int)flags & (int)SimFlags::BackWall)
                        grid->collideBack(PHOTONS.at(i).x, PHOTONS.at(i).mu, ds, reflect);
                }
                
                //-Check sides
                if (!((int)flags & (int)SimFlags::PeriodicXY) && ((int)flags & (int)SimFlags::RadialWall))
                    grid->collideSide(PHOTONS.at(i).x, PHOTONS.at(i).mu, ds, reflect);
                
                //-Check time boundary
                if ((int)flags & (int)SimFlags::TimeResolved) {
                    ds2 = (tsim - PHOTONS.at(i).t) * CONST_C;
                    if (ds2 < ds and ds2 >= 0) {
                        ds = ds2;
                        reflect = 0;
                    }
                }
                
                //Move photon as necessary and update the incremental shift
                PHOTONS.at(i).S -= ds * k;
                PHOTONS.at(i).x += PHOTONS.at(i).mu * ds * k / medium.ke(PHOTONS.at(i).v);
                PHOTONS.at(i).t += ds / CONST_C * n;
                
                //Enforce periodic BC's
                if ((int)flags & (int)SimFlags::PeriodicXY)
                    grid->boundXY(PHOTONS.at(i).x);
                if ((int)flags & (int)SimFlags::PeriodicZ)
                    grid->boundZ(PHOTONS.at(i).x);
                
                //Get normal vector, m, and reflected direction
                ///NOTE: May need to clip photons to just inside (e.g., L-CONST_EPS) the boundary to make sure nothing weird happens
                norm = grid->getNormal(reflect, PHOTONS.at(i).x, PHOTONS.at(i).mu, u);
                if (reflect == 1)
                    m = n/nx;
                else if (reflect == 2)
                    m = n/n0;
                else if (reflect > 3)
                    m = n/nr;
                
                //TODO: Offload the following onto grid, and have grid own stats.
                // Instead of imgT and imgR, grid should hold two additional 2D arrays
                // to store reflect/transmit images (for sphere, e.g., they should map to 
                // the hemispherical front or back surface).
                //
                // Then we call reflect or transmit function with "reflect" argument to
                // determine side, along with given normal, u, and m, and it should 
                // add the contribution to imgT/imgR (now members of grid) and update the
                // actual stats at the same time.
                //
                // We can also offload the transmitted photon generation/propagation to 
                // the camera class so we just call, e.g., cam->transmitAndImage(PHOTON.at(i), norm, m, 1-tempR);
                // where tempR = 0 for single photon mode
                
                //And implement the reflections
                if (reflect) {
                
                    //Calculate reflection coefficient and sint
                    tempR = stats.getR(PHOTONS.at(i).mu, norm, m, sint);
                    c = norm.dot(PHOTONS.at(i).mu);
                    
                    //Single-photon mode or photon packet mode
                    if ((int)flags & (int)SimFlags::SinglePhoton) {
                    
                        //Figure out if we reflect or transmit
                        eps = util::roll();
                        if (eps < tempR) { //Reflect
                            //If reflecting off higher index, flip phase
                            if (m < 1)
                                PHOTONS.at(i).flipPhase();
                            
                            //Reflect the packet
                            PHOTONS.at(i).mu = u;
                            PHOTONS.at(i).Reflect();
                            
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
                            break;
                        }
                    } else {
                        //Increment images
                        if (reflect == 1)
                            imgT.image(PHOTONS.at(i).x.X, PHOTONS.at(i).x.Y, PHOTONS.at(i).W * (1-tempR), PHOTONS.at(i).phi());
                        else if (reflect == 2)
                            imgR.image(PHOTONS.at(i).x.X, PHOTONS.at(i).x.Y, PHOTONS.at(i).W * (1-tempR), PHOTONS.at(i).phi());
                            
                        //Now image transmitted photon
                        if (PHOTONS.at(i).W > Emin) {
                            pt = Photon(PHOTONS.at(i));
                            pt.mu = pt.mu * m + norm*(m*c - sqrt(1.0-m*m*(1-c*c)));
                            pt.W = pt.W * (1-tempR);
                            cam.image(pt);
                        }
                        
                        //Update stats and photon
                        stats.reflect(PHOTONS.at(i), reflect, tempR, sint);
                        PHOTONS.at(i).W *= tempR;
                        if (m < 1)
                            PHOTONS.at(i).flipPhase();
                            
                        //And reflect the remainder of the photon packet
                        PHOTONS.at(i).mu = u;
                        PHOTONS.at(i).Reflect();
                    }
                }
                
                //Break here if we're dynamic and the photon time reaches tsim
                if (((int)flags & (int)SimFlags::TimeResolved) and (PHOTONS.at(i).t >= tsim))
                    break;
                    
                //... or if we're static and we reach a small enough S value  
                if (!((int)flags & (int)SimFlags::TimeResolved) and (PHOTONS.at(i).S <= CONST_EPS)) {
                    PHOTONS.at(i).S = 0;
                    break;
                }
                
                //Check if we should terminate this packet
                //PHOTONS.at(i).roulette();
                if ((PHOTONS.at(i).W <= Emin) && (PHOTONS.at(i).W > -CONST_EPS)
                         && !((int)flags & (int)SimFlags::SinglePhoton))
                    PHOTONS.at(i).roulette();
            }
        }
        
        //Add new photon packets as necessary
        if ((int)flags & (int)SimFlags::Fluorescence) {
            
            //-Time-resolved version 
            if ((int)flags & (int)SimFlags::TimeResolved) {
                //Branch based on single- or multi-photon packets
                //For single photons, copy into simulation when we reach the emission time
                if ((int)flags & (int)SimFlags::SinglePhoton) {
                    if (NEWPHOTONS.size() > 0) {
                        auto upper = NEWPHOTONS.lower_bound(tsim+dT);
                        auto lower = NEWPHOTONS.upper_bound(0);
                        std::transform( lower, upper, back_inserter( PHOTONS ), [](auto &kv){ return kv.second;});
                        NEWPHOTONS.erase(lower, upper);
                    }
                }
                        
                //Otherwise, generate packets based on fraction of light emitted within this time step
                else if (tsim - tfluor >= dTf)
                    genFluor();
                
            //-Static version: immediately add photons in single-photon mode, add when out of primary photons otherwise
            } else {
                if ((int)flags & (int)SimFlags::SinglePhoton) {
                    std::transform( NEWPHOTONS.begin(), NEWPHOTONS.end(), back_inserter( PHOTONS ), [](auto &kv){ return kv.second;});
                    NEWPHOTONS.clear();
                } else if (PHOTONS.size() == 0)
                    genFluor();
            }
        }
        
        //Print status
        t2 = chrono::system_clock::now();
        if (STEP == 0)
            printstatusheader();
        printstatus(STEP);
        t1 = t2;
        
        //Increment T-step
        STEP += 1;
        
        //Check END conditions
        if (PHOTONS.size() == 0) {
            
            //No fluorescence so just end it
            if ( !((int)flags & (int)SimFlags::Fluorescence) )
                END = true;
            
            //If time-resolved, check if we have to wait until later to generate photons and increment, or exit
            else if ((int)flags & (int)SimFlags::TimeResolved) {
                
                //Time-resolved and single-photon, so check if we have NEWPHOTONS
                if ( (int)flags & (int)SimFlags::SinglePhoton ) {
                    
                    //Exit if no new photons here
                    if (NEWPHOTONS.size() == 0)
                        END = true;
                    //Otherwise incrememnt to the next photon emission time
                    else
                        tsim = NEWPHOTONS.begin()->first;
                }
                //Time-resolved and packet mode, so check if we have particles left in the excited state
                else {
                    
                    //Exit if grid is empty
                    if (grid->less(4, Emin))
                        END = true;
                    //Otherwise increment by the fluorescence time
                    else
                        tsim += dTf;
                }
            }
            //If not time-resolved, we just continue until NEWPHOTONS is also empty
            else if (NEWPHOTONS.size() == 0)
                END = true;
            else
                tsim += dT;
        }
        //Otherwise, we just stay the course and increment time. NOTE: don't need to increment if not time-resolved so it's OK that one non-time-resolved case in branch above doesn't increment dT if not end
        else
            tsim += dT;
        
        //End if we set the condition
        if (END) {
            cerr<<"Finished calculating energy deposition (in "<<STEP<<" steps)"<<endl;
            break;
        }
        
        //Print the grids if requested
        if ((printSteps) and ((STEP % printSteps) == 0)) {
            oout << "Step: ";
            if ((int)flags & (int)SimFlags::TimeResolved)
                oout << tsim << "  ps" << endl;
            else
                oout << STEP << endl;
            grid->print(oout);
            oout << "--------------------" << endl;
        }
    }
    
    //Print final grid only if we don't ask for steps
    if (!printSteps) {
        oout << "Final state: ";
        if ((int)flags & (int)SimFlags::TimeResolved)
            oout << tsim << "  ps" << endl;
        else
            oout << STEP << endl;
        grid->print(oout, true);
    }
    
    //Report status
    oout<<"Calculation complete in "<<STEP-1<<" steps ("<<maxstep<<" allowed)"<<endl<<endl;
}

bool Simulation::set(const string &key, const vector<string>& val) {
    //Temporary variables
    bool tmpbool = false;
    unsigned long int tmpflags = (unsigned long int)flags;

    //Run through the giant list of parameters to check - return true if we set one, otherwise false
    
    //Test Medium subclass first
    tmpbool = medium.set(key, val);
    if (tmpbool)
        return true;
        
    //And the beam subclass
    tmpbool = beam.set(key, val);
    if (tmpbool)
        return true;
    
    //If we're here, the key either belongs to Simulation or doesn't exist
    if (!key.compare("index-front"))
        n0 = stod(val.at(0));
    else if (!key.compare("index-back"))
        nx = stod(val.at(0));
    else if (!key.compare("index-side"))
        nr = stod(val.at(0));
    else if (!key.compare("simulation-packets"))
        N0 = stoul(val.at(0));
    else if (!key.compare("simulation-length"))
        L = stod(val.at(0));
    else if (!key.compare("simulation-radius"))
        R = stod(val.at(0));
    else if (!key.compare("simulation-timestep"))
        dT = stod(val.at(0));
    else if (!key.compare("simulation-frequency-range")) {
        Fmin = stod(val.at(0));
        Fmax = stod(val.at(1));
    } else if (!key.compare("simulation-timestep-fluorescence"))
        dTf = stod(val.at(0));
    else if (!key.compare("roulette-newweight"))
        Photon::Wm = stod(val.at(0));
    else if (!key.compare("roulette-minweight"))
        Wmin = stod(val.at(0));
    else if (!key.compare("max-steps"))
        maxstep = stoul(val.at(0));
    else if (!key.compare("database-filename"))
        dbfile = val.at(0);
    else if (!key.compare("base-filename")) {
        basename = val.at(0);
        outputfile = true;
    }
    else if (!key.compare("print-steps"))
        printSteps = stoul(val.at(0));
    else if (!key.compare("track-photon"))
        trackPhoton = stoul(val.at(0));
    else if (!key.compare("export-paths")) {
        storepaths = stoul(val.at(0));
        RayPath::ofbasename = val.at(1);
    }
    else if (!key.compare("grid-z-points"))
        Zres = stoul(val.at(0));
    else if (!key.compare("grid-xy-points"))
        Rres = stoul(val.at(0));
    else if (!key.compare("grid-t-points"))
        Tres = stoul(val.at(0));
    else if (!key.compare("grid-theta-points"))
        THres = stoul(val.at(0));
    else if (!key.compare("grid-freq-points"))
        Fres = stoul(val.at(0));
    else if (!key.compare("moments"))
        momentlvl = stoul(val.at(0));
    else if (!key.compare("enable-backwall")) {
        tmpbool = stoi(val.at(0));
        if (tmpbool)
            tmpflags |= (int)SimFlags::BackWall;
        else
            tmpflags &= ~((int)SimFlags::BackWall);
        flags = (SimFlags) tmpflags;
    } else if (!key.compare("enable-frontwall")) {
        tmpbool = stoi(val.at(0));
        if (tmpbool)
            tmpflags |= (int)SimFlags::FrontWall;
        else
            tmpflags &= ~((int)SimFlags::FrontWall);
        flags = (SimFlags) tmpflags;
    } else if (!key.compare("enable-sidewall")) {
        tmpbool = stoi(val.at(0));
        if (tmpbool)
            tmpflags |= (int)SimFlags::RadialWall;
        else
            tmpflags &= ~((int)SimFlags::RadialWall);
        flags = (SimFlags) tmpflags;
    } else if (!key.compare("enable-interference")) {
        tmpbool = stoi(val.at(0));
        if (tmpbool)
            tmpflags |= (int)SimFlags::Interference;
        else
            tmpflags &= ~((int)SimFlags::Interference);
        flags = (SimFlags) tmpflags;
    } else if (!key.compare("enable-saturation")) {
        tmpbool = stoi(val.at(0));
        if (tmpbool)
            tmpflags |= (int)SimFlags::Saturation;
        else
            tmpflags &= ~((int)SimFlags::Saturation);
        flags = (SimFlags) tmpflags;
    } else if (!key.compare("enable-single-photon")) {
        tmpbool = stoi(val.at(0));
        if (tmpbool)
            tmpflags |= (int)SimFlags::SinglePhoton;
        else
            tmpflags &= ~((int)SimFlags::SinglePhoton);
        flags = (SimFlags) tmpflags;
    } else if (!key.compare("enable-time-dependent")) {
        tmpbool = stoi(val.at(0));
        if (tmpbool)
            tmpflags |= (int)SimFlags::TimeResolved;
        else
            tmpflags &= ~((int)SimFlags::TimeResolved);
        flags = (SimFlags) tmpflags;
    } else if (!key.compare("enable-fluorescence")) {
        tmpbool = stoi(val.at(0));
        if (tmpbool)
            tmpflags |= (int)SimFlags::Fluorescence;
        else
            tmpflags &= ~((int)SimFlags::Fluorescence);
        flags = (SimFlags) tmpflags;
    } else if (!key.compare("enable-fluorescence-trapping")) {
        tmpbool = stoi(val.at(0));
        if (tmpbool)
            tmpflags |= (int)SimFlags::FluorescenceTrapping;
        else
            tmpflags &= ~((int)SimFlags::FluorescenceTrapping);
        flags = (SimFlags) tmpflags;
    } else if (!key.compare("enable-cartesian")) {
        tmpbool = stoi(val.at(0));
        if (tmpbool)
            tmpflags |= (int)SimFlags::Cartesian;
        else
            tmpflags &= ~((int)SimFlags::Cartesian);
        flags = (SimFlags) tmpflags;
    } else if (!key.compare("enable-periodic-xy")) {
        tmpbool = stoi(val.at(0));
        if (tmpbool)
            tmpflags |= (int)SimFlags::PeriodicXY;
        else
            tmpflags &= ~((int)SimFlags::PeriodicXY);
        flags = (SimFlags) tmpflags;
    } else if (!key.compare("enable-periodic-z")) {
        tmpbool = stoi(val.at(0));
        if (tmpbool)
            tmpflags |= (int)SimFlags::PeriodicZ;
        else
            tmpflags &= ~((int)SimFlags::PeriodicZ);
        flags = (SimFlags) tmpflags;
    } else
        return false;
    return true;
}

void Simulation::load(stringstream& settings) {
        //Allocate some memory and open the file
        string cmd, key, tmp;
        vector<string> vals;
        stringstream cmdstr;
        bool done;
        
        //Process line by line
        while (getline(settings,cmd)) {
            
            //Skip empty strings
            if (cmd.empty())
                continue;
                
            //Make into a stringstream and split into key/value pairs
            cmdstr.str(cmd);
            cmdstr >> key;
            while (cmdstr >> tmp)
                vals.push_back(tmp);
                
            //Check values to set
            done = set(key, vals);
            
            //Commands to run the program
            if (!done)
                cerr << "Unknown input option: " << cmd << endl;
                
            //Clear values
            key.clear();
            cmd.clear();
            vals.clear();
            cmdstr.clear();
        }
        return;
}

void Simulation::writedbheader(ostream& OF) const {
    OF << "N0,L,R,dT,dTf,flags,n0,nx,nr";
}

void Simulation::writedb(ostream& OF) const {
    OF << setprecision(8) << N0 << ",";
    OF << setprecision(8) << L << ",";
    OF << setprecision(8) << R << ",";
    OF << setprecision(8) << dT << ",";
    OF << setprecision(8) << dTf << ",";
    OF << setprecision(8) << util::bitstring((unsigned short)flags) << ",";
    OF << setprecision(8) << n0 << ",";
    OF << setprecision(8) << nx << ",";
    OF << setprecision(8) << nr;
}


string Simulation::getdbline() const {
    //Create a string stream to store the text in
    stringstream line;
    
    //Generate it chunk by chunk
    writedb(line);
    line << ",";
    beam.writedb(line);
    line << ",";
    medium.writedb(line);
    line << ",";
    stats.writedb(line);
    line << endl;
    
    return line.str();
}

void Simulation::write(string s, string line) const {
    //Check if string is empty
    if (s.empty()) {
        //Don't use ID in filename since we want to write all to the same file. Write lock issues?
        s = util::makefilename(dbfile, "csv", -1);
        if (s.empty())
            return;
    }
    
    //Open file in append mode
    ofstream OF(s, ios_base::app);
    
    //If file is empty, write a header
    if (OF.tellp() == 0) {
        writedbheader(OF);
        OF << ",";
        beam.writedbheader(OF);
        OF << ",";
        medium.writedbheader(OF);
        OF << ",";
        stats.writedbheader(OF);
        OF << endl;
    }
    
    //Now write the row
    OF << line;
    
    //Close file
    OF.close();
}

void Simulation::readArgs(const vector<string>& args) {
    //Loop through args and parse them
    string cmd;
    vector<string> arg;
    unsigned long int i = 0;
    while(i < args.size()) {
        //Get command
        cmd = args.at(i);
        i ++;
        
        //Get value if we didn't get the end of the string
        if (i < args.size()) {
            arg = util::splitcomma(args.at(i));
            i ++;
        } else
            arg = {""};
            
        //Now interpret them
        set(cmd, arg);
    }
    return;
}

void Simulation::exec() {
    //Now generate beam, setup, etc., and run the program 
    genBeam();
    setup();
    run();
    print();
}

#endif
