#include <vector>
#include <iostream>

#include "vec.h"
#include "photon.h"

#ifndef _STATS_H_
#define _STATS_H_

struct Stats {
    //Settings
    int Tres, THETAres, momentlvl, Fres;
    double dT, dTf, dTHETA, T, PHI, dF, F0;
    double Rdiffuse, Tdiffuse, Tballistic, Rspec, Ldiffuse;
    
    //Output arrays for energy deposition, scatter, etc.
    std::vector<double> Rtheta, Rstheta, Ttheta;
    
    //Incident spectra and time
    std::vector<double> It, Iw;
    
    //Initialize moments and other outputs
    std::vector<double> MOMENTS, SMOMENTS, THMOMENT, TMOMENT, STMOMENT;
    double RAD2, RAD2F, WRAD, WRADF;
    
    //Time-dependence
    std::vector<double> Tt, Rt, Lt, At;
    
    //Wavelength/frequency-dependence
    std::vector<double> Tw, Rw, Lw, Aw;
    
    //Diffusion coefficients
    double DX, DY, DZ, ND;
    
    //Transmission/reflection scattering order
    double Nr, Nt, N2r, N2t;
    double Na, ta, N2a, t2a;
    
    //Transmission/reflection times
    double tRd, tTd, t2Rd, t2Td;
    
    //Fluorescence properties
    double Ff, Fb, Fl, Ffballistic, Fbballistic, Flballistic;   //Front, Back, Side transmission (and ballistic)
    std::vector<double> Fftheta, Fbtheta;                            //Front and back angle-dependence
    std::vector<double> Fft, Fbt, Flt, Fat, Ft;                      //Front, back, side, reabsorption, emission time-dependence
    std::vector<double> Ffw, Fbw, Flw, Faw, Fw;                      //Front, back, side, reabsorption, emission time-dependence
    std::vector<double> FMOMENTS, FSMOMENTS, FTHMOMENT, FTMOMENT, FSTMOMENT;
    
    //Cache
    double ADEP, ASCAT, FDEP, FSCAT, FGEN;                      //Incident absorption/scatter, fluorescence absorption/scatter/generation
    
    //Functions
    //--Constructors
    Stats();
    Stats(int n, double dT, double Fmin, double Fmax, int LVL=4);
    Stats(int nT, int nTH, int nF, double dT, double dTf, double Fmin, double Fmax, int LVL=4);
    
    //--Accessors
    int getBinT(double t) const;
    int getBinTf(double t) const;
    int getBinTh(double sint) const;
    int getBinF(double F) const;
    
    //--Print results functions
    void setup();
    void print(std::ostream& oout) const;
    void printTime(std::ostream& oout) const;
    void printFreq(std::ostream& oout) const;
    
    void writedb(std::ostream& OF) const;
    void writedbheader(std::ostream& OF) const;
    
    //--Functions to aggregate data
    void scatter(const Photon& p, double w);
    void absorb(const Photon& p, double w);
    void emit(const Photon& p);
    void diffusion(const vec& x, double w, double t);
    double initialize(const Photon& x, double m);
    void reflect(const Photon& p, int surf, double R, double sint);
    double getR(const vec& mu, const vec& n, double m, double& sint) const;
};


#endif
