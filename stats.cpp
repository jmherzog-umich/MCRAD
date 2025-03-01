#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>

#include "vec.h"
#include "photon.h"
#include "stats.h"

#define CONST_PI            3.1415926535897932
#define CONST_C             299.792458
#define CONST_HBAR          1.054571817679489e-22
#define CONST_HK            47.992430733662212
#define CONST_APERY         1.2020569031595942
#define CONST_WIEN          1.5936242600400401
#define CONST_EPS           1e-10
#define CONST_PLANCKMAX     0.6476102378919149

Stats::Stats() {
    Stats(20, 1.0, 0, 1000, 4);
}

Stats::Stats(int n, double dT, double Fmin, double Fmax, int LVL) {
    Stats(n, n, n, dT, dT, Fmin, Fmax, LVL);
}

Stats::Stats(int nT, int nTH, int nF, double dT, double dTf, double Fmin, double Fmax, int LVL) {
    //Set default resolution
    Tres = nT;             //Number of time samples
    Fres = nF;             //Number of Frequency samples
    THETAres = nTH;        //Number of cos(theta) values for reflection
    momentlvl = LVL;       //Order of moments to calculate (only up to 4th order is implemented)
    this->dT = dT;
    this->dTf = dTf;
    this->dF = (Fmax-Fmin) / nF;
    F0 = Fmin;
    PHI = 0;               //Total number of photons
    
    //Reflection and transmission
    Rdiffuse = 0;
    Tdiffuse = 0;
    Ldiffuse = 0;
    Tballistic = 0;
    Rspec = 0;
    PHI = 0;
    
    //Diffusion coefficients
    DX = DY = DZ = 0;
    ND = 0;
    
    //Scattering order and time stats
    Nr = Nt = N2r = N2t = 0;
    tRd = tTd = t2Rd = t2Td = 0;
    Na = ta = N2a = t2a = 0;
    
    //Cache values
    ADEP = 0;
    ASCAT = 0;

    //Fluorescence values
    Ff = Fb = Fl = Ffballistic = Fbballistic = Flballistic = 0;
    FDEP = FGEN = FSCAT = 0;
    
    //Grid resolution
    T = dT * Tres;
    dTHETA = CONST_PI/THETAres/2;
    
    //Initialize moments and other outputs
    MOMENTS = std::vector<double>(14, 0.0);
    THMOMENT = std::vector<double>(8, 0.0);
    TMOMENT = std::vector<double>(4, 0.0);
    SMOMENTS = std::vector<double>(14, 0.0);
    STMOMENT = std::vector<double>(4, 0.0);
    
    //Fluorescence moments
    FMOMENTS = std::vector<double>(14, 0.0);
    FTHMOMENT = std::vector<double>(8, 0.0);
    FTMOMENT = std::vector<double>(4, 0.0);
    FSMOMENTS = std::vector<double>(14, 0.0);
    FSTMOMENT = std::vector<double>(4, 0.0);
}

void Stats::setup() {
    //Grid resolution - note; recalculate here just in case
    T = dT * Tres;
    dTHETA = CONST_PI/THETAres/2;
    
    //Time arrays
    Tt = std::vector<double>(Tres, 0.0);
    Rt = std::vector<double>(Tres, 0.0);
    Lt = std::vector<double>(Tres, 0.0);
    At = std::vector<double>(Tres, 0.0);
    It = std::vector<double>(Tres, 0.0);
    
    //Frequency arrays
    Tw = std::vector<double>(Fres, 0.0);
    Rw = std::vector<double>(Fres, 0.0);
    Lw = std::vector<double>(Fres, 0.0);
    Aw = std::vector<double>(Fres, 0.0);
    Iw = std::vector<double>(Fres, 0.0);
    
    //Reflection and transmission
    Rtheta = std::vector<double>(THETAres, 0.0);
    Rstheta = std::vector<double>(THETAres, 0.0);
    Ttheta = std::vector<double>(THETAres, 0.0);
    
    //Fluorescence values
    Fftheta = std::vector<double>(THETAres, 0.0);
    Fbtheta = std::vector<double>(THETAres, 0.0);
    
    Fft = std::vector<double>(Tres, 0.0);
    Fbt = std::vector<double>(Tres, 0.0);
    Flt = std::vector<double>(Tres, 0.0);
    Fat = std::vector<double>(Tres, 0.0);
    Ft = std::vector<double>(Tres, 0.0);
    
    Ffw = std::vector<double>(Fres, 0.0);
    Fbw = std::vector<double>(Fres, 0.0);
    Flw = std::vector<double>(Fres, 0.0);
    Faw = std::vector<double>(Fres, 0.0);
    Fw = std::vector<double>(Fres, 0.0);
    
    //Initialize moments and other outputs
    MOMENTS = std::vector<double>(14, 0.0);
    THMOMENT = std::vector<double>(8, 0.0);
    TMOMENT = std::vector<double>(4, 0.0);
    SMOMENTS = std::vector<double>(14, 0.0);
    STMOMENT = std::vector<double>(4, 0.0);
    
    //Fluorescence moments
    FMOMENTS = std::vector<double>(14, 0.0);
    FTHMOMENT = std::vector<double>(8, 0.0);
    FTMOMENT = std::vector<double>(4, 0.0);
    FSMOMENTS = std::vector<double>(14, 0.0);
    FSTMOMENT = std::vector<double>(4, 0.0);
}

int Stats::getBinT(double t) const {
    int binT = (int) std::floor(t/dT);
    if ((binT >= Tres) or (binT < 0))
        return -1;
    return binT;
}

int Stats::getBinTf(double t) const {
    int binT = (int) std::floor((t-dTf)/dTf);
    if ((binT >= Tres) or (binT < 0))
        return -1;
    return binT;
}

int Stats::getBinTh(double sint) const {
    int binTheta = (int) std::floor(sint/dTHETA);
    if (binTheta >= THETAres)
        binTheta = THETAres-1;
    else if (binTheta < 0)
        return -1;
    return binTheta;
}

int Stats::getBinF(double w) const {
    int binW = (int) std::floor((w-F0)/dF);
    if ((binW >= Fres) or (binW < 0))
        return -1;
    return binW;
}

void Stats::emit(const Photon& p) {
    FGEN += p.W;
    int binT = getBinTf(p.t);
    if (binT >= 0)
        Ft.at(binT) += p.W;
    int binW = getBinF(p.v);
    if (binW >= 0)
        Fw.at(binW) += p.W;
}

void Stats::absorb(const Photon& p, double W) {
    //Initialize
    double t = p.t;
    const vec& x = p.x;
    int binT;
    
    //Calculate moments
    if (p.isFluorescence()) {
        //Temporal/frequency absorption data
        FDEP += W;
        binT = getBinTf(t);
        if (binT >= 0)
            Fat.at(binT) += W;
        binT = getBinF(p.v);
        if (binT >= 0)
            Faw.at(binT) += W;
        
        //Absorption moments for fluorescence
        if (momentlvl > 0) {
            FMOMENTS.at(0) += W * x.Z;
            FMOMENTS.at(1) += W * x.r();
            FTMOMENT.at(0) += W * t;
        } if (momentlvl > 1) {
            FMOMENTS.at(2) += W * std::pow(x.Z,2);
            FMOMENTS.at(3) += W * x.r() * x.Z;
            FMOMENTS.at(4) += W * std::pow(x.r(),2);
            FTMOMENT.at(1) += W * std::pow(t,2);
        } if (momentlvl > 2) {
            FMOMENTS.at(5) += W * std::pow(x.Z,3);
            FMOMENTS.at(6) += W * x.r() * std::pow(x.Z,2);
            FMOMENTS.at(7) += W * x.Z * std::pow(x.r(),2);
            FMOMENTS.at(8) += W * std::pow(x.r(),3);
            FTMOMENT.at(2) += W * std::pow(t,3);
        } if (momentlvl > 3) {
            FMOMENTS.at(9) += W * std::pow(x.Z,4);
            FMOMENTS.at(10) += W * x.r() * std::pow(x.Z,3);
            FMOMENTS.at(11) += W * std::pow(x.Z * x.r(),2);
            FMOMENTS.at(12) += W * x.Z * std::pow(x.r(),3);
            FMOMENTS.at(13) += W * std::pow(x.r(),4);
            FTMOMENT.at(3) += W * std::pow(t,4);
        }
        
        //And we're done
        return;
    }
    
    //Normal photons
    binT = getBinT(t);
    if (momentlvl > 0) {
        MOMENTS.at(0) += W * x.Z;
        MOMENTS.at(1) += W * x.r();
        TMOMENT.at(0) += W * t;
    } if (momentlvl > 1) {
        MOMENTS.at(2) += W * std::pow(x.Z,2);
        MOMENTS.at(3) += W * x.r() * x.Z;
        MOMENTS.at(4) += W * std::pow(x.r(),2);
        TMOMENT.at(1) += W * std::pow(t,2);
    } if (momentlvl > 2) {
        MOMENTS.at(5) += W * std::pow(x.Z,3);
        MOMENTS.at(6) += W * x.r() * std::pow(x.Z,2);
        MOMENTS.at(7) += W * x.Z * std::pow(x.r(),2);
        MOMENTS.at(8) += W * std::pow(x.r(),3);
        TMOMENT.at(2) += W * std::pow(t,3);
    } if (momentlvl > 3) {
        MOMENTS.at(9) += W * std::pow(x.Z,4);
        MOMENTS.at(10) += W * x.r() * std::pow(x.Z,3);
        MOMENTS.at(11) += W * std::pow(x.Z * x.r(),2);
        MOMENTS.at(12) += W * x.Z * std::pow(x.r(),3);
        MOMENTS.at(13) += W * std::pow(x.r(),4);
        TMOMENT.at(3) += W * std::pow(t,4);
    }
    
    //Increment deposition value
    ADEP += W;
    
    //Increment absorption order
    Na += (double)p.n * W;
    N2a += (double)(p.n*p.n) * W;
    ta += p.t * W;
    t2a += p.t*p.t * W;
    
    //Time vector
    if (binT >= 0)
        At.at(binT) += W;
    binT = getBinF(p.v);
    if (binT >= 0)
        Aw.at(binT) += W;
}

void Stats::scatter(const Photon& p, double W) {
    //Initialize
    double t = p.t;
    const vec x = p.x;
    
    //Calculate moments
    if (p.isFluorescence()) {
        //Scattering data
        FSCAT += W;
        
        //Scattering fluorescence moments
        if (momentlvl > 0) {
            FSMOMENTS.at(0) += W * x.Z;
            FSMOMENTS.at(1) += W * x.r();
            FSTMOMENT.at(0) += W * t;
        } if (momentlvl > 1) {
            FSMOMENTS.at(2) += W * std::pow(x.Z,2);
            FSMOMENTS.at(3) += W * x.r() * x.Z;
            FSMOMENTS.at(4) += W * std::pow(x.r(),2);
            FSTMOMENT.at(1) += W * std::pow(t,2);
        } if (momentlvl > 2) {
            FSMOMENTS.at(5) += W * std::pow(x.Z,3);
            FSMOMENTS.at(6) += W * x.r() * std::pow(x.Z,2);
            FSMOMENTS.at(7) += W * x.Z * std::pow(x.r(),2);
            FSMOMENTS.at(8) += W * std::pow(x.r(),3);
            FSTMOMENT.at(2) += W * std::pow(t,3);
        } if (momentlvl > 3) {
            FSMOMENTS.at(9) += W * std::pow(x.Z,4);
            FSMOMENTS.at(10) += W * x.r() * std::pow(x.Z,3);
            FSMOMENTS.at(11) += W * std::pow(x.Z * x.r(),2);
            FSMOMENTS.at(12) += W * x.Z * std::pow(x.r(),3);
            FSMOMENTS.at(13) += W * std::pow(x.r(),4);
            FSTMOMENT.at(3) += W * std::pow(t,4);
        }
        return;
    }

    //Calculate scattering moments
    if (momentlvl > 0) {
        SMOMENTS.at(0) += W * x.Z;
        SMOMENTS.at(1) += W * x.r();
        STMOMENT.at(0) += W * t;
    } if (momentlvl > 1) {
        SMOMENTS.at(2) += W * std::pow(x.Z,2);
        SMOMENTS.at(3) += W * x.r() * x.Z;
        SMOMENTS.at(4) += W * std::pow(x.r(),2);
        STMOMENT.at(1) += W * std::pow(t,2);
    } if (momentlvl > 2) {
        SMOMENTS.at(5) += W * std::pow(x.Z,3);
        SMOMENTS.at(6) += W * x.r() * std::pow(x.Z,2);
        SMOMENTS.at(7) += W * x.Z * std::pow(x.r(),2);
        SMOMENTS.at(8) += W * std::pow(x.r(),3);
        STMOMENT.at(2) += W * std::pow(t,3);
    } if (momentlvl > 3) {
        SMOMENTS.at(9) += W * std::pow(x.Z,4);
        SMOMENTS.at(10) += W * x.r() * std::pow(x.Z,3);
        SMOMENTS.at(11) += W * std::pow(x.Z * x.r(),2);
        SMOMENTS.at(12) += W * x.Z * std::pow(x.r(),3);
        SMOMENTS.at(13) += W * std::pow(x.r(),4);
        STMOMENT.at(3) += W * std::pow(t,4);
    }
    
    //Increment scatter and deposition values
    ASCAT += W;
    
    //Calculate diffusion properties
    diffusion(x, W, t);
}

void Stats::diffusion(const vec& x, double w, double t) {
    DX += w * x.X*x.X / t; //Add to average diffusion coefficient
    DY += w * x.Y*x.Y / t; //Y component
    DZ += w * x.Z*x.Z / t; //Z component
    ND += w;
}

double Stats::getR(const vec& mu, const vec& n, double m, double& sint) const {    
    //Otherwise perform calculation
    double R, cost;
    cost = -n.dot(mu);
    if (1-cost < CONST_EPS) {
        R = (1.0-m)*(1.0-m)/(1.0+m)/(1.0+m);
        sint = 0;
    } else {
        cost = std::acos(cost); //theta_i
        if (1-std::abs(m*std::sin(cost)) < CONST_EPS) {
            R = 1;
            sint = 0;
        } else {
            sint = std::asin(m*std::sin(cost));       //theta_t
            R = 0.5*(std::pow(std::sin(cost-sint)/std::sin(cost+sint),2) 
                    + std::pow(std::tan(cost-sint)/std::tan(cost+sint),2));
        }
    }
    return R;
}

void Stats::reflect(const Photon& p, int surf, double R, double sint) {
    //Calculate theta from Snell's law
    double W;
    int binTHETA, binT, binW;
    
    //Calculate the intensity to add, and the T scale if we're static
    W = p.W;
    if (T == 0)
        T = 1;
    
    //Do moments first since they depend only on the power (not amplitude)
    //Bin in exit angle
    binTHETA = getBinTh(sint);
    
    //Do fluorescence first then quit
    if (p.isFluorescence()) {
        binT = getBinTf(p.t);
        binW = getBinF(p.v);
        switch (surf) {
            case 1:
                Fb += (1-R) * W;
                if (binTHETA+1)
                    Fbtheta.at(binTHETA) += W * (1-R);
                if (binT >= 0)
                    Fbt.at(binT) += (1-R) * W;
                if (binW >= 0)
                    Fbw.at(binW) += (1-R) * W;
                if (p.isBallistic())
                    Fbballistic += W * (1-R);
                break;
            case 2:
                Ff += (1-R) * W;
                if (binTHETA+1)
                    Fftheta.at(binTHETA) += W * (1-R);
                if (binT >= 0)
                    Fft.at(binT) += (1-R) * W;
                if (binW >= 0)
                    Ffw.at(binW) += (1-R) * W;
                if (p.isBallistic())
                    Ffballistic += W * (1-R);
                break;
            case 3:
                Fl += (1-R) * W;
                if (binT >= 0)
                    Flt.at(binT) += (1-R) * W;
                if (binW >= 0)
                    Flw.at(binW) += (1-R) * W;
                if (p.isBallistic())
                    Flballistic += W * (1-R);
                break;
        }
        
        //Transmission angle moments
        if (momentlvl > 0)
            FTHMOMENT.at(0+((surf == 1) ? 4 : 0)) += W * (1-R) * sint;
        if (momentlvl > 1)
            FTHMOMENT.at(1+((surf == 1) ? 4 : 0)) += W * (1-R) * std::pow(sint,2);
        if (momentlvl > 2)
            FTHMOMENT.at(2+((surf == 1) ? 4 : 0)) += W * (1-R) * std::pow(sint,3);
        if (momentlvl > 3)
            FTHMOMENT.at(3+((surf == 1) ? 4 : 0)) += W * (1-R) * std::pow(sint,4);
        
        return;
    }
    
    //Get moments
    binT = getBinT(p.t);
    binW = getBinF(p.v);
    if (p.n > 0) {
        if (momentlvl > 0)
            THMOMENT.at(0+((surf == 1) ? 4 : 0)) += W * (1-R) * sint;
        if (momentlvl > 1)
            THMOMENT.at(1+((surf == 1) ? 4 : 0)) += W * (1-R) * std::pow(sint,2);
        if (momentlvl > 2)
            THMOMENT.at(2+((surf == 1) ? 4 : 0)) += W * (1-R) * std::pow(sint,3);
        if (momentlvl > 3)
            THMOMENT.at(3+((surf == 1) ? 4 : 0)) += W * (1-R) * std::pow(sint,4);
    }    
        
    //Increment reflection/transmission times and orders but only for non-specular reflection
    if (surf == 1) {
        Nt += (double)p.n * (1-R) * W;
        N2t += (double)(p.n*p.n) * (1-R) * W;
        tTd += p.t * (1-R) * W;
        t2Td += p.t*p.t * (1-R) * W;
    } else if ((surf == 2) and (p.n > 0)) {
        Nr += (double)p.n * (1-R) * W;
        N2r += (double)(p.n*p.n) * (1-R) * W;
        tRd += p.t * (1-R) * W;
        t2Rd += p.t*p.t * (1-R) * W;
    }  
    
    //Increment transmission
    if (surf == 1) {
        Tdiffuse += (1-R) * W;
        if (binTHETA+1)
            Ttheta.at(binTHETA) += W * (1-R);
        if (binT >= 0)
            Tt.at(binT) += (1-R) * W;
        if (binW >= 0)
            Tw.at(binW) += (1-R) * W;
        if (p.isBallistic())
            Tballistic += W * (1-R);
    } else if (surf == 2) {
        Rdiffuse += (1-R) * W;
        if (binTHETA+1)
            Rtheta.at(binTHETA) += W * (1-R);
        if (binT >= 0)
            Rt.at(binT) += (1-R) * W;
        if (binW >= 0)
            Rw.at(binW) += (1-R) * W;
    } else if (surf == 3) {
        Ldiffuse += (1-R) * W;
        if (binT >= 0)
            Lt.at(binT) += (1-R) * W;
        if (binW >= 0)
            Lw.at(binW) += (1-R) * W;
    }
    
    //Increment beam radii
    if (surf == 1) {
        WRAD += W;
        RAD2 += W * p.x.r2();
    } else if (surf == 2) {
        WRAD += W;
        RAD2 += W * p.x.r2();
    }
}

double Stats::initialize(const Photon& p, double m) {
    //Calculate reflection coefficient and increment specular reflection
    double R, cost, sint;
    int binTHETA, binF, binT;
    if (1-std::abs(p.mu.Z) < CONST_EPS) {
        R = (1.0-m)*(1.0-m)/(1.0+m)/(1.0+m);
        sint = 0;
    } else {
        cost = std::acos(std::abs(p.mu.Z)); //theta_i
        if (1-std::abs(m*std::sin(cost)) < CONST_EPS) {
            R = 1;
            sint = 0;
        } else {
            sint = std::asin(m*std::sin(cost));       //theta_t
            R = 0.5*(std::pow(std::sin(cost-sint)/std::sin(cost+sint),2) 
                    + std::pow(std::tan(cost-sint)/std::tan(cost+sint),2));
        }
    }
    
    //Increment tracked photon energy
    PHI += p.W;
    
    //Calculate binTHETA, store theta resolved Rspec and Rspec
    binTHETA = getBinTh(sint);
    if (binTHETA+1)
        Rstheta.at(binTHETA) += p.W*R;
    Rspec += p.W*R;
    
    //Frequency spectrum
    binF = getBinF(p.v);
    if (binF + 1)
        Iw.at(binF) += p.W;
    
    //Time-domain
    binT = getBinT(p.t);
    if (binT + 1)
        It.at(binT) += p.W;
        
    //Return reflection coefficient
    return R;
}

void Stats::printTime(std::ostream& oout) const {
    //Print header
    oout << std::endl << "Temporal profiles [photons/time step]:" << std::endl;
    oout << "  T [ps]: ";
    for (int i = 0; i < Tres; i ++)
        oout << std::scientific << std::setw(18) << dT * i;
        
    //Now print the reflection/transmission/absorption data
    oout << std::endl << "Incident: ";
    for (int i = 0; i < Tres; i ++)
        oout << std::scientific << std::setw(18) << It.at(i);
    oout << std::endl << "Rdiffuse: ";
    for (int i = 0; i < Tres; i ++)
        oout << std::scientific << std::setw(18) << Rt.at(i);
    oout << std::endl << "Tdiffuse: ";
    for (int i = 0; i < Tres; i ++)
        oout << std::scientific << std::setw(18) << Tt.at(i);
    oout << std::endl << "Ldiffuse: ";
    for (int i = 0; i < Tres; i ++)
        oout << std::scientific << std::setw(18) << Lt.at(i);
    oout << std::endl << "Adep:     ";
    for (int i = 0; i < Tres; i ++)
        oout << std::scientific << std::setw(18) << At.at(i);
    
    //If no fluorescence, exit here
    if (FGEN < CONST_EPS)
        return;
        
    //Otherwise print fluorescence
    oout << std::endl << std::endl;
    oout << " Tf [ps]: ";
    for (int i = 0; i < Tres; i ++)
        oout << std::scientific << std::setw(18) << dTf * (i+1);
    oout << std::endl << "Femitted: ";
    for (int i = 0; i < Tres; i ++)
        oout << std::scientific << std::setw(18) << Ft.at(i);
    oout << std::endl << "Ffront:   ";
    for (int i = 0; i < Tres; i ++)
        oout << std::scientific << std::setw(18) << Fft.at(i);
    oout << std::endl << "Fback:    ";
    for (int i = 0; i < Tres; i ++)
        oout << std::scientific << std::setw(18) << Fbt.at(i);
    oout << std::endl << "Fside:    ";
    for (int i = 0; i < Tres; i ++)
        oout << std::scientific << std::setw(18) << Flt.at(i);
    oout << std::endl << "Freabs:   ";
    for (int i = 0; i < Tres; i ++)
        oout << std::scientific << std::setw(18) << Fat.at(i);
}

void Stats::printFreq(std::ostream& oout) const {
    //Print header
    oout << std::endl << "Radiation frequency profiles [photons/freq. bin]:" << std::endl;
    oout << " W [THz]: ";
    for (int i = 0; i < Fres; i ++)
        oout << std::scientific << std::setw(18) << dF * i + F0;
        
    //Now print the reflection/transmission/absorption data
    oout << std::endl << "Incident: ";
    for (int i = 0; i < Fres; i ++)
        oout << std::scientific << std::setw(18) << Iw.at(i);
    oout << std::endl << "Rdiffuse: ";
    for (int i = 0; i < Fres; i ++)
        oout << std::scientific << std::setw(18) << Rw.at(i);
    oout << std::endl << "Tdiffuse: ";
    for (int i = 0; i < Fres; i ++)
        oout << std::scientific << std::setw(18) << Tw.at(i);
    oout << std::endl << "Ldiffuse: ";
    for (int i = 0; i < Fres; i ++)
        oout << std::scientific << std::setw(18) << Lw.at(i);
    oout << std::endl << "Adep:     ";
    for (int i = 0; i < Fres; i ++)
        oout << std::scientific << std::setw(18) << Aw.at(i);
    
    //If no fluorescence, exit here
    if (FGEN < CONST_EPS)
        return;
        
    //Otherwise print fluorescence
    oout << std::endl << "Femitted: ";
    for (int i = 0; i < Fres; i ++)
        oout << std::scientific << std::setw(18) << Fw.at(i);
    oout << std::endl << "Ffront:   ";
    for (int i = 0; i < Fres; i ++)
        oout << std::scientific << std::setw(18) << Ffw.at(i);
    oout << std::endl << "Fback:    ";
    for (int i = 0; i < Fres; i ++)
        oout << std::scientific << std::setw(18) << Fbw.at(i);
    oout << std::endl << "Fside:    ";
    for (int i = 0; i < Fres; i ++)
        oout << std::scientific << std::setw(18) << Flw.at(i);
    oout << std::endl << "Freabs:   ";
    for (int i = 0; i < Fres; i ++)
        oout << std::scientific << std::setw(18) << Faw.at(i);
}

void Stats::print(std::ostream& oout) const {
    
    //Theta-dependent reflection/transmission/fluorescence
    oout << std::endl;
    oout << "theta:      ";
    for (int i = 0; i < THETAres; i ++)
        oout << std::fixed << std::setw(18) << i*dTHETA;
    oout << std::endl;
    oout << "R(theta):       ";
    for (int i = 0; i < THETAres; i ++)
        oout << std::scientific << std::setw(18) << Rtheta.at(i)/PHI;
    oout << std::endl;
    oout << "Rs(theta):      ";
    for (int i = 0; i < THETAres; i ++)
        oout << std::scientific << std::setw(18) << Rstheta.at(i)/PHI;
    oout << std::endl;
    oout << "T(theta):       ";
    for (int i = 0; i < THETAres; i ++)
        oout << std::scientific << std::setw(18) << Ttheta.at(i)/PHI;
    oout << std::endl;
    
    //Print fluorescence part
    if (FGEN > CONST_EPS) {
        oout << "Ffront(theta):  ";
        for (int i = 0; i < THETAres; i ++)
            oout << std::scientific << std::setw(18) << Fftheta.at(i)/PHI;
        oout << std::endl;
        oout << "Fback(theta):   ";
        for (int i = 0; i < THETAres; i ++)
            oout << std::scientific << std::setw(18) << Fbtheta.at(i)/PHI;
        oout << std::endl;
    }
    
    //Deposition energy moments (time)
    if (momentlvl > 0) {
        oout << std::endl;
        oout << "Incident scatter";
        if (ADEP > CONST_EPS)
            oout << ", Incident absorption";
        if (FGEN > CONST_EPS)
            oout << ", Fluorescence scatter, Fluorescence reabsorption";
        oout << " moments (time) [ps**n]:" << std::endl;
        oout << "  E[T]     = " << STMOMENT[0]/ASCAT;
        if (ADEP > CONST_EPS)
            oout << "     " << TMOMENT[0]/ADEP;
        if (FGEN > CONST_EPS)
            oout << "     " << FSTMOMENT[0]/FSCAT << "     " << FTMOMENT[0]/FDEP;
            
    } if (momentlvl > 1) {
        oout << std::endl << "  E[T2]    = " << STMOMENT[1]/ASCAT;
        if (ADEP > CONST_EPS)
            oout << "     " << TMOMENT[1]/ADEP;
        if (FGEN > CONST_EPS)
            oout << "     " << FSTMOMENT[1]/FSCAT << "     " << FTMOMENT[1]/FDEP;
        
    } if (momentlvl > 2) {
        oout << std::endl << "  E[T3]    = " << STMOMENT[2]/ASCAT;
        if (ADEP > CONST_EPS)
            oout << "     " << TMOMENT[2]/ADEP;
        if (FGEN > CONST_EPS)
            oout << "     " << FSTMOMENT[2]/FSCAT << "     " << FTMOMENT[2]/FDEP;
            
    } if (momentlvl > 3) {
        oout << std::endl << "  E[T4]    = " << STMOMENT[3]/ASCAT;
        if (ADEP > CONST_EPS)
            oout << "     " << TMOMENT[3]/ADEP;
        if (FGEN > CONST_EPS)
            oout << "     " << FSTMOMENT[3]/FSCAT << "     " << FTMOMENT[3]/FDEP;
    }
    oout << std::endl;
    
    //Deposition energy moments (space)
    if (momentlvl > 0) {
    
        //Header
        oout << std::endl;
        oout << "Incident scatter";
        if (ADEP > CONST_EPS)
            oout << ", Incident absorption";
        if (FGEN > CONST_EPS)
            oout << ", Fluorescence scatter, Fluorescence reabsorption";
        oout << " moments (space) [um**n]:" << std::endl;
        
        //Moment Z
        oout << "  E[Z]    = " << SMOMENTS[0]/ASCAT;
        if (ADEP > CONST_EPS)
            oout << "     " << MOMENTS[0]/ADEP;
        if (FGEN > CONST_EPS)
            oout << "     " << FSMOMENTS[0]/FSCAT << "     " << FMOMENTS[0]/FDEP;
        
        //Moment R
        oout << std::endl << "  E[R]    = " << SMOMENTS[1]/ASCAT;
        if (ADEP > CONST_EPS)
            oout << "     " << MOMENTS[1]/ADEP;
        if (FGEN > CONST_EPS)
            oout << "     " << FSMOMENTS[1]/FSCAT << "     " << FMOMENTS[1]/FDEP;
            
    } if (momentlvl > 1) {
        //Moment Z*Z
        oout << std::endl << "  E[Z2]   = " << SMOMENTS[2]/ASCAT;
        if (ADEP > CONST_EPS)
            oout << "     " << MOMENTS[2]/ADEP;
        if (FGEN > CONST_EPS)
            oout << "     " << FSMOMENTS[2]/FSCAT << "     " << FMOMENTS[2]/FDEP;
        
        //Moment R*Z
        oout << std::endl << "  E[RZ]   = " << SMOMENTS[3]/ASCAT;
        if (ADEP > CONST_EPS)
            oout << "     " << MOMENTS[3]/ADEP;
        if (FGEN > CONST_EPS)
            oout << "     " << FSMOMENTS[3]/FSCAT << "     " << FMOMENTS[3]/FDEP;
        
        //Moment R*R
        oout << std::endl << "  E[R2]   = " << SMOMENTS[4]/ASCAT;
        if (ADEP > CONST_EPS)
            oout << "     " << MOMENTS[4]/ADEP;
        if (FGEN > CONST_EPS)
            oout << "     " << FSMOMENTS[4]/FSCAT << "     " << FMOMENTS[4]/FDEP;
        
    } if (momentlvl > 2) {
        //Moment Z*Z*Z
        oout << std::endl << "  E[Z3]   = " << SMOMENTS[5]/ASCAT;
        if (ADEP > CONST_EPS)
            oout << "     " << MOMENTS[5]/ADEP;
        if (FGEN > CONST_EPS)
            oout << "     " << FSMOMENTS[5]/FSCAT << "     " << FMOMENTS[5]/FDEP;

        //Moment Z*Z*R
        oout << std::endl << "  E[Z2R]  = " << SMOMENTS[6]/ASCAT;
        if (ADEP > CONST_EPS)
            oout << "     " << MOMENTS[6]/ADEP;
        if (FGEN > CONST_EPS)
            oout << "     " << FSMOMENTS[6]/FSCAT << "     " << FMOMENTS[6]/FDEP;
        
        //Moment Z*R*R
        oout << std::endl << "  E[ZR2]  = " << SMOMENTS[7]/ASCAT;
        if (ADEP > CONST_EPS)
            oout << "     " << MOMENTS[7]/ADEP;
        if (FGEN > CONST_EPS)
            oout << "     " << FSMOMENTS[7]/FSCAT << "     " << FMOMENTS[7]/FDEP;
        
        //Moment R*R*R
        oout << std::endl << "  E[R3]   = " << SMOMENTS[8]/ASCAT;
        if (ADEP > CONST_EPS)
            oout << "     " << MOMENTS[8]/ADEP;
        if (FGEN > CONST_EPS)
            oout << "     " << FSMOMENTS[8]/FSCAT << "     " << FMOMENTS[8]/FDEP;
            
    } if (momentlvl > 3) {
        //Moment Z*Z*Z*Z
        oout << std::endl << "  E[Z4]   = " << SMOMENTS[9]/ASCAT;
        if (ADEP > CONST_EPS)
            oout << "     " << MOMENTS[9]/ADEP;
        if (FGEN > CONST_EPS)
            oout << "     " << FSMOMENTS[9]/FSCAT << "     " << FMOMENTS[9]/FDEP;
        
        //Moment Z*Z*Z*R
        oout << std::endl << "  E[Z3R]  = " << SMOMENTS[10]/ASCAT;
        if (ADEP > CONST_EPS)
            oout << "     " << MOMENTS[10]/ADEP;
        if (FGEN > CONST_EPS)
            oout << "     " << FSMOMENTS[10]/FSCAT << "     " << FMOMENTS[10]/FDEP;
        
        //Moment Z*Z*R*R
        oout << std::endl << "  E[Z2R2] = " << SMOMENTS[11]/ASCAT;
        if (ADEP > CONST_EPS)
            oout << "     " << MOMENTS[11]/ADEP;
        if (FGEN > CONST_EPS)
            oout << "     " << FSMOMENTS[11]/FSCAT << "     " << FMOMENTS[11]/FDEP;
        
        //Moment Z*R*R
        oout << std::endl << "  E[ZR3]  = " << SMOMENTS[12]/ASCAT;
        if (ADEP > CONST_EPS)
            oout << "     " << MOMENTS[12]/ADEP;
        if (FGEN > CONST_EPS)
            oout << "     " << FSMOMENTS[12]/FSCAT << "     " << FMOMENTS[12]/FDEP;
        
        //Moment R*R*R*R
        oout << std::endl << "  E[R4]   = " << SMOMENTS[13]/ASCAT;
        if (ADEP > CONST_EPS)
            oout << "     " << MOMENTS[13]/ADEP;
        if (FGEN > CONST_EPS)
            oout << "     " << FSMOMENTS[13]/FSCAT << "     " << FMOMENTS[13]/FDEP;
    }
    
    //Transmission angle moments
    if (momentlvl > 0) {
        oout << std::endl << std::endl;
        oout << "Transmission";
        if (FGEN > CONST_EPS)
            oout << "/Emission";
        oout << " angle moments (R/T";
        if (FGEN > CONST_EPS)
            oout << "/Ff/Fb";
        oout << " [-], E[THETA**n] [rad**n]):";
        oout << std::endl << "  R     = " << std::fixed << Rdiffuse/PHI << ",";
        for (int i = 0; i < 4; i ++)
            oout << " " << std::fixed << std::setw(10) << THMOMENT[i]/Rdiffuse;
        oout << std::endl << "  T     = " << std::fixed << Tdiffuse/PHI << ",";
        for (int i = 4; i < 8; i ++)
            oout << " " << std::fixed << std::setw(10) << THMOMENT[i]/Tdiffuse;
        if (FGEN > CONST_EPS) {
            oout << std::endl << "  Ff    = " << std::fixed << Ff/PHI << ",";
            for (int i = 0; i < 4; i ++)
                oout << " " << std::fixed << std::setw(10) << FTHMOMENT[i]/Ff;
            oout << std::endl << "  Fb    = " << std::fixed << Fb/PHI << ",";
            for (int i = 4; i < 8; i ++)
                oout << " " << std::fixed << std::setw(10) << FTHMOMENT[i]/Fb;
        
        }
    }
    oout << std::endl;
    
    //Print time dependence
    printTime(oout);
    oout << std::endl;
    
    //Print frequency dependence
    printFreq(oout);
    oout << std::endl << std::endl;
    
    //Reflection coefficients (total)
    oout << "Reflection/Transmission/Absorption";
    if (FGEN > CONST_EPS)
        oout << "/Emission";
    oout << " coefficients [photon/source photon]:" << std::endl;
    oout << "  Rdiffuse   = " << std::fixed << std::setw(10) << Rdiffuse/PHI;
    if (FGEN > CONST_EPS)
        oout << "     " << "  Fgen       = " << std::fixed << std::setw(10) << FGEN/PHI;
    oout << std::endl << "  Rspecular  = " << std::fixed << std::setw(10) << Rspec/PHI;
    if (FGEN > CONST_EPS)
        oout << "     " << "  Ffront     = " << std::fixed << std::setw(10) << Ff/PHI;
    oout << std::endl << "  Tdiffuse   = " << std::fixed << std::setw(10) << (Tdiffuse-Tballistic)/PHI;
    if (FGEN > CONST_EPS)
        oout << "     " << "  Fback      = " << std::fixed << std::setw(10) << Fb/PHI;
    oout << std::endl << "  Tballistic = " << std::fixed << std::setw(10) << Tballistic/PHI;
    if (FGEN > CONST_EPS)
        oout << "     " << "  Fside      = " << std::fixed << std::setw(10) << Fl/PHI;
    oout << std::endl << "  Ldiffuse   = " << std::fixed << std::setw(10) << Ldiffuse/PHI;
    if (FGEN > CONST_EPS)
        oout << "     " << "  Fabs       = " << std::fixed << std::setw(10) << FDEP/PHI;
    oout << std::endl << "  A          = " << std::fixed << std::setw(10) << ADEP/PHI << std::endl << std::endl;
    oout << std::endl << "Checking for math errors:" << std::endl;
    oout << "     1 - (Rd + Rs + Td + Tb + A + Ld)   ?=  0  =  " << 1-(Rdiffuse + Rspec + Tdiffuse + ADEP + Ldiffuse) / PHI << std::endl;
    if (FGEN > CONST_EPS)
        oout << "     Fg - (Ff + Fb + Fl + Fa)      ?=  0  =  " << (FGEN - Ff - Fb - Fl - FDEP) / PHI << std::endl;
    if (FGEN > CONST_EPS)
        oout << "     Fg / A   =  " << FGEN / ADEP << std::endl;
    
    //Reflection time/order
    oout << std::endl << std::endl;
    oout << "Diffuse Reflection/Transmission scattering order (mean, std) [-] and duration (mean, std) [ps]:" << std::endl;
    if (Rdiffuse > CONST_EPS)
        oout << "  N[Rdiffuse]      = " << std::fixed << std::setw(10) << Nr / Rdiffuse << std::endl;
    if (Tdiffuse > CONST_EPS)
        oout << "  N[Tdiffuse]      = " << std::fixed << std::setw(10) << Nt / Tdiffuse << std::endl;
    if (ADEP > CONST_EPS)
        oout << "  N[A]             = " << std::fixed << std::setw(10) << Na / ADEP << std::endl;
    if (Rdiffuse > CONST_EPS)
        oout << "  dN[Rdiffuse]     = " << std::fixed << std::setw(10) << std::sqrt(N2r*Rdiffuse - Nr*Nr) / Rdiffuse << std::endl;
    if (Tdiffuse > CONST_EPS)
        oout << "  dN[Tdiffuse]     = " << std::fixed << std::setw(10) << std::sqrt(N2t*Tdiffuse - Nt*Nt) / Tdiffuse << std::endl;
    if (ADEP > CONST_EPS)
        oout << "  dN[A]            = " << std::fixed << std::setw(10) << std::sqrt(N2a*ADEP - Na*Na) / ADEP << std::endl;
    if (Rdiffuse > CONST_EPS)
        oout << "  t[Rdiffuse]      = " << std::fixed << std::setw(10) << tRd / Rdiffuse << std::endl;
    if (Tdiffuse > CONST_EPS)
        oout << "  t[Tdiffuse]      = " << std::fixed << std::setw(10) << tTd / Tdiffuse << std::endl;
    if (ADEP > CONST_EPS)
        oout << "  t[A]             = " << std::fixed << std::setw(10) << ta / ADEP << std::endl;
    if (Rdiffuse > CONST_EPS)
        oout << "  dt[Rdiffuse]     = " << std::fixed << std::setw(10) << std::sqrt(t2Rd*Rdiffuse - tRd*tRd) / Rdiffuse << std::endl;
    if (Tdiffuse > CONST_EPS)
        oout << "  dt[Tdiffuse]     = " << std::fixed << std::setw(10) << std::sqrt(t2Td*Tdiffuse - tTd*tTd) / Tdiffuse << std::endl;
    if (ADEP > CONST_EPS)
        oout << "  dt[A]            = " << std::fixed << std::setw(10) << std::sqrt(t2a*ADEP - ta*ta) / ADEP << std::endl;

    //Diffusion coefficient
    oout << std::endl << std::endl;
    oout << "Diffusion coefficients [um2/ps]:" << std::endl;
    oout << "  Dx = " << DX/ND/2.0 << std::endl;
    oout << "  Dy = " << DY/ND/2.0 << std::endl;
    oout << "  Dz = " << DZ/ND/2.0 << std::endl;
}

void Stats::writedbheader(std::ostream& OF) const {
    OF << "Trans,Ref0,RefD,TransS,TransB,Abs,FGen,FAbs,Ff,Fb,Fs,Ffb,Fbb,Fsb,Pen,RadL,Rad0";
}

void Stats::writedb(std::ostream& OF) const {
    OF << std::scientific << std::setprecision(8) << Tdiffuse/PHI << ",";
    OF << std::scientific << std::setprecision(8) << Rspec/PHI << ",";
    OF << std::scientific << std::setprecision(8) << Rdiffuse/PHI << ",";
    OF << std::scientific << std::setprecision(8) << Ldiffuse/PHI << ",";
    OF << std::scientific << std::setprecision(8) << Tballistic/PHI << ",";
    OF << std::scientific << std::setprecision(8) << ADEP/PHI << ",";
    OF << std::scientific << std::setprecision(8) << FGEN/PHI << ",";
    OF << std::scientific << std::setprecision(8) << FDEP/PHI << ",";
    OF << std::scientific << std::setprecision(8) << Ff/PHI << ",";
    OF << std::scientific << std::setprecision(8) << Fb/PHI << ",";
    OF << std::scientific << std::setprecision(8) << Fl/PHI << ",";
    OF << std::scientific << std::setprecision(8) << Ffballistic/PHI << ",";
    OF << std::scientific << std::setprecision(8) << Fbballistic/PHI << ",";
    OF << std::scientific << std::setprecision(8) << Flballistic/PHI << ",";
    OF << std::scientific << std::setprecision(8) << SMOMENTS[0]/ASCAT << ",";
    OF << std::scientific << std::setprecision(8) << std::sqrt(RAD2/WRAD) << ",";
    OF << std::scientific << std::setprecision(8) << std::sqrt(RAD2F/WRADF);
}

