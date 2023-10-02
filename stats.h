#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>

#include "vec.h"
#include "photon.h"

#ifndef _STATS_H_
#define _STATS_H_

#define CONST_PI 3.1415926535
#define CONST_EPS 1e-10

using namespace std;

struct Stats {
    //Settings
    int Tres, THETAres, momentlvl;
    double dT, dTf, dTHETA, T, PHI;
    double Rdiffuse, Tdiffuse, Tballistic, Rspec, Ldiffuse;
    
    //Output arrays for energy deposition, scatter, etc.
    vector<double> Rtheta, Rstheta, Ttheta;
    
    //Initialize moments and other outputs
    vector<double> MOMENTS, SMOMENTS, THMOMENT, TMOMENT, STMOMENT;
    double RAD2, RAD2F, WRAD, WRADF;
    
    //Time-dependence
    vector<double> Tt, Rt, Lt, At;
    
    //Diffusion coefficients
    double DX, DY, DZ, ND;
    
    //Transmission/reflection scattering order
    double Nr, Nt, N2r, N2t;
    double Na, ta, N2a, t2a;
    
    //Transmission/reflection times
    double tRd, tTd, t2Rd, t2Td;
    
    //Fluorescence properties
    double Ff, Fb, Fl, Ffballistic, Fbballistic, Flballistic;   //Front, Back, Side transmission (and ballistic)
    vector<double> Fftheta, Fbtheta;                            //Front and back angle-dependence
    vector<double> Fft, Fbt, Flt, Fat, Ft;                      //Front, back, side, reabsorption, emission time-dependence
    vector<double> FMOMENTS, FSMOMENTS, FTHMOMENT, FTMOMENT, FSTMOMENT;
    
    //Cache
    double ADEP, ASCAT, FDEP, FSCAT, FGEN;                      //Incident absorption/scatter, fluorescence absorption/scatter/generation
    
    //Functions
    //--Constructors
    Stats();
    Stats(int N, double dT, int LVL=4);
    Stats(int nT, int nTH, double dT, double dTf, int LVL=4);
    
    //--Accessors
    int getBinT(double t) const;
    int getBinTf(double t) const;
    int getBinTh(double sint) const;
    
    //--Print results functions
    void setup();
    void print() const;
    void printTime() const;
    
    void writedb(ofstream& OF) const;
    void writedbheader(ofstream& OF) const;
    
    //--Functions to aggregate data
    void scatter(const Photon& p, double w);
    void absorb(const Photon& p, double w);
    void emit(const Photon& p, double w);
    void diffusion(const vec& x, double w, double t);
    double initialize(const Photon& x, double m);
    void reflect(const Photon& p, int surf, double R, double sint);
    double getR(const vec& mu, const vec& n, double m, double& sint) const;
};

Stats::Stats() {
    Stats(20, 1.0, 4);
}

Stats::Stats(int n, double dT, int LVL) {
    Stats(n, n, dT, LVL);
}

Stats::Stats(int nT, int nTH, double dT, double dTf, int LVL) {
    //Set default resolution
    Tres = nT;             //Number of time samples
    THETAres = nTH;        //Number of cos(theta) values for reflection
    momentlvl = LVL;       //Order of moments to calculate (only up to 4th order is implemented)
    this->dT = dT;
    this->dTf = dTf;
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
    MOMENTS = vector<double>(14, 0.0);
    THMOMENT = vector<double>(8, 0.0);
    TMOMENT = vector<double>(4, 0.0);
    SMOMENTS = vector<double>(14, 0.0);
    STMOMENT = vector<double>(4, 0.0);
    
    //Fluorescence moments
    FMOMENTS = vector<double>(14, 0.0);
    FTHMOMENT = vector<double>(8, 0.0);
    FTMOMENT = vector<double>(4, 0.0);
    FSMOMENTS = vector<double>(14, 0.0);
    FSTMOMENT = vector<double>(4, 0.0);
}

void Stats::setup() {
    //Grid resolution - note; recalculate here just in case
    T = dT * Tres;
    dTHETA = CONST_PI/THETAres/2;
    
    //Time arrays
    Tt = vector<double>(Tres, 0.0);
    Rt = vector<double>(Tres, 0.0);
    Lt = vector<double>(Tres, 0.0);
    At = vector<double>(Tres, 0.0);
    
    //Reflection and transmission
    Rtheta = vector<double>(THETAres, 0.0);
    Rstheta = vector<double>(THETAres, 0.0);
    Ttheta = vector<double>(THETAres, 0.0);
    
    //Fluorescence values
    Fftheta = vector<double>(THETAres, 0.0);
    Fbtheta = vector<double>(THETAres, 0.0);
    Fft = vector<double>(Tres, 0.0);
    Fbt = vector<double>(Tres, 0.0);
    Flt = vector<double>(Tres, 0.0);
    Fat = vector<double>(Tres, 0.0);
    Ft = vector<double>(Tres, 0.0);
    
    //Initialize moments and other outputs
    MOMENTS = vector<double>(14, 0.0);
    THMOMENT = vector<double>(8, 0.0);
    TMOMENT = vector<double>(4, 0.0);
    SMOMENTS = vector<double>(14, 0.0);
    STMOMENT = vector<double>(4, 0.0);
    
    //Fluorescence moments
    FMOMENTS = vector<double>(14, 0.0);
    FTHMOMENT = vector<double>(8, 0.0);
    FTMOMENT = vector<double>(4, 0.0);
    FSMOMENTS = vector<double>(14, 0.0);
    FSTMOMENT = vector<double>(4, 0.0);
}

int Stats::getBinT(double t) const {
    int binT = (int)(t/dT);
    if ((binT >= Tres) or (binT < 0))
        return -1;
    return binT;
}

int Stats::getBinTf(double t) const {
    int binT = (int)((t-dTf)/dTf);
    if ((binT >= Tres) or (binT < 0))
        return -1;
    return binT;
}

int Stats::getBinTh(double sint) const {
    int binTheta = (int) (sint/dTHETA);
    if (binTheta >= THETAres)
        binTheta = THETAres-1;
    else if (binTheta < 0)
        return -1;
    return binTheta;
}

void Stats::emit(const Photon& p, double W) {
    FGEN += W;
    int binT = getBinTf(p.t);
    if (binT >= 0)
        Ft.at(binT) += W;
}

void Stats::absorb(const Photon& p, double W) {
    //Initialize
    double t = p.t;
    const vec& x = p.x;
    int binT;
    
    //Calculate moments
    if (p.isFluorescence()) {
        //Temporal absorption data
        FDEP += W;
        binT = getBinTf(t);
        if (binT >= 0)
            Fat.at(binT) += W;
        
        //Absorption moments for fluorescence
        if (momentlvl > 0) {
            FMOMENTS.at(0) += W * x.Z;
            FMOMENTS.at(1) += W * x.r();
            FTMOMENT.at(0) += W * t;
        } if (momentlvl > 1) {
            FMOMENTS.at(2) += W * pow(x.Z,2);
            FMOMENTS.at(3) += W * x.r() * x.Z;
            FMOMENTS.at(4) += W * pow(x.r(),2);
            FTMOMENT.at(1) += W * pow(t,2);
        } if (momentlvl > 2) {
            FMOMENTS.at(5) += W * pow(x.Z,3);
            FMOMENTS.at(6) += W * x.r() * pow(x.Z,2);
            FMOMENTS.at(7) += W * x.Z * pow(x.r(),2);
            FMOMENTS.at(8) += W * pow(x.r(),3);
            FTMOMENT.at(2) += W * pow(t,3);
        } if (momentlvl > 3) {
            FMOMENTS.at(9) += W * pow(x.Z,4);
            FMOMENTS.at(10) += W * x.r() * pow(x.Z,3);
            FMOMENTS.at(11) += W * pow(x.Z * x.r(),2);
            FMOMENTS.at(12) += W * x.Z * pow(x.r(),3);
            FMOMENTS.at(13) += W * pow(x.r(),4);
            FTMOMENT.at(3) += W * pow(t,4);
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
        MOMENTS.at(2) += W * pow(x.Z,2);
        MOMENTS.at(3) += W * x.r() * x.Z;
        MOMENTS.at(4) += W * pow(x.r(),2);
        TMOMENT.at(1) += W * pow(t,2);
    } if (momentlvl > 2) {
        MOMENTS.at(5) += W * pow(x.Z,3);
        MOMENTS.at(6) += W * x.r() * pow(x.Z,2);
        MOMENTS.at(7) += W * x.Z * pow(x.r(),2);
        MOMENTS.at(8) += W * pow(x.r(),3);
        TMOMENT.at(2) += W * pow(t,3);
    } if (momentlvl > 3) {
        MOMENTS.at(9) += W * pow(x.Z,4);
        MOMENTS.at(10) += W * x.r() * pow(x.Z,3);
        MOMENTS.at(11) += W * pow(x.Z * x.r(),2);
        MOMENTS.at(12) += W * x.Z * pow(x.r(),3);
        MOMENTS.at(13) += W * pow(x.r(),4);
        TMOMENT.at(3) += W * pow(t,4);
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
            FSMOMENTS.at(2) += W * pow(x.Z,2);
            FSMOMENTS.at(3) += W * x.r() * x.Z;
            FSMOMENTS.at(4) += W * pow(x.r(),2);
            FSTMOMENT.at(1) += W * pow(t,2);
        } if (momentlvl > 2) {
            FSMOMENTS.at(5) += W * pow(x.Z,3);
            FSMOMENTS.at(6) += W * x.r() * pow(x.Z,2);
            FSMOMENTS.at(7) += W * x.Z * pow(x.r(),2);
            FSMOMENTS.at(8) += W * pow(x.r(),3);
            FSTMOMENT.at(2) += W * pow(t,3);
        } if (momentlvl > 3) {
            FSMOMENTS.at(9) += W * pow(x.Z,4);
            FSMOMENTS.at(10) += W * x.r() * pow(x.Z,3);
            FSMOMENTS.at(11) += W * pow(x.Z * x.r(),2);
            FSMOMENTS.at(12) += W * x.Z * pow(x.r(),3);
            FSMOMENTS.at(13) += W * pow(x.r(),4);
            FSTMOMENT.at(3) += W * pow(t,4);
        }
        return;
    }

    //Calculate scattering moments
    if (momentlvl > 0) {
        SMOMENTS.at(0) += W * x.Z;
        SMOMENTS.at(1) += W * x.r();
        STMOMENT.at(0) += W * t;
    } if (momentlvl > 1) {
        SMOMENTS.at(2) += W * pow(x.Z,2);
        SMOMENTS.at(3) += W * x.r() * x.Z;
        SMOMENTS.at(4) += W * pow(x.r(),2);
        STMOMENT.at(1) += W * pow(t,2);
    } if (momentlvl > 2) {
        SMOMENTS.at(5) += W * pow(x.Z,3);
        SMOMENTS.at(6) += W * x.r() * pow(x.Z,2);
        SMOMENTS.at(7) += W * x.Z * pow(x.r(),2);
        SMOMENTS.at(8) += W * pow(x.r(),3);
        STMOMENT.at(2) += W * pow(t,3);
    } if (momentlvl > 3) {
        SMOMENTS.at(9) += W * pow(x.Z,4);
        SMOMENTS.at(10) += W * x.r() * pow(x.Z,3);
        SMOMENTS.at(11) += W * pow(x.Z * x.r(),2);
        SMOMENTS.at(12) += W * x.Z * pow(x.r(),3);
        SMOMENTS.at(13) += W * pow(x.r(),4);
        STMOMENT.at(3) += W * pow(t,4);
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
    double R, cost;
    cost = -n.dot(mu);
    if (1-cost < CONST_EPS) {
        R = (1.0-m)*(1.0-m)/(1.0+m)/(1.0+m);
        sint = 0;
    } else {
        cost = acos(cost); //theta_i
        if (1-abs(m*sin(cost)) < CONST_EPS) {
            R = 1;
            sint = 0;
        } else {
            sint = asin(m*sin(cost));       //theta_t
            R = 0.5*(pow(sin(cost-sint)/sin(cost+sint),2) 
                    + pow(tan(cost-sint)/tan(cost+sint),2));
        }
    }
    return R;
}

void Stats::reflect(const Photon& p, int surf, double R, double sint) {
    //Calculate theta from Snell's law
    double W;
    int binTHETA, binT;
    
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
        switch (surf) {
            case 1:
                Fb += (1-R) * W;
                if (binTHETA+1)
                    Fbtheta.at(binTHETA) += W * (1-R);
                if (binT >= 0)
                    Fbt.at(binT) += (1-R) * W;
                if (p.isBallistic())
                    Fbballistic += W * (1-R);
                break;
            case 2:
                Ff += (1-R) * W;
                if (binTHETA+1)
                    Fftheta.at(binTHETA) += W * (1-R);
                if (binT >= 0)
                    Fft.at(binT) += (1-R) * W;
                if (p.isBallistic())
                    Ffballistic += W * (1-R);
                break;
            case 3:
                Fl += (1-R) * W;
                if (binT >= 0)
                    Flt.at(binT) += (1-R) * W;
                if (p.isBallistic())
                    Flballistic += W * (1-R);
                break;
        }
        
        //Transmission angle moments
        if (momentlvl > 0)
            FTHMOMENT.at(0+((surf == 1) ? 4 : 0)) += W * (1-R) * sint;
        if (momentlvl > 1)
            FTHMOMENT.at(1+((surf == 1) ? 4 : 0)) += W * (1-R) * pow(sint,2);
        if (momentlvl > 2)
            FTHMOMENT.at(2+((surf == 1) ? 4 : 0)) += W * (1-R) * pow(sint,3);
        if (momentlvl > 3)
            FTHMOMENT.at(3+((surf == 1) ? 4 : 0)) += W * (1-R) * pow(sint,4);
        
        return;
    }
    
    //Get moments
    binT = getBinT(p.t);
    if (p.n > 0) {
        if (momentlvl > 0)
            THMOMENT.at(0+((surf == 1) ? 4 : 0)) += W * (1-R) * sint;
        if (momentlvl > 1)
            THMOMENT.at(1+((surf == 1) ? 4 : 0)) += W * (1-R) * pow(sint,2);
        if (momentlvl > 2)
            THMOMENT.at(2+((surf == 1) ? 4 : 0)) += W * (1-R) * pow(sint,3);
        if (momentlvl > 3)
            THMOMENT.at(3+((surf == 1) ? 4 : 0)) += W * (1-R) * pow(sint,4);
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
        if (p.isBallistic())
            Tballistic += W * (1-R);
    } else if (surf == 2) {
        Rdiffuse += (1-R) * W;
        if (binTHETA+1)
            Rtheta.at(binTHETA) += W * (1-R);
        if (binT >= 0)
            Rt.at(binT) += (1-R) * W;
    } else if (surf == 3) {
        Ldiffuse += (1-R) * W;
        if (binT >= 0)
            Lt.at(binT) += (1-R) * W;
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
    int binTHETA;
    if (1-abs(p.mu.Z) < CONST_EPS) {
        R = (1.0-m)*(1.0-m)/(1.0+m)/(1.0+m);
        sint = 0;
    } else {
        cost = acos(abs(p.mu.Z)); //theta_i
        if (1-abs(m*sin(cost)) < CONST_EPS) {
            R = 1;
            sint = 0;
        } else {
            sint = asin(m*sin(cost));       //theta_t
            R = 0.5*(pow(sin(cost-sint)/sin(cost+sint),2) 
                    + pow(tan(cost-sint)/tan(cost+sint),2));
        }
    }
    
    //Increment tracked photon energy
    PHI += p.W;
    
    //Calculate binTHETA, store theta resolved Rspec and Rspec
    binTHETA = getBinTh(sint);
    
    //Increment Rspec
    if (binTHETA+1)
        Rstheta.at(binTHETA) += p.W*R;
    Rspec += p.W*R;
        
    //Return reflection coefficient
    return R;
}

void Stats::printTime() const {
    //Print header
    cout << endl << "Temporal profiles [photons/time step]:" << endl;
    cout << "  T [ps]: ";
    for (int i = 0; i < Tres; i ++)
        cout << scientific << setw(18) << dT * i;
        
    //Now print the reflection/transmission/absorption data
    cout << endl << "Rdiffuse: ";
    for (int i = 0; i < Tres; i ++)
        cout << scientific << setw(18) << Rt.at(i);
    cout << endl << "Tdiffuse: ";
    for (int i = 0; i < Tres; i ++)
        cout << scientific << setw(18) << Tt.at(i);
    cout << endl << "Ldiffuse: ";
    for (int i = 0; i < Tres; i ++)
        cout << scientific << setw(18) << Lt.at(i);
    cout << endl << "Adep:     ";
    for (int i = 0; i < Tres; i ++)
        cout << scientific << setw(18) << At.at(i);
    
    //If no fluorescence, exit here
    if (FGEN < CONST_EPS)
        return;
        
    //Otherwise print fluorescence
    cout << endl << endl;
    cout << " Tf [ps]: ";
    for (int i = 0; i < Tres; i ++)
        cout << scientific << setw(18) << dTf * (i+1);
    cout << endl << "Femitted: ";
    for (int i = 0; i < Tres; i ++)
        cout << scientific << setw(18) << Ft.at(i);
    cout << endl << "Ffront:   ";
    for (int i = 0; i < Tres; i ++)
        cout << scientific << setw(18) << Fft.at(i);
    cout << endl << "Fback:    ";
    for (int i = 0; i < Tres; i ++)
        cout << scientific << setw(18) << Fbt.at(i);
    cout << endl << "Fside:    ";
    for (int i = 0; i < Tres; i ++)
        cout << scientific << setw(18) << Flt.at(i);
    cout << endl << "Freabs:   ";
    for (int i = 0; i < Tres; i ++)
        cout << scientific << setw(18) << Fat.at(i);
}

void Stats::print() const {
    
    //Theta-dependent reflection/transmission/fluorescence
    cout<<endl;
    cout<<"cos(theta):     ";
    for (int i = 0; i < THETAres; i ++)
        cout << fixed << setw(18) << i*dTHETA;
    cout<<endl;
    cout<<"R(theta):       ";
    for (int i = 0; i < THETAres; i ++)
        cout << scientific << setw(18) << Rtheta.at(i)/PHI;
    cout<<endl;
    cout<<"Rs(theta):      ";
    for (int i = 0; i < THETAres; i ++)
        cout << scientific << setw(18) << Rstheta.at(i)/PHI;
    cout<<endl;
    cout<<"T(theta):       ";
    for (int i = 0; i < THETAres; i ++)
        cout << scientific << setw(18) << Ttheta.at(i)/PHI;
    cout<<endl;
    
    //Print fluorescence part
    if (FGEN > CONST_EPS) {
        cout<<"Ffront(theta):  ";
        for (int i = 0; i < THETAres; i ++)
            cout << scientific << setw(18) << Fftheta.at(i)/PHI;
        cout<<endl;
        cout<<"Fback(theta):   ";
        for (int i = 0; i < THETAres; i ++)
            cout << scientific << setw(18) << Fbtheta.at(i)/PHI;
        cout<<endl;
    }
    
    //Deposition energy moments (time)
    if (momentlvl > 0) {
        cout <<endl;
        cout << "Incident scatter";
        if (ADEP > CONST_EPS)
            cout << ", Incident absorption";
        if (FGEN > CONST_EPS)
            cout << ", Fluorescence scatter, Fluorescence reabsorption";
        cout << " moments (time) [ps**n]:" << endl;
        cout << "  E[T]     = " << STMOMENT[0]/ASCAT;
        if (ADEP > CONST_EPS)
            cout << "     " << TMOMENT[0]/ADEP;
        if (FGEN > CONST_EPS)
            cout << "     " << FSTMOMENT[0]/FSCAT << "     " << FTMOMENT[0]/FDEP;
            
    } if (momentlvl > 1) {
        cout << endl << "  E[T2]    = " << STMOMENT[1]/ASCAT;
        if (ADEP > CONST_EPS)
            cout << "     " << TMOMENT[1]/ADEP;
        if (FGEN > CONST_EPS)
            cout << "     " << FSTMOMENT[1]/FSCAT << "     " << FTMOMENT[1]/FDEP;
        
    } if (momentlvl > 2) {
        cout << endl << "  E[T3]    = " << STMOMENT[2]/ASCAT;
        if (ADEP > CONST_EPS)
            cout << "     " << TMOMENT[2]/ADEP;
        if (FGEN > CONST_EPS)
            cout << "     " << FSTMOMENT[2]/FSCAT << "     " << FTMOMENT[2]/FDEP;
            
    } if (momentlvl > 3) {
        cout << endl << "  E[T4]    = " << STMOMENT[3]/ASCAT;
        if (ADEP > CONST_EPS)
            cout << "     " << TMOMENT[3]/ADEP;
        if (FGEN > CONST_EPS)
            cout << "     " << FSTMOMENT[3]/FSCAT << "     " << FTMOMENT[3]/FDEP;
    }
    cout << endl;
    
    //Deposition energy moments (space)
    if (momentlvl > 0) {
    
        //Header
        cout <<endl;
        cout << "Incident scatter";
        if (ADEP > CONST_EPS)
            cout << ", Incident absorption";
        if (FGEN > CONST_EPS)
            cout << ", Fluorescence scatter, Fluorescence reabsorption";
        cout << " moments (space) [um**n]:" << endl;
        
        //Moment Z
        cout << "  E[Z]    = " << SMOMENTS[0]/ASCAT;
        if (ADEP > CONST_EPS)
            cout << "     " << MOMENTS[0]/ADEP;
        if (FGEN > CONST_EPS)
            cout << "     " << FSMOMENTS[0]/FSCAT << "     " << FMOMENTS[0]/FDEP;
        
        //Moment R
        cout << endl << "  E[R]    = " << SMOMENTS[1]/ASCAT;
        if (ADEP > CONST_EPS)
            cout << "     " << MOMENTS[1]/ADEP;
        if (FGEN > CONST_EPS)
            cout << "     " << FSMOMENTS[1]/FSCAT << "     " << FMOMENTS[1]/FDEP;
            
    } if (momentlvl > 1) {
        //Moment Z*Z
        cout << endl << "  E[Z2]   = " << SMOMENTS[2]/ASCAT;
        if (ADEP > CONST_EPS)
            cout << "     " << MOMENTS[2]/ADEP;
        if (FGEN > CONST_EPS)
            cout << "     " << FSMOMENTS[2]/FSCAT << "     " << FMOMENTS[2]/FDEP;
        
        //Moment R*Z
        cout << endl << "  E[RZ]   = " << SMOMENTS[3]/ASCAT;
        if (ADEP > CONST_EPS)
            cout << "     " << MOMENTS[3]/ADEP;
        if (FGEN > CONST_EPS)
            cout << "     " << FSMOMENTS[3]/FSCAT << "     " << FMOMENTS[3]/FDEP;
        
        //Moment R*R
        cout << endl << "  E[R2]   = " << SMOMENTS[4]/ASCAT;
        if (ADEP > CONST_EPS)
            cout << "     " << MOMENTS[4]/ADEP;
        if (FGEN > CONST_EPS)
            cout << "     " << FSMOMENTS[4]/FSCAT << "     " << FMOMENTS[4]/FDEP;
        
    } if (momentlvl > 2) {
        //Moment Z*Z*Z
        cout << endl << "  E[Z3]   = " << SMOMENTS[5]/ASCAT;
        if (ADEP > CONST_EPS)
            cout << "     " << MOMENTS[5]/ADEP;
        if (FGEN > CONST_EPS)
            cout << "     " << FSMOMENTS[5]/FSCAT << "     " << FMOMENTS[5]/FDEP;

        //Moment Z*Z*R
        cout << endl << "  E[Z2R]  = " << SMOMENTS[6]/ASCAT;
        if (ADEP > CONST_EPS)
            cout << "     " << MOMENTS[6]/ADEP;
        if (FGEN > CONST_EPS)
            cout << "     " << FSMOMENTS[6]/FSCAT << "     " << FMOMENTS[6]/FDEP;
        
        //Moment Z*R*R
        cout << endl << "  E[ZR2]  = " << SMOMENTS[7]/ASCAT;
        if (ADEP > CONST_EPS)
            cout << "     " << MOMENTS[7]/ADEP;
        if (FGEN > CONST_EPS)
            cout << "     " << FSMOMENTS[7]/FSCAT << "     " << FMOMENTS[7]/FDEP;
        
        //Moment R*R*R
        cout << endl << "  E[R3]   = " << SMOMENTS[8]/ASCAT;
        if (ADEP > CONST_EPS)
            cout << "     " << MOMENTS[8]/ADEP;
        if (FGEN > CONST_EPS)
            cout << "     " << FSMOMENTS[8]/FSCAT << "     " << FMOMENTS[8]/FDEP;
            
    } if (momentlvl > 3) {
        //Moment Z*Z*Z*Z
        cout << endl << "  E[Z4]   = " << SMOMENTS[9]/ASCAT;
        if (ADEP > CONST_EPS)
            cout << "     " << MOMENTS[9]/ADEP;
        if (FGEN > CONST_EPS)
            cout << "     " << FSMOMENTS[9]/FSCAT << "     " << FMOMENTS[9]/FDEP;
        
        //Moment Z*Z*Z*R
        cout << endl << "  E[Z3R]  = " << SMOMENTS[10]/ASCAT;
        if (ADEP > CONST_EPS)
            cout << "     " << MOMENTS[10]/ADEP;
        if (FGEN > CONST_EPS)
            cout << "     " << FSMOMENTS[10]/FSCAT << "     " << FMOMENTS[10]/FDEP;
        
        //Moment Z*Z*R*R
        cout << endl << "  E[Z2R2] = " << SMOMENTS[11]/ASCAT;
        if (ADEP > CONST_EPS)
            cout << "     " << MOMENTS[11]/ADEP;
        if (FGEN > CONST_EPS)
            cout << "     " << FSMOMENTS[11]/FSCAT << "     " << FMOMENTS[11]/FDEP;
        
        //Moment Z*R*R
        cout << endl << "  E[ZR3]  = " << SMOMENTS[12]/ASCAT;
        if (ADEP > CONST_EPS)
            cout << "     " << MOMENTS[12]/ADEP;
        if (FGEN > CONST_EPS)
            cout << "     " << FSMOMENTS[12]/FSCAT << "     " << FMOMENTS[12]/FDEP;
        
        //Moment R*R*R*R
        cout << endl << "  E[R4]   = " << SMOMENTS[13]/ASCAT;
        if (ADEP > CONST_EPS)
            cout << "     " << MOMENTS[13]/ADEP;
        if (FGEN > CONST_EPS)
            cout << "     " << FSMOMENTS[13]/FSCAT << "     " << FMOMENTS[13]/FDEP;
            
    }
    
    //Transmission angle moments
    if (momentlvl > 0) {
        cout << endl << endl;
        cout << "Transmission";
        if (FGEN > CONST_EPS)
            cout << "/Emission";
        cout << " angle moments (R/T";
        if (FGEN > CONST_EPS)
            cout << "/Ff/Fb";
        cout << " [-], E[THETA**n] [rad**n]):";
        cout << endl << "  R     = " << fixed << Rdiffuse/PHI << ",";
        for (int i = 0; i < 4; i ++)
            cout << " " << fixed << setw(10) << THMOMENT[i]/Rdiffuse;
        cout << endl << "  T     = " << fixed << Tdiffuse/PHI << ",";
        for (int i = 4; i < 8; i ++)
            cout << " " << fixed << setw(10) << THMOMENT[i]/Tdiffuse;
        if (FGEN > CONST_EPS) {
            cout << endl << "  Ff    = " << fixed << Ff/PHI << ",";
            for (int i = 0; i < 4; i ++)
                cout << " " << fixed << setw(10) << FTHMOMENT[i]/Ff;
            cout << endl << "  Fb    = " << fixed << Fb/PHI << ",";
            for (int i = 4; i < 8; i ++)
                cout << " " << fixed << setw(10) << FTHMOMENT[i]/Fb;
        
        }
    }
    cout<<endl;
    
    //Print time dependence
    printTime();
    
    //Reflection coefficients (total)
    cout << endl << endl;
    cout << "Reflection/Transmission/Absorption";
    if (FGEN > CONST_EPS)
        cout << "/Emission";
    cout << " coefficients [photon/source photon]:" << endl;
    cout << "  Rdiffuse   = " << fixed << setw(10) << Rdiffuse/PHI;
    if (FGEN > CONST_EPS)
        cout << "     " << "  Fgen       = " << fixed << setw(10) << FGEN/PHI;
    cout << endl << "  Rspecular  = " << fixed << setw(10) << Rspec/PHI;
    if (FGEN > CONST_EPS)
        cout << "     " << "  Ffront     = " << fixed << setw(10) << Ff/PHI;
    cout << endl << "  Tdiffuse   = " << fixed << setw(10) << Tdiffuse/PHI;
    if (FGEN > CONST_EPS)
        cout << "     " << "  Fback      = " << fixed << setw(10) << Fb/PHI;
    cout << endl << "  Tballistic = " << fixed << setw(10) << Tballistic/PHI;
    if (FGEN > CONST_EPS)
        cout << "     " << "  Fside      = " << fixed << setw(10) << Fl/PHI;
    cout << endl << "  Ldiffuse   = " << fixed << setw(10) << Ldiffuse/PHI;
    if (FGEN > CONST_EPS)
        cout << "     " << "  Fabs       = " << fixed << setw(10) << FDEP/PHI;
    cout << endl << "  A          = " << fixed << setw(10) << ADEP/PHI << endl << endl;
    cout << "1 - (Rd + Rs + Td + A + Ld)   ?=  0  =  " << 1-(Rdiffuse + Rspec + Tdiffuse + ADEP + Ldiffuse) / PHI << endl;
    if (FGEN > CONST_EPS)
        cout << "Fg - (Ff + Fb + Fl + Fa)      ?=  0  =  " << (FGEN - Ff - Fb - Fl - FDEP) / PHI << endl;
    if (FGEN > CONST_EPS)
        cout << "Fg / Fabs   =  " << FGEN / ADEP << endl;
    
    //Reflection time/order
    cout << endl << endl;
    cout << "Diffuse Reflection/Transmission scattering order (mean, std) [-] and duration (mean, std) [ps]:" << endl;
    if (Rdiffuse > CONST_EPS)
        cout << "  N[Rdiffuse]      = " << fixed << setw(10) << Nr / Rdiffuse << endl;
    if (Tdiffuse > CONST_EPS)
        cout << "  N[Tdiffuse]      = " << fixed << setw(10) << Nt / Tdiffuse << endl;
    if (ADEP > CONST_EPS)
        cout << "  N[A]             = " << fixed << setw(10) << Na / ADEP << endl;
    if (Rdiffuse > CONST_EPS)
        cout << "  dN[Rdiffuse]     = " << fixed << setw(10) << sqrt(N2r*Rdiffuse - Nr*Nr) / Rdiffuse << endl;
    if (Tdiffuse > CONST_EPS)
        cout << "  dN[Tdiffuse]     = " << fixed << setw(10) << sqrt(N2t*Tdiffuse - Nt*Nt) / Tdiffuse << endl;
    if (ADEP > CONST_EPS)
        cout << "  dN[A]            = " << fixed << setw(10) << sqrt(N2a*ADEP - Na*Na) / ADEP << endl;
    if (Rdiffuse > CONST_EPS)
        cout << "  t[Rdiffuse]      = " << fixed << setw(10) << tRd / Rdiffuse << endl;
    if (Tdiffuse > CONST_EPS)
        cout << "  t[Tdiffuse]      = " << fixed << setw(10) << tTd / Tdiffuse << endl;
    if (ADEP > CONST_EPS)
        cout << "  t[A]             = " << fixed << setw(10) << ta / ADEP << endl;
    if (Rdiffuse > CONST_EPS)
        cout << "  dt[Rdiffuse]     = " << fixed << setw(10) << sqrt(t2Rd*Rdiffuse - tRd*tRd) / Rdiffuse << endl;
    if (Tdiffuse > CONST_EPS)
        cout << "  dt[Tdiffuse]     = " << fixed << setw(10) << sqrt(t2Td*Tdiffuse - tTd*tTd) / Tdiffuse << endl;
    if (ADEP > CONST_EPS)
        cout << "  dt[A]            = " << fixed << setw(10) << sqrt(t2a*ADEP - ta*ta) / ADEP << endl;

    //Diffusion coefficient
    cout << endl << endl;
    cout << "Diffusion coefficients [um2/ps]:" << endl;
    cout << "  Dx = " << DX/ND/2.0 << endl;
    cout << "  Dy = " << DY/ND/2.0 << endl;
    cout << "  Dz = " << DZ/ND/2.0 << endl;
}

void Stats::writedbheader(ofstream& OF) const {
    OF << "Trans,Ref0,RefD,TransS,TransB,Abs,FGen,FAbs,Ff,Fb,Fs,Ffb,Fbb,Fsb,Pen,RadL,Rad0";
}

void Stats::writedb(ofstream& OF) const {
    OF << scientific << setprecision(8) << Tdiffuse/PHI << ",";
    OF << scientific << setprecision(8) << Rspec/PHI << ",";
    OF << scientific << setprecision(8) << Rdiffuse/PHI << ",";
    OF << scientific << setprecision(8) << Ldiffuse/PHI << ",";
    OF << scientific << setprecision(8) << Tballistic/PHI << ",";
    OF << scientific << setprecision(8) << ADEP/PHI << ",";
    OF << scientific << setprecision(8) << FGEN/PHI << ",";
    OF << scientific << setprecision(8) << FDEP/PHI << ",";
    OF << scientific << setprecision(8) << Ff/PHI << ",";
    OF << scientific << setprecision(8) << Fb/PHI << ",";
    OF << scientific << setprecision(8) << Fl/PHI << ",";
    OF << scientific << setprecision(8) << Ffballistic/PHI << ",";
    OF << scientific << setprecision(8) << Fbballistic/PHI << ",";
    OF << scientific << setprecision(8) << Flballistic/PHI << ",";
    OF << scientific << setprecision(8) << SMOMENTS[0]/ASCAT << ",";
    OF << scientific << setprecision(8) << sqrt(RAD2/WRAD) << ",";
    OF << scientific << setprecision(8) << sqrt(RAD2F/WRADF);
}

#endif
