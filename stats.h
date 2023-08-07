#include <vector>
#include <cmath>
#include <iostream>
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
    double dT, dTHETA, T, PHI;
    double Rdiffuse, Tdiffuse, Tballistic, Rspec, Ldiffuse;
    
    //Output arrays for energy deposition, scatter, etc.
    vector<double> Rtheta, Rstheta, Ttheta;
    
    //Initialize moments and other outputs
    vector<double> MOMENTS, SMOMENTS, THMOMENT, TMOMENT, STMOMENT;
    
    //Time-dependence
    vector<double> Tt, Rt, Lt, At;
    
    //Diffusion coefficients
    double DX, DY, DZ, ND;
    
    //Transmission/reflection scattering order
    double Nr, Nt, N2r, N2t;
    double Na, ta, N2a, t2a;
    
    //Transmission/reflection times
    double tRd, tTd, t2Rd, t2Td;
    
    //Cache
    double ADEP, ASCAT;
    
    //Functions
    //--Constructors
    Stats();
    Stats(int N, double dT, int LVL=4);
    Stats(int nT, int nTH, double dT, int LVL=4);
    
    //--Accessors
    int getBinT(double t) const;
    int getBinTh(double sint) const;
    
    //--Print results functions
    void setup();
    void print() const;
    void printTime() const;
    
    //--Functions to aggregate data
    void scatter(const Photon& p, double w);
    void absorb(const Photon& p, double w);
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

Stats::Stats(int nT, int nTH, double dT, int LVL) {
    //Set default resolution
    Tres = nT;             //Number of time samples
    THETAres = nTH;        //Number of cos(theta) values for reflection
    momentlvl = LVL;       //Order of moments to calculate (only up to 4th order is implemented)
    this->dT = dT;
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
    
    //Grid resolution
    T = dT * Tres;
    dTHETA = CONST_PI/THETAres/2;
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
    
    //Initialize moments and other outputs
    MOMENTS = vector<double>(14, 0.0);
    THMOMENT = vector<double>(8, 0.0);
    TMOMENT = vector<double>(4, 0.0);
    SMOMENTS = vector<double>(14, 0.0);
    STMOMENT = vector<double>(4, 0.0);
}

int Stats::getBinT(double t) const {
    int binT = (int)(t/dT);
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

void Stats::absorb(const Photon& p, double W) {
    //Initialize
    double t = p.t;
    const vec x = p.x;
    
    //Calculate moments
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
    int binT = getBinT(t);
    if (binT >= 0)
        At.at(binT) += W;
}

void Stats::scatter(const Photon& p, double W) {
    //Initialize
    double t = p.t;
    const vec x = p.x;

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
    binT = getBinT(p.t);
    
    //Get moments
    if (p.n > 0) {
        if (momentlvl > 3)
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
        if (p.isBallistic)
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
    cout << endl << "Temporal profiles [photons/time step]:" << endl;
    cout << "  T [ps]: ";
    for (int i = 0; i < Tres; i ++)
        cout << scientific << setw(18) << dT * i;
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
    
    //Deposition energy moments (time)
    if (momentlvl > 0) {
        cout <<endl;
        cout << "Energy scatter, absorption moments (time) [ps**n]:" << endl;
        cout << "  E[T]     = " << STMOMENT[0]/ASCAT << "     " << TMOMENT[0]/ADEP << endl;
    } if (momentlvl > 1)
        cout << "  E[T2]    = " << STMOMENT[1]/ASCAT << "     " << TMOMENT[1]/ADEP << endl;
    if (momentlvl > 2)
        cout << "  E[T3]    = " << STMOMENT[2]/ASCAT << "     " << TMOMENT[2]/ADEP << endl;
    if (momentlvl > 3)
        cout << "  E[T4]    = " << STMOMENT[3]/ASCAT << "     " << TMOMENT[3]/ADEP << endl;
    
    //Deposition energy moments (space)
    if (momentlvl > 0) {
        cout <<endl;
        cout << "Energy scatter, absorption moments (space) [um**n]:" << endl;
        cout << "  E[Z]    = " << SMOMENTS[0]/ASCAT << "     " << MOMENTS[0]/ADEP << endl;
        cout << "  E[R]    = " << SMOMENTS[1]/ASCAT << "     " << MOMENTS[1]/ADEP << endl;
    } if (momentlvl > 1) {
        cout << "  E[Z2]   = " << SMOMENTS[2]/ASCAT << "     " << MOMENTS[2]/ADEP << endl;
        cout << "  E[RZ]   = " << SMOMENTS[3]/ASCAT << "     " << MOMENTS[3]/ADEP << endl;
        cout << "  E[R2]   = " << SMOMENTS[4]/ASCAT << "     " << MOMENTS[4]/ADEP << endl;
    } if (momentlvl > 2) {
        cout << "  E[Z3]   = " << SMOMENTS[5]/ASCAT << "     " << MOMENTS[5]/ADEP << endl;
        cout << "  E[Z2R]  = " << SMOMENTS[6]/ASCAT << "     " << MOMENTS[6]/ADEP << endl;
        cout << "  E[ZR2]  = " << SMOMENTS[7]/ASCAT << "     " << MOMENTS[7]/ADEP << endl;
        cout << "  E[R3]   = " << SMOMENTS[8]/ASCAT << "     " << MOMENTS[8]/ADEP << endl;
    } if (momentlvl > 3) {
        cout << "  E[Z4]   = " << SMOMENTS[9]/ASCAT << "     " << MOMENTS[9]/ADEP << endl;
        cout << "  E[Z3R]  = " << SMOMENTS[10]/ASCAT << "     " << MOMENTS[10]/ADEP << endl;
        cout << "  E[Z2R2] = " << SMOMENTS[11]/ASCAT << "     " << MOMENTS[11]/ADEP << endl;
        cout << "  E[ZR3]  = " << SMOMENTS[12]/ASCAT << "     " << MOMENTS[12]/ADEP << endl;
        cout << "  E[R4]   = " << SMOMENTS[13]/ASCAT << "     " << MOMENTS[13]/ADEP << endl;
    }
    
    //Transmission angle moments
    if (momentlvl > 0) {
        cout << endl << endl;
        cout << "Transmission angle moments (R/T [-], E[THETA**n] [rad**n]):";
        cout << endl << "  R     = " << Rdiffuse/PHI << ",";
        for (int i = 0; i < 4; i ++)
            cout << " " << fixed << setw(10) << THMOMENT[i]/Rdiffuse;
    } if (momentlvl > 1) {
        cout << endl << "  T     = " << Tdiffuse/PHI << ",";
        for (int i = 4; i < 8; i ++)
            cout << " " << fixed << setw(10) << THMOMENT[i]/Tdiffuse;
    }
    cout<<endl;
    
    //Print time dependence
    printTime();
    
    //Reflection coefficients (total)
    cout << endl << endl;
    cout << "Reflection/Transmission/Absorption coefficients [-]:" << endl;
    cout << "  Rdiffuse   = " << fixed << setw(10) << Rdiffuse/PHI << endl;
    cout << "  Rspecular  = " << fixed << setw(10) << Rspec/PHI << endl;
    cout << "  Tdiffuse   = " << fixed << setw(10) << Tdiffuse/PHI << endl;
    cout << "  Tballistic = " << fixed << setw(10) << Tballistic/PHI << endl;
    cout << "  Ldiffuse   = " << fixed << setw(10) << Ldiffuse/PHI << endl;
    cout << "  A          = " << fixed << setw(10) << ADEP/PHI << endl << endl;
    cout << "Rd + Rs + Td + A + Ld = " << (Rdiffuse + Rspec + Tdiffuse + ADEP + Ldiffuse) / PHI << endl;
    
    //Reflection time/order
    cout << endl << endl;
    cout << "Diffuse Reflection/Transmission scattering order (mean, std) [-] and duration (mean, std) [ps]:" << endl;
    cout << "  N[Rdiffuse]      = " << fixed << setw(10) << Nr / Rdiffuse << endl;
    cout << "  N[Tdiffuse]      = " << fixed << setw(10) << Nt / Tdiffuse << endl;
    cout << "  N[A]             = " << fixed << setw(10) << Na / ADEP << endl;
    cout << "  dN[Rdiffuse]     = " << fixed << setw(10) << sqrt(N2r*Rdiffuse - Nr*Nr) / Rdiffuse << endl;
    cout << "  dN[Tdiffuse]     = " << fixed << setw(10) << sqrt(N2t*Tdiffuse - Nt*Nt) / Tdiffuse << endl;
    cout << "  dN[A]            = " << fixed << setw(10) << sqrt(N2a*ADEP - Na*Na) / ADEP << endl;
    cout << "  t[Rdiffuse]      = " << fixed << setw(10) << tRd / Rdiffuse << endl;
    cout << "  t[Tdiffuse]      = " << fixed << setw(10) << tTd / Tdiffuse << endl;
    cout << "  t[A]             = " << fixed << setw(10) << ta / ADEP << endl;
    cout << "  dt[Rdiffuse]     = " << fixed << setw(10) << sqrt(t2Rd*Rdiffuse - tRd*tRd) / Rdiffuse << endl;
    cout << "  dt[Tdiffuse]     = " << fixed << setw(10) << sqrt(t2Td*Tdiffuse - tTd*tTd) / Tdiffuse << endl;
    cout << "  dt[A]            = " << fixed << setw(10) << sqrt(t2a*ADEP - ta*ta) / ADEP << endl;

    //Diffusion coefficient
    cout << endl << endl;
    cout << "Diffusion coefficients [um2/ps]:" << endl;
    cout << "  Dx = " << DX/ND/2.0 << endl;
    cout << "  Dy = " << DY/ND/2.0 << endl;
    cout << "  Dz = " << DZ/ND/2.0 << endl;
}
    
#endif
