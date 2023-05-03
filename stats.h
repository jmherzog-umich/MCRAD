#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>

#include "vec.h"

#ifndef _STATS_H_
#define _STATS_H_

#define CONST_PI 3.1415926535

using namespace std;

struct Stats {
    //Settings
    int N0, Zres, Rres, Tres, THETAres, momentlvl;
    double dR, dZ, dT, dTHETA, depLifetime;
    double Rdiffuse, Tdiffuse;
    
    //Output arrays for energy deposition, scatter, etc.
    vector<double> DEP, SCAT, PENDEPTH, RADIUS, Rtheta, Ttheta;
    vector<bool> NONZERO;        //true if a particle interacted here
    
    //Initialize moments and other outputs
    vector<double> MOMENTS, SMOMENTS, THMOMENT, TMOMENT, STMOMENT;
    
    //Diffusion coefficients
    double DX, DY, DZ;
    int ND;
    
    //Energy deposited cache
    double _ADEP;
    
    //Functions
    Stats();
    Stats(int N0, int Z, int R, int T, int TH, int LVL);
    void print();
    void scatter(const vec& x, double t, double w0, double wf);
    void terminate(const vec& x, double t);
    double reflect(const vec& mu, double W, double m, bool BACK);
};

Stats::Stats() {
    Stats(1e6, 100, 100, 1, 50, 4);
}

Stats::Stats(int N0, int Z=100, int R=100, int T=1, int TH=50, int LVL=4) {
    //Set default resolution
    this->N0 = N0;
    Zres = Z;              //Number of samples to grid in Z direction
    Rres = R;              //Number of radial samples to grid
    Tres = T;              //Number of time samples
    THETAres = TH;         //Number of cos(theta) values for reflection
    momentlvl = LVL;       //Order of moments to calculate (only up to 4th order is implemented)
    
    //Grid resolution
    dR = 1.0/Rres;
    dZ = 1.0/Zres;
    dT = 1.0/Tres;
    dTHETA = CONST_PI/THETAres/2;
    depLifetime = 0;
    _ADEP = 0;
    
    //Reflection and transmission
    Rdiffuse = 0;
    Tdiffuse = 0;
    
    //Diffusion coefficients
    DX = DY = DZ = 0;
    ND = 0;
    
    //Output arrays for energy deposition, scatter, etc.
    DEP = vector<double>(Zres*Rres, 0.0);
    SCAT = vector<double>(Zres*Rres, 0.0);
    NONZERO = vector<bool>(Zres*Rres, false);
    
    //Penetration depth and beam radius
    PENDEPTH = vector<double>(Zres,0.0);
    RADIUS = vector<double>(Zres,0.0);
    
    //Reflection and transmission
    Rtheta = vector<double>(THETAres, 0.0);
    Ttheta = vector<double>(THETAres, 0.0);
    
    //Initialize moments and other outputs
    MOMENTS = vector<double>(14, 0.0);
    SMOMENTS = vector<double>(14, 0.0);
    THMOMENT = vector<double>(8, 0.0);
    TMOMENT = vector<double>(4, 0.0);
    STMOMENT = vector<double>(4, 0.0);
    
    //Other output of initialized variables
    cout << "Grid resolution: " << Rres << " x " << Zres << endl;
}

void Stats::scatter(const vec& x, double t, double w0, double wf) {
    //Initialize
    int binR, binZ;

    //Calculate bins
    binR = (x.r()/dR >= Rres) ? Rres-1 : (int)(x.r()/dR);
    binZ = (x.Z/dZ >= Zres) ? Zres-1 : (int)(x.Z/dZ);
    
    //Basic data
    PENDEPTH.at(binZ) = PENDEPTH.at(binZ) + w0;
    RADIUS.at(binZ) = RADIUS.at(binZ) + w0 * pow(x.r(), 2);
    DEP.at(binZ*Rres + binR) = DEP.at(binZ*Rres + binR) + wf;
    SCAT.at(binZ*Rres + binR) = SCAT.at(binZ*Rres + binR) + w0;
    NONZERO.at(binZ*Rres + binR) = true;
    
    //Calculate moments
    if (momentlvl > 0) {
        MOMENTS.at(0) += wf*x.Z;
        MOMENTS.at(1) += wf*x.r();
        TMOMENT.at(0) += wf*t;
    } if (momentlvl > 1) {
        MOMENTS.at(2) += wf*pow(x.Z,2);
        MOMENTS.at(3) += wf*x.r()*x.Z;
        MOMENTS.at(4) += wf*pow(x.r(),2);
        TMOMENT.at(1) += wf*pow(t,2);
    } if (momentlvl > 2) {
        MOMENTS.at(5) += wf*pow(x.Z,3);
        MOMENTS.at(6) += wf*x.r()*pow(x.Z,2);
        MOMENTS.at(7) += wf*x.Z*pow(x.r(),2);
        MOMENTS.at(8) += wf*pow(x.r(),3);
        TMOMENT.at(2) += wf*pow(t,3);
    } if (momentlvl > 3) {
        MOMENTS.at(9) += wf*pow(x.Z,4);
        MOMENTS.at(10) += wf*x.r()*pow(x.Z,3);
        MOMENTS.at(11) += wf*pow(x.Z*x.r(),2);
        MOMENTS.at(12) += wf*x.Z*pow(x.r(),3);
        MOMENTS.at(13) += wf*pow(x.r(),4);
        TMOMENT.at(3) += wf*pow(t,4);
    }
    
    //Calculate scattering moments
    if (momentlvl > 0) {
        SMOMENTS.at(0) += w0*x.Z;
        SMOMENTS.at(1) += w0*x.r();
        STMOMENT.at(0) += w0*t;
    } if (momentlvl > 1) {
        SMOMENTS.at(2) += w0*pow(x.Z,2);
        SMOMENTS.at(3) += w0*x.r()*x.Z;
        SMOMENTS.at(4) += w0*pow(x.r(),2);
        STMOMENT.at(1) += w0*pow(t,2);
    } if (momentlvl > 2) {
        SMOMENTS.at(5) += w0*pow(x.Z,3);
        SMOMENTS.at(6) += w0*x.r()*pow(x.Z,2);
        SMOMENTS.at(7) += w0*x.Z*pow(x.r(),2);
        SMOMENTS.at(8) += w0*pow(x.r(),3);
        STMOMENT.at(2) += w0*pow(t,3);
    } if (momentlvl > 3) {
        SMOMENTS.at(9) += w0*pow(x.Z,4);
        SMOMENTS.at(10) += w0*x.r()*pow(x.Z,3);
        SMOMENTS.at(11) += w0*pow(x.Z*x.r(),2);
        SMOMENTS.at(12) += w0*x.Z*pow(x.r(),3);
        SMOMENTS.at(13) += w0*pow(x.r(),4);
        STMOMENT.at(3) += w0*pow(t,4);
    }
}

void Stats::terminate(const vec& x, double t) {
    DX += pow(x.X,2)/t; //Add to average diffusion coefficient
    DY += pow(x.Y,2)/t; //Y component
    DZ += pow(x.Z,2)/t; //Z component
    depLifetime += t;
    ND ++;
}

double Stats::reflect(const vec& mu, double W, double m, bool BACK) {
    //Calculate theta from Snell's law
    double R, cost, sint;
    double R0 = (1.0-m)*(1.0-m)/(1.0+m)/(1.0+m);
    int binTHETA;
    if (1-abs(mu.Z) < 1e-6) {
        R = R0;
    } else {
        cost = acos(abs(mu.Z)); //theta_i
        if (1-abs(m*sin(cost)) < 1e-6) {
            R = 1;
            sint = 0;
        } else {
            sint = asin(m*sin(cost));       //theta_t
            R = 0.5*(pow(sin(cost-sint)/sin(cost+sint),2) 
                    + pow(tan(cost-sint)/tan(cost+sint),2));
        }
    }
    
    //Bin in exit angle
    binTHETA = (sint/dTHETA >= THETAres) ? THETAres-1 : (int) (sint/dTHETA);
    
    //Increment transmission
    if (BACK) {
        Tdiffuse += (1-R) * W;
        Ttheta.at(binTHETA) += W * (1-R);
    }else{
        Rdiffuse += (1-R) * W;
        Rtheta.at(binTHETA) += W * (1-R);
    }
    
    //Get moments
    if (momentlvl > 3)
        THMOMENT.at(0+(BACK ? 4 : 0)) += W * (1-R) * sint;
    if (momentlvl > 1)
        THMOMENT.at(1+(BACK ? 4 : 0)) += W * (1-R) * pow(sint,2);
    if (momentlvl > 2)
        THMOMENT.at(2+(BACK ? 4 : 0)) += W * (1-R) * pow(sint,3);
    if (momentlvl > 3)
        THMOMENT.at(3+(BACK ? 4 : 0)) += W * (1-R) * pow(sint,4);
        
    return R;
}

void Stats::print() {
    _ADEP = 0;
    cout<<"ENERGY DEPOSITED"<<endl;
    for (int j = 0; j < Zres; j ++) {
        for (int k = 0; k < Rres; k ++) {
            cout << scientific << setw(14) << DEP.at(j*Rres + k)/N0;
            _ADEP += DEP.at(j*Rres + k);
        }
        cout<<endl;
    }
    cout<<endl;
    
    //energy scattered (in mJ/mJ-incident)
    double ASCAT = 0;
    cout<<"ENERGY SCATTERED"<<endl;
    for (int j = 0; j < Zres; j ++) {
        for (int k = 0; k < Rres; k ++) {
            cout << scientific << setw(14) << SCAT.at(j*Rres + k)/N0;
            ASCAT += SCAT.at(j*Rres + k);
        }
        cout<<endl;
    }
    cout<<endl;
    
    //Energy scattered from depth Z
    cout<<"ENERGY SCATTERED FROM DEPTH Z"<<endl;
    for (int j = 0; j < Zres; j ++) {
        cout << scientific << setw(14) << PENDEPTH.at(j)/N0;
    }
    cout<<endl<<endl;
    
    //Beam radius at depth Z
    cout<<"BEAM RADIUS (STANDARD DEVIATION) AT DEPTH Z"<<endl;
    for (int j = 0; j < Zres; j ++) {
        RADIUS.at(j) = sqrt(RADIUS.at(j)/PENDEPTH.at(j));
        cout << scientific << setw(14) << RADIUS.at(j); 
    }
    cout<<endl;
    
    //Theta-dependent reflection/transmission/fluorescence
    cout<<endl<<endl;
    cout<<"cos(theta):     ";
    for (int i = 0; i < THETAres; i ++)
        cout << fixed << setw(14) << i*dTHETA;
    cout<<endl;
    cout<<"R(theta):       ";
    for (int i = 0; i < THETAres; i ++)
        cout << scientific << setw(14) << Rtheta.at(i)/N0;
    cout<<endl;
    cout<<"T(theta):       ";
    for (int i = 0; i < THETAres; i ++)
        cout << scientific << setw(14) << Ttheta.at(i)/N0;
    cout<<endl;
    
    //Deposition energy moments (time)
    if (momentlvl > 0) {
        cout <<endl <<endl;
        cout << "Energy deposition moments (time) [ns**n]:" << endl;
        cout << "  E[T]     = " << TMOMENT[0]/_ADEP << endl;
    } if (momentlvl > 1)
        cout << "  E[T2]    = " << TMOMENT[1]/_ADEP << endl;
    if (momentlvl > 2)
        cout << "  E[T3]    = " << TMOMENT[2]/_ADEP << endl;
    if (momentlvl > 3)
        cout << "  E[T4]    = " << TMOMENT[3]/_ADEP << endl;
    
    //Deposition energy moments (space)
    if (momentlvl > 0) {
        cout <<endl <<endl;
        cout << "Energy deposition moments (space) [nm**n]:" << endl;
        cout << "  E[Z]    = " << MOMENTS[0]/_ADEP << endl;
        cout << "  E[R]    = " << MOMENTS[1]/_ADEP << endl;
    } if (momentlvl > 1) {
        cout << "  E[Z2]   = " << MOMENTS[2]/_ADEP << endl;
        cout << "  E[RZ]   = " << MOMENTS[3]/_ADEP << endl;
        cout << "  E[R2]   = " << MOMENTS[4]/_ADEP << endl;
    } if (momentlvl > 2) {
        cout << "  E[Z3]   = " << MOMENTS[5]/_ADEP << endl;
        cout << "  E[Z2R]  = " << MOMENTS[6]/_ADEP << endl;
        cout << "  E[ZR2]  = " << MOMENTS[7]/_ADEP << endl;
        cout << "  E[R3]   = " << MOMENTS[8]/_ADEP << endl;
    } if (momentlvl > 3) {
        cout << "  E[Z4]   = " << MOMENTS[9]/_ADEP << endl;
        cout << "  E[Z3R]  = " << MOMENTS[10]/_ADEP << endl;
        cout << "  E[Z2R2] = " << MOMENTS[11]/_ADEP << endl;
        cout << "  E[ZR3]  = " << MOMENTS[12]/_ADEP << endl;
        cout << "  E[R4]   = " << MOMENTS[13]/_ADEP << endl;
    }
    
    //Scattered energy moments (time)
    if (momentlvl > 0) {
        cout <<endl <<endl;
        cout << "Energy scatter moments (time) [ns**n]:" << endl;
        cout << "  E[T]     = " << STMOMENT[0]/ASCAT << endl;
    } if (momentlvl > 1)
        cout << "  E[T2]    = " << STMOMENT[1]/ASCAT << endl;
    if (momentlvl > 2)
        cout << "  E[T3]    = " << STMOMENT[2]/ASCAT << endl;
    if (momentlvl > 3)
        cout << "  E[T4]    = " << STMOMENT[3]/ASCAT << endl;
    
    //Scattered energy moments (space)
    if (momentlvl > 0) {
        cout <<endl <<endl;
        cout << "Energy scatter moments (space) [nm**n]:" << endl;
        cout << "  E[Z]    = " << SMOMENTS[0]/ASCAT << endl;
        cout << "  E[R]    = " << SMOMENTS[1]/ASCAT << endl;
    } if (momentlvl > 1) {
        cout << "  E[Z2]   = " << SMOMENTS[2]/ASCAT << endl;
        cout << "  E[RZ]   = " << SMOMENTS[3]/ASCAT << endl;
        cout << "  E[R2]   = " << SMOMENTS[4]/ASCAT << endl;
    } if (momentlvl > 2) {
        cout << "  E[Z3]   = " << SMOMENTS[5]/ASCAT << endl;
        cout << "  E[Z2R]  = " << SMOMENTS[6]/ASCAT << endl;
        cout << "  E[ZR2]  = " << SMOMENTS[7]/ASCAT << endl;
        cout << "  E[R3]   = " << SMOMENTS[8]/ASCAT << endl;
    } if (momentlvl > 3) {
        cout << "  E[Z4]   = " << SMOMENTS[9]/ASCAT << endl;
        cout << "  E[Z3R]  = " << SMOMENTS[10]/ASCAT << endl;
        cout << "  E[Z2R2] = " << SMOMENTS[11]/ASCAT << endl;
        cout << "  E[ZR3]  = " << SMOMENTS[12]/ASCAT << endl;
        cout << "  E[R4]   = " << SMOMENTS[13]/ASCAT << endl;
    }
    
    //Transmission angle moments
    if (momentlvl > 0) {
        cout << endl << endl;
        cout << "Transmission angle moments (R/T [-], E[THETA**n] [rad**n]):";
        cout << endl << "  R     = " << Rdiffuse/N0 << ",";
        for (int i = 0; i < 4; i ++)
            cout << " " << fixed << setw(10) << THMOMENT[i]/Rdiffuse;
    } if (momentlvl > 1) {
        cout << endl << "  T     = " << Tdiffuse/N0 << ",";
        for (int i = 4; i < 8; i ++)
            cout << " " << fixed << setw(10) << THMOMENT[i]/Tdiffuse;
    }
    
    //Reflection coefficients (total)
    cout << endl << endl;
    cout << "Reflection/Transmission/Absorption coefficients [-] (relative [-]):" << endl;
    cout << "  Rdiffuse = " << fixed << setw(10) << Rdiffuse/N0 << endl;
    cout << "  Tdiffuse = " << fixed << setw(10) << Tdiffuse/N0 << endl;
    cout << "  A        = " << fixed << setw(10) << _ADEP/N0 << endl;

    //Diffusion coefficient
    cout << endl << endl;
    cout << "Diffusion coefficients [nm2/ns]:" << endl;
    cout << "  Dx = " << DX/ND/2 << endl;
    cout << "  Dy = " << DY/ND/2 << endl;
    cout << "  Dz = " << DZ/ND/2 << endl;
}
    
#endif
