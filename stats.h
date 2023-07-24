#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>

#include "vec.h"
#include "photon.h"

#ifndef _STATS_H_
#define _STATS_H_

#define CONST_PI 3.1415926535
#define CONST_EPS 1e-6

using namespace std;

struct Stats {
    //Settings
    int N0, Zres, Rres, Tres, THETAres, momentlvl;
    double dR, dZ, dT, dTHETA;
    double Rdiffuse, Tdiffuse, Tballistic;
    bool saturate, interfere;
    
    //Output arrays for energy deposition, scatter, etc.
    vector<double> DEP, SCAT, PENDEPTH, RADIUS, Rtheta, Ttheta;
    vector<double> DEPOF, SCATOF;
    vector<bool> NONZERO;        //true if a particle interacted here
    
    //Input and output beam arrays
    vector<double> I0, It, Ir;
    
    //Initialize moments and other outputs
    vector<double> MOMENTS, SMOMENTS, THMOMENT, TMOMENT, STMOMENT;
    
    //Diffusion coefficients
    double DX, DY, DZ;
    int ND;
    
    //Transmission/reflection scattering order
    double Nr, Nt, N2r, N2t;
    
    //Transmission/reflection times
    double tRd, tTd, t2Rd, t2Td;
    
    //Cache
    double _ADEP, _ATOT, _ASCAT;
    
    //Functions
    //--Constructors
    Stats();
    Stats(int N0, int Z, int R, int T, int TH, int LVL);
    
    //--Accessors
    double E(const vec& x) const;
    double E(int binZ, int binR) const;
    double Volume(int binR) const;
    int getBinR(const vec& x) const;
    int getBinZ(const vec& x) const;
    
    //--Print results functions
    void finalize();
    void print() const;
    void printIC() const;
    void printProfiles() const;
    
    //--Functions to aggregate data
    void scatter(const Photon& p, const vec& norm, double T, double wscale, double a);
    void diffusion(const vec& x, double w, double t);
    void initialize(const Photon& x, const vec& norm);
    double reflect(const Photon& p, const vec& norm, double T, double wscale, double m, int surf, unsigned int N);
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
    saturate = false;      //Whether we include saturation behavior or not
    
    //Grid resolution
    dR = 1.0/Rres;
    dZ = 1.0/Zres;
    dT = 1.0/Tres;
    dTHETA = CONST_PI/THETAres/2;
    
    //Reflection and transmission
    Rdiffuse = 0;
    Tdiffuse = 0;
    Tballistic = 0;
    
    //Diffusion coefficients
    DX = DY = DZ = 0;
    ND = 0;
    
    //Scattering order and time stats
    Nr = Nt = N2r = N2t = 0;
    tRd = tTd = t2Rd = t2Td = 0;
    
    //Cache values
    _ADEP = 0;
    _ATOT = 0;
    _ASCAT = 0;
    
    //Output arrays for energy deposition, scatter, etc.
    if (saturate)
        SCAT = vector<double>(Zres*Rres, 0.0);
    DEP = vector<double>(Zres*Rres, 0.0);
    NONZERO = vector<bool>(Zres*Rres, false);
    
    //Beam array images
    I0 = vector<double>((2*Rres+1)*(2*Rres+1), 0.0);
    It = vector<double>((2*Rres+1)*(2*Rres+1), 0.0);
    Ir = vector<double>((2*Rres+1)*(2*Rres+1), 0.0);
    
    //Overflow arrays
    if (saturate)
        SCATOF = vector<double>(2*Rres+Zres+2, 0.0);
    DEPOF = vector<double>(2*Rres+Zres+2, 0.0);
    
    //Penetration depth and beam radius
    PENDEPTH = vector<double>(Zres,0.0);
    RADIUS = vector<double>(Zres,0.0);
    
    //Reflection and transmission
    Rtheta = vector<double>(THETAres, 0.0);
    Ttheta = vector<double>(THETAres, 0.0);
    
    //Initialize moments and other outputs
    MOMENTS = vector<double>(14, 0.0);
    THMOMENT = vector<double>(8, 0.0);
    TMOMENT = vector<double>(4, 0.0);
    if (saturate) {
        SMOMENTS = vector<double>(14, 0.0);
        STMOMENT = vector<double>(4, 0.0);
    }
}

int Stats::getBinR(const vec& x) const {
    int binR = (int)(x.r()/dR);
    if (binR > Rres)
        binR = Rres;
    return binR;
}

int Stats::getBinZ(const vec& x) const {
    int binZ = (int)(x.Z/dZ);
    if (binZ < -1)
        binZ = -1;
    else if (binZ > Zres)
        binZ = Zres;
    return binZ;
}

void Stats::scatter(const Photon& p, const vec& norm, double T, double wscale, double a) {
    //Initialize
    int binR, binZ;
    double w0 = p.W * wscale;
    double wf, t = p.t/T;
    const vec x = p.x / norm;

    //Calculate bins and cap at +/-1 from 0 and Res
    binR = getBinR(x);
    binZ = getBinZ(x);
    
    //If saturation is off, incident light and absorbed light are identical to
    // within a factor dependent on albedo, so use incident light here to avoid
    // roundoff from very small a values. 
    if (!saturate)
        wf = w0;
    else
        wf = w0 * a;

    //Calculate moments first, because these should be done in the absence of interference
    if (momentlvl > 0) {
        MOMENTS.at(0) += wf * x.Z;
        MOMENTS.at(1) += wf * x.r();
        TMOMENT.at(0) += wf * t;
    } if (momentlvl > 1) {
        MOMENTS.at(2) += wf * pow(x.Z,2);
        MOMENTS.at(3) += wf * x.r() * x.Z;
        MOMENTS.at(4) += wf * pow(x.r(),2);
        TMOMENT.at(1) += wf * pow(t,2);
    } if (momentlvl > 2) {
        MOMENTS.at(5) += wf * pow(x.Z,3);
        MOMENTS.at(6) += wf * x.r() * pow(x.Z,2);
        MOMENTS.at(7) += wf * x.Z * pow(x.r(),2);
        MOMENTS.at(8) += wf * pow(x.r(),3);
        TMOMENT.at(2) += wf * pow(t,3);
    } if (momentlvl > 3) {
        MOMENTS.at(9) += wf * pow(x.Z,4);
        MOMENTS.at(10) += wf * x.r() * pow(x.Z,3);
        MOMENTS.at(11) += wf * pow(x.Z * x.r(),2);
        MOMENTS.at(12) += wf * x.Z * pow(x.r(),3);
        MOMENTS.at(13) += wf * pow(x.r(),4);
        TMOMENT.at(3) += wf * pow(t,4);
    }
    
    //Calculate scattering moments
    if (saturate) {
        if (momentlvl > 0) {
            SMOMENTS.at(0) += w0 * x.Z;
            SMOMENTS.at(1) += w0 * x.r();
            STMOMENT.at(0) += w0 * t;
        } if (momentlvl > 1) {
            SMOMENTS.at(2) += w0 * pow(x.Z,2);
            SMOMENTS.at(3) += w0 * x.r() * x.Z;
            SMOMENTS.at(4) += w0 * pow(x.r(),2);
            STMOMENT.at(1) += w0 * pow(t,2);
        } if (momentlvl > 2) {
            SMOMENTS.at(5) += w0 * pow(x.Z,3);
            SMOMENTS.at(6) += w0 * x.r() * pow(x.Z,2);
            SMOMENTS.at(7) += w0 * x.Z * pow(x.r(),2);
            SMOMENTS.at(8) += w0 * pow(x.r(),3);
            STMOMENT.at(2) += w0 * pow(t,3);
        } if (momentlvl > 3) {
            SMOMENTS.at(9) += w0 * pow(x.Z,4);
            SMOMENTS.at(10) += w0 * x.r() * pow(x.Z,3);
            SMOMENTS.at(11) += w0 * pow(x.Z * x.r(),2);
            SMOMENTS.at(12) += w0 * x.Z * pow(x.r(),3);
            SMOMENTS.at(13) += w0 * pow(x.r(),4);
            STMOMENT.at(3) += w0 * pow(t,4);
        }
    }

    //Penetration depth and beam radius
    if (binZ >= 0 and binZ < Zres) {
        PENDEPTH.at(binZ) = PENDEPTH.at(binZ) + w0;
        RADIUS.at(binZ) = RADIUS.at(binZ) + w0 * pow(x.r(), 2);
    }
    
    //2D Arrays
    if (binZ < 0) {
        if (saturate)
            SCATOF.at(binR) = SCATOF.at(binR) + w0;
        DEPOF.at(binR) = DEPOF.at(binR) + wf;
    } else if (binZ >= Zres)  {
        if (saturate)
            SCATOF.at(1+Rres+Zres+(Rres-binR)) = SCATOF.at(1+Rres+Zres+(Rres-binR)) + w0;
        DEPOF.at(1+Rres+Zres+(Rres-binR)) = DEPOF.at(1+Rres+Zres+(Rres-binR)) + wf;
    } else if (binR >= Rres) {
        if (saturate)
            SCATOF.at(1+Rres+binZ) = SCATOF.at(1+Rres+binZ) + w0;
        DEPOF.at(1+Rres+binZ) = DEPOF.at(1+Rres+binZ) + wf;
    } else {
        if (saturate)
            SCAT.at(binZ*Rres + binR) = SCAT.at(binZ*Rres + binR) + w0;
        DEP.at(binZ*Rres + binR) = DEP.at(binZ*Rres + binR) + wf;
        NONZERO.at(binZ*Rres + binR) = true;
    }
    
    //Calculate diffusion properties
    diffusion(x, w0, t);
}

void Stats::diffusion(const vec& x, double w, double t) {
    DX += w * x.X*x.X / t; //Add to average diffusion coefficient
    DY += w * x.Y*x.Y / t; //Y component
    DZ += w * x.Z*x.Z / t; //Z component
    ND += w;
}

double Stats::reflect(const Photon& p, const vec& norm, double T, double wscale, double m, int surf, unsigned int N) {
    //Calculate theta from Snell's law
    double R, cost, sint, W;
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
    
    //Calculate the intensity to add
    W = p.W * wscale;
    
    //Do moments first since they depend only on the power (not amplitude)
    //Bin in exit angle
    binTHETA = (sint/dTHETA >= THETAres) ? THETAres-1 : (int) (sint/dTHETA);
    
    //Get moments
    if (momentlvl > 3)
        THMOMENT.at(0+((surf == 1) ? 4 : 0)) += W * (1-R) * sint;
    if (momentlvl > 1)
        THMOMENT.at(1+((surf == 1) ? 4 : 0)) += W * (1-R) * pow(sint,2);
    if (momentlvl > 2)
        THMOMENT.at(2+((surf == 1) ? 4 : 0)) += W * (1-R) * pow(sint,3);
    if (momentlvl > 3)
        THMOMENT.at(3+((surf == 1) ? 4 : 0)) += W * (1-R) * pow(sint,4);
        
    //Increment reflection/transmission times and orders
    if ((surf == 1)) {
        Nt += N * (1-R) * W;
        N2t += N*N * (1-R) * W;
        tTd += p.t/T * (1-R) * W;
        t2Td += p.t*p.t/T/T * (1-R) * W;
    } else if (surf == 2) {
        Nr += N * (1-R) * W;
        N2r += N*N * (1-R) * W;
        tRd += p.t/T * (1-R) * W;
        t2Rd += p.t*p.t/T/T * (1-R) * W;
    }  
    
    //Increment transmission
    if (surf == 1) {
        Tdiffuse += (1-R) * W;
        Ttheta.at(binTHETA) += W * (1-R);
        if (p.isBallistic)
            Tballistic += W * (1-R);
    } else if (surf == 2) {
        Rdiffuse += (1-R) * W;
        Rtheta.at(binTHETA) += W * (1-R);
    }
    
    //Account for interference in reflection/transmission (but not weighted quantities)
    if (interfere)
        W = sqrt(W) * p.phi();
    
    //Get x and y bins, and increment transmission/reflection beam images
    int binX = (int)(p.x.X/norm.X/dR) + Rres;
    int binY = (int)(p.x.Y/norm.Y/dR) + Rres;
    if ((binX >= 0) && (binX <= 2*Rres) && (binY >= 0) && (binY <= 2*Rres)) {
        if (surf == 1) {
            It.at(binY*(2*Rres+1)+binX) += (1-R) * W;
        } else if (surf == 2) {
            Ir.at(binY*(2*Rres+1)+binX) += (1-R) * W;
        }
    }
    
    //Return reflection coeff
    return R;
}

void Stats::initialize(const Photon& p, const vec& norm) {
    int binX = (int)(p.x.X/norm.X/dR) + Rres;
    int binY = (int)(p.x.Y/norm.Y/dR) + Rres;
    double w = p.W;
    if ((binX >= 0) && (binX <= 2*Rres) && (binY >= 0) && (binY <= 2*Rres))
        I0.at(binY*(2*Rres+1) + binX) += w;
}

void Stats::printIC() const {
    double total = 0;
    for (unsigned long int j = 0; j < I0.size(); j ++) total += I0.at(j);
    cout << "INITIAL CONDITION:" << endl;
    for (int j = 0; j <= 2*Rres; j ++) {
        for (int k = 0; k <= 2*Rres; k ++) {
            cout << scientific << setw(14) << I0.at(j*(2*Rres+1) + k) / total;
        }
        cout << endl;
    }
    cout << endl;
}

void Stats::printProfiles() const {
    double total1 = 0;
    double total2 = 0;
    for (unsigned long int j = 0; j < Ir.size(); j ++) {
        total1 += Ir.at(j);
        total2 += It.at(j);
    }
    cout << "TRANSMITTED BEAM PROFILE:" << endl;
    for (int j = 0; j <= 2*Rres; j ++) {
        for (int k = 0; k <= 2*Rres; k ++) {
            cout << scientific << setw(14) << It.at(j*(2*Rres+1) + k) / total2;
        }
        cout << endl;
    }
    cout << endl;
    
    cout << "REFLECTED BEAM PROFILE:" << endl;
    for (int j = 0; j <= 2*Rres; j ++) {
        for (int k = 0; k <= 2*Rres; k ++) {
            cout << scientific << setw(14) << Ir.at(j*(2*Rres+1) + k) / total1;
        }
        cout << endl;
    }
    cout << endl;
}

void Stats::print() const {
    //Print initial condition first
    printIC();
    printProfiles();
    
    //Now print energy deposition
    if (saturate)
        cout<<"ENERGY ABSORBED"<<endl;
    else
        cout << "ENERGY INCIDENT" << endl;
    for (int j = 0; j < Zres; j ++) {
        for (int k = 0; k < Rres; k ++)
            cout << scientific << setw(14) << DEP.at(j*Rres + k)/N0;
        cout<<endl;
    }
    cout<<endl;
    
    //energy scattered (in mJ/mJ-incident)
    if (saturate) {
        cout<<"ENERGY SCATTERED"<<endl;
        for (int j = 0; j < Zres; j ++) {
            for (int k = 0; k < Rres; k ++)
                cout << scientific << setw(14) << SCAT.at(j*Rres + k)/N0;
            cout<<endl;
        }
        cout<<endl;
    }
    
    //Energy scattered from depth Z
    if (saturate)
        cout<<"ENERGY INCIDENT AT DEPTH Z"<<endl;
    else
        cout<<"ENERGY ABSORBED AT DEPTH Z"<<endl;
    for (int j = 0; j < Zres; j ++) {
        cout << scientific << setw(14) << PENDEPTH.at(j)/N0;
    }
    cout<<endl<<endl;
    
    //Overflow matrices
    if (saturate)
        cout<<"ENERGY ABSORBED OUTSIDE GRID"<<endl;
    else
        cout<<"ENERGY INCIDENT OUTSIDE GRID"<<endl;
    cout<<"A(r; z < 0):        ";
    for (int j = 0; j < Rres; j ++)
        cout << scientific << setw(14) << DEPOF.at(j)/N0;
    cout<<endl<<"A(r; z > L):        ";
    for (int j = Rres; j > 0; j --)
        cout << scientific << setw(14) << DEPOF.at(1+Rres+Zres + j)/N0;
    cout<<endl<<"A(z; r > R):        ";
    for (int j = 0; j < Zres; j ++)
        cout << scientific << setw(14) << DEPOF.at(1+Rres+ j)/N0;
    cout<<endl<<"A(z < 0, r > R):    " << scientific << setw(14) << DEPOF.at(Rres)/N0;
    cout<<endl<<"A(z > L, r > R):    " << scientific << setw(14) << DEPOF.at(1+Rres+Zres)/N0;
    cout<<endl<<endl;
    
    if (saturate) {
        cout<<"ENERGY SCATTERED OUTSIDE GRID"<<endl;
        cout<<"S(r; z < 0):        ";
        for (int j = 0; j < Rres; j ++)
            cout << scientific << setw(14) << SCATOF.at(j)/N0;
        cout<<endl<<"S(r; z > L):        ";
        for (int j = Rres; j > 0; j --)
            cout << scientific << setw(14) << SCATOF.at(1+Rres+Zres + j)/N0;
        cout<<endl<<"S(z; r > R):        ";
        for (int j = 0; j < Zres; j ++)
            cout << scientific << setw(14) << SCATOF.at(1+Rres+ j)/N0;
        cout<<endl<<"S(z < 0, r > R):    " << scientific << setw(14) << SCATOF.at(Rres)/N0;
        cout<<endl<<"S(z > L, r > R):    " << scientific << setw(14) << SCATOF.at(1+Rres+Zres)/N0;
        cout<<endl<<endl;
    }
    
    //Beam radius at depth Z
    cout<<"BEAM RADIUS (STANDARD DEVIATION) AT DEPTH Z [nm/nm]"<<endl;
    for (int j = 0; j < Zres; j ++) {
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
        if (saturate)
            cout << "Energy absorption moments (time) [(ns/ns)**n]:" << endl;
        else
            cout << "Incident energy moments (time) [(ns/ns)**n]:" << endl;
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
        if (saturate)
            cout << "Energy absorption moments (space) [(nm/nm)**n]:" << endl;
        else
            cout << "Incident energy moments (space) [(nm/nm)**n]:" << endl;
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
    if (saturate) {
        if (momentlvl > 0) {
            cout <<endl <<endl;
            cout << "Energy scatter moments (time) [(ns/ns)**n]:" << endl;
            cout << "  E[T]     = " << STMOMENT[0]/_ASCAT << endl;
        } if (momentlvl > 1)
            cout << "  E[T2]    = " << STMOMENT[1]/_ASCAT << endl;
        if (momentlvl > 2)
            cout << "  E[T3]    = " << STMOMENT[2]/_ASCAT << endl;
        if (momentlvl > 3)
            cout << "  E[T4]    = " << STMOMENT[3]/_ASCAT << endl;
    }
    
    //Scattered energy moments (space)
    if (saturate) {
        if (momentlvl > 0) {
            cout <<endl <<endl;
            cout << "Energy scatter moments (space) [(nm/nm)**n]:" << endl;
            cout << "  E[Z]    = " << SMOMENTS[0]/_ASCAT << endl;
            cout << "  E[R]    = " << SMOMENTS[1]/_ASCAT << endl;
        } if (momentlvl > 1) {
            cout << "  E[Z2]   = " << SMOMENTS[2]/_ASCAT << endl;
            cout << "  E[RZ]   = " << SMOMENTS[3]/_ASCAT << endl;
            cout << "  E[R2]   = " << SMOMENTS[4]/_ASCAT << endl;
        } if (momentlvl > 2) {
            cout << "  E[Z3]   = " << SMOMENTS[5]/_ASCAT << endl;
            cout << "  E[Z2R]  = " << SMOMENTS[6]/_ASCAT << endl;
            cout << "  E[ZR2]  = " << SMOMENTS[7]/_ASCAT << endl;
            cout << "  E[R3]   = " << SMOMENTS[8]/_ASCAT << endl;
        } if (momentlvl > 3) {
            cout << "  E[Z4]   = " << SMOMENTS[9]/_ASCAT << endl;
            cout << "  E[Z3R]  = " << SMOMENTS[10]/_ASCAT << endl;
            cout << "  E[Z2R2] = " << SMOMENTS[11]/_ASCAT << endl;
            cout << "  E[ZR3]  = " << SMOMENTS[12]/_ASCAT << endl;
            cout << "  E[R4]   = " << SMOMENTS[13]/_ASCAT << endl;
        }
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
    cout << "Reflection/Transmission/Absorption coefficients [-]:" << endl;
    cout << "  Rdiffuse   = " << fixed << setw(10) << Rdiffuse/N0 << endl;
    cout << "  Tdiffuse   = " << fixed << setw(10) << Tdiffuse/N0 << endl;
    cout << "  Tballistic = " << fixed << setw(10) << Tballistic/N0 << endl;
    cout << "  A          = " << fixed << setw(10) << _ADEP/N0 << endl;
    cout << "  Atotal     = " << fixed << setw(10) << _ATOT/N0 << endl;
    
    //Reflection time/order
    cout << endl << endl;
    cout << "Diffuse Reflection/Transmission scattering order (mean, std) [-] and duration (mean, std) [ns/ns]:" << endl;
    cout << "  N[Rdiffuse]      = " << fixed << setw(10) << Nr / Rdiffuse << endl;
    cout << "  N[Tdiffuse]      = " << fixed << setw(10) << Nt / Tdiffuse << endl;
    cout << "  dN[Rdiffuse]     = " << fixed << setw(10) << sqrt(N2r*Rdiffuse - Nr*Nr) / Rdiffuse << endl;
    cout << "  dN[Tdiffuse]     = " << fixed << setw(10) << sqrt(N2t*Tdiffuse - Nt*Nt) / Tdiffuse << endl;
    cout << "  t[Rdiffuse]      = " << fixed << setw(10) << tRd / Rdiffuse << endl;
    cout << "  t[Tdiffuse]      = " << fixed << setw(10) << tTd / Tdiffuse << endl;
    cout << "  dt[Rdiffuse]     = " << fixed << setw(10) << sqrt(t2Rd*Rdiffuse - tRd*tRd) / Rdiffuse << endl;
    cout << "  dt[Tdiffuse]     = " << fixed << setw(10) << sqrt(t2Td*Tdiffuse - tTd*tTd) / Tdiffuse << endl;

    //Diffusion coefficient
    cout << endl << endl;
    cout << "Diffusion coefficients [(nm2/ns)/(nm2/ns)]:" << endl;
    cout << "  Dx = " << DX/ND/2.0 << endl;
    cout << "  Dy = " << DY/ND/2.0 << endl;
    cout << "  Dz = " << DZ/ND/2.0 << endl;
}

void Stats::finalize() {
    //Go through and correct intensities in case of coherence
    if (interfere)
        for (long unsigned int j = 0; j < It.size(); j ++) {
            It.at(j) *= It.at(j);
            Ir.at(j) *= Ir.at(j);
            I0.at(j) *= I0.at(j);
        }

    //Calculate total absorbed light inside volume
    _ADEP = 0;
    for (long unsigned int j = 0; j < DEP.size(); j ++)
        _ADEP += DEP.at(j);
        
    //Calculate total absorbed even outside volume
    _ATOT = _ADEP;
    for (long unsigned int j = 0; j < DEPOF.size(); j ++)
      _ATOT += DEPOF.at(j);
      
    //And scattered light if its different
    _ASCAT = 0;
    if (saturate)
        for (long unsigned int j = 0; j < SCAT.size(); j ++)
            _ASCAT += SCAT.at(j);
            
    //Calculate beam radius
    for (int j = 0; j < Zres; j ++)
        RADIUS.at(j) = sqrt(RADIUS.at(j)/PENDEPTH.at(j));
}

double Stats::E(const vec& x) const {
    int binR = getBinR(x);
    int binZ = getBinZ(x);
    return E(binZ, binR);
}

double Stats::Volume(int binR) const {
    if (binR == 0)
        return CONST_PI / Rres / Rres / Zres;
    else
        return CONST_PI / Rres / Rres / Zres * (1 + 2.0*binR);
}

double Stats::E(int binZ, int binR) const {
    if ( ((binZ) < 0) or (binZ >= Zres) or (binR >= Rres) )
        return 0;
    else
        return DEP.at(binZ*Rres + binR) / Volume(binR);
}
    
#endif
