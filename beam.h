#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>

#include "utility.h"
#include "photon.h"

#ifndef _BEAM_H_
#define _BEAM_H_

#define CONST_PI 3.1415926535
#define CONST_EPS 1e-10
#define CONST_HBAR 1.0545718176e-22

using namespace std;

struct Beam {

    //Definition for beam geometries and random spreading functions
    enum struct BeamSpread { Collimated = 0, Gaussian = 1, Lambertian = 2, Isotropic = 3 };
    enum struct BeamType { Uniform = 0, Gaussian = 1, UniformEllipse = 2, GaussianEllipse = 3, UniformAnnulus = 4, GaussianAnnulus = 5 };
    enum struct BeamDuration { Uniform = 0, Gaussian = 1, Cauchy = 2 };
    enum struct BeamSpectrum { Uniform = 0, Gaussian = 1, Cauchy = 2 };


    //Data members
    double E=1e5;                   //Photons per packet
    double sin0=0.0;                //Incident sin(theta) = 0 for normal incidence (-1 to 1)
    double Rb=5e2;                  //Beam radius [um]
    double Pb=0;                    //Beam spread parameter: standard deviation (gauss)/width (uniform)
    double Sb=2.5e2;                //Beam profile parameter: linewidth (UniformEllipse, Guass1D), Inner radius (UniformAnnulus, GaussianAnnulus)
    double Zb = 0;                  //Focus location of beam (0 for infinity/unfocused)
    double Tb = 0;                  //Laser pulse duration
    double wb=5.5e2;                //Angular frequency in THz (5.5e5 ~ 545 nm light)
    double dwb=0;                   //Angular frequency spread parameter (THz);
    double N0 = 0;
    double Rmax = 0;                //Maximum allowed radius 
    
    BeamType beamprofile = BeamType::Uniform;
    BeamSpread spreadfxn = BeamSpread::Collimated;
    BeamDuration beamdur = BeamDuration::Uniform;
    BeamSpectrum beamspec = BeamSpectrum::Uniform;

    Beam();
    void print() const;
    
    vector<Photon> sampleBeam(unsigned long int N);
};

Beam::Beam() {
    spreadfxn = BeamSpread::Gaussian;
    beamprofile = BeamType::Uniform;
    beamdur = BeamDuration::Uniform;
    beamspec = BeamSpectrum::Uniform;
    sin0=0.0; Rb=5e2; wb=5.5e2; dwb=0;
    Zb = 0; Tb = 0; Pb = 0; N0 = 0; E = 1e5;
}

void Beam::print() const {
    double Ab;
    cout << "Beam profile: ";
    switch (beamprofile) {
        case BeamType::Uniform : cout << "Uniform (r = " << Rb << " um)" << endl; Ab = CONST_PI * Rb * Rb; break;
        case BeamType::Gaussian : cout << "Gaussian (sigma = " << Rb << " um)" << endl; Ab = CONST_PI * Rb * Rb; break;
        case BeamType::UniformEllipse : cout << "Elliptical Uniform (a = " << Rb << " um, b = " << Sb << " um)" << endl;  Ab = CONST_PI*Rb*Sb; break;
        case BeamType::GaussianEllipse : cout << "Elliptical Gaussian (sigmax = " << Rb << " um, sigmay = " << Sb << " um)" << endl;  Ab = CONST_PI*Rb*Sb; break;
        case BeamType::UniformAnnulus : cout << "Annular Uniform (ri = " << Sb << " um, ro = " << Rb << " um)" << endl; Ab = CONST_PI*(Rb*Rb-Sb*Sb); break;
        case BeamType::GaussianAnnulus : cout << "Annular Gaussian (mu = " << Rb << " um, sigma = " << Sb << " um)" << endl; Ab = CONST_PI*(Rb*Rb-Sb*Sb); break;
        default : cout << "Other" << endl; Ab = 0; break;
    }
    
    cout << "Beam spread function: ";
    switch (spreadfxn) {
        case BeamSpread::Collimated : cout << "Collimated" << endl; break;
        case BeamSpread::Gaussian : cout << "Gaussian (sigma = " << Pb << ")" << endl; break;
        case BeamSpread::Lambertian : cout << "Lambertian [0, " << Pb << ")" << endl; break;
        case BeamSpread::Isotropic : cout << "Isotropic [0, " << Pb << ")" << endl; break;
        default: cout << "Other" << endl; break;
    }
    
    cout << "Beam temporal profile: ";
    switch (beamdur) {
        case BeamDuration::Uniform : cout << "Uniform (" << Tb << " ps)" << endl; break;
        case BeamDuration::Gaussian : cout << "Gaussian (sigma = " << Tb << " ps)" << endl; break;
        case BeamDuration::Cauchy : cout << "Gaussian (gamma = " << Tb << " ps)" << endl; break;
        default: cout << "Other" << endl; break;
    }
    
    cout << "Beam frequency spectrum: ";
    switch (beamspec) {
        case BeamSpectrum::Uniform : cout << "Uniform (" << wb << " +/- " << dwb << " THz)" << endl; break;
        case BeamSpectrum::Gaussian : cout << "Gaussian (sigma = " << wb << " +/- " << dwb << " THz)" << endl; break;
        case BeamSpectrum::Cauchy : cout << "Gaussian (gamma = " << wb << " +/- " << dwb << " THz)" << endl; break;
        default: cout << "Other" << endl; break;
    }
    
    cout << "Focus at depth: " << ((Zb > 0) ? Zb : INFINITY) << ((Zb > 0) ? " um" : "") << endl;
    cout << "Incident sin(theta): " << sin0 << endl;
    cout << "Photons per packet: " << E << endl;
    cout << "Total photons: " << E * N0 << endl;
    if (Ab > 0)
        cout << "Average fluence: " << E * N0 / Ab << " photons/um2     " << E * N0 * CONST_HBAR * wb / Ab * 1e8 << " J/cm2" << endl;
    if ((Ab > 0) and (Tb > 0))
        cout << "Average fluence rate: " << E * N0 / Ab / Tb << " photons/um2/ps" << E * N0 * CONST_HBAR * wb / Ab / Tb * 1e20 << " W/cm2" << endl;
    cout << endl;

}

vector<Photon> Beam::sampleBeam(unsigned long int N0) {
    //Initialize some values
    vector<Photon> PHOTONS(N0, Photon(vec(),vec(sin0, 0, sqrt(1-sin0*sin0)), E));
    this->N0 += N0;
    double eps, eps2;
    
    //Loop through and generate beam if needed
    double x,y;
    double r2;
    double Xf = Zb*tan(asin(sin0));
    double Sp = sin(Pb);
    for (long unsigned int i = 0; i < N0; i ++) {
        
        ///
        //  Calculate beam profile
        ///
        if (Rb > 1) {      
        
            //Rejection sampling to make sure we're less than Rmax
            r2 = 2*Rmax*Rmax;
            do {
                //Random numbers
                eps2 = roll();
                eps = roll();
            
                //Use the randos to generate x,y values based on model
                switch (beamprofile) {
                    case BeamType::Uniform :
                        x = Rb * sqrt(eps) * cos(2.0 * CONST_PI * eps2);
                        y = Rb * sqrt(eps) * sin(2.0 * CONST_PI * eps2);
                        break;
                        
                    case BeamType::Gaussian :
                        x = Rb * erfinvf((float)eps);
                        y = x * sin(2.0 * CONST_PI * eps2);
                        x = x * cos(2.0 * CONST_PI * eps2);
                        break;
                        
                    case BeamType::UniformEllipse :
                        eps2 = atan(Sb/Rb*tan(2*CONST_PI*eps2));
                        eps = Sb*Rb/sqrt(pow(Sb*cos(eps2),2)+pow(Rb*sin(eps2),2)) * sqrt(eps);
                        x = eps*cos(eps2);
                        y = eps*sin(eps2);
                        break;
                        
                    case BeamType::GaussianEllipse :
                        x = Rb * erfinvf((float)eps);
                        y = Sb * erfinvf((float)eps2);
                        break;
                        
                    case BeamType::UniformAnnulus :
                        x = sqrt(Sb*Sb + (Rb*Rb-Sb*Sb)*eps) * cos(2.0 * CONST_PI * eps2);
                        y = sqrt(Sb*Sb + (Rb*Rb-Sb*Sb)*eps) * sin(2.0 * CONST_PI * eps2);
                        break;
                        
                    case BeamType::GaussianAnnulus :
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
                PHOTONS.at(i).x.X = x;
                PHOTONS.at(i).x.Y = y;
                r2 = PHOTONS.at(i).x.r2();
                
            } while ((Rmax>0) and (r2>Rmax*Rmax));
        }
        
        ///
        //  Calculate beam frequency
        ///
        PHOTONS.at(i).v = wb;
        if (dwb > 0) {
            eps = roll();
            switch (beamspec) {
                case BeamSpectrum::Uniform :
                    PHOTONS.at(i).v += dwb * (eps-0.5);
                    break;
                case BeamSpectrum::Gaussian :
                    PHOTONS.at(i).v += dwb * erfinvf((float)eps);
                    break;
                case BeamSpectrum::Cauchy :
                    PHOTONS.at(i).v += - dwb / tan(CONST_PI * eps);
                    break;
            }
        }
        
        ///
        //  Calculate photon packet time
        ///
        if (Tb > 0) {
            eps = roll();
            switch (beamdur) {
                case BeamDuration::Uniform :
                    PHOTONS.at(i).t = eps * Tb;
                    break;
                case BeamDuration::Gaussian :
                    PHOTONS.at(i).t = Tb * erfinvf((float)eps);
                    break;
                    
                case BeamDuration::Cauchy :
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
    
        //Implement diffraction here? Maybe need to have entire set of photons first though...
        ///
        // TODO: DIFFRACTION CALCULATION
        ///
    
        //Exit now if we aren't further spreading the beam
        if (spreadfxn == BeamSpread::Collimated)
            continue;
            
        //Roll new random numbers
        eps2 = roll()*2*CONST_PI; //Azimuthal angle
        eps = roll();
        
        //Sample distributions for dispersion on top of focus
        x = 1;
        switch (spreadfxn) {
            case BeamSpread::Collimated :
                x = 1.0;
                break;
            case BeamSpread::Gaussian :
                x = 1.0 - Sp * erfinvf((float)eps);
                break;
            case BeamSpread::Lambertian :
                if (Pb > 0)
                    x = sqrt(1.0 - pow(eps*Sp,2));
                else
                    x = sqrt(1.0 - pow(eps,2));
                break;
            case BeamSpread::Isotropic :
                if (Pb > 0)
                    x = 1.0 - eps*Sp;
                else
                    x = 1.0 - eps;
                break;
        }
        PHOTONS.at(i).mu = PHOTONS.at(i).mu * x + PHOTONS.at(i).mu.perp(eps2) * sqrt(1.0-x*x);
    }
    return PHOTONS;
}

#endif