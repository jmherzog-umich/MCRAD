#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>

#include "utility.h"
#include "photon.h"
#include "spectrum.h"

#include "beam.h"

#define CONST_PI            3.1415926535897932
#define CONST_C             299.792458
#define CONST_HBAR          1.054571817679489e-22
#define CONST_HK            47.992430733662212
#define CONST_APERY         1.2020569031595942
#define CONST_WIEN          1.5936242600400401
#define CONST_EPS           1e-10
#define CONST_PLANCKMAX     0.6476102378919149

Beam::Beam() {
    spreadfxn = BeamSpread::Collimated;
    beamprofile = BeamType::Uniform;
    beamdur = BeamDuration::Uniform;
    beamspec = Spectrum(Spectrum::SpectrumModel::Impulse, {800});
    sin0=0.0; Rb=5e2;
    Zb = 0; Tb = 0; Pb = 0; N0 = 0; E = 1e5;
}

void Beam::print(std::ostream& oout) const {
    double Ab;
    oout << "Beam profile: ";
    switch (beamprofile) {
        case BeamType::Uniform : oout << "Uniform (r = " << Rb << " um)" << std::endl; Ab = CONST_PI * Rb * Rb; break;
        case BeamType::Gaussian : oout << "Gaussian (sigma = " << Rb << " um)" << std::endl; Ab = CONST_PI * Rb * Rb; break;
        case BeamType::UniformEllipse : oout << "Elliptical Uniform (a = " << Rb << " um, b = " << Sb << " um)" << std::endl;  Ab = CONST_PI*Rb*Sb; break;
        case BeamType::GaussianEllipse : oout << "Elliptical Gaussian (sigmax = " << Rb << " um, sigmay = " << Sb << " um)" << std::endl;  Ab = CONST_PI*Rb*Sb; break;
        case BeamType::UniformAnnulus : oout << "Annular Uniform (ri = " << Sb << " um, ro = " << Rb << " um)" << std::endl; Ab = CONST_PI*(Rb*Rb-Sb*Sb); break;
        case BeamType::GaussianAnnulus : oout << "Annular Gaussian (mu = " << Rb << " um, sigma = " << Sb << " um)" << std::endl; Ab = CONST_PI*(Rb*Rb-Sb*Sb); break;
        default : oout << "Other" << std::endl; Ab = 0; break;
    }
    
    oout << "Beam spread function: ";
    switch (spreadfxn) {
        case BeamSpread::Collimated : oout << "Collimated" << std::endl; break;
        case BeamSpread::Gaussian : oout << "Gaussian (sigma = " << Pb << ")" << std::endl; break;
        case BeamSpread::Lambertian : oout << "Lambertian [0, " << Pb << ")" << std::endl; break;
        case BeamSpread::Isotropic : oout << "Isotropic [0, " << Pb << ")" << std::endl; break;
        default: oout << "Other" << std::endl; break;
    }
    
    oout << "Beam temporal profile: ";
    switch (beamdur) {
        case BeamDuration::Uniform : oout << "Uniform (" << Tb << " ps)" << std::endl; break;
        case BeamDuration::Gaussian : oout << "Gaussian (sigma = " << Tb << " ps)" << std::endl; break;
        case BeamDuration::Cauchy : oout << "Gaussian (gamma = " << Tb << " ps)" << std::endl; break;
        default: oout << "Other" << std::endl; break;
    }
    
    oout << "Beam frequency spectrum: ";
    oout << beamspec.name() << std::endl << "   ";
    oout << beamspec.paramstring() << std::endl;
    
    oout << "Focus at depth: " << ((Zb > 0) ? Zb : INFINITY) << ((Zb > 0) ? " um" : "") << std::endl;
    oout << "Incident sin(theta): " << sin0 << std::endl;
    oout << "Photons per packet: " << E << std::endl;
    oout << "Total photons: " << E * N0 << std::endl;
    if (Ab > 0)
        oout << "Average fluence: " << E * N0 / Ab << " photons/um2     ~" << E * N0 * CONST_HBAR * beamspec.peak() / Ab * 1e8 << " J/cm2" << std::endl;
    if ((Ab > 0) and (Tb > 0))
        oout << "Average fluence rate: " << E * N0 / Ab / Tb << " photons/um2/ps     ~" << E * N0 * CONST_HBAR * beamspec.peak() / Ab / Tb * 1e20 << " W/cm2" << std::endl;
    oout << std::endl;

}

std::vector<Photon> Beam::sampleBeam(unsigned long int N0) {
    //Initialize some values
    std::vector<Photon> PHOTONS(N0, Photon(vec(0,0,CONST_EPS),vec(sin0, 0, sqrt(1-sin0*sin0)), E));
    this->N0 += N0;
    double eps, eps2;
    
    //Loop through and generate beam if needed
    double x,y;
    double r2;
    double Xf = Zb*std::tan(std::asin(sin0));
    double Sp = sin(Pb);
    for (long unsigned int i = 0; i < N0; i ++) {
        
        ///
        //  Calculate beam profile
        ///
        if (Rb > 1) {      
        
            //Rejection sampling to make sure we're less than Rmax
            r2 = 2*Rmax*Rmax;
            do {
                //Use the randos to generate x,y values based on model
                switch (beamprofile) {
                    case BeamType::Uniform :
                        //Random numbers
                        eps2 = util::roll();
                        eps = util::roll();
                        
                        //Calculation
                        x = Rb * std::sqrt(eps) * std::cos(2.0 * CONST_PI * eps2);
                        y = Rb * std::sqrt(eps) * std::sin(2.0 * CONST_PI * eps2);
                        break;
                        
                    case BeamType::Gaussian :
                        //Random numbers
                        eps2 = util::roll();
                        
                        //Calculation
                        x = Rb * util::erfinvf(util::rollf());
                        y = x * std::sin(2.0 * CONST_PI * eps2);
                        x = x * std::cos(2.0 * CONST_PI * eps2);
                        break;
                        
                    case BeamType::UniformEllipse :
                        eps2 = std::atan(Sb/Rb*tan(2*CONST_PI*util::roll()));
                        eps = Sb*Rb/std::hypot(Sb*std::cos(eps2), Rb*std::sin(eps2)) * std::sqrt(util::roll());
                        x = eps*std::cos(eps2);
                        y = eps*std::sin(eps2);
                        break;
                        
                    case BeamType::GaussianEllipse :
                        x = Rb * util::erfinvf(util::rollf());
                        y = Sb * util::erfinvf(util::rollf());
                        break;
                        
                    case BeamType::UniformAnnulus :
                        //Random numbers
                        eps2 = util::roll();
                        eps = util::roll();
                        
                        //Calculation
                        x = std::sqrt(Sb*Sb + (Rb*Rb-Sb*Sb)*eps) * std::cos(2.0 * CONST_PI * eps2);
                        y = std::sqrt(Sb*Sb + (Rb*Rb-Sb*Sb)*eps) * std::sin(2.0 * CONST_PI * eps2);
                        break;
                        
                    case BeamType::GaussianAnnulus :
                        //Random numbers
                        eps2 = util::roll();
                        
                        //Calculation
                        x = Rb + (Sb - Rb) * util::erfinvf(util::rollf());
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
        PHOTONS.at(i).v = beamspec.sample();
        
        ///
        //  Calculate photon packet time
        ///
        if (Tb > 0) {
            switch (beamdur) {
                case BeamDuration::Uniform :
                    PHOTONS.at(i).t = util::roll() * Tb;
                    break;
                case BeamDuration::Gaussian :
                    PHOTONS.at(i).t = Tb * util::erfinvf(util::rollf());
                    break;
                    
                case BeamDuration::Cauchy :
                    PHOTONS.at(i).t = -Tb / std::tan(CONST_PI * util::roll());
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
            
        //Sample distributions for dispersion on top of focus
        x = 1;
        switch (spreadfxn) {
            case BeamSpread::Collimated :
                x = 1.0;
                break;
            case BeamSpread::Gaussian :
                x = 1.0 - Sp * util::erfinvf(util::rollf());
                break;
            case BeamSpread::Lambertian :
                if (Pb > 0)
                    x = std::sqrt(1.0 - std::pow(util::roll()*Sp,2));
                else
                    x = std::sqrt(1.0 - std::pow(util::roll(),2));
                break;
            case BeamSpread::Isotropic :
                if (Pb > 0)
                    x = 1.0 - util::roll()*Sp;
                else
                    x = 1.0 - util::roll();
                break;
        }
        PHOTONS.at(i).mu = PHOTONS.at(i).mu * x + PHOTONS.at(i).mu.perp(util::roll()*2*CONST_PI) * sqrt(1.0-x*x);
    }
    return PHOTONS;
}

void Beam::writedbheader(std::ostream& OF) const {
    OF << "E";
}

void Beam::writedb(std::ostream& OF) const {
    OF << std::scientific << std::setprecision(8) << E;
}

bool Beam::set(const std::string& key, const std::vector<std::string>& val) {
    if (!key.compare("beam-photons-per-packet"))
        E = std::stod(val.at(0));
    else if (!key.compare("beam-initial-angle"))
        sin0 = std::stod(val.at(0));
    else if (!key.compare("beam-radius-cutoff"))
        Rmax = std::stod(val.at(0));
    else if (!key.compare("beam-focal-depth"))
        Zb = std::stod(val.at(0));
    else if (!key.compare("beam-profile")) {
        Rb = std::max(1.0, std::stod(val.at(0)));
        beamprofile = (val.size()>1) ? (Beam::BeamType)std::stoul(val.at(1)) : Beam::BeamType::Uniform;
        if (val.size() > 2) Sb = std::stod(val.at(2));
    } else if (!key.compare("beam-divergence")) {
        Pb = std::stod(val.at(0));
        spreadfxn = (val.size()>1) ? (Beam::BeamSpread)std::stoul(val.at(1)) : Beam::BeamSpread::Gaussian;
    } else if (!key.compare("beam-pulse-duration")) {
        Tb = std::stod(val.at(0));
        beamdur = (val.size()>1) ? (Beam::BeamDuration)std::stoul(val.at(1)) : Beam::BeamDuration::Uniform;
    } else if (!key.compare("beam-spectrum")) {
        std::vector<double> vals;
        for (unsigned long int i = 1; i < val.size(); i ++)
            vals.push_back(std::stod(val.at(i)));
        beamspec = Spectrum((Spectrum::SpectrumModel)std::stoul(val.at(0)), vals);
    } else
        return false;
    return true;
}

