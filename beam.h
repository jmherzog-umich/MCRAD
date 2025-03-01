#include <iostream>
#include <vector>

#include "utility.h"
#include "photon.h"
#include "spectrum.h"

#ifndef _BEAM_H_
#define _BEAM_H_

struct Beam {

    //Definition for beam geometries and random spreading functions
    enum struct BeamSpread { Collimated = 0, Gaussian = 1, Lambertian = 2, Isotropic = 3 };
    enum struct BeamType { Uniform = 0, Gaussian = 1, UniformEllipse = 2, GaussianEllipse = 3, UniformAnnulus = 4, GaussianAnnulus = 5 };
    enum struct BeamDuration { Uniform = 0, Gaussian = 1, Cauchy = 2 };

    //Data members
    double E=1e5;                   //Photons per packet
    double sin0=0.0;                //Incident sin(theta) = 0 for normal incidence (-1 to 1)
    double Rb=5e2;                  //Beam radius [um]
    double Pb=0;                    //Beam spread parameter: standard deviation (gauss)/width (uniform)
    double Sb=2.5e2;                //Beam profile parameter: linewidth (UniformEllipse, Guass1D), Inner radius (UniformAnnulus, GaussianAnnulus)
    double Zb = 0;                  //Focus location of beam (0 for infinity/unfocused)
    double Tb = 0;                  //Laser pulse duration
    double N0 = 0;
    double Rmax = 0;                //Maximum allowed radius 
    
    BeamType beamprofile = BeamType::Uniform;
    BeamSpread spreadfxn = BeamSpread::Collimated;
    BeamDuration beamdur = BeamDuration::Uniform;
    Spectrum beamspec;

    Beam();
    void print(std::ostream& oout) const;
    
    bool set(const std::string& key, const std::vector<std::string>& val);
    
    void writedbheader(std::ostream& OF) const;
    void writedb(std::ostream& OF) const;
    
    std::vector<Photon> sampleBeam(unsigned long int N);
};

#endif
