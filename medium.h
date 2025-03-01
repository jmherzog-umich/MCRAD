#include <iostream>
#include <vector>

#include "spectrum.h"
#include "utility.h"

#ifndef _MEDIUM_H_
#define _MEDIUM_H_

struct Medium {

    //Function definitions for wavelength dependence in cross-sections and g, and 
    enum struct AnisotropyModel { Constant=0 };
    enum struct PhaseFunction { HenyeyGreenstein=0, Rayleigh=1, Draine=2 };

    double _dens;                       //Particle density um-3
    double _n;                          //Refractive index
    std::vector<double> _g;             //Phase function parameters
    double _Ss;                         //Scattering cross-section parameters
    double _Sa;                         //Absorption cross-section parameters
    double _tau;                        //Excited state lifetime
    double _FQY;                        //Fluorescence quantum yield
    
    double Fmax;                        //Maximum frequency of simulation
    double Fmin;                        //Minimum frequency of simulation
     
    PhaseFunction phase;                //Phase function to use for calculation
    Spectrum xca, xcs, xf;              //Cross-section models for absorption and scattering
    
    void print(std::ostream& oout) const;
    void print_at_f(std::ostream& oout, double f) const;
    void writedbheader(std::ostream& OF) const;
    void writedb(std::ostream& OF) const;
    
    bool set(const std::string& key, const std::vector<std::string>& val);
    
    double dens() const;
    double FQY() const;
    double g(double w = 0) const;
    double g2(double w = 0) const;
    double tau() const;
    
    double Ss(double w) const;
    double Sa(double w) const;
    double Se(double w) const;
    
    double Ss() const;
    double Sa() const;
    double Se() const;
    
    double ks(double w = 0) const;
    double ka(double w = 0) const;
    double ke(double w = 0) const;
    
    double ws(double w = 0) const;
    double wa(double w = 0) const;
    double we(double w = 0) const;
    
    double ls(double w = 0) const;
    double la(double w = 0) const;
    double le(double w = 0) const;
    
    double albedo(double w = 0) const;
    double n(double w = 0) const;
    
    double scatter(double eps, double w = 0) const;
    double pscatter(double cost, double w = 0) const;

    double peak_v() const;
    double emit_tau(double eps, double w = 0) const;
    double emit_v(double eps, double w = 0) const;
    
    Medium(double Fmin=0, double Fmax=1000);
    
};

#endif
