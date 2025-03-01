#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>

#include "spectrum.h"
#include "utility.h"

#include "medium.h"

#define CONST_PI            3.1415926535897932
#define CONST_C             299.792458
#define CONST_HBAR          1.054571817679489e-22
#define CONST_HK            47.992430733662212
#define CONST_APERY         1.2020569031595942
#define CONST_WIEN          1.5936242600400401
#define CONST_EPS           1e-10
#define CONST_PLANCKMAX     0.6476102378919149

std::string _phasemodelname(Medium::PhaseFunction f) {
    switch(f) {
        case Medium::PhaseFunction::HenyeyGreenstein: return "Henyey-Greenstein"; break;
        case Medium::PhaseFunction::Rayleigh: return "Rayleigh"; break;
        case Medium::PhaseFunction::Draine: return "Draine"; break;
        default: return "Other"; break;
    }
}

void Medium::print_at_f(std::ostream& oout, double f) const {
    oout << "   Medium refractive index: " << n(f) << std::endl;
    oout << "   Scattering anisotropy <cos theta>: " << g(f) << std::endl;
    oout << "   Scattering cross-section: " << Ss(f) << " um2" << std::endl;
    oout << "   Absorption cross-section: " << Sa(f) << " um2" << std::endl;
    oout << "   Albedo: " << albedo(f) << std::endl;
    
    //print MFPs
    oout << "   Scattering mean free path: " << ls(f) << " um" << std::endl;
    oout << "   Scattering rate: " << ws(f) << " THz" << std::endl;
    oout << "   Absorption mean free path: " << la(f) << " um" << std::endl;
    oout << "   Absorption rate: " << wa(f) << " THz" << std::endl;
    oout << "   Collision mean free path: " << le(f) << " um" << std::endl;
    oout << "   Collision rate: " << we(f) << " THz" << std::endl << std::endl;
}

void Medium::print(std::ostream& oout) const {
    oout << "Medium density: " << dens() << " um-3" << std::endl;
    oout << "Phase function: " << _phasemodelname(phase) << std::endl;

    oout << "Absorption cross-section model: " << xca.name() << std::endl;
    oout << "   Integrated cross-section: " << xca.norm()*_Sa << " um2-THz" << std::endl;
    oout << "   " << (xca.peakLabel() ? "Peak" : "Specified") << " value: " << xca.max()*_Sa;
    oout << " um2 at " << xca.peak() << " THz" << std::endl;
    oout << "   Width: " << xca.width() << " THz" << std::endl;

    oout << "Scattering cross-section model: " << xcs.name() << std::endl;
    oout << "   Integrated cross-section: " << xcs.norm()*_Ss << " um2-THz" << std::endl;
    oout << "   " << (xcs.peakLabel() ? "Peak" : "Specified") << " value: " << xcs.max()*_Ss;
    oout << " um2 at " << xcs.peak() << " THz" << std::endl;
    oout << "   Width: " << xcs.width() << " THz" << std::endl;
    
    oout << "Fluorescence spectrum model: " << xf.name() << std::endl;
    oout << "   Quantum yield: " << _FQY << std::endl;
    oout << "   " << (xf.peakLabel() ? "Peak" : "Specified") << " value: " << xf.max()*_FQY/xf.norm();
    oout << " THz-1 at " << xf.peak() << " THz" << std::endl;
    oout << "   Width: " << xf.width() << " THz" << std::endl;
    
    oout << std::endl;
}

Medium::Medium(double Wmin, double Wmax) {
    Fmax = Wmax; Fmin = Wmin;
    Spectrum::setFreqLimits(Fmin, Fmax);
    _g = {0.98, 1};
    phase = PhaseFunction::HenyeyGreenstein;
    _FQY = 1; _Ss = 0.5; _Sa = 0.1;
    _dens = 1e-4; _n = 1.33;
    _tau = 1000;
}

double Medium::peak_v() const {
    return xf.peak();
}

double Medium::g(double w) const {
    return _g.at(0);
}
double Medium::g2(double w) const {
    if (_g.size() > 1)
        return _g.at(1);
    else
        return 0;
}
double Medium::Ss(double w) const {
    return _Ss * xcs.evaluate(w);
}
double Medium::Sa(double w) const {
    return _Sa * xca.evaluate(w);
}
double Medium::Ss() const {
    return _Ss * xcs.norm();
}
double Medium::Sa() const {
    return _Sa * xca.norm();
}
double Medium::dens() const {
    return _dens;
}
double Medium::n(double w) const {
    return _n;
}
double Medium::FQY() const {
    return _FQY;
}
double Medium::tau() const {
    return _tau;
}

//These are trivial
double Medium::Se(double w) const {
    return Sa(w) + Ss(w);
}
double Medium::Se() const {
    return Sa() + Ss();
}
double Medium::ks(double w) const {
    return Ss(w)*_dens;
}
double Medium::ka(double w) const {
    return Sa(w)*_dens;
}
double Medium::ke(double w) const {
    return Se(w)*_dens;
}
double Medium::ls(double w) const {
    return 1/ks(w);
}
double Medium::la(double w) const {
    return 1/ka(w);
}
double Medium::le(double w) const {
    return 1/ke(w);
}
double Medium::ws(double w) const {
    return ks(w)*CONST_C/n(w);
}
double Medium::wa(double w) const {
    return ka(w)*CONST_C/n(w);
}
double Medium::we(double w) const {
    return ke(w)*CONST_C/n(w);
}
double Medium::albedo(double w) const {
    return Ss(w) / (Ss(w) + Sa(w));
}

double Medium::scatter(double eps, double w) const {
    double gg = g(w); double gg2 = g2(w);
    double t0, t1, t2, t3, t4, t5, t6, t7, t8, t9;
    switch (phase) {
        case PhaseFunction::HenyeyGreenstein :
            if (abs(gg)>CONST_EPS)
                return (1+gg*gg - std::pow((1 - gg * gg) / (1-gg + 2.0 * gg * eps), 2))/2.0/gg;
            else
                return 1-2*eps;
                
        case PhaseFunction::Rayleigh :
            eps = 2*eps - 1;
            t0 = std::sqrt(4*eps*eps+1)+2*eps;
            return (std::pow(t0, 2.0/3.0)-1)/std::pow(t0, 1.0/3.0);
            
        case PhaseFunction::Draine :
            t0 = gg2*(1-gg*gg);
            t1 = gg2*(std::pow(gg,4)-1);
            t2 = -3*(4*(std::pow(gg,4) - std::pow(gg,2)) + t1 * (1+gg*gg));
            t3 = gg*(2*eps - 1);
            t4 = 3*gg*gg * (1 + t3) * gg2 * (2 + gg*gg*(1 + (1 + 2*gg*gg)*t3));
            t5 = t0 * (t1 * t2 + t4*t4) + pow(t1, 3);
            t6 = t0 * 4 * (std::pow(gg,4)-gg*gg);
            t7 = std::pow(t5 + std::sqrt(t5*t5 - std::pow(t6,3)), 1.0/3.0);
            t8 = 2*(t1 + t6/t7 + t7)/t0;
            t9 = sqrt(6*(1+gg*gg)+t8);
            return gg/2 + (1 - 0.25*std::pow(sqrt(6*(1+gg*gg)-t8+8*t4/t0/t9)-t9 ,2))/2/gg;
            
        default:
            return 1-2*eps;
    }
}

double Medium::pscatter(double cost, double w) const {
    double gg = g(w); double gg2 = g2(w);
    switch (phase) {
        case PhaseFunction::HenyeyGreenstein :
            return 0.25 / CONST_PI * (1-gg*gg) / std::pow(1 + gg*gg - 2*gg*cost, 1.5);
        
        case PhaseFunction::Draine :
            return 0.25 / CONST_PI * (1-gg*gg) / std::pow(1 + gg*gg - 2*gg*cost, 1.5) 
                    * (1 + gg2*cost*cost) / (1 + gg2 * (1 + 2 * gg*gg)/3.0);
                    
        case PhaseFunction::Rayleigh :
            return 0.75 * (1 - cost*cost);
            
        default :
            return 0;
    }
}

double Medium::emit_tau(double eps, double w) const {
    return _tau * eps;
}

double Medium::emit_v(double eps, double w) const {
    return xf.sample();
}

void Medium::writedbheader(std::ostream& OF) const {
    OF << "n,dens,phaseFxn,xca,xcs,xf,g,Ss,Sa,a,tau,FQY";
}

void Medium::writedb(std::ostream& OF) const {
    OF << std::setprecision(8) << n() << ",";
    OF << std::setprecision(8) << dens() << ",";
    OF << _phasemodelname(phase) << ",";
    OF << xca.name() << ",";
    OF << xcs.name() << ",";
    OF << xf.name() << ",";
    OF << std::scientific << std::setprecision(8) << g() << ",";
    OF << std::scientific << std::setprecision(8) << Ss() << ",";
    OF << std::scientific << std::setprecision(8) << Sa() << ",";
    OF << std::scientific << std::setprecision(8) << albedo() << ",";
    OF << std::scientific << std::setprecision(8) << tau() << ",";
    OF << std::scientific << std::setprecision(8) << FQY();
}

bool Medium::set(const std::string& key, const std::vector<std::string>& val) {
    if (!key.compare("spectrum-scattering")) {
        _Ss = std::stod(val.at(0));
        if (val.size() > 1) {
            std::vector<double> vals;
            for (unsigned long int i = 2; i < val.size(); i ++)
                vals.push_back(stod(val.at(i)));
            xcs = Spectrum((Spectrum::SpectrumModel)std::stoul(val.at(1)), vals);
        }
    } else if (!key.compare("spectrum-absorption")) {
        _Sa = std::stod(val.at(0));
        if (val.size() > 1) {
            std::vector<double> vals;
            for (unsigned long int i = 2; i < val.size(); i ++)
                vals.push_back(stod(val.at(i)));
            xca = Spectrum((Spectrum::SpectrumModel)std::stoul(val.at(1)), vals);
        }
    } else if (!key.compare("spectrum-fluorescence")) {
        _FQY = std::stod(val.at(0));
        if (val.size() > 1) {
            std::vector<double> vals;
            for (unsigned long int i = 2; i < val.size(); i ++)
                vals.push_back(stod(val.at(i)));
            xf = Spectrum((Spectrum::SpectrumModel)std::stoul(val.at(1)), vals);
        }
    } else if (!key.compare("scattering-anisotropy")) {
        _g = std::vector<double>();
        phase = (Medium::PhaseFunction)std::stoul(val.at(0));
        for (unsigned long int i = 1; i < val.size(); i ++)
            _g.push_back(std::stod(val.at(i)));
    } else if (!key.compare("lifetime-fluorescence"))
        _tau = std::stod(val.at(0));
    else if (!key.compare("index-medium"))
        _n = std::stod(val.at(0));
    else if (!key.compare("density"))
        _dens = std::stod(val.at(0));
    else
        return false;
    return true;
}

