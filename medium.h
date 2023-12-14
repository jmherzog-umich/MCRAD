#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>

#include "utility.h"

#ifndef _MEDIUM_H_
#define _MEDIUM_H_

#define CONST_PI            3.1415926535897932
#define CONST_C             299.792458
#define CONST_HBAR          1.054571817679489e-22
#define CONST_HK            7.6382325822577381
#define CONST_APERY         1.2020569031595942
#define CONST_WIEN          1.5936242600400401
#define CONST_EPS           1e-10
#define CONST_PLANCKMAX     0.6476102378919149
#define CONST_BMAX          35      //at x=35, Planck distribution CDF > 1-CONST_EPS

using namespace std;

struct Medium {

    //Function definitions for wavelength dependence in cross-sections and g, and 
    enum struct SpectrumModel { Constant=0, Gaussian=1, Cauchy=2, Laplace=3, Impulse=4, Rayleigh=5, Blackbody=6, Binary=12, Uniform=13 };
    enum struct AnisotropyModel { Constant=0 };
    enum struct PhaseFunction { HenyeyGreenstein=0, Rayleigh=1, Draine=2 };

    double _dens;                       //Particle density um-3
    double _n;                          //Refractive index
    vector<double> _g;                  //Phase function parameters
    double _Ss;                         //Scattering cross-section parameters
    double _Sa;                         //Absorption cross-section parameters
    double _tau;                        //Excited state lifetime
    double _FQY;                        //Fluorescence quantum yield
    vector<double> _cs;                 //Scattering cross-section distribution parameters
    vector<double> _ca;                 //Absorption cross-section distribution parameters
    vector<double> _cf;                 //Fluorescence spectrum parameters
    
    double Fmax;                        //Maximum frequency of simulation
    double Fmin;                        //Minimum frequency of simulation
     
    PhaseFunction phase;                //Phase function to use for calculation
    SpectrumModel xca, xcs, xf;         //Cross-section models for absorption and scattering
    
    void print() const;
    void print_at_f(double f) const;
    void writedbheader(ofstream& OF) const;
    void writedb(ofstream& OF) const;
    
    bool set(const string& key, const vector<string>& val);
    
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
    
    private:
        double _evalSpec(SpectrumModel x, double w, const vector<double>& S) const;
        double _meanSpec(SpectrumModel x, const vector<double>& S) const;
        double _sampleSpec(SpectrumModel x, double eps, const vector<double>& S) const;
    
};

string _xcmodelname(Medium::SpectrumModel x) {
    switch (x) {
        case Medium::SpectrumModel::Constant: return "Constant"; break;
        case Medium::SpectrumModel::Gaussian: return "Gaussian"; break;
        case Medium::SpectrumModel::Cauchy: return "Cauchy"; break;
        case Medium::SpectrumModel::Laplace: return "Laplace"; break;
        case Medium::SpectrumModel::Binary: return "Binary"; break;
        case Medium::SpectrumModel::Uniform: return "Uniform"; break;
        case Medium::SpectrumModel::Impulse: return "Impulse"; break;
        case Medium::SpectrumModel::Blackbody: return "Blackbody"; break;
        case Medium::SpectrumModel::Rayleigh: return "Rayleigh"; break;
        default: return "Other"; break;
    }
}

string _phasemodelname(Medium::PhaseFunction f) {
    switch(f) {
        case Medium::PhaseFunction::HenyeyGreenstein: return "Henyey-Greenstein"; break;
        case Medium::PhaseFunction::Rayleigh: return "Rayleigh"; break;
        case Medium::PhaseFunction::Draine: return "Draine"; break;
        default: return "Other"; break;
    }
}

void Medium::print_at_f(double f) const {
    cout << "   Medium refractive index: " << n(f) << endl;
    cout << "   Scattering anisotropy <cos theta>: " << g(f) << endl;
    cout << "   Scattering cross-section: " << Ss(f) << " um2" << endl;
    cout << "   Absorption cross-section: " << Sa(f) << " um2" << endl;
    cout << "   Albedo: " << albedo(f) << endl;
    
    //Print MFPs
    cout << "   Scattering mean free path: " << ls(f) << " um" << endl;
    cout << "   Scattering rate: " << ws(f) << " THz" << endl;
    cout << "   Absorption mean free path: " << la(f) << " um" << endl;
    cout << "   Absorption rate: " << wa(f) << " THz" << endl;
    cout << "   Collision mean free path: " << le(f) << " um" << endl;
    cout << "   Collision rate: " << we(f) << " THz" << endl << endl;
}

void Medium::print() const {
    cout << "Phase function: " << _phasemodelname(phase) << endl;
    cout << "Absorption cross-section model: " << _xcmodelname(xca) << endl;
    cout << "   Integrated absorption cross-section: " << Sa() << " um2-THz" << endl;
    cout << "Scattering cross-section model: " << _xcmodelname(xcs) << endl;
    cout << "   Integrated scattering cross-section: " << Ss() << " um2-THz" << endl;
    cout << "Medium density: " << dens() << " um-3" << endl;
    cout << "Fluorescence spectrum model: " << _xcmodelname(xf) << endl;
    cout << "   Fluorescence quantum yield: " << FQY() << endl;
    cout << "   Peak fluorescence frequency: " << peak_v() << " THz" << endl << endl;
}

Medium::Medium(double Wmin, double Wmax) {
    Fmax = Wmax; Fmin = Wmin;
    _g = {0.98, 1};
    phase = PhaseFunction::HenyeyGreenstein;
    xcs = SpectrumModel::Constant;
    xca = SpectrumModel::Constant;
    xf = SpectrumModel::Impulse;
    
    _FQY = 1; _Ss = 0.5; _Sa = 0.1;
    _dens = 1e-4; _n = 1.33;
    _tau = 1000;
}

//Evaluates the spectral distribution function (just the PDF!)
double Medium::_evalSpec(Medium::SpectrumModel x, double w, const vector<double>& S) const {
    switch (x) {
        default:
        case Medium::SpectrumModel::Constant :
            return 1.0;
            break;
        case Medium::SpectrumModel::Gaussian :
            return 1.0 / S.at(1) / sqrt(2.*CONST_PI) * exp(-pow((w-S.at(0))/S.at(1),2)/2);
            break;
        case Medium::SpectrumModel::Cauchy :
            return 1.0 / CONST_PI / S.at(1) / (1 + pow((w-S.at(0))/S.at(1),2));
            break;
        case Medium::SpectrumModel::Laplace :
            return 0.5 / S.at(1) * exp(-abs(w-S.at(0))/S.at(1));
            break;
        case Medium::SpectrumModel::Binary :
            return (w>S.at(0)) ? S.at(2) : S.at(1);
            break;
        case Medium::SpectrumModel::Impulse :
            return (abs(w-S.at(1)) < CONST_EPS) ? HUGE_VAL : 0.0;
            break;
        case Medium::SpectrumModel::Blackbody :
            return pow(CONST_HK/S.at(0),3)/2.0/CONST_APERY * pow(w,2)/(exp(CONST_HK*w/S.at(0))-1);
            break;
        case Medium::SpectrumModel::Uniform : 
            if (w <= S.at(1) and w >= S.at(0))
                return 1.0/(S.at(1)-S.at(0));
            else
                return 0;
            break;
        case Medium::SpectrumModel::Rayleigh :
            return pow(w/S.at(0), 4);
    }
    return S.at(0);
}

//Calculates the mean frequency of the spectrum model, or for non-convergent models the specified/target wavelength
double Medium::_meanSpec(Medium::SpectrumModel x, const vector<double>& S) const {
    switch (x) {
        case Medium::SpectrumModel::Constant :
            return (Fmax+Fmin)/2.0;
            break;
        case Medium::SpectrumModel::Binary :
            return S.at(0) - S.at(0)*(S.at(1) - 0.5);
            break;
        case Medium::SpectrumModel::Uniform :
            return (S.at(0) + S.at(1))/2.0;
            break;
        case Medium::SpectrumModel::Blackbody :
            return S.at(0) / CONST_HK * CONST_WIEN;
            break;
        case Medium::SpectrumModel::Rayleigh :
        case Medium::SpectrumModel::Gaussian :
        case Medium::SpectrumModel::Cauchy :
        case Medium::SpectrumModel::Laplace :
        case Medium::SpectrumModel::Impulse :
            return S.at(0);
            break;
        default:
            cerr << "Error: no mean available for spectrum model """ << _xcmodelname(x) << """" << endl;
            exit(-1);
            return 0;
    }
}

//Randomly samples the spectrum model
double Medium::_sampleSpec(Medium::SpectrumModel x, double eps, const vector<double>& S) const {
    double tmp;
    switch (x) {
        case Medium::SpectrumModel::Constant :
            return Fmin + (Fmax-Fmin)*eps;
            break;
        case Medium::SpectrumModel::Binary :
            return ( eps <= S.at(1) ) ? (S.at(0)*eps/S.at(1)) : (S.at(0)*(1.0 + (eps-S.at(1))/(1-S.at(1))));
            break;
        case Medium::SpectrumModel::Uniform : 
            return S.at(0) + (S.at(1)-S.at(0))*eps;
            break;
        case Medium::SpectrumModel::Gaussian :
            return S.at(0) + S.at(1) * sqrt(2) * erfinvf(2*(float)eps-1);
            break;
        case Medium::SpectrumModel::Cauchy :
            return S.at(0) + S.at(0) * tan(CONST_PI * (eps-0.5));
            break;
        case Medium::SpectrumModel::Laplace :
            return S.at(0) + S.at(1) * ( (eps <= 0.5) ? log(2*eps) : -log(2-2*eps) );
            break;
        case Medium::SpectrumModel::Impulse :
            return S.at(0);
            break;
        case Medium::SpectrumModel::Blackbody :
            //Use rejection sampling
            while(true) {
                tmp = roll() * CONST_BMAX;
                if (roll() <= (tmp*tmp/(exp(tmp)-1)/CONST_PLANCKMAX))
                    return tmp/CONST_HK*S.at(0);
            } break;
        default:
        case Medium::SpectrumModel::Rayleigh :
            //Can't sample these
            cerr << "Error sampling spectrum model """ << _xcmodelname(x) << """" << endl;
            exit(-1);
            return 0;
    }
}

double Medium::peak_v() const {
    return _meanSpec(xf, _cf);
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
    return _Ss * _evalSpec(xcs, w, _cs);
}
double Medium::Sa(double w) const {
    return _Sa * _evalSpec(xca, w, _ca);
}
double Medium::Ss() const {
    if (xcs == Medium::SpectrumModel::Constant)
        return HUGE_VAL;
    return _Ss;
}
double Medium::Sa() const {
    if (xca == Medium::SpectrumModel::Constant)
        return HUGE_VAL;
    return _Sa;
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
                return (1+gg*gg - pow((1 - gg * gg) / (1-gg + 2.0 * gg * eps), 2))/2.0/gg;
            else
                return 1-2*eps;
                
        case PhaseFunction::Rayleigh :
            eps = 2*eps - 1;
            t0 = sqrt(4*eps*eps+1)+2*eps;
            return (pow(t0, 2.0/3.0)-1)/pow(t0, 1.0/3.0);
            
        case PhaseFunction::Draine :
            t0 = gg2*(1-gg*gg);
            t1 = gg2*(pow(gg,4)-1);
            t2 = -3*(4*(pow(gg,4) - pow(gg,2)) + t1 * (1+gg*gg));
            t3 = gg*(2*eps - 1);
            t4 = 3*gg*gg * (1 + t3) * gg2 * (2 + gg*gg*(1 + (1 + 2*gg*gg)*t3));
            t5 = t0 * (t1 * t2 + t4*t4) + pow(t1, 3);
            t6 = t0 * 4 * (pow(gg,4)-gg*gg);
            t7 = pow(t5 + sqrt(t5*t5 - pow(t6,3)), 1.0/3.0);
            t8 = 2*(t1 + t6/t7 + t7)/t0;
            t9 = sqrt(6*(1+gg*gg)+t8);
            return gg/2 + (1 - 0.25*pow(sqrt(6*(1+gg*gg)-t8+8*t4/t0/t9)-t9 ,2))/2/gg;
            
        default:
            return 1-2*eps;
    }
}

double Medium::pscatter(double cost, double w) const {
    double gg = g(w); double gg2 = g2(w);
    switch (phase) {
        case PhaseFunction::HenyeyGreenstein :
            return 0.25 / CONST_PI * (1-gg*gg) / pow(1 + gg*gg - 2*gg*cost, 1.5);
        
        case PhaseFunction::Draine :
            return 0.25 / CONST_PI * (1-gg*gg) / pow(1 + gg*gg - 2*gg*cost, 1.5) 
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
    return _sampleSpec(xf, eps, _cf);
}

void Medium::writedbheader(ofstream& OF) const {
    OF << "n,dens,phaseFxn,xca,xcs,xf,g,Ss,Sa,a,tau,FQY";
}

void Medium::writedb(ofstream& OF) const {
    OF << setprecision(8) << n() << ",";
    OF << setprecision(8) << dens() << ",";
    OF << _phasemodelname(phase) << ",";
    OF << _xcmodelname(xca) << ",";
    OF << _xcmodelname(xcs) << ",";
    OF << _xcmodelname(xf) << ",";
    OF << scientific << setprecision(8) << g() << ",";
    OF << scientific << setprecision(8) << Ss() << ",";
    OF << scientific << setprecision(8) << Sa() << ",";
    OF << scientific << setprecision(8) << albedo();
    OF << scientific << setprecision(8) << tau();
    OF << scientific << setprecision(8) << FQY();
}

bool Medium::set(const string& key, const vector<string>& val) {
    if (!key.compare("spectrum-scattering")) {
        _cs = vector<double>();
        _Ss = stod(val.at(0));
        xcs = (val.size()>1) ? (Medium::SpectrumModel)stoul(val.at(1)) : Medium::SpectrumModel::Constant;
        for (unsigned long int i = 2; i < val.size(); i ++)
            _cs.push_back(stod(val.at(i)));
    } else if (!key.compare("spectrum-absorption")) {
        _ca = vector<double>();
        _Sa = stod(val.at(0));
        xca = (val.size()>1) ? (Medium::SpectrumModel)stoul(val.at(1)) : Medium::SpectrumModel::Constant;
        for (unsigned long int i = 2; i < val.size(); i ++)
            _ca.push_back(stod(val.at(i)));
    } else if (!key.compare("spectrum-fluorescence")) {
        _cf = vector<double>();
        _FQY = stod(val.at(0));
        xf = (val.size()>1) ? (Medium::SpectrumModel)stoul(val.at(1)) : Medium::SpectrumModel::Impulse;
        for (unsigned long int i = 2; i < val.size(); i ++)
            _cf.push_back(stod(val.at(i)));
    } else if (!key.compare("scattering-anisotropy")) {
        _g = vector<double>();
        phase = (Medium::PhaseFunction)stoul(val.at(0));
        for (unsigned long int i = 1; i < val.size(); i ++)
            _g.push_back(stod(val.at(i)));
    } else if (!key.compare("lifetime-fluorescence"))
        _tau = stod(val.at(0));
    else if (!key.compare("index-medium"))
        _n = stod(val.at(0));
    else if (!key.compare("density"))
        _dens = stod(val.at(0));
    else
        return false;
    return true;
}

#endif
