#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>

#ifndef _MEDIUM_H_
#define _MEDIUM_H_

#define CONST_PI 3.1415926535
#define CONST_C  299.8
#define CONST_H  1
#define CONST_EPS 1e-10

using namespace std;

struct Medium {

    //Function definitions for wavelength dependence in cross-sections and g, and 
    enum struct CrossSectionModel { Constant=0, Gaussian=1, Cauchy=2, Exponential=3,
            Polynomial=4, Rayleigh=5, PowerLaw=6, ADTsphere=7, ADTcyl=8, ADTpath=9,
            Blackbody=10, CCM1D=11, Binary=12 };
    enum struct AnisotropyModel { Constant=0 };
    enum struct PhaseFunction { HenyeyGreenstein=0, Rayleigh=1, Draine=2 };

    double _FQY;                //FQY of material
    double _dens;               //Particle density um-3
    double _n;                  //Refractive index
    double _f;                  //Frequency of emitted photons
    vector<double> _g;          //Phase function parameters
    vector<double> _Ss;         //Scattering cross-section parameters
    vector<double> _Sa;         //Absorption cross-section parameters
    double _tau;                //Excited state lifetime
     
    PhaseFunction phase;        //Phase function to use for calculation
    CrossSectionModel xca, xcs; //Cross-section models for absorption and scattering
    
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
    
    Medium();
    
};

string _xcmodelname(Medium::CrossSectionModel x) {
    switch (x) {
        case Medium::CrossSectionModel::Constant: return "Constant"; break;
        case Medium::CrossSectionModel::Gaussian: return "Gaussian"; break;
        case Medium::CrossSectionModel::Cauchy: return "Cauchy"; break;
        case Medium::CrossSectionModel::Exponential: return "Exponential"; break;
        case Medium::CrossSectionModel::Polynomial: return "Polynomial"; break;
        case Medium::CrossSectionModel::Rayleigh: return "Rayleigh"; break;
        case Medium::CrossSectionModel::PowerLaw: return "PowerLaw"; break;
        case Medium::CrossSectionModel::ADTsphere: return "ADTsphere"; break;
        case Medium::CrossSectionModel::ADTcyl: return "ADTcyl"; break;
        case Medium::CrossSectionModel::ADTpath: return "ADTpath"; break;
        case Medium::CrossSectionModel::Blackbody: return "Blackbody"; break;
        case Medium::CrossSectionModel::CCM1D: return "CCM1D"; break;
        case Medium::CrossSectionModel::Binary: return "Binary"; break;
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

double Medium::peak_v() const {
    return _f;
}

void Medium::print_at_f(double f) const {
    cout << "   Medium refractive index: " << n(f) << endl;
    cout << "   Scattering anisotropy <cos theta>: " << g(f) << endl;
    cout << "   Integrated scattering cross-section: " << Ss(f) << " um2" << endl;
    cout << "   Integrated absorption cross-section: " << Sa(f) << " um2" << endl;
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
    cout << "Scattering cross-section model: " << _xcmodelname(xcs) << endl;
    cout << "Medium density: " << dens() << " um-3" << endl;
    cout << "Fluorescence quantum yield: " << FQY() << endl << endl;
}

Medium::Medium() {
    _g = {0.98, 1}; _FQY = 0.5;
    phase = PhaseFunction::HenyeyGreenstein;
    xcs = CrossSectionModel::Constant;
    xca = CrossSectionModel::Binary;
    _f = 550; _Ss = {0.5}; _Sa = {0.1, 700, 0.5};
    _dens = 1e-4; _n = 1.33;
    _tau = 1000;
}

//Actually implement these
double _evalXc(Medium::CrossSectionModel x, double w, const vector<double>& S) {
    switch (x) {
        default:
        case Medium::CrossSectionModel::Constant :
            return S.at(0);
            break;
        case Medium::CrossSectionModel::Gaussian :
            return S.at(0) / S.at(2) / sqrt(2.*CONST_PI) * exp(-pow((w-S.at(1))/S.at(2),2)/2);
            break;
        case Medium::CrossSectionModel::Cauchy :
            return S.at(0) / CONST_PI / S.at(2) / (1 + pow((w-S.at(1))/S.at(2),2));
            break;
        case Medium::CrossSectionModel::Exponential :
            return S.at(0)/2.0/S.at(2) * exp(-abs(w-S.at(1))/S.at(2));
            break;
        case Medium::CrossSectionModel::Binary :
            return (w>S.at(1)) ? S.at(2) : S.at(0);
        case Medium::CrossSectionModel::Polynomial :
        case Medium::CrossSectionModel::Rayleigh :
        case Medium::CrossSectionModel::PowerLaw :
        case Medium::CrossSectionModel::ADTsphere :
        case Medium::CrossSectionModel::ADTcyl :
        case Medium::CrossSectionModel::ADTpath :
        case Medium::CrossSectionModel::Blackbody :
        case Medium::CrossSectionModel::CCM1D :
            return 0;
            break;
    }
    return S.at(0);
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
    return _evalXc(xcs, w, _Ss);
}
double Medium::Sa(double w) const {
    return _evalXc(xca, w, _Sa);
}
double Medium::Ss() const {
    return _Ss.at(0);
}
double Medium::Sa() const {
    return _Sa.at(0);
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
    return Se()*_dens;
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
    return _f;
}

void Medium::writedbheader(ofstream& OF) const {
    OF << "n,dens,phaseFxn,xca,xcs,g,Ss,Sa,a";
}

void Medium::writedb(ofstream& OF) const {
    OF << setprecision(8) << n() << ",";
    OF << setprecision(8) << dens() << ",";
    OF << _phasemodelname(phase) << ",";
    OF << _xcmodelname(xca) << ",";
    OF << _xcmodelname(xcs) << ",";
    OF << scientific << setprecision(8) << g() << ",";
    OF << scientific << setprecision(8) << Ss() << ",";
    OF << scientific << setprecision(8) << Sa() << ",";
    OF << scientific << setprecision(8) << albedo();
}

bool Medium::set(const string& key, const vector<string>& val) {
    if (!key.compare("Ss")) {
        _Ss = vector<double>();
        for (unsigned long int i = 0; i < val.size(); i ++)
            _Ss.push_back(stod(val.at(i)));
    } else if (!key.compare("Sa")) {
        _Sa = vector<double>();
        for (unsigned long int i = 0; i < val.size(); i ++)
            _Sa.push_back(stod(val.at(i)));
    } else if (!key.compare("g")) {
        _g = vector<double>();
        for (unsigned long int i = 0; i < val.size(); i ++)
            _g.push_back(stod(val.at(i)));
    } else if (!key.compare("n"))
        _n = stod(val.at(0));
    else if (!key.compare("dens"))
        _dens = stod(val.at(0));
    else if (!key.compare("FQY"))
        _FQY = stod(val.at(0));
    else if (!key.compare("femit"))
        _f = stod(val.at(0));
    else if (!key.compare("phase"))
        phase = (Medium::PhaseFunction)stoul(val.at(0));
    else if (!key.compare("xca"))
        xca = (Medium::CrossSectionModel)stoul(val.at(0));
    else if (!key.compare("xcs"))
        xcs = (Medium::CrossSectionModel)stoul(val.at(0));
    else
        return false;
    return true;
}

#endif
