#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>

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
            Polynomial=4, Rayleigh=5, PowerLaw=6, ADTsphere=6, ADTcyl=7, ADTpath=8,
            Blackbody=9, CCM1D=10 };
    enum struct AnisotropyModel { Constant=0 };
    enum struct PhaseFunction { HenyeyGreenstein=0, Rayleigh=1, Draine=2 };

    double _FQY;                //FQY of material
    double _dens;               //Particle density um-3
    double _n;                  //Refractive index
    double _f;                  //Frequency at which parameters are given
    double _g;                  //Average scattering cos(a)
    double _g2;                 //2nd scattering parameter
    double _Ss;                 //Scattering cross-section (peak)
    double _Sa;                 //Absorption cross-section (peak)
    double _cs1, _cs2, _cs3;    //Scattering cross-section model parameters
    double _ca1, _ca2, _ca3;    //Absorption cross-section model parameters
    double _tau;                //Excited state lifetime
     
    PhaseFunction phase;        //Phase function to use for calculation
    CrossSectionModel xc;       //Phase function to use for calculation
    
    void print() const;
    void writedbheader(ofstream& OF) const;
    void writedb(ofstream& OF) const;
    
    bool set(const string& key, const string& val);
    
    double dens() const;
    double FQY() const;
    double g(double l = 0) const;
    double g2(double l = 0) const;
    double tau() const;
    
    double Ss(double l = 0) const;
    double Sa(double l = 0) const;
    double Se(double l = 0) const;
    
    double ks(double l = 0) const;
    double ka(double l = 0) const;
    double ke(double l = 0) const;
    
    double ws(double l = 0) const;
    double wa(double l = 0) const;
    double we(double l = 0) const;
    
    double ls(double l = 0) const;
    double la(double l = 0) const;
    double le(double l = 0) const;
    
    double albedo(double l = 0) const;
    double n(double l = 0) const;
    
    double scatter(double eps, double l = 0) const;
    double pscatter(double cost, double l = 0) const;

    double emit_tau(double eps, double l = 0) const;
    double emit_v(double eps, double l = 0) const;
    
    Medium();
    
};

void Medium::print() const {
    cout << "Phase function: ";
    switch (phase) {
        case Medium::PhaseFunction::HenyeyGreenstein: cout << "Henyey-Greenstein" << endl; break;
        case Medium::PhaseFunction::Rayleigh: cout << "Rayleigh" << endl; break;
        case Medium::PhaseFunction::Draine: cout << "Draine" << endl; break;
        default: cout << "Other" << endl; break;
    }
    
    cout << "Medium refractive index: " << n() << endl;
    cout << "Medium density: " << dens() << " um-3" << endl << endl;
    
    cout << "Scattering anisotropy <cos theta>: " << g() << endl;
    cout << "Secondary scattering parameter (alpha): " << g2() << endl;
    cout << "Scattering cross-section: " << Ss() << " um2" << endl;
    cout << "Absorption cross-section: " << Sa() << " um2" << endl;
    cout << "Albedo: " << albedo() << endl;
    cout << "Fluorescence quantum yield: " << FQY() << endl << endl;
    
    //Print MFPs
    cout << "Scattering mean free path: " << ls() << " um" << endl;
    cout << "Scattering rate: " << ws() << " THz" << endl;
    cout << "Absorption mean free path: " << la() << " um" << endl;
    cout << "Absorption rate: " << wa() << " THz" << endl;
    cout << "Collision mean free path: " << le() << " um" << endl;
    cout << "Collision rate: " << we() << " THz" << endl << endl;
}

Medium::Medium() {
    _g = 0.98; _g2 = 1; _FQY = 0.5;
    phase = PhaseFunction::HenyeyGreenstein;
    xc = CrossSectionModel::Constant;
    _f = 550; _Ss = 0.5; _Sa = 0.5;
    _dens = 1e-4; _n = 1.33;
    _cs1 = 0; _cs2 = 0; _cs3 = 0;
    _ca1 = 0; _ca2 = 0; _ca3 = 0;
    _tau = 1000;
}

//Actually implement these
double Medium::g(double l) const {
    return _g;
}
double Medium::g2(double l) const {
    return _g2;
}
double Medium::Ss(double l) const {
    return _Ss;
}
double Medium::Sa(double l) const {
    return _Sa;
}
double Medium::dens() const {
    return _dens;
}
double Medium::n(double l) const {
    return _n;
}
double Medium::FQY() const {
    return _FQY;
}
double Medium::tau() const {
    return _tau;
}

//These are trivial
double Medium::Se(double l) const {
    return Sa(l) + Ss(l);
}
double Medium::ks(double l) const {
    return Ss(l)*_dens;
}
double Medium::ka(double l) const {
    return Sa(l)*_dens;
}
double Medium::ke(double l) const {
    return Se()*_dens;
}
double Medium::ls(double l) const {
    return 1/ks(l);
}
double Medium::la(double l) const {
    return 1/ka(l);
}
double Medium::le(double l) const {
    return 1/ke(l);
}
double Medium::ws(double l) const {
    return ka(l)*CONST_C/n(l);
}
double Medium::wa(double l) const {
    return ka(l)*CONST_C/n(l);
}
double Medium::we(double l) const {
    return ke(l)*CONST_C/n(l);
}
double Medium::albedo(double l) const {
    return Ss(l) / (Ss(l) + Sa(l));
}

double Medium::scatter(double eps, double l) const {
    double gg = g(l); double gg2 = g2(l);
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

double Medium::pscatter(double cost, double l) const {
    double gg = g(l); double gg2 = g2(l);
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

double Medium::emit_tau(double eps, double l) const {
    return _tau * eps;
}

double Medium::emit_v(double eps, double l) const {
    return _f;
}

void Medium::writedbheader(ofstream& OF) const {
    OF << "n,dens,phaseFxn,g,Ss,Sa,a";
}

void Medium::writedb(ofstream& OF) const {
    OF << setprecision(8) << n() << ",";
    OF << setprecision(8) << dens() << ",";
    switch (phase) {
        case PhaseFunction::HenyeyGreenstein : 
            OF << "Henyey-Greenstein,";
            break;
        case PhaseFunction::Rayleigh :
            OF << "Rayleigh,";
            break;
        case PhaseFunction::Draine :
            OF << "Draine,";
            break;
    }
    OF << scientific << setprecision(8) << g() << ",";
    OF << scientific << setprecision(8) << Ss() << ",";
    OF << scientific << setprecision(8) << Sa() << ",";
    OF << scientific << setprecision(8) << albedo();
}

bool Medium::set(const string& key, const string& val) {
    if (!key.compare("Ss"))
        _Ss = stod(val);
    else if (!key.compare("Sa"))
        _Sa = stod(val);
    else if (!key.compare("g"))
        _g = stod(val);
    else if (!key.compare("g2"))
        _g2 = stod(val);
    else if (!key.compare("n"))
        _n = stod(val);
    else if (!key.compare("dens"))
        _dens = stod(val);
    else if (!key.compare("FQY"))
        _FQY = stod(val);
    else if (!key.compare("phase"))
        phase = (Medium::PhaseFunction)stoul(val);
    else
        return false;
    return true;
}

#endif
