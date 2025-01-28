#include <string>
#include <cmath>
#include <iostream>
#include <vector>

#include "utility.h"

#ifndef _SPECTRUM_H_
#define _SPECTRUM_H_

#define CONST_PI            3.1415926535897932
#define CONST_C             299.792458
#define CONST_HBAR          1.054571817679489e-22
#define CONST_HK            47.992430733662212
#define CONST_APERY         1.2020569031595942
#define CONST_WIEN          1.5936242600400401
#define CONST_EPS           1e-10
#define CONST_PLANCKMAX     0.6476102378919149

using namespace std;

class Spectrum {
    
    public:
        enum struct SpectrumModel { Constant=0, Gaussian=1, Cauchy=2, Laplace=3, Impulse=4, Rayleigh=5, Blackbody=6, DoubleGauss=9, Binary=12, Uniform=13 };
    
        Spectrum();
        Spectrum(SpectrumModel m, vector<double> vals);
        
        double sample() const;
        double evaluate(double x) const;
        double peak() const;
        double max() const;
        double norm() const;
        double width() const;
        bool peakLabel() const;
        
        string name() const;
        string paramstring() const;
        
        static void setFreqLimits(double f1, double f2);
        
    private:
        vector<double> S;
        SpectrumModel type;
        static double Fmin, Fmax;
};

double Spectrum::Fmin = 0;
double Spectrum::Fmax = 1000;

void Spectrum::setFreqLimits(double f1, double f2) {
    Fmin = f1;
    Fmax = f2;
}

Spectrum::Spectrum() {
    type = SpectrumModel::Constant;
}

Spectrum::Spectrum(SpectrumModel m, vector<double> vals) {
    type = m;
    S = vals;
}

bool Spectrum::peakLabel() const {
    switch (type) {
        default:
        case SpectrumModel::Gaussian :
        case SpectrumModel::Cauchy :
        case SpectrumModel::Laplace :
        case SpectrumModel::Blackbody :
        case SpectrumModel::Impulse :
            return true;
            break;
        case SpectrumModel::Uniform :
        case SpectrumModel::Constant : 
        case SpectrumModel::Rayleigh :
        case SpectrumModel::Binary :
        case SpectrumModel::DoubleGauss :
            return false;
            break;
    }
}

string Spectrum::name() const {
switch(type) {
        case SpectrumModel::Constant: return "Constant"; break;
        case SpectrumModel::Gaussian: return "Gaussian"; break;
        case SpectrumModel::DoubleGauss: return "Two-Gaussian"; break;
        case SpectrumModel::Cauchy: return "Cauchy"; break;
        case SpectrumModel::Laplace: return "Laplace"; break;
        case SpectrumModel::Binary: return "Binary"; break;
        case SpectrumModel::Uniform: return "Uniform"; break;
        case SpectrumModel::Impulse: return "Impulse"; break;
        case SpectrumModel::Blackbody: return "Blackbody"; break;
        case SpectrumModel::Rayleigh: return "Rayleigh"; break;
        default: return "Other"; break;
    }
}

string Spectrum::paramstring() const {
    return "Frequency: " + to_string(peak()) + "   +/-   " + to_string(width()) + " THz";
}

double Spectrum::norm() const {
    switch(type) {
        case SpectrumModel::Constant :
            return Fmax-Fmin;
            break;
        case SpectrumModel::Uniform :
            return S.at(1)-S.at(0);
            break;
        case SpectrumModel::Rayleigh :
            return (pow(Fmax,5)-pow(Fmin,5))/(5*pow(S.at(0), 4));
            break;
        case SpectrumModel::Gaussian :
            return  S.at(1) * sqrt(2.*CONST_PI);
            break;
        case SpectrumModel::Cauchy :
            return CONST_PI * S.at(1);
            break;
        case SpectrumModel::Laplace :
            return 2 * S.at(1);
            break;
        case SpectrumModel::DoubleGauss :
            return  sqrt(2.*CONST_PI) * (S.at(1) * S.at(2) + (1-S.at(2))* S.at(4));
            break;
        case SpectrumModel::Binary :
            return (S.at(0)-Fmin)*S.at(1) + (Fmax-S.at(0))*S.at(2);
        case SpectrumModel::Impulse :
        case SpectrumModel::Blackbody :
            return 2 * CONST_APERY * S.at(0) / CONST_HK / CONST_PLANCKMAX;
        default:
            return 1;
    }
}

double Spectrum::evaluate(double w) const {
    switch (type) {
        default:
        case SpectrumModel::Constant :
            return 1.0;
            break;
        case SpectrumModel::Gaussian :
            return exp(-pow((w-S.at(0))/S.at(1),2)/2);
            break;
        case SpectrumModel::Cauchy :
            return 1.0 / (1 + pow((w-S.at(0))/S.at(1),2));
            break;
        case SpectrumModel::Laplace :
            return exp(-abs(w-S.at(0))/S.at(1));
            break;
        case SpectrumModel::Binary :
            return (w>S.at(0)) ? S.at(2) : S.at(1);
            break;
        case SpectrumModel::Impulse :
            return (abs(w-S.at(1)) < CONST_EPS) ? HUGE_VAL : 0.0;
            break;
        case SpectrumModel::Blackbody :
            return pow(CONST_HK*w/S.at(0),2) / (exp(CONST_HK*w/S.at(0))-1) / CONST_PLANCKMAX;
            break;
        case SpectrumModel::DoubleGauss :
            return exp(-pow((w-S.at(0))/S.at(1),2)/2)*S.at(2) + (1-S.at(2))*exp(-pow((w-S.at(3))/S.at(4),2)/2);
            break;
        case SpectrumModel::Uniform : 
            if (w <= S.at(1) and w >= S.at(0))
                return 1.0;
            else
                return 0;
            break;
        case SpectrumModel::Rayleigh :
            return pow(w/S.at(0), 4);
    }
    return S.at(0);
}

double Spectrum::max() const {
    switch (type) {
        default:
        case SpectrumModel::Constant :
        case SpectrumModel::Gaussian :
        case SpectrumModel::Cauchy :
        case SpectrumModel::Laplace :
        case SpectrumModel::Uniform : 
        case SpectrumModel::Rayleigh :
        case SpectrumModel::Blackbody :
            return 1.0;
            break;
        case SpectrumModel::Binary :
            return (S.at(2)>S.at(1)) ? S.at(2) : S.at(1);
            break;
        case SpectrumModel::Impulse :
            return HUGE_VAL;
            break;
        case SpectrumModel::DoubleGauss :
            return evaluate(peak());
            break;
    }
    return S.at(0);
}

//Calculates the mean frequency of the spectrum model, or for non-convergent models the specified/target wavelength
double Spectrum::peak() const {
    switch (type) {
        case SpectrumModel::Constant :
            return (Fmax+Fmin)/2.0;
            break;
        case SpectrumModel::Binary :
            return (S.at(1) > S.at(2)) ? (Fmin+S.at(0))/2 : (Fmax+S.at(0))/2;
            break;
        case SpectrumModel::Uniform :
            return (S.at(0) + S.at(1))/2.0;
            break;
        case SpectrumModel::Blackbody :
            return S.at(0) / CONST_HK * CONST_WIEN;
            break;
        case SpectrumModel::DoubleGauss :
            return S.at(2)>0.5 ? S.at(0) : S.at(3);
            break;
        case SpectrumModel::Rayleigh :
        case SpectrumModel::Gaussian :
        case SpectrumModel::Cauchy :
        case SpectrumModel::Laplace :
        case SpectrumModel::Impulse :
            return S.at(0);
            break;
        default:
            cerr << "Error: no mean available for spectrum model """ << name() << """" << endl;
            exit(-1);
            return 0;
    }
}

//Calculates the mean frequency of the spectrum model, or for non-convergent models the specified/target wavelength
double Spectrum::width() const {
    switch (type) {
        case SpectrumModel::Rayleigh :
        case SpectrumModel::Blackbody :
        case SpectrumModel::Binary :
        case SpectrumModel::Constant :
            return HUGE_VAL;
            break;
        case SpectrumModel::Uniform :
            return (S.at(1) - S.at(0))/2.0;
            break;
        case SpectrumModel::Gaussian :
        case SpectrumModel::Cauchy :
        case SpectrumModel::Laplace :
            return S.at(1);
            break;
        case SpectrumModel::Impulse :
            return 0;
            break;
        case SpectrumModel::DoubleGauss :
            return sqrt(S.at(1)*S.at(1)*S.at(2) + S.at(4)*S.at(4)*(1-S.at(2)));
            break;
        default:
            cerr << "Error: no width available for spectrum model """ << name() << """" << endl;
            exit(-1);
            return 0;
    }
}

//Randomly samples the spectrum model
double Spectrum::sample() const {
    double tmp, eps = util::roll();
    switch (type) {
        case SpectrumModel::Constant :
            return Fmin + (Fmax-Fmin)*eps;
            break;
        case SpectrumModel::Binary :
            tmp = S.at(1)/(S.at(1)+S.at(2));
            return (util::roll() <= tmp) ? (Fmin + (S.at(0)-Fmin)*eps) : (S.at(0) + (Fmax-S.at(0))*eps);
            break;
        case SpectrumModel::Uniform : 
            return S.at(0) + (S.at(1)-S.at(0))*eps;
            break;
        case SpectrumModel::Gaussian :
            return S.at(0) + S.at(1) * sqrt(2) * util::erfinvf(2*(float)eps-1);
            break;
        case SpectrumModel::Cauchy :
            return S.at(0) + S.at(0) * tan(CONST_PI * (eps-0.5));
            break;
        case SpectrumModel::Laplace :
            return S.at(0) + S.at(1) * ( (eps <= 0.5) ? log(2*eps) : -log(2-2*eps) );
            break;
        case SpectrumModel::Impulse :
            return S.at(0);
            break;
        case SpectrumModel::DoubleGauss :
            if (util::roll() <= S.at(2))
                return S.at(0) + S.at(1) * sqrt(2) * util::erfinvf(2*(float)eps-1);
            else
                return S.at(3) + S.at(4) * sqrt(2) * util::erfinvf(2*(float)eps-1);
            break;
        case SpectrumModel::Blackbody :
            //Use rejection sampling
            while(true) {
                tmp = (Fmin + (Fmax-Fmin)*util::roll()) * CONST_HK / S.at(0);
                if (util::roll() <= (tmp*tmp/(exp(tmp)-1)/CONST_PLANCKMAX))
                    return tmp/CONST_HK*S.at(0);
            } break;
        default:
        case SpectrumModel::Rayleigh :
            //Can't sample these
            cerr << "Error sampling spectrum model """ << name() << """" << endl;
            exit(-1);
            return 0;
    }
}

#endif
