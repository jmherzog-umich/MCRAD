#include <string>
#include <vector>

#ifndef _SPECTRUM_H_
#define _SPECTRUM_H_

class Spectrum {
    
    public:
        enum struct SpectrumModel { Constant=0, Gaussian=1, Cauchy=2, Laplace=3, Impulse=4, Rayleigh=5, Blackbody=6, DoubleGauss=9, Binary=12, Uniform=13 };
    
        Spectrum();
        Spectrum(SpectrumModel m, std::vector<double> vals);
        
        double sample() const;
        double evaluate(double x) const;
        double peak() const;
        double max() const;
        double norm() const;
        double width() const;
        bool peakLabel() const;
        
        std::string name() const;
        std::string paramstring() const;
        
        static void setFreqLimits(double f1, double f2);
        
    private:
        std::vector<double> S;
        SpectrumModel type;
        static double Fmin, Fmax;
};

#endif
