#include <cmath>
#include "vec.h"

#define CONST_PI 3.1415926535
#define CONST_EPS 1e-6

#ifndef _PHOTON_H_
#define _PHOTON_H_

using namespace std;

///
//Henyey-Greenstein phase function
///
double p_henyeyGreenstein(double g, const vec& mu1, const vec& mu2) {
    return 0.5 * (1-g*g)/pow(1 + g*g - 2*g* mu1.dot(mu2), 1.5);
}
double cos_henyeyGreenstein(double g, double eps) {
    if (abs(g)>CONST_EPS)
        return (1+g*g-pow((1.0-g*g)/(1-g+2.0*g*eps), 2))/2.0/g;
    else
        return 1-2*eps;
}

///
//Rayleigh phase function
///
double p_rayleigh(const vec& mu1, const vec& mu2) {
    return 0.75 * (1 - pow(mu1.dot(mu2),2));
}
double cos_rayleigh(double eps) {
    eps = 2*eps - 1;      //Eps is cos(theta) so allow negative values too
    double tmp = sqrt(4*eps*eps+1)+2*eps;
    return (pow(tmp, 2.0/3.0)-1)/pow(tmp, 1.0/3.0);
}


///
//Photon class
///
struct Photon {
    public:
        //Phase function definitions
        enum struct PhaseFunction {HenyeyGreenstein=0, Rayleigh=1};
    
        //Data members
        vec x;     //Position
        vec mu;    //Orientation
        double W;  //Weight remaining
        double t;  //Total elapsed lifetime
        double g;  //Average scattering cos(a)
        double f;  //Angular frequency (GHz)
        double S;  //Distance left to travel
        bool isBallistic;   //Has the photon scattered yet?
        bool flipPhase;
        
        static PhaseFunction phase;  //Phase function to use for calculation
        
        //Functions
        //Default (empty) constructor
        Photon() {
            x.X = 0; x.Y = 0; x.Z = 0;
            mu.X = 0; mu.Y = 0; mu.Z = 1;
            W = 1; t = 0; g = 0; f = 0;
            isBallistic = true;
            phase = PhaseFunction::HenyeyGreenstein;
        };
        
        //Constructor for specific photon direction and position
        Photon(const vec& x, const vec& mu, double g=0, double w=1, double t0=0) {
            this->x = x; this->mu = mu; this->g = g;
            W = w; t = t0; isBallistic = true;
        };
        
        //Scatter the photon
        void Scatter(double eps1, double eps2) {
            //Calculate new cos and sin
            double cost;
            switch (phase) {
                case PhaseFunction::HenyeyGreenstein:
                    cost = cos_henyeyGreenstein(g, eps1);
                    break;
                case PhaseFunction::Rayleigh:
                    cost = cos_rayleigh(eps1);
                    break;
                default:
                    cost = 0;
                    break;
            }
            double sint = sqrt(1.0-cost*cost);
            
            //If we're moving very close to +/- Z, the solution simplifies
            double newX, newY, newZ;
            if (abs(mu.Z-1.0) < CONST_EPS) { //Z=1
                newZ = -cost;
                newX = sint*cos(eps2*2.0*CONST_PI);
                newY = sint*sin(eps2*2.0*CONST_PI);
            } else if (abs(mu.Z+1.0) < CONST_EPS) { //Z=-1
                newX = sint*cos(eps2*2.0*CONST_PI);
                newY = sint*sin(eps2*2.0*CONST_PI);
                newZ = cost;
            } else {
                newZ = -sqrt(1-pow(mu.Z,2)) * sint 
                        * cos(eps2*2.0*CONST_PI) + mu.Z*cost;
                newY = sint*(mu.Y*mu.Z*cos(eps2*2.0*CONST_PI)
                        + mu.X*sin(eps2*2.0*CONST_PI))
                        / sqrt(1-pow(mu.Z,2)) + mu.Y*cost;
                newX = sint*(mu.X*mu.Z*cos(eps2*2.0*CONST_PI)
                        - mu.Y*sin(eps2*2.0*CONST_PI)) 
                        / sqrt(1-pow(mu.Z,2)) + mu.X*cost;
            };
            mu = vec(newX, newY, newZ);
            isBallistic = false;
        }
        
        double phi() const {
            return cos(f * t) * (flipPhase ? -1 : 1);
        }
        
        //Math functions
        double muR() const {
            vec xx = x;
            xx.Z = 0;
            return mu.dot(xx / xx.norm());
        }
        
        double intersectR(double R) const {
            //Some edge cases
            if (abs(mu.r2()) < CONST_EPS)
                return -1;
        
            //Some preliminary variables
            double a = mu.X*mu.X + mu.Y*mu.Y;
            double b = (mu.X*x.X + mu.Y*x.Y)/a;
            double c = x.r2() - R*R;
            if (c/a > b*b)
                return -1;
                    
            //Calculate the roots
            double r1 = sqrt(b*b - c/a) - b;
            double r2 = -sqrt(b*b - c/a) - b;
            
            //Return the right root
            if (r2 > 0)
                return r2;
            else
                return r1;
        }
        
        bool operator<(const Photon& rhs) const {
            return this->t < rhs.t;
        }
};

Photon::PhaseFunction Photon::phase = Photon::PhaseFunction::HenyeyGreenstein;

#endif
