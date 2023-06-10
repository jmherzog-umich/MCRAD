#include <cmath>
#include "vec.h"

#define CONST_PI 3.1415926535

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
    if (abs(g)>1e-5)
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
        bool isBallistic;   //Has the photon scattered yet?
        
        static PhaseFunction phase;  //Phase function to use for calculation
        
        //Functions
        //Default (empty) constructor
        Photon() {
            x.X = 0; x.Y = 0; x.Z = 0;
            mu.X = 0; mu.Y = 0; mu.Z = 1;
            W = 1; t = 0; g = 0;
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
            }
            double sint = sqrt(1.0-cost*cost);
            
            //If we're moving very close to +/- Z, the solution simplifies
            double newX, newY, newZ;
            if (abs(mu.Z-1.0) < 1e-6) { //Z=1
                newZ = -cost;
                newX = sint*cos(eps2*2.0*CONST_PI);
                newY = sint*sin(eps2*2.0*CONST_PI);
            } else if (abs(mu.Z+1.0) < 1e-6) { //Z=-1
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
};

Photon::PhaseFunction Photon::phase = Photon::PhaseFunction::HenyeyGreenstein;

#endif
