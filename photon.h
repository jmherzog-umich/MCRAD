#include <cmath>

#include "vec.h"
#include "medium.h"
#include "utility.h"

#define CONST_PI 3.1415926535
#define CONST_EPS 1e-10

#ifndef _PHOTON_H_
#define _PHOTON_H_

using namespace std;

///
//Photon class
///
struct Photon {
    
    public:
    
        //Data members
        vec x;     //Position
        vec mu;    //Orientation
        double W;  //Weight remaining
        double t;  //Total elapsed lifetime
        double v;  //Angular frequency (GHz)
        double S;  //Distance left to travel
        int n;     //Scattering order
        bool isBallistic;   //Has the photon scattered yet?
        bool flipPhase;
        
        Photon();
        Photon(const Photon& p);
        Photon(const vec& x, const vec& mu, double w, double t0);
        
        //Helper functions
        void Scatter(double eps1, double eps2, const Medium& p);
        
        double phi() const;
        double muR() const;
        double intersectR(double R) const;
        double intersectXY(double R, bool& xy) const;
        bool operator<(const Photon& rhs) const;  
};

//Functions
//Default (empty) constructor
Photon::Photon() {
    x.X = 0; x.Y = 0; x.Z = 0;
    mu.X = 0; mu.Y = 0; mu.Z = 1;
    W = 1; t = 0; v = 550; n = 0;
    isBallistic = true; flipPhase = false;
    S = 0;
};

//Copy constructor
Photon::Photon(const Photon& p) {
    this->x = p.x;
    this->mu = p.mu;
    this->W = p.W;
    this->t = p.t;
    this->isBallistic = p.isBallistic;
    this->v = p.v;
    this->n = p.n;
    this->S = p.S;
    this->flipPhase = p.flipPhase;
}

//Constructor for specific photon direction and position
Photon::Photon(const vec& x, const vec& mu, double w=1, double t0=0) {
    this->x = x; this->mu = mu;
    W = w; t = t0; isBallistic = true;
    v = 0; n = 0; flipPhase = false;
    S = 0;
};

//Scatter the photon
void Photon::Scatter(double eps1, double eps2, const Medium& p) {
    //Increment scattering order
    n++;
    
    //Calculate new cos and sin
    double cost = p.scatter(eps1, v);
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

double Photon::phi() const {
    return cos(v * t) * (flipPhase ? -1 : 1);
}

//Math functions
double Photon::muR() const {
    vec xx = x;
    xx.Z = 0;
    return mu.dot(xx / xx.norm());
}

double Photon::intersectXY(double R, bool& xy) const {
    //Allocate variables
    double sx, sy;
    
    //X direction
    if (abs(mu.X) < CONST_EPS)
        sx = -1;
    else if (mu.X > 0)
        sx = (R-x.X)/mu.X;
    else
        sx = (-R-x.X)/mu.X;
    
    //Y direction
    if (abs(mu.Y) < CONST_EPS)
        sy = -1;
    else if (mu.Y > 0)
        sy = (R-x.Y)/mu.Y;
    else
        sy = (-R-x.Y)/mu.Y;
        
    //Report the smallest value greater than zero
    if ((sx > 0) and (sy > 0))
        if (sx < sy)
            xy = true;
        else
            xy = false;
    else if (sx > 0)
        xy = true;
    else if (sy > 0)
        xy = false;
    else
        xy = false;
    if (xy)
        return sx;
    else
        return sy;
}

double Photon::intersectR(double R) const {
    //Some edge cases
    if (abs(mu.r2()) < CONST_EPS)
        return -1;

    //Some preliminary variables
    double a = mu.X*mu.X + mu.Y*mu.Y;
    double b = (mu.X*x.X + mu.Y*x.Y)/a;
    double c = (x.r2() - R*R)/a;
    if (c/a > b*b)
        return -1;
            
    //Calculate the roots
    double r1 = -sqrt(b*b - c) - b;
    double r2 = sqrt(b*b - c) - b;
    
    //If r1 is negative, return r2 regardless; or return r2 if its smaller than r1
    if ((r1 < 0) or (r2 < r1))
        return r2;
    else
        return r1;
}

bool Photon::operator<(const Photon& rhs) const {
    return this->t < rhs.t;
}

#endif
