#include <cmath>

#include "vec.h"
#include "medium.h"
#include "utility.h"
#include "raypath.h"

#include "photon.h"

#define CONST_PI            3.1415926535897932
#define CONST_C             299.792458
#define CONST_HBAR          1.054571817679489e-22
#define CONST_HK            47.992430733662212
#define CONST_APERY         1.2020569031595942
#define CONST_WIEN          1.5936242600400401
#define CONST_EPS           1e-10
#define CONST_PLANCKMAX     0.6476102378919149

double Photon::Wm = 0.1;

//Functions
//Default (empty) constructor
Photon::Photon() {
    x.X = 0; x.Y = 0; x.Z = 0;
    mu.X = 0; mu.Y = 0; mu.Z = 1;
    W = 1; t = 0; v = 550; n = 0;
    flags = PhotonFlags::isBallistic;
    S = 0;
    pth = nullptr;
};

//Copy constructor
Photon::Photon(const Photon& p) {
    this->x = p.x;
    this->mu = p.mu;
    this->W = p.W;
    this->t = p.t;
    this->flags = p.flags;
    this->v = p.v;
    this->n = p.n;
    this->S = p.S;
    this->pth = p.pth;
}

//Constructor for specific photon direction and position
Photon::Photon(const vec& x, const vec& mu, double w, double t0) {
    this->x = x; this->mu = mu;
    W = w; t = t0; flags = PhotonFlags::isBallistic;
    v = 0; n = 0; S = 0;
    pth = nullptr;
};

void Photon::Reflect() {
    if (pth)
        pth->collide(x, t, W);
}

//Scatter the photon
void Photon::Scatter(double eps1, double eps2, const Medium& p, double eps3) {
    //Increment scattering order
    n++;
    
    //Calculate new cos and sin
    double cost = p.scatter(eps1, v);
    double sint = std::sqrt(1.0-cost*cost);
    
    //If we're moving very close to +/- Z, the solution simplifies
    double newX, newY, newZ;
    if (abs(mu.Z-1.0) < CONST_EPS) { //Z=1
        newZ = cost;
        newX = sint*cos(eps2*2.0*CONST_PI);
        newY = sint*sin(eps2*2.0*CONST_PI);
    } else if (abs(mu.Z+1.0) < CONST_EPS) { //Z=-1
        newX = sint*cos(eps2*2.0*CONST_PI);
        newY = sint*sin(eps2*2.0*CONST_PI);
        newZ = -cost;
    } else {
        newZ = -std::sqrt(1-pow(mu.Z,2)) * sint 
                * std::cos(eps2*2.0*CONST_PI) + mu.Z*cost;
        newY = sint*(mu.Y*mu.Z*std::cos(eps2*2.0*CONST_PI)
                + mu.X*std::sin(eps2*2.0*CONST_PI))
                / std::sqrt(1-std::pow(mu.Z,2)) + mu.Y*cost;
        newX = sint*(mu.X*mu.Z*std::cos(eps2*2.0*CONST_PI)
                - mu.Y*std::sin(eps2*2.0*CONST_PI)) 
                / std::sqrt(1-std::pow(mu.Z,2)) + mu.X*cost;
    };
    mu = vec(newX, newY, newZ);
    flags = (PhotonFlags) (~(int)PhotonFlags::isBallistic & (int)flags);
    
    //Store raypath
    if (pth)
        pth->collide(x, t, W);
        
    //Roll new S
    S = eps3;
}

void Photon::Kill() {
    W = -1;
    S = 0;
}

double Photon::phi() const {
    return std::cos(v * t) * (((int)flags & (flippedPhase())) ? -1 : 1);
}

//Math functions
double Photon::muR() const {
    vec xx = x;
    xx.Z = 0;
    return mu.dot(xx / xx.norm());
}

bool Photon::operator<(const Photon& rhs) const {
    return this->t < rhs.t;
}

bool Photon::operator<(double t) const {
    return this->t < t;
}

Photon& Photon::operator=(const Photon& rhs) {
    x = rhs.x;
    mu = rhs.mu;
    W = rhs.W;
    t = rhs.t;
    flags = rhs.flags;
    v = rhs.v;
    n = rhs.n;
    S = rhs.S;
    pth = rhs.pth;
    return *this;
}

void Photon::storeRayPath(RayPath& path) {
    pth = &path;
    pth->f = v;
    pth->collide(x, t, W);
}

void Photon::printstatus() const {
    std::cerr << "Position: " << x << "   Direction: " << mu << "   Weight: " << W << "   RayPath: " << pth << "   n=" << n << "   t=" << t << "   Flags: " << (int)flags << std::endl;
}

bool Photon::hasPath() const {
    return (pth != nullptr);
}

bool Photon::isFluorescence() const {
    return (int)flags & (int)PhotonFlags::isFluoresce;
}

bool Photon::isBallistic() const {
    return (int)flags & (int)PhotonFlags::isBallistic;
}

bool Photon::flippedPhase() const {
    return (int)flags & (int)PhotonFlags::flippedPhase;
}

void Photon::flipPhase() {
    flags = (PhotonFlags) ( (int)flags ^ (int)PhotonFlags::flippedPhase );
}

void Photon::clearRayPath() {
    pth = nullptr;
}

void Photon::roulette() {
    double eps = util::roll();
    if (eps <= Photon::Wm)
        W /= Photon::Wm;
    else {
        W = 0;
    }
}


