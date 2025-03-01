#include "vec.h"
#include "medium.h"
#include "utility.h"
#include "raypath.h"

#ifndef _PHOTON_H_
#define _PHOTON_H_

///
//Photon class
///
struct Photon {
    
    public:
        
        //Static members
        static double Wm;
        
        enum struct PhotonFlags {
            isBallistic      = 0b00000001,      //Whether the photon has collided yet
            flippedPhase     = 0b00000010,      //Whether the phase should be flipped because of reflection from a higher index medium
            isFluoresce      = 0b00000100,      //Whether this is a fluorescence emission packet
        };
    
        //Data members
        vec x;     //Position
        vec mu;    //Orientation
        double W;  //Weight remaining
        double t;  //Total elapsed lifetime
        double v;  //Angular frequency (THz)
        double S;  //Distance left to travel
        int n;     //Scattering order
        
        //Flags
        PhotonFlags flags = PhotonFlags::isBallistic;
        
        Photon();
        Photon(const Photon& p);
        Photon(const vec& x, const vec& mu, double w=1, double t0=0);
        
        //Helper functions
        void Scatter(double eps1, double eps2, const Medium& p, double eps3);
        void Reflect();
        
        void storeRayPath(RayPath& path);
        void clearRayPath();
        
        void printstatus() const;
        bool hasPath() const;
        
        bool isFluorescence() const;
        bool isBallistic() const;
        bool flippedPhase() const;
        
        void flipPhase();
        
        void roulette();
        void Kill();
        
        double phi() const;
        double muR() const;
        bool operator<(const Photon& rhs) const;
        bool operator<(double rhs) const;
        Photon& operator=(const Photon& rhs);
        
        private:
            RayPath* pth;
};

#endif
