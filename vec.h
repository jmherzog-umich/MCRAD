#include <cmath>

using namespace std;

#ifndef _VEC_H_
#define _VEC_H_

struct vec {
    public:
        //DATA MEMBERS
        double X;
        double Y;
        double Z;
        
        //Constructor
        vec() {
            X = 0; Y = 0; Z = 0;
        };
        vec(double x, double y=0, double z=0) {
            X = x; Y = y; Z = z;        
        };
        
        //FUNCTIONS
        vec operator+(const vec& rhs) const {
            return vec(this->X+rhs.X, this->Y+rhs.Y, this->Z+rhs.Z);
        };
        vec operator-(const vec& rhs) const {
            return vec(this->X-rhs.X, this->Y-rhs.Y, this->Z-rhs.Z);
        };
        void operator=(const vec& rhs) {
            this->X = rhs.X;
            this->Y = rhs.Y;
            this->Z = rhs.Z;
        };
        vec operator*(double a) const {
            return vec(this->X*a, this->Y*a, this->Z*a);
        }
        vec operator/(const vec& rhs) const {
            return vec(this->X/rhs.X, this->Y/rhs.Y, this->Z/rhs.Z);
        }        
		
		double dot(const vec& rhs) const {
			return this->X * rhs.X + this->Y * rhs.Y + this->Z * rhs.Z;
		}

        //Math
        double norm() const {
            return sqrt(X*X + Y*Y + Z*Z);
        };
        double r() const {
            return sqrt(X*X + Y*Y);
        };
};

#endif
