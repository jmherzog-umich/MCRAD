#include <cmath>
#include <iostream>

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
        vec operator-() const {
            return vec(-X, -Y, -Z);
        };
        void operator=(const vec& rhs) {
            this->X = rhs.X;
            this->Y = rhs.Y;
            this->Z = rhs.Z;
        };
        vec operator*(double a) const {
            return vec(this->X*a, this->Y*a, this->Z*a);
        }
        vec operator/(double a) const {
            return vec(this->X/a, this->Y/a, this->Z/a);
        }  
        vec operator/(const vec& rhs) const {
            return vec(this->X/rhs.X, this->Y/rhs.Y, this->Z/rhs.Z);
        }        
		
		//Other operators
		void operator+=(const vec& rhs) {
            this->X += rhs.X;
            this->Y += rhs.Y;
            this->Z += rhs.Z;
        };
        void operator-=(const vec& rhs) {
            this->X -= rhs.X;
            this->Y -= rhs.Y;
            this->Z -= rhs.Z;
        };
        void operator/=(double r) {
            this->X /= r;
            this->Y /= r;
            this->Z /= r;
        };
        void operator*=(double r) {
            this->X *= r;
            this->Y *= r;
            this->Z *= r;
        };
		
		double dot(const vec& rhs) const {
			return this->X * rhs.X + this->Y * rhs.Y + this->Z * rhs.Z;
		}

        //Math
        double norm() const {
            return hypot(X,Y,Z);
        };
        
        double r() const {
            return hypot(X,Y);
        };
        
        double r2() const {
            return X*X + Y*Y;
        };
        
        vec perp(double theta) const {
            double cost = cos(theta);
            double sint = pow(1.0 - cost*cost, 0.5);
            double d, den;
            if ((abs(Z) >= abs(X)) and (abs(Z) >= abs(Y))) {
                d = X * cost + Y * sint;
                den = hypot(Z,d);
                return vec(Z*cost / den, Z*sint / den, -d / den);
            } else if ((abs(Y) >= abs(X)) and (abs(Y) >= abs(Z))) {
                d = X * cost + Z * sint;
                den = hypot(Y,d);
                return vec(Y*cost / den, -d / den, Y*sint / den);
            } else {
                d = Z * cost + Y * sint;
                den = hypot(X,d);
                return vec(-d / den, X*sint / den, X*cost / den);
            }
        };
        
        double costheta() const {
            return X/hypot(X,Y);
        };
        
        double cosphi() const {
            double r2 = X*X+Y*Y;
            return sqrt(r2/(r2 + Z*Z));
        };
        
        double costheta2() const {
            double x2 = X*X;
            return x2/(x2+Y*Y);
        };
        
        double cosphi2() const {
            double r2 = X*X+Y*Y;
            return r2/(r2 + Z*Z);
        };
        
        double theta() const {
            return atan2(Y, X);
        };
        
        double phi() const {
            return atan2(Z, r());
        };
        
};

ostream& operator<<(ostream& os, const vec& in) {
   os << "<" << in.X << ", " << in.Y << ", " << in.Z << ">";
   return os;
}

#endif
