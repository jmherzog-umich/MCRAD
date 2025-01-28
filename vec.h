#include <iostream>

#ifndef _VEC_H_
#define _VEC_H_

struct vec {
    public:
        //DATA MEMBERS
        double X;
        double Y;
        double Z;
        
        //Constructor
        vec();
        vec(double x, double y=0, double z=0);
        
        //FUNCTIONS
        vec operator+(const vec& rhs) const;
        vec operator-(const vec& rhs) const;
        vec operator-() const;
        void operator=(const vec& rhs);
        vec operator*(double a) const;
        vec operator/(double a) const;
        vec operator/(const vec& rhs) const;
		
		//Other operators
		void operator+=(const vec& rhs);
        void operator-=(const vec& rhs);
        void operator/=(double r);
        void operator*=(double r);
		
		double dot(const vec& rhs) const;

        //Math
        double norm() const;
        double r() const;
        double r2() const;
        vec perp(double theta) const;
        
        double costheta() const;
        double cosphi() const;
        double costheta2() const;
        double cosphi2() const;
        double theta() const;
        double phi() const;
        
};

std::ostream& operator<<(std::ostream& os, const vec& in);

#endif
