#include <cmath>
#include <iostream>

#include "vec.h"

vec::vec() {
    X = 0; Y = 0; Z = 0;
}

vec::vec(double x, double y, double z) {
    X = x; Y = y; Z = z;        
}

vec vec::operator+(const vec& rhs) const {
    return vec(this->X+rhs.X, this->Y+rhs.Y, this->Z+rhs.Z);
}

vec vec::operator-(const vec& rhs) const {
    return vec(this->X-rhs.X, this->Y-rhs.Y, this->Z-rhs.Z);
}

vec vec::operator-() const {
    return vec(-X, -Y, -Z);
}

void vec::operator=(const vec& rhs) {
    this->X = rhs.X;
    this->Y = rhs.Y;
    this->Z = rhs.Z;
}

vec vec::operator*(double a) const {
    return vec(this->X*a, this->Y*a, this->Z*a);
}

vec vec::operator/(double a) const {
    return vec(this->X/a, this->Y/a, this->Z/a);
}  

vec vec::operator/(const vec& rhs) const {
    return vec(this->X/rhs.X, this->Y/rhs.Y, this->Z/rhs.Z);
}        

void vec::operator+=(const vec& rhs) {
    this->X += rhs.X;
    this->Y += rhs.Y;
    this->Z += rhs.Z;
}

void vec::operator-=(const vec& rhs) {
    this->X -= rhs.X;
    this->Y -= rhs.Y;
    this->Z -= rhs.Z;
}

void vec::operator/=(double r) {
    this->X /= r;
    this->Y /= r;
    this->Z /= r;
}

void vec::operator*=(double r) {
    this->X *= r;
    this->Y *= r;
    this->Z *= r;
}

double vec::dot(const vec& rhs) const {
	return this->X * rhs.X + this->Y * rhs.Y + this->Z * rhs.Z;
}

double vec::norm() const {
    return std::hypot(X,Y,Z);
}

double vec::r() const {
    return std::hypot(X,Y);
}

double vec::r2() const {
    return X*X + Y*Y;
}

vec vec::perp(double theta) const {
    double cost = std::cos(theta);
    double sint = std::pow(1.0 - cost*cost, 0.5);
    double d, den;
    if ((std::abs(Z) >= std::abs(X)) and (std::abs(Z) >= std::abs(Y))) {
        d = X * cost + Y * sint;
        den = std::hypot(Z,d);
        return vec(Z*cost / den, Z*sint / den, -d / den);
    } else if ((std::abs(Y) >= std::abs(X)) and (std::abs(Y) >= std::abs(Z))) {
        d = X * cost + Z * sint;
        den = std::hypot(Y,d);
        return vec(Y*cost / den, -d / den, Y*sint / den);
    } else {
        d = Z * cost + Y * sint;
        den = std::hypot(X,d);
        return vec(-d / den, X*sint / den, X*cost / den);
    }
}

double vec::costheta() const {
    return X/std::hypot(X,Y);
}

double vec::cosphi() const {
    double r2 = X*X+Y*Y;
    return std::sqrt(r2/(r2 + Z*Z));
}

double vec::costheta2() const {
    double x2 = X*X;
    return x2/(x2+Y*Y);
}

double vec::cosphi2() const {
    double r2 = X*X+Y*Y;
    return r2/(r2 + Z*Z);
}

double vec::theta() const {
    return std::atan2(Y, X);
}

double vec::phi() const {
    return std::atan2(Z, r());
}

std::ostream& operator<<(std::ostream& os, const vec& in) {
    os << "<" << in.X << ", " << in.Y << ", " << in.Z << ">";
    return os;
}

