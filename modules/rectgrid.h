#include <vector>
#include <string>
#include <iostream>

#include "../vec.h"
#include "../grid.h"

#ifndef __RECTGRID_H
#define __RECTGRID_H

struct RectGrid : public Grid {

    //Constructors
    RectGrid();
    RectGrid(double x, double y, double z, unsigned int nx=0, unsigned int ny=0, unsigned int nz=0);

    //Get first index of each dimension that corresponds to inside the grid
    virtual unsigned long int geti0() const;
    virtual unsigned long int getj0() const;
    virtual unsigned long int getk0() const;
    virtual long unsigned int ntotal() const;

    //Data accessors
    virtual long unsigned int ind1(const vec& x) const;
    virtual long unsigned int ind2(const vec& x) const;
    virtual long unsigned int ind3(const vec& x) const;
    
    virtual long unsigned int ind2sub1(long unsigned int i) const;
    virtual long unsigned int ind2sub2(long unsigned int i) const;
    virtual long unsigned int ind2sub3(long unsigned int i) const;
    
    virtual long unsigned int sub2ind(long unsigned int i, long unsigned int j, long unsigned int k) const;

    //Geometry operations on cells
    virtual double Volume(long unsigned int i, long unsigned int j, long unsigned int k) const;
    virtual vec rand(long unsigned int i, long unsigned int j, long unsigned int k) const;
    virtual bool contains(const vec& x) const;
    virtual double intersect(const vec& x, const vec& mu) const;
    
    //Print grid specific info
    virtual void printOutside(std::ostream& oout, int id) const;
    virtual void printGrid(std::ostream& oout) const;
    virtual void printKeys(std::ostream& oout) const;

    //Collision checks
    //virtual void collideFront(const vec& x, const vec& mu, double& ds, int& reflect) const;
    //virtual void collideBack(const vec& x, const vec& mu, double& ds, int& reflect) const;
    //virtual void collideSide(const vec& x, const vec& mu, double& ds, int& reflect) const;
    //virtual void boundXY(vec& x) const;
    //virtual void boundZ(vec& x) const;
    virtual vec getNormal(int reflect, const vec& x, const vec& mu, vec& mur);
    
};

#endif
