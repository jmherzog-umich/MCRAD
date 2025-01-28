#include <vector>
#include <string>
#include <iostream>

#include "vec.h"

#ifndef __GRID_H
#define __GRID_H

struct Grid {

    //Grid iterator
    struct iterator {
        //private:
            unsigned long int ix, iy, iz;
            const Grid& g;
            bool DONE;
        public:
            iterator(const Grid& p);
            iterator& operator++(int);
            unsigned long int index() const;
            bool end() const;
    };
  
    //Data  members
    std::vector<std::vector<double> > m;    //Data members, use a map to store multiple things
    std::vector<std::string> labels;        //Data member labels
    long unsigned int nx, ny, nz;           //Number of grid elements in each dimension
    double Lx, Ly, Lz;                      //Dimensions of gridded region

    //Get first index of each dimension that corresponds to inside the grid
    virtual unsigned long int geti0() const = 0;
    virtual unsigned long int getj0() const = 0;
    virtual unsigned long int getk0() const = 0;
    virtual long unsigned int ntotal() const = 0;

    //Data accessors
    virtual long unsigned int ind1(const vec& x) const = 0;
    virtual long unsigned int ind2(const vec& x) const = 0;
    virtual long unsigned int ind3(const vec& x) const = 0;
    
    virtual long unsigned int ind2sub1(long unsigned int i) const = 0;
    virtual long unsigned int ind2sub2(long unsigned int i) const = 0;
    virtual long unsigned int ind2sub3(long unsigned int i) const = 0;
    
    virtual long unsigned int sub2ind(long unsigned int i, long unsigned int j, long unsigned int k) const = 0;

    //Geometry operations on cells
    virtual double Volume(long unsigned int i, long unsigned int j, long unsigned int k) const = 0;
    virtual vec rand(long unsigned int i, long unsigned int j, long unsigned int k) const = 0;
    virtual bool contains(const vec& x) const = 0;
    virtual double intersect(const vec& x, const vec& mu) const = 0;
    
    //Print grid specific info
    virtual void printOutside(std::ostream& oout, int id) const = 0;
    virtual void printGrid(std::ostream& oout) const = 0;
    virtual void printKeys(std::ostream& oout) const = 0;

    //Collision checks
    void collideCell(const vec& x, const vec& mu, double& ds, int& reflect) const;
    virtual void collideFront(const vec& x, const vec& mu, double& ds, int& reflect) const;
    virtual void collideBack(const vec& x, const vec& mu, double& ds, int& reflect) const;
    virtual void collideSide(const vec& x, const vec& mu, double& ds, int& reflect) const;
    virtual void boundXY(vec& x) const;
    virtual void boundZ(vec& x) const;
    virtual vec getNormal(int reflect, const vec& x, const vec& mu, vec& mur) = 0;

    //Simple iterator alternative
    long unsigned int ncell() const;
    iterator begin() const;

    double& at(int id, long unsigned int i, long unsigned int j, long unsigned int k);
    double& at(int id, long unsigned int i);
    double& at(int id, const vec& x);
    double& at(const vec& x);
    
    double dens(int id, const vec& x) const;
    
    double Volume(long unsigned int i) const;
    double Volume(const vec& x) const;
    
    vec rand(long unsigned int i) const;
    
    bool isempty(int id) const;
    double sum(int id) const;
    bool less(int id, double val) const;
    
    std::vector<double> cum(int id, double& sum) const;
    std::vector<double> cum(int id) const;
    
    void newval(std::string label);
    void clear();
    
    void print(std::ostream& oout, bool ALL = false) const;
    
};

#endif
