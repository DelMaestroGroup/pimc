/**
 * @file container.h 
 * @author Adrian Del Maestro
 * @date 11.24.2008
 *
 * @brief The simulation cell.
 */

#ifndef CONTAINER_H 
#define CONTAINER_H

#include "common.h"

// ========================================================================  
// Container Class
// ========================================================================  
/** 
 * The base class which holds details on the generalized box that our 
 * system will be simulated inside of.
 *
 * This class implements boundary conditions, and holds details on the
 * physical dimensions of the simulation cell. 
 */
class Container {
    public:
        Container();
        virtual ~Container();

        TinyVector <unsigned int, NDIM> periodic;   ///< Determines which dimensions have periodic bc

        dVec side;                          ///< The linear dimensions of the box
        dVec sideInv;                       ///< The inverse box dimensions
        dVec sideInv2;                      ///< 2 times the inverse box dimensions

        double volume;                      ///< The volume of the container in A^3
        double rcut2;                       ///< The smallest separation squared
        double maxSep;                      ///< The maximum possible separation for 2 beads on the same timeslice

        string name;                        ///< The name of the container

        int numGrid;                        ///< The number of grid boxes for the position grid
        bool fullyPeriodic;                 ///< Is the prism fully periodic?

        dVec gridSize;                      ///< The grid size in each dimension

        /** Place a vector in boundary conditions. */
        /** Algorithm C4 from 
         * @see: Z. Phys. Chem. 227 (2013) 345–352
         */
        void putInBC(dVec & r) const {
            r -= pSide*blitz::floor(r*sideInv + 0.5);
            /* int k; */
            /* for (int i = 0; i < NDIM; ++i) { */
            /*     k = int(r[i]*sideInv[i] + ((r[i]>=0.0)?0.5:-0.5)); */
            /*     r[i] -= k*pSide[i]; */
            /* } */
        }

        /* An old version */
        /* void putInBC1(dVec & r) const { */
        /*     for (int i = 0; i < NDIM; ++i) { */
        /*         r[i] -= (r[i] >= 0.5*side[i])*pSide[i]; */
        /*         r[i] += (r[i] < -0.5*side[i])*pSide[i]; */
        /*     } */
        /* } */

        /* Place a vector in boundary conditions. */
        void putInBC1(dVec & r) const {
            for (int i = 0; i < NDIM; ++i) {
                while (r[i] >= 0.5*side[i] && periodic[i])
                    r[i] -= pSide[i];
                while (r[i] < -0.5*side[i] && periodic[i])
                    r[i] += pSide[i];
            }
        }

        /** Place a vector inside the simulation cell */
        virtual void putInside(dVec &) const = 0;   
        
        /** Random position inside a box. */
        virtual dVec randPosition(MTRand &) const = 0;                  

        /** Random updated position inside a box. */
        virtual dVec randUpdate(MTRand &, const dVec &) const = 0;

        /** Map a position into a grid index */
        virtual int gridIndex(const dVec &) const = 0;

        /** The physical size of a NDIM-dimensional grid box */
        virtual double gridBoxVolume(const int) const = 0;

        /** The radius of a grid box */
        double gridRadius2(const int) const;


    protected:
        dVec pSide;         ///< Periodic * side
};

// ========================================================================  
// Prism Class
// ========================================================================  
/** 
 * A NDIM-dimensional hyperprism with periodic boundary conditions.
 */
class Prism: public Container {
    public:
        Prism(const double, const int);
        Prism(const dVec &, const iVec &_periodic=1);
        ~Prism();

        /** For PBC, this is identical to putInBC */
        void putInside(dVec &r) const {
            putInBC(r);

            /* This is slow!  Try function pointers? Or create a separate hard
             * top box? */
            if (!fullyPeriodic) {
                for (int i = 0; i < NDIM; i++) {
                    if (!periodic[i]) {
                        if (r[i] >= 0.5*side[i])
                            r[i] = 0.5*side[i] - 2*EPS;
                        if (r[i] < -0.5*side[i]) 
                            r[i] = -0.5*side[i] + 2*EPS;
                    }
                }
            }
        }

        dVec randPosition(MTRand &) const;                  
        dVec randUpdate(MTRand &, const dVec &) const;
        int gridIndex(const dVec &) const;
        double gridBoxVolume(const int) const;
};

// ========================================================================  
// Cylinder Class
// ========================================================================  
/**
 * A three dimensional cylinder with fixed boundary conditions in the 
 * x and y directions and periodic boundary conditions in the z direction.
 */
class Cylinder: public Container {
    public:
        Cylinder(const double, const double, const int);
        Cylinder(const double, const double);
        ~Cylinder();

        /* Place a vector inside the cylinder */
        void putInside(dVec &r) const;

        /* The various types of random positions inside the cylinder */
        dVec randPosition(MTRand &) const;                  
        dVec randUpdate(MTRand &, const dVec &) const;

        dVec randUpdateJumpShell(MTRand &, const dVec &) const;
        dVec randUpdateSmall(MTRand &, const dVec &) const;

        int gridIndex(const dVec &) const;
        double gridBoxVolume(const int) const;
};
#endif
