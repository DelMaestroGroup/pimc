/**
 * @file potential.h
 * @author Adrian Del Maestro
 * @date 10.14.2008
 *
 * @brief All possible potential classes.
 */

#ifndef POTENTIAL_H
#define POTENTIAL_H

#include "common.h"
#include "constants.h"
#include "gpkernel.h"
#include <boost/program_options.hpp>

#ifdef USE_GPU
#include "common_gpu.h"
#endif

class Path;
class LookupTable;
class Container;

namespace po = boost::program_options;

// ========================================================================  
// PotentialBase Class
// ========================================================================  
/** 
 * The base class from which all specific potentials are derived from.
 *
 * This class contains methods which return the actual value of 
 * the potential, an effective potential related to the pair product
 * approximation,  one which returns the gradient of the potential and a final
 * one which generates a sensible initial particle configuration.
 * We require knowledge of both the interaction as well as external potential
 * to run the simulation.
 */
class PotentialBase {

    public:
        PotentialBase ();
        virtual ~PotentialBase();
    
#if _cplusplus >= 202302L
        /** The potential */
        template <typename Self>
        inline __attribute__((always_inline)) double V(this Self&& s, const dVec &pos) {
            return std::forward<Self>(s).V(pos);
        }

        /** The effective potential for the pair product approximation */
        template <typename Self>
        inline __attribute__((always_inline)) double V(this Self&& s, const dVec &sep1, const dVec &sep2) {
            return std::forward<Self>(s).V(sep1, sep2);
        }

        /** The gradient of the potential*/
        template <typename Self>
        inline __attribute__((always_inline)) dVec gradV(this Self&& s, const dVec &pos) {
            return std::forward<Self>(s).gradV(pos);
        }

        /** Grad^2 of the potential*/
        template <typename Self>
        inline __attribute__((always_inline)) dVec grad2V(this Self&& s, const dVec &r)) {
            return std::forward<Self>(s).grad2V(r);
        }

        /** The derivative of the effective potential with respect to lambda
         *  and tau */
        template <typename Self>
        inline __attribute__((always_inline)) double dVdlambda(this Self&& s, const dVec &sep1, const dVec &sep2) {
            return std::forward<Self>(s).dVdlambda(sep1, sep2);
        }
        template <typename Self>
        inline __attribute__((always_inline)) double dVdtau(this Self&& s, const dVec &sep1, const dVec &sep2) {
            return std::forward<Self>(s).dVdtau(sep1, sep2);
        }
#else
        /** The potential */
        virtual double V(const dVec &) { return 0.0; }

        /** Batched potential values for external-potential style calls. */
        virtual void V(const dVec*, double*, int);

        /** The effective potential for the pair product approximation */
        virtual double V(const dVec &, const dVec &) { return 0.0; }

        /** The gradient of the potential*/
        virtual dVec gradV(const dVec &) { return dVec{}; }

        /** Grad^2 of the potential*/
        virtual double grad2V(const dVec &) { return 0.0; }

        /** The derivative of the effective potential with respect to lambda
         *  and tau */
        virtual double dVdlambda(const dVec &, const dVec &) {return 0.0;}
        virtual double dVdtau(const dVec &, const dVec &) {return 0.0;}
#endif
        
        /** Default Initial configuration of particles*/
        virtual DynamicArray<dVec,1> initialConfig(const Container*, MTRand &, const int); 

        /** A debug method that output's the potential based on user-supplied testPositions.  */
        void output();

        double tailV;       ///< Tail correction factor.

        /** Array to hold data elements*/
        virtual DynamicArray<double,1> getExcLen();

    protected:
        inline __attribute__((always_inline)) double deltaSeparation(double sep1,double sep2) const;
        DynamicArray<dVec,1> testPositions;     // An array that can be used for debugging a potential.
};

// ========================================================================  
// TabulatedPotential Class
// ========================================================================  
/** 
 * Pre-tabulated potential for complicated functions.
 *
 * In order to speed up the evaluation of complicated potentials, we
 * use a 2-point Newton-Gregory spline fit to perform the actual 
 * interpolation.
 */
template <typename T>
class TabulatedPotential {
    public:
        /**************************************************************************//**
         * Constructor.
        ******************************************************************************/
        TabulatedPotential() {
            extV.fill(0.0);
            extdVdr.fill(0.0);
            extd2Vdr2.fill(0.0);
        }

        /**************************************************************************//**
         * Destructor. 
        ******************************************************************************/
        ~TabulatedPotential() {
        }

    protected:
        DynamicArray <double,1> lookupV;           ///< A potential lookup table
        DynamicArray <double,1> lookupdVdr;        ///< A lookup table for dVint/dr
        DynamicArray <double,1> lookupd2Vdr2;      ///< A lookup table for d2Vint/dr2

        double dr;                          ///< The discretization for the lookup table
        int tableLength;                    ///< The number of elements in the lookup table

        std::array<double,2> extV;          ///< Extremal value of V
        std::array<double,2> extdVdr;       ///< Extremal value of dV/dr
        std::array<double,2> extd2Vdr2;     ///< Extremal value of d2V/dr2

        /**************************************************************************//**
         *  Given a discretization factor and the system size, create and fill
         *  the lookup tables for the potential and its derivative.
        ******************************************************************************/
        void initLookupTable(const double _dr, const double maxSep) {

            /* We now calculate the lookup tables for the interaction potential and 
             * its first and second derivatives. */
            dr = _dr;
            tableLength = int(maxSep/dr);
            lookupV.resize(tableLength);
            lookupdVdr.resize(tableLength);
            lookupd2Vdr2.resize(tableLength);
            lookupV.fill(0.0);
            lookupdVdr.fill(0.0);
            lookupd2Vdr2.fill(0.0);

            double r = 0;

            for (int n = 0; n < tableLength; n++) {
                lookupV(n)    = valueV(r);
                lookupdVdr(n) = valuedVdr(r);
                lookupd2Vdr2(n) = valued2Vdr2(r);
                r += dr;
            }

            /* r = 0.0; */
            /* for (int n = 0; n < tableLength; n++) { */
            /*     communicate()->file("debug")->stream() << format("%24.16e %24.16e\n") */
            /*         % r % lookupV(n); */
            /*     r += dr; */
            /* }; */

            /* exit(-1); */


            //      std::cout << format("%16.8E%16.8E%16.8E%16.8E%16.8E%16.8E%16.8E\n") % r % lookupV(n) % valueV(r) % 

            //  double rc = constants()->rc();
            //  for (int n = 0; n < tableLength; n++) {
            //      r += dr;
            //      if (r <= rc) {
            //          lookupV(n) = valueV(r) - valueV(rc) - valuedVdr(rc)*(r-rc);
            //          lookupdVdr(n) = valuedVdr(r) - valuedVdr(rc);
            //      }
            //      else {
            //          lookupV(n) = 0.0;
            //          lookupdVdr(n) = 0.0;
            //      }
            //      std::cout << format("%16.8E%16.8E%16.8E%16.8E%16.8E%16.8E%16.8E\n") % r % lookupV(n) % valueV(r) % 
            //          lookupdVdr(n) % valuedVdr(r) % (lookupV(n) - valueV(r)) % (lookupdVdr(n) - valuedVdr(r));
            //  }
        }

        /**************************************************************************//**
         *  Use the Newton-Gregory forward difference method to do a 2-point lookup
         *  on the potential table.  
         *
         *  @see M.P. Allen and D.J. Tildesley, "Computer Simulation of Liquids" 
         *  (Oxford Press, London, England) p 144 (2004).
        ******************************************************************************/
        inline __attribute__((always_inline)) double newtonGregory(const DynamicArray<double,1> &VTable, 
                const std::array<double,2> &extVal, const double r) {

            double rdr = r/dr;
            int k = int(rdr);

            if (k <= 0) 
                return extVal[0];

            if (k >= tableLength)
                return extVal[1];

            double xi = rdr - 1.0*k;
            double vkm1 = VTable(k-1);
            double vk = VTable(k);
            double vkp1 = VTable(k+1);

            double T1 = vkm1 + (vk - vkm1) * xi;
            double T2 = vk + (vkp1 - vk) * (xi - 1.0);

            return (T1 + 0.5 * (T2 - T1) * xi);
        }

        /**************************************************************************//**
         *  Use a direct lookup for the potential table.
         *
         *  This is faster thant Newton-Gregory and may give similar results for a fine
         *  enough mesh.
        ******************************************************************************/
        inline __attribute__((always_inline)) double direct(const DynamicArray<double,1> &VTable, 
                const std::array<double,2> &extVal, const double r) {

            int k = int(r/dr);
            if (k <= 0) 
                return extVal[0];

            if (k >= tableLength)
                return extVal[1];

            return VTable(k);
        }


        /** The functional value of V */
        inline __attribute__((always_inline)) double valueV (const double r) {
            return static_cast<T*>(this)->valueV(r);
        }

        /** The functional value of dV/dr */
        inline __attribute__((always_inline)) double valuedVdr (const double r) {
            return static_cast<T*>(this)->valuedVdr(r);
        }

        /** The functional value of d2V/dr2 */
        inline __attribute__((always_inline)) double valued2Vdr2 (const double r) {
            return static_cast<T*>(this)->valued2Vdr2(r);
        }
};

// ========================================================================  
// GaussianProcessPotential Class
// ========================================================================  
/** 
 * A generic Gaussian Process potential.
 *
 * Given a kernel, implements GP inference. 
 */
class GaussianProcessPotential {

    public:
        GaussianProcessPotential(const Container*, const po::variables_map &);
        virtual ~GaussianProcessPotential();

        double GP(const dVec&);
        void GP(const dVec*, double*, int);

    protected:
        std::unique_ptr<GaussianProcessKernelBase> kernelPtr; // the kernel
        DynamicArray<dVec,1> trainX;                // Normalized training data points
        DynamicArray<double,1> KinvY;               // The right hand vector K^{-1}(Y -μ) where Y is standardised 

    private:
        std::string kernelType;         // the type of kernel (e.g. matern)
        std::string meanType;           // the type of mean (e.g. constant)
        dVec ℓ;                         // the length scale hyperparameter
        uint32 numTrainingPoints;       // Number of training points

        dVec normOffset;                // normalization offset             
        dVec normScale;                 // normalization scale

        double dataStandardMean;        // data standardization mean
        double dataStandardStd;         // data standardization standard deviation
        double sigma2;                  // the scale of the kernel σ² 
                            
        double μ;                       // the GP mean
        double maternNu;
#ifdef USE_GPU
        gpu_stream_t gpStream;
        double *d_gpData = nullptr;
        double *d_positions = nullptr;
        double *d_values = nullptr;
        int gpuBufferCapacity = 0;
#endif
};

// ========================================================================  
// MultiFidelityGaussianProcessPotential Class
// ========================================================================  
// !!NOTE!! We decided we don't need this for now. 
/** 
 * A Multifidelity Gaussian Process potential.
 *
 * Given kernel, implements MFGP inference. 
 */
// class MultiFidelityGaussianProcessPotential {

//     public:
//         MultiFidelityGaussianProcessPotential(const Container*, const po::variables_map &);
//         virtual ~MultiFidelityGaussianProcessPotential() = default;

//         double GP(const dVec&);

//     protected:
//         std::vector<std::unique_ptr<GaussianProcessKernelBase>> kernelPtrs; // pointers to the needed kernels
//         DynamicArray<dVec,1> trainX;                // Normalized training data points
//         DynamicArray<double,1> KinvY;               // The right hand vector K^{-1}(Y -μ) where Y is standardised 
//         DynamicArray<double,1> fidelity;            // The fidelity associated with each point

//     private:
//         std::string kernelType;         // the type of kernel (e.g. matern)
//         std::string meanType;           // the type of mean (e.g. constant)
//         std::vector<dVec> ℓ;            // the length scale hyperparameter (for each fidelity)
//         uint32 numTrainingPoints;       // Number of training points

//         dVec normOffset;                // normalization offset             
//         dVec normScale;                 // normalization scale

//         double dataStandardMean;        // data standardization mean
//         double dataStandardStd;         // data standardization standard deviation
//         double MFPower;                 // relative power in multifidelity kernel
                            
//         double μ;                       // the GP mean
                                        
//         uint32 numFidelity;             // The number of kernels
// };

// ========================================================================  
// FreePotential Class
// ========================================================================  
/** 
 * Free potential.
 */
class FreePotential: public PotentialBase {
    public:
        FreePotential();
        ~FreePotential();

        /** The potential. */
        inline __attribute__((always_inline)) double V(const dVec &sep) { return 0.0*sep[0]; };

        /** The gradient of the potential. */
        inline __attribute__((always_inline)) dVec gradV(const dVec &pos) {
            return (0.0*pos);
        }
};

// ========================================================================  
// HarmonicPotential Class
// ========================================================================  
/** 
 * Computes the potential energy for an external harmonic potential.  
 *
 * We work in generalized units where hbar omega / k_B = 1.
 */
class HarmonicPotential : public PotentialBase {
    public:
        HarmonicPotential ();
        HarmonicPotential (double);
        ~HarmonicPotential ();
    
        double omega2;       //The SHO frequency in units of hbar

        /** The potential. */
        inline __attribute__((always_inline)) double V(const dVec &r) { 
            return (omega2*dot(r,r)/(4.0*constants()->lambda()));
        }

        /** The gradient of the potential. */
        inline __attribute__((always_inline)) dVec gradV(const dVec &r) {
            dVec tempr;
            tempr = r;
            return (omega2*tempr/(2.0*constants()->lambda()));
        }

        /** Initial configuration corresponding to Harmonic potential */
	DynamicArray<dVec,1> initialConfig(const Container*, MTRand &, const int); 
};

// ========================================================================  
// SingleWellPotential Class
// ========================================================================  
/**
 * Computes the potential energy for an external single well potential.
 */
class SingleWellPotential : public PotentialBase {
    public:
        SingleWellPotential ();
        ~SingleWellPotential ();

        /** The potential */
        inline __attribute__((always_inline)) double V(const dVec &r) { 
            double r2 = dot(r,r);
            return ( 0.5*r2 + r2*r2 );
        }

        /** The gradient of the potential. */
        inline __attribute__((always_inline)) dVec gradV(const dVec &r) {
            double r2 = dot(r,r);
            dVec tempr;
            tempr = r;
            return ((1.0 + 4.0*r2)*tempr);
        }
};

// ========================================================================  
// HarmonicCylinderPotential Class
// ========================================================================  
/**
 * Computes the potential energy for an external harmonic potential with
 * axial symmetry.
 */
class HarmonicCylinderPotential : public PotentialBase {
    public:
        HarmonicCylinderPotential (const double);
        ~HarmonicCylinderPotential ();

        /** The potential. */
        inline __attribute__((always_inline)) double V(const dVec &r) { 
            double r2 = 0.0;
            for (int i=0; i < NDIM-1; i++)
                r2 += r[i]*r[i];
            return ( 0.5 * c * constants()->m() * w * w * r2 ); 
        }

        /** The gradient of the potential. */
        inline __attribute__((always_inline)) dVec gradV(const dVec &r) {
            dVec tempr{};
            for (int i=0; i < NDIM-1; i++)
                tempr[i] = r[i];
            return ( c * constants()->m() * w * w * tempr );
        }

    private:
        double w;               // The confining frequency
        double c;               // A dimension-full constant
};

// ========================================================================  
// DeltaPotential Class
// ========================================================================  
/** 
 * Computes the potential energy for delta function interaction potential,
 * approximated here as the limit of a Cauchy distribution.
 *
 * NB: Need to add an explanation to compare with the previous version of the
 *     Gaussian.
 */
class DeltaPotential : public PotentialBase  {
    public:
        DeltaPotential (double,double);
        ~DeltaPotential ();

        /**
         * Return the delta function potential with strength g approximated as
         * the limit of a Gaussian distribution.
         * Tested and working in Mathematica.
         */
        inline __attribute__((always_inline)) double V(const dVec &r) {
            return (norm*exp(-dot(r,r)*inv2sigma2));
        }

        /**
         * Return the gradient of the delta function potential with strength 
         * g approximated as the limit of a Gaussian distribution.
         * Tested and working in Mathematica.
         */
        inline __attribute__((always_inline)) dVec gradV(const dVec &r) {
            return (-2.0*r*norm*inv2sigma2*exp(-dot(r,r)*inv2sigma2));
        }

    private:
        double norm;            // A normalization constant for fixed strength
        double inv2sigma2;      // 1/(2\sigma^2) where \sigma^2 is the variance
};


// ========================================================================  
// LorentzianPotential Class
// ========================================================================  
/**
 * Computes the potential energy for delta function interaction potential,
 * approximated here as the limit of a Cauchy distribution.
 */
class LorentzianPotential : public PotentialBase  {
    public:
        LorentzianPotential (double,double);
        ~LorentzianPotential ();

        /** 
         * Return the delta function potential with strength 2c approximated as
         * the limit of a Lorentzian distribution.
         * Tested and working in Mathematica.
         */
        inline __attribute__((always_inline)) double V(const dVec &r) {
            return (norm / (a*a + dot(r,r)));
        }

        /** 
         * Return the gradient of the delta function potential with strength 
         * 2c approximated as the limit of a Lorentzian distribution.
         * Tested and working in Mathematica.
         */
        inline __attribute__((always_inline)) dVec gradV(const dVec &r) {
            double b = a*a + dot(r,r);
            return ((-(2.0*norm*a)/(b*b))*r);
        }

    private:
        double c;               // The strength of the delta function
        double norm;            // A normalization constant for fixed strength
        double a;               // The order of the limit
};

// ========================================================================  
// SutherlandPotential Class
// ========================================================================  
/** 
 * Computes the potential energy for the periodic Sutherland model which
 * approximates long-range 1/r^2 interactions on a ring.
 * @see:  B. Sutherland, Phys Rev A 4, 2019 (1971) 
 */
class SutherlandPotential : public PotentialBase  {
    public:
        SutherlandPotential (double);
        ~SutherlandPotential ();

        /**
         * Return the Sutherland potential g/r^2.
         */
        inline __attribute__((always_inline)) double V(const dVec &r) {
            double x = pioL*(sqrt(dot(r,r)) + EPS);
            return g * pioL * pioL / (sin(x)*sin(x));
        }

        /**
         * Return the gradient of the Sutherland potential.
         */
        inline __attribute__((always_inline)) dVec gradV(const dVec &r) {
            double rnorm = sqrt(dot(r,r)) + EPS;
            double x = pioL*rnorm;
            double s = sin(x);
            return (-2.0* g * pioL * pioL * pioL * cos(x) / (s*s*s*rnorm)) * r;
        }

        /**
         * Return the Laplacian of the Sutherland potential.
         */
        inline __attribute__((always_inline)) double grad2V(const dVec &r) {
            double x = pioL*(sqrt(dot(r,r))+EPS);
            double s = sin(x);
            return 2.0* g * pioL * pioL * pioL * pioL * (2.0+cos(2*x)) / 
                (s*s*s*s);
        }

    private:
        double g;               // The interaction constant g
        double pioL;            // \pi/L
};

// ========================================================================  
// Dipolar Potential Class
// ========================================================================  
/** 
 * Computes the potential energy for polarized electrical dipoles with strength
 * D in reduced units where lengths are measured in units of a = m D / \hbar^2 and
 * energies in units of \hbar^2 / m a^2.
 * @see:  http://link.aps.org/doi/10.1103/PhysRevA.87.063604
 */
class DipolePotential : public PotentialBase  {
    public:
        DipolePotential ();
        ~DipolePotential ();

        /**
         * Return the dipole potential 1/r^3.
         */
        inline __attribute__((always_inline)) double V(const dVec &r) {
            double x = sqrt(dot(r,r));
            if (x < EPS)
                return LBIG;
            return 1.0/(x*x*x);
        }

        /**
         * Return the gradient of the dipole potential.
         * \partial (1/r^3) \partial r \hat{r} = -3r^{-5} \vec{r}
         */
        inline __attribute__((always_inline)) dVec gradV(const dVec &r) {
            double x = sqrt(dot(r,r));
            if (x < EPS)
                return dVec{};
            return (-3.0/(x*x*x*x*x)) * r;
        }

        /**
         * Return the Laplacian of the dipolar potential.
         */
        inline __attribute__((always_inline)) double grad2V(const dVec &r) {
            double x = sqrt(dot(r,r));
            if (x < EPS)
                return 0.0;
            return 6.0/(x*x*x*x*x);
        }
};

#if NDIM > 2
// ========================================================================  
// Hard Cylinder Potential Class
// ========================================================================  
/** 
 * Computes the value of the external wall potential for a hard-walled 
 * cylindrical cavity.
 */
class HardCylinderPotential : public PotentialBase {
    public:
        HardCylinderPotential (const double);
        ~HardCylinderPotential ();

        /** A step function at rho=R. */
        inline __attribute__((always_inline)) double V(const dVec &r) {
            if (sqrt(r[0]*r[0]+r[1]*r[1]) >= R)
                return LBIG;
            else
                return 0.0;
        }

        /** A delta function at rho=R. */
        inline __attribute__((always_inline)) dVec gradV(const dVec &r) {
            dVec tempr;
            tempr = r;
            tempr[2] = 0.0;
            if (abs((r[0]*r[0]+r[1]*r[1])-R*R)<1.0E-3)
                return LBIG*tempr;
            else
                return 0.0*tempr;
        }

    private:
        double R;       // Radius of the tube
};
#endif

#if NDIM > 2
// ========================================================================  
// Plated LJ Cylinder Potential Class
// ========================================================================  
/** 
 * Computes the value of the external wall potential for a plated cylindrical 
 * cavity. 
 */
class PlatedLJCylinderPotential : public PotentialBase, public TabulatedPotential<PlatedLJCylinderPotential> {
    public:
        PlatedLJCylinderPotential (const double, const double, const double, const double, const double);
        ~PlatedLJCylinderPotential ();

        /** The integrated LJ Wall potential. */
        inline __attribute__((always_inline)) double V(const dVec &r) {
            int k = static_cast<int>(sqrt(r[0]*r[0] + r[1]*r[1])/dR);
            if (k >= tableLength)
                return extV[1];
            else
                return lookupV(k);
        }

        /* The gradient of the LJ Wall potential */
        inline __attribute__((always_inline)) dVec gradV(const dVec &);

        /* Laplacian of the LJ Wall potential */
        inline __attribute__((always_inline)) double grad2V(const dVec &);

        /** Initial configuration corresponding to the LJ cylinder potential */
	DynamicArray<dVec,1> initialConfig(const Container*, MTRand &, const int); 

        /* Used to construct the lookup tables */
        inline __attribute__((always_inline)) double valueV (const double);               
        inline __attribute__((always_inline)) double valuedVdr (const double);                
        inline __attribute__((always_inline)) double valued2Vdr2 (const double);

    private:
        /** local methods for computing the potential of a LJ cylinder */
        double V_ (const double, const double, const double, const double, const double);
        double dVdr_ (const double, const double, const double, const double, const double);

        /* All the parameters needed for the LJ wall potential */
        double densityPlated;
        double sigmaPlated;
        double epsilonPlated;
        double density;
        double sigma;
        double epsilon;

        double Ri;      // Inner radius of the tube
        double Ro;      // Outer radius of the tube
        double dR;      // Discretization for the lookup table

        double minV;    // The minimum value of the potential
};

// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// INLINE FUNCTION DEFINITIONS
// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

/** 
 * Return the gradient of aziz potential for separation r using a 
 * lookup table. 
 */
inline __attribute__((always_inline)) dVec PlatedLJCylinderPotential::gradV(const dVec &r) {
    double rnorm = sqrt(r[0]*r[0] + r[1]*r[1]);
    dVec tempr;
    tempr = r;
    tempr[2] = 0.0;
    int k = static_cast<int>(rnorm/dR);
    dVec gV;
    if (k >= tableLength)
        gV = (extdVdr[1]/rnorm)*tempr;
    else
        gV = (lookupdVdr(k)/rnorm)*tempr;
    return gV;
}

/** 
 * Return the Laplacian of aziz potential for separation r using a 
 * lookup table. 
 */
inline __attribute__((always_inline)) double PlatedLJCylinderPotential::grad2V(const dVec &r) {
    double rnorm = sqrt(r[0]*r[0] + r[1]*r[1]);
    dVec tempr;
    tempr = r;
    tempr[2] = 0.0; // PBC in z-direction
    int k = static_cast<int>(rnorm/dR);
    double g2V;
    if (k >= tableLength)
        g2V = extd2Vdr2[1];
    else
        g2V = lookupd2Vdr2(k);
    return g2V;
}
#endif

#if NDIM > 2
// ========================================================================  
// LJ Cylinder Potential Class
// ========================================================================  
/** 
 * Computes the value of the external wall potential for a cylindrical 
 * cavity. 
 */
class LJCylinderPotential : public PotentialBase, public TabulatedPotential<LJCylinderPotential> {
    public:
        LJCylinderPotential (const double, const double, const double, const double);
        ~LJCylinderPotential ();

        /** The integrated LJ Wall potential. */
        inline __attribute__((always_inline)) double V(const dVec &r) {
            int k = static_cast<int>(sqrt(r[0]*r[0] + r[1]*r[1])/dR);
            if (k >= tableLength)
                return extV[1];
            else
                return lookupV(k);
        }

        /* The gradient of the LJ Wall potential */
        inline __attribute__((always_inline)) dVec gradV(const dVec &);

        /* Laplacian of the LJ Wall potential */
        inline __attribute__((always_inline)) double grad2V(const dVec &);

        /** Initial configuration corresponding to the LJ cylinder potential */
	DynamicArray<dVec,1> initialConfig(const Container*, MTRand &, const int); 

        /* Used to construct the lookup tables */
        inline __attribute__((always_inline)) double valueV (const double);               
        inline __attribute__((always_inline)) double valuedVdr (const double);                
        inline __attribute__((always_inline)) double valued2Vdr2 (const double);

    private:
        /* All the parameters needed for the LJ wall potential */
        double density;
        double sigma;
        double epsilon;

        double R;       // Radius of the tube
        double dR;      // Discretization for the lookup table

        double minV;    // The minimum value of the potential
};

// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// INLINE FUNCTION DEFINITIONS
// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

/** 
 * Return the gradient of aziz potential for separation r using a 
 * lookup table. 
 */
inline __attribute__((always_inline)) dVec LJCylinderPotential::gradV(const dVec &r) {
    double rnorm = sqrt(r[0]*r[0] + r[1]*r[1]);
    dVec tempr;
    tempr = r;
    tempr[2] = 0.0;
    int k = static_cast<int>(rnorm/dR);
    dVec gV;
    if (k >= tableLength)
        gV = (extdVdr[1]/rnorm)*tempr;
    else
        gV = (lookupdVdr(k)/rnorm)*tempr;
    return gV;
}

/** 
 * Return the Laplacian of aziz potential for separation r using a 
 * lookup table. 
 */
inline __attribute__((always_inline)) double LJCylinderPotential::grad2V(const dVec &r) {
    double rnorm = sqrt(r[0]*r[0] + r[1]*r[1]);
    dVec tempr;
    tempr = r;
    tempr[2] = 0.0; // PBC in z-direction
    int k = static_cast<int>(rnorm/dR);
    double g2V;
    if (k >= tableLength)
        g2V = extd2Vdr2[1];
    else
        g2V = lookupd2Vdr2(k);
    return g2V;
}
#endif

#if NDIM > 2
// ========================================================================  
// LJ Hour Glass Potential Class
// ========================================================================  
/** 
 * Computes the value of the external wall potential for an hour-glass shaped
 * cavity. 
 */
class LJHourGlassPotential : public PotentialBase { 

    public:
        LJHourGlassPotential (const double, const double, const double);
        ~LJHourGlassPotential ();

        /** The integrated LJ hour glass potential. */
        inline __attribute__((always_inline)) double V(const dVec &); 

        /** The gradient of the potential. Use the prrimitive approx. */
        inline __attribute__((always_inline)) dVec gradV(const dVec &pos) {
            return (0.0*pos);
        }

        /* Laplacian of the LJ Wall potential. Use primitive approx. */
        inline __attribute__((always_inline)) double grad2V(const dVec & pos) {
            return 0.0;
        }

        /** Initial configuration corresponding to the LJ cylinder potential */
	DynamicArray<dVec,1> initialConfig(const Container*, MTRand &, const int); 
	DynamicArray<dVec,1> initialConfig1(const Container*, MTRand &, const int); 

    private:
        /* All the parameters needed for the LJ wall potential */
        double density;
        double sigma;
        double epsilon;

        double R;       // Radius of the tube at +/- L/2
        double dR;      // variation of radius

        /* double minV; // The minimum value of the potential */

        double L;       // The legnth of the pore
        double invd;    // The inverse of the variation length scale
        double R0;      // A constant used in the tanh radius function

        /* Different functional forms for the hourglass radius */
        double Rlinear(double z) {
            return (2.0*dR*abs(z)/L + R - dR);
        }

        double Rtanh(double z) {
            double t = R0*tanh(z*invd);
            return (dR*t*t + R - dR);
        }
};
#endif

// ========================================================================  
// Aziz Potential Class
// ========================================================================  
/** 
 * Computes the value of the semi-empircal Aziz potential that is known
 * to be accurate for He-4.
 */
class AzizPotential : public PotentialBase, public TabulatedPotential<AzizPotential> {
    public:
        AzizPotential (const int, const Container *);
        ~AzizPotential ();

        /* The Aziz HFDHE2 Potential */
        inline __attribute__((always_inline)) double V(const dVec &);

        /* The gradient of the Aziz potential */
        inline __attribute__((always_inline)) dVec gradV(const dVec &);

        /* The Laplacian of the Aziz potential */
        inline __attribute__((always_inline)) double grad2V(const dVec &);

        /* Used to construct the lookup tables */
        inline __attribute__((always_inline)) double valueV (const double);               
        inline __attribute__((always_inline)) double valuedVdr (const double);                    
        inline __attribute__((always_inline)) double valued2Vdr2 (const double);

    private:
        /* All the parameters of the Aziz potential */
        double rm, A, epsilon, alpha, beta, D, C6, C8, C10;

        /* The F-function needed for the Aziz potential */
        double F(const double x) {
            return (x < D ? exp(-(D/x - 1.0)*(D/x - 1.0)) : 1.0 );
        }

        /* The derivative of the F-function needed for the Aziz potential */
        double dF(const double x) {
            double ix = 1.0/x;
            double r = 2.0*D*ix*ix*(D*ix-1.0)*exp(-(D*ix - 1.0)*(D*ix - 1.0));
            return (x < D ? r : 0.0 );
        }   
     
        /* The 2nd derivative of the F-function needed for the Aziz potential
         * Double checked with Mathematica --MTG */
        double d2F(const double x) {
            double ix = 1.0/x;
            double r = 2.0*D*ix*ix*ix*( 2.0*D*D*D*ix*ix*ix - 4.0*D*D*ix*ix 
                    - D*ix + 2.0) * exp(-(D*ix - 1.0)*(D*ix - 1.0));
            return (x < D ? r : 0.0 );
        }
};

// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// INLINE FUNCTION DEFINITIONS
// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

/** 
 * Return the aziz potential for separation r using a lookup table. 
 */
inline __attribute__((always_inline)) double AzizPotential::V(const dVec &r) {
    //double rnorm = sqrt(dot(r,r));
    //return newtonGregory(lookupV,extV,rnorm);
    return direct(lookupV,extV,sqrt(dot(r,r)));
}

// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

/** 
 * Return the gradient of aziz potential for separation r using a 
 * lookup table. 
 */
inline __attribute__((always_inline)) dVec AzizPotential::gradV(const dVec &r) {
    double rnorm = sqrt(dot(r,r));
    dVec gV;
    //gV = (newtonGregory(lookupdVdr,extdVdr,rnorm)/rnorm)*r;
    gV = (direct(lookupdVdr,extdVdr,rnorm)/rnorm)*r;
    return gV;
}

// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

/** 
 * Return the Laplacian of aziz potential for separation r using a 
 * lookup table. 
 */

inline __attribute__((always_inline)) double AzizPotential::grad2V(const dVec &r) {
    double rnorm = sqrt(dot(r,r));
    double g2V;
    //g2V = (newtonGregory(lookupd2Vdr2,extd2Vdr2,rnorm)/rnorm)*r;
    g2V = direct(lookupd2Vdr2,extd2Vdr2,rnorm);
    return g2V;
}

// ========================================================================
// Hydrogen Potential Class
// ========================================================================

class H2LJ : public PotentialBase, public TabulatedPotential<H2LJ>
{
     public:
         H2LJ (const Container *);
         ~H2LJ ();

         /* The Lennard-Jones potential of H2 */
         inline __attribute__((always_inline)) double V(const dVec &);

         /* The gradient of the Lennard-Jones potential */
         inline __attribute__((always_inline)) dVec gradV(const dVec &);

         /* The Laplacian of the Lennard-Jones potential */
         inline __attribute__((always_inline)) double grad2V(const dVec &);

         /* Used to construct the lookup tables */
         inline __attribute__((always_inline)) double valueV (const double);
         inline __attribute__((always_inline)) double valuedVdr (const double);
         inline __attribute__((always_inline)) double valued2Vdr2 (const double);

     private:
         /* All the parameters of the Lennard Jones potential */
         double EPSILON, SIGMA;
 };


 inline __attribute__((always_inline)) double H2LJ::V(const dVec &r)
 {
     return direct(lookupV,extV,sqrt(dot(r,r)));
 }

 inline __attribute__((always_inline)) dVec H2LJ::gradV(const dVec &r)
 {
     double rnorm = sqrt(dot(r,r));
     dVec gV;
     gV = (direct(lookupdVdr,extdVdr,rnorm)/rnorm)*r;
     return gV;
 }

 inline __attribute__((always_inline)) double H2LJ::grad2V(const dVec &r)
 {
     double rnorm = sqrt(dot(r,r));
     double g2V;
     g2V = direct(lookupd2Vdr2,extd2Vdr2,rnorm);
     return g2V;
 }


 // ========================================================================
 // Silvera-Goldman Potential Class
 // ========================================================================

 // See: I.F. Silvera and V.V. Goldman, J. Chem. Phys. 69, 4209 (1978).
 // This implementation includes the effective Axilrod-Teller interactions.

 class SilveraPotential : public PotentialBase, public TabulatedPotential<SilveraPotential> {
     public:
         SilveraPotential (const Container *);
         ~SilveraPotential ();

         /* The Silvera-Goldman Potential */
         inline __attribute__((always_inline)) double V(const dVec &);

         /* The gradient of the Silvera-Goldman potential */
         inline __attribute__((always_inline)) dVec gradV(const dVec &);

         /* The Laplacian of the Aziz potential */
         inline __attribute__((always_inline)) double grad2V(const dVec &);

         /* Used to construct the lookup tables */
         inline __attribute__((always_inline)) double valueV (const double);
         inline __attribute__((always_inline)) double valuedVdr (const double);
         inline __attribute__((always_inline)) double valued2Vdr2 (const double);

     private:
         /* All the parameters of the Silvera-Goldman potential */
         double ALPHA, BETA, GAMMA, C6, C8, C9, C10, Rc;
         double BohrPerAngstrom, KelvinPerHartree;

         // The damping function needed for the Silvera-Golman potential.
         // The argument must be in units of Bohr radii.
         double F(const double r)
         {
             double term = -1.0*(Rc/r - 1.0)*(Rc/r - 1.0);
             return (r < Rc ? exp(term) : 1.0);
         }

         // The 1st derivative of damping function.
         // Calculated with Mathematica.  03/05/25 TRP.
         double dF(const double r)
         {
             double ir = 1.0/r;
             double ir3 = ir*ir*ir;
             double term = -1.0*(Rc/r - 1.0)*(Rc/r - 1.0);
             double express = 2.0*Rc*(Rc-r)*ir3*exp(term);
             return (r < Rc ? express : 0.0);
         }

         // The 2nd derivative of the damping function.
         // Calculated witwh Mathematica.  03/05/25 TRP.
           double d2F(const double r)
           {
               double ir = 1.0/r;
               double ir6 = (ir*ir*ir)*(ir*ir*ir);
               double t1 = 2.0*Rc*Rc*Rc;
               double t2 = -4.0*Rc*Rc*r;
               double t3 = -1.0*Rc*r*r;
               double t4 = 2.0*r*r*r;
               double poly = t1 + t2 + t3 + t4;
               double term = -1.0*(Rc/r - 1.0)*(Rc/r - 1.0);
               double express = 2.0*Rc*exp(term)*poly*ir6;
               return (r < Rc ? express : 0.0);
         }

         // Repulsive part of potential, along with first and second derivatives.
         // Atomic units (r in Bohr; V in Hartree).
         double Vrep(const double r)
         {
             return exp(ALPHA - BETA*r - GAMMA*r*r);
         }

         double dVrep(const double r)
         {
             return (BETA + 2.0*GAMMA*r)*(-1.0*Vrep(r));
         }

         double d2Vrep(const double r)
         {
             double term1 = 2.0*GAMMA*r*r - 1;
             double term2 = BETA*BETA + 4.0*BETA*GAMMA*r + 2.0*GAMMA*term1;
             return term2*Vrep(r);
         }

         // Attractive part of potential (plus three-body interaction).
         // Atomic units (r in Bohr; V in Hartree.)
         double Vatt(const double r)
         {
             double ir = 1.0/r;
             double ir6 = (ir*ir*ir)*(ir*ir*ir);
             double ir8 = ir6*ir*ir;
             double ir10 = ir8*ir*ir;
             double ir9 = ir8*ir;
             double attract = -(C6*ir6 + C8*ir8 + C10*ir10);
             double threebody = C9*ir9;

             return attract + threebody;
         }

         double dVatt(const double r)
         {
             double ir = 1.0/r;
             double ir7 = (ir*ir*ir)*(ir*ir*ir*ir);
             double ir9 = ir7*(ir*ir);
             double ir10 = ir9*ir;
             double ir11 = ir10*ir;

             return 6.0*C6*ir7 + 8.0*C8*ir9 - 9.0*C9*ir10 + 10.0*C10*ir11;
         }

         double d2Vatt(const double r)
         {
             double term1 = 7.0*C6*r*r*r + 12.0*C8*r - 15.0*C9;
             double numerator = 3.0*r*term1 + 55*C10;
             double ir = 1.0/r;
             double ir4 = ir*ir*ir*ir;
             double ir12 = ir4*ir4*ir4;
             return -2.0*numerator*ir12;
         }
 };


 inline __attribute__((always_inline)) double SilveraPotential::V(const dVec &r)
 {
     return direct(lookupV,extV,sqrt(dot(r,r)));
 }

 inline __attribute__((always_inline)) dVec SilveraPotential::gradV(const dVec &r)
 {
     double rnorm = sqrt(dot(r,r));
     dVec gV;
     gV = (direct(lookupdVdr,extdVdr,rnorm)/rnorm)*r;
     return gV;
 }

 inline __attribute__((always_inline)) double SilveraPotential::grad2V(const dVec &r)
 {
     double rnorm = sqrt(dot(r,r));
     double g2V;
     g2V = direct(lookupd2Vdr2,extd2Vdr2,rnorm);
     return g2V;
 }

// ========================================================================
// PatkowskiPotential Class
// ========================================================================
//
//  Spherically-averaged (isotropic) H2-H2 pair potential of
//  K. Patkowski, W. Cencek, P. Jankowski, K. Szalewicz, J. B. Mehl,
//  G. Garberoglio, A. H. Harvey,  J. Chem. Phys. 129, 094304 (2008).
//
//  The supplementary material contains Fortran source code for the isotropic 
// potential, which is used here to implement the PatkowskiPotential class. 

class PatkowskiPotential : public PotentialBase,
                           public TabulatedPotential<PatkowskiPotential> {
    public:
        PatkowskiPotential (const Container *);
        ~PatkowskiPotential ();

        inline __attribute__((always_inline)) double V(const dVec &);
        inline __attribute__((always_inline)) dVec gradV(const dVec &);

        // Radial 2nd derivative d^2V/dr^2 (NOT the 3D Laplacian), K/A^2.
        inline __attribute__((always_inline)) double grad2V(const dVec &);

        inline __attribute__((always_inline)) double valueV     (const double);
        inline __attribute__((always_inline)) double valuedVdr  (const double);
        inline __attribute__((always_inline)) double valued2Vdr2(const double);

    private:
        // All parameters preserved as they appear in the EPAPS
        // Fortran block-data (isoH2H2.f / fit_param.txt).
        double cex[2];      // short-range exp prefactor params
        double csp[4];      // short-range polynomial coeffs (K, K/bohr, ...)
        double cdata[3];    // long-range C_n^000 for n=6,8,10  (a.u., negative)

        // Inner hard-wall plateau r < 0.75 A
        double r_inner_A;     // Angstroms
        double V_inner_K;     // Kelvin

        // Unit conversions.
        double BohrPerAngstrom;     // R[bohr] = r[A] * BohrPerAngstrom
        double xK2au;               // 1 K = xK2au atomic units

        // Tang-Toennies damping function.
        double TT_f (const int n, const double x) const 
        {
            double term = 1.0, sum = 1.0;
            for (int k = 1; k <= n; ++k) 
            {
                term *= x / static_cast<double>(k);
                sum  += term;
            }
            return 1.0 - exp(-x) * sum;
        }

        double TT_df(const int n, const double x) const 
        {
            double xn_over_nfact = 1.0;
            for (int k = 1; k <= n; ++k)
                xn_over_nfact *= x / static_cast<double>(k);
            return exp(-x) * xn_over_nfact;
        }

        double TT_d2f(const int n, const double x) const 
        {
            if (x == 0.0) return 0.0;
            return TT_df(n, x) * (static_cast<double>(n) - x) / x;
        }
};


// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// INLINE FUNCTION DEFINITIONS
// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

inline __attribute__((always_inline))
double PatkowskiPotential::V(const dVec &r)
{
    return direct(lookupV, extV, sqrt(dot(r,r)));
}

inline __attribute__((always_inline))
dVec PatkowskiPotential::gradV(const dVec &r)
{
    double rnorm = sqrt(dot(r,r));
    dVec gV;
    gV = (direct(lookupdVdr, extdVdr, rnorm) / rnorm) * r;
    return gV;
}

inline __attribute__((always_inline))
double PatkowskiPotential::grad2V(const dVec &r)
{
    double rnorm = sqrt(dot(r,r));
    return direct(lookupd2Vdr2, extd2Vdr2, rnorm);
}

// ========================================================================  
// Szalewicz Potential Class
// ========================================================================  
/** 
 * Computes the value of the semi-empircal Szalewicz potential that is known
 * to be accurate for He-4.
 */
class SzalewiczPotential : public PotentialBase, public TabulatedPotential<SzalewiczPotential> {
    public:
        SzalewiczPotential (const Container *);
        ~SzalewiczPotential ();

        /* The Szalewicz HFDHE2 Potential */
        inline __attribute__((always_inline)) double V(const dVec &);

        /* The gradient of the Szalewicz potential */
        inline __attribute__((always_inline)) dVec gradV(const dVec &);

        /* The Laplacian of the Szalewicz potential */
        inline __attribute__((always_inline)) double grad2V(const dVec &);

        /* Used to construct the lookup tables */
        inline __attribute__((always_inline)) double valueV (const double);               
        inline __attribute__((always_inline)) double valuedVdr (const double);                    
        inline __attribute__((always_inline)) double valued2Vdr2 (const double);

    private:
        /* All the parameters of the Szalewicz potential */
        double rm = 2.9262186279335958;
        double Cn[17] = {   0.0,
                            0.0,
                            0.0,
                            0.000000577235,
                            -0.000035322,
                            0.000001377841,
                            1.461830,
                            0.0,
                            14.12350,
                            0.0,
                            183.7497,
                            -0.7674e2,
                            0.3372e4,
                            -0.3806e4,
                            0.8534e5,
                            -0.1707e6,
                            0.286e7
                            };
        double a = 3.64890303652830;
        double b = 2.36824871743591;
        double eta = 4.09423805117871;
        double Pn[3] = {-25.4701669416621, 269.244425630616, -56.3879970402079};
        double Qn[2] = {38.7957487310071, -2.76577136772754};
        long int factorials[17] = {     1,
                                        1,
                                        2,
                                        6,
                                        24,
                                        120,
                                        720,
                                        5040,
                                        40320,
                                        362880,
                                        3628800,
                                        39916800,
                                        479001600,
                                        6227020800,
                                        87178291200,
                                        1307674368000,
                                        20922789888000};

        /* The Tang-Toennies damping function needed for the Szalewicz potential */
        /* Can be described as upper incomplete gamma function.*/
        double fn(const double x, const int n) {
            double s1 = 0.0;
            for (int i = 0; i < n + 1; i++) {
            s1 += pow(x,i)/factorials[i];           
            }
            return 1.0 - (exp(-x)*s1);
        }

        /* The derivative of the Tang-Toennies damping function needed for the Szalewicz potential */
        double dfn(const double x, const int n) {
            return (exp(-x)*pow(x,n)/factorials[n]);
        }   
     
        /* The 2nd derivative of the Tang-Toennies needed for the Szalewicz potential*/
        double d2fn(const double x, const int n) {
            return (exp(-x)*(n-x)*pow(x,(n-1))/factorials[n]);
        }

        /* The Tang-Toennies damping function expanded for small r */
        /* needed for the Szalewicz potential */
        /* Can be described as upper incomplete gamma function. */
        double fn2(const double x, const int n) {
        
            double s1 = 0.0;
            for (int i = 0; i < 17; i++) {
                s1 += (pow(-1,i)*pow(x,(i+1)))/(factorials[i]*(n+i+1));         
            }
            s1 *= pow(x,n)/factorials[n];
            return s1;
        }

};

// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// INLINE FUNCTION DEFINITIONS
// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

/** 
 * Return the Szalewicz potential for separation r using a lookup table. 
 */
inline __attribute__((always_inline)) double SzalewiczPotential::V(const dVec &r) {
    //double rnorm = sqrt(dot(r,r));
    //return newtonGregory(lookupV,extV,rnorm);
    return direct(lookupV,extV,sqrt(dot(r,r)));
}

// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

/** 
 * Return the gradient of Szalewicz potential for separation r using a 
 * lookup table. 
 */
inline __attribute__((always_inline)) dVec SzalewiczPotential::gradV(const dVec &r) {
    double rnorm = sqrt(dot(r,r));
    dVec gV;
    //gV = (newtonGregory(lookupdVdr,extdVdr,rnorm)/rnorm)*r;
    gV = (direct(lookupdVdr,extdVdr,rnorm)/rnorm)*r;
    return gV;
}

// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

/** 
 * Return the Laplacian of Szalewicz potential for separation r using a 
 * lookup table. 
 */

inline __attribute__((always_inline)) double SzalewiczPotential::grad2V(const dVec &r) {
    double rnorm = sqrt(dot(r,r));
    double g2V;
    //g2V = (newtonGregory(lookupd2Vdr2,extd2Vdr2,rnorm)/rnorm)*r;
    //g2V = (direct(lookupd2Vdr2,extd2Vdr2,rnorm)/rnorm)*r;
    g2V = direct(lookupd2Vdr2,extd2Vdr2,rnorm);
    return g2V;
}

// ========================================================================  
// FixedAzizPotential Class
// ========================================================================  
/** 
 * Computes the potential energy resulting from a series of fixed helium
 * atoms that are not updated and provide a static 'external' potential.
 *
 * We require a Aziz potential interaction pointer to properly compute
 * the interaction with the static particles.
 */
class FixedAzizPotential : public PotentialBase  {

    public:
        FixedAzizPotential(const Container *);
        ~FixedAzizPotential();

        /* Return the sum of the Aziz 'interaction energy' between the supplied
         * particle and all fixed particles. */
        inline __attribute__((always_inline)) double V(const dVec &r);

        /* Return the gradient of the sum of the Aziz 'interaction energy' */
        inline __attribute__((always_inline)) dVec gradV(const dVec &r);

        /** Initial configuration corresponding to FixedAziz potential */
	DynamicArray<dVec,1> initialConfig(const Container*, MTRand &, const int); 

    private:
        AzizPotential aziz;                 // A copy of the aziz potential
	DynamicArray <dVec,1> fixedParticles;      // The location of the fixed particles
	DynamicArray <int,2> fixedBeadsInGrid;     // The local grid hash table
	DynamicArray <int,1> numFixedBeadsInGrid;  // The number of fixed particles in each grid box
        int numFixedParticles;              // The total number of fixed particles
        LookupTable *lookupPtr;             // A lookup table pointer
        double rc2;                         // A local copy of the potential cutoff squared

};

#if NDIM > 2
// ========================================================================  
// FixedPositionLJPotential Class
// ========================================================================  
/** 
 * @brief Returns Lennard-Jones potential between adatoms and fixed postions in FILENAME.
 *
 * Author: Nathan Nichols & Sutirtha Paul
 * Returns the potential energy resulting from interaction between adatom and fixed positions
 * given in FILENAME.
 */
class FixedPositionLJPotential: public PotentialBase  {

    public:
        FixedPositionLJPotential(const double, const double, const Container*);
        ~FixedPositionLJPotential();

        /* Return the sum of the Lennard-Jones potential between the supplied
         * particle and the fixed positions found in FILENAME. */
        inline __attribute__((always_inline)) double V(const dVec &r);

    private:
        const Container *boxPtr;
        // double Lz;
        // double Lx;
        // double Ly;
        // double Wallcz;
        // double Wallcy; 
        // double Wallcx;
        // double invWallWidth; 

        DynamicArray <dVec,1> fixedAtoms;     // The location of the fixed particles
        DynamicArray <int,1> atomType;        // The atom types 
        DynamicArray <double,1> sigma;        // The Lennard-Jones σ in Å
        DynamicArray <double,1> epsilon;      // The Lennard-Jones ϵ in K
                                                              
        int numFixedAtoms;               // The total number of fixed particles
        int numAtomTypes;                // The number of different atom types
};      
#endif

#if NDIM > 2
// ========================================================================  
// Excluded Volume Class (volume excluded w/ large potential)
// ========================================================================  
/**
 * Computes potential energy for Gasparini potential.
 * This is the case of gradV = 0 for all positions in cell.
 */
class Gasparini_1_Potential : public PotentialBase {
    public:
        Gasparini_1_Potential (double, double, const Container*);
        ~Gasparini_1_Potential ();

        /** The potential */
        inline __attribute__((always_inline)) double V(const dVec &r){ 
            double r2 = 0.0;
            if ((r[2] >= -excZ) && (r[2] <= excZ) && (r[1] >= -excY) && (r[1] <= excY))
                r2 = 1.0*V0;
            return r2;
        }

        /** The gradient of the potential. */
        inline __attribute__((always_inline)) dVec gradV(const dVec &) { return dVec{}; }

        /** Laplacian of the potential. */
        inline __attribute__((always_inline)) double grad2V(const dVec &r) { return 0.0; }

    
        /** Initial configuration corresponding to FixedAziz potential */
	DynamicArray<dVec,1> initialConfig(const Container*, MTRand &, const int); 

        /** get the exclusion lengths ay and az */
	DynamicArray<double,1> getExcLen(); 

        /* parameters needed for Gasp Potential_1 */
        const double excZ;      //half amt. of exclusion (z)
        const double excY;      //half amt. of exclusion (y)
        const double V0;        //scales the potential step
};
#endif

// ========================================================================  
// Hard Sphere Potential Class
// ========================================================================  
/** 
 * Computes the effective potential from the exact two-body density matrix
 * for hard spheres in 3D.  
 *
 * @see: S. Pilati, K. Sakkos, J. Boronat, J. Casulleras, and 
 *       S. Giorgini, Phys Rev A 74, 043621 (2006).
 */
class HardSpherePotential : public PotentialBase  {
    public:
        HardSpherePotential (double);
        ~HardSpherePotential ();

        /** The classical potential */
        inline __attribute__((always_inline)) double V(const dVec &r) { 
            return ((sqrt(dot(r,r)) <= a) ? LBIG : 0.0);
        }

        /** The effective potential */
        inline __attribute__((always_inline)) double V(const dVec &, const dVec &);
        inline __attribute__((always_inline)) double dVdlambda(const dVec &, const dVec &);
        inline __attribute__((always_inline)) double dVdtau(const dVec &, const dVec &);

    private:
        double a;               // The strength of the delta function
};


// ========================================================================
// 1D Delta Potential Class
// ========================================================================
/**
 * Computes the effective potential from the exact two-body density matrix
 * for delta interactions in 1D
 *       
 */
class Delta1DPotential : public PotentialBase  {
public:
    Delta1DPotential (double);
    ~Delta1DPotential ();
    
    /** The classical potential */
    inline __attribute__((always_inline)) double V(const dVec &r) {
        return (0.0);
    }
    
    /** The effective potential */
    inline __attribute__((always_inline)) double V(const dVec &, const dVec &);
    inline __attribute__((always_inline)) double dVdlambda(const dVec &, const dVec &);
    inline __attribute__((always_inline)) double dVdtau(const dVec &, const dVec &);
    
private:

    double erfCO;   //cutoff of erfc for asymptotic form
    double g;       // The strength of the delta function
    double l0;      // The COM kinetic length scale: 2\sqrt{\lambda\tau}
    double xi;      // The ratio l0/li
    double li;      // The COM interaction length scale: 2\lambda/g
    
    double xiSqOver2; // (xi^2)/2.0
    double xiSqrtPIOver2; // xi*\sqrt{\pi/2.0}
    
    /** Interaction weight and derivatives **/
    double Wint(double yt, double dxt);
    double dWdyt(double yt, double dxt);
    double dWdxi(double yt, double dxt);
    double dWddxt(double yt, double dxt);
};

// ========================================================================  
// Hard Rod Potential Class
// ========================================================================  
/** 
 * Computes the effective potential from the exact two-body density matrix
 * for hard rods in 1D.  
 *
 */
class HardRodPotential : public PotentialBase  {
    public:
        HardRodPotential (double);
        ~HardRodPotential ();

        /** The classical potential */
        inline __attribute__((always_inline)) double V(const dVec &r) { 
            return ((sqrt(dot(r,r)) <= a) ? LBIG : 0.0);
        }

        /** The effective potential */
        inline __attribute__((always_inline)) double V(const dVec &, const dVec &);
        inline __attribute__((always_inline)) double dVdlambda(const dVec &, const dVec &);
        inline __attribute__((always_inline)) double dVdtau(const dVec &, const dVec &);

    private:
        double a;               // The strength of the delta function
};

// ========================================================================  
// Carbon Nanotube Potential Class
// ========================================================================  
/** 
 * The smooth non-corregated version of the helium-carbon nanotube potential.
 * @see http://prb.aps.org/abstract/PRB/v62/i3/p2173_1
 */
//class CarbonNanotubePotential : public PotentialBase, public TabulatedPotential {
//  public:
//      CarbonNanotubePotential(const double);
//      ~CarbonNanotubePotential();
//
//      /** The cylindrically symmetric potential. */
//      double V(const dVec &r) {
//          int k = int(sqrt(r[0]*r[0] + r[1]*r[1])/dR);
//          if (k >= tableLength)
//              return extV[1];
//          else
//              return lookupV(k);
//      }
//
//      /* The gradient of the CNT potential */
//      dVec gradV(const dVec &);
//
//      /** Initial configuration corresponding to the CNT potential */
//      Array<dVec,1> initialConfig(const Container*, MTRand &, const int); 
//
//  private:
//      /* All the parameters needed for the LJ wall potential */
//      double density;
//      double sigmaHe,sigmaC;
//      double epsilonHe,epsilonC;
//
//      double R;       // Radius of the tube
//      double dR;      // Discretization for the lookup table
//
//      double minV;    // The minimum value of the potential
//
//      /* Used to construct the lookup tables */
//      double valueV (const double);               
//      double valuedVdr (const double);                    
//};

#if NDIM > 2
// ========================================================================  
// GraphenePotential Class
// ========================================================================  
/** 
 * \brief Returns van der Waals' potential between a helium adatom and a graphene sheet using summation in reciprocal space.
 *
 * Author: Nathan Nichols
 * Returns the potential energy resulting from a van der Waals' interaction
 * between a helium adatom and a fixed infinite graphene lattice
 */
class GraphenePotential: public PotentialBase  {

    public:
        GraphenePotential(const double, const double, const double, const double, const double);
        ~GraphenePotential();

        /* Return the sum of the van der Waals' interaction energy between the supplied
         * particle and the fixed graphene lattice. */
        inline __attribute__((always_inline)) double V(const dVec &r);

        /** Initial configuration corresponding to graphene-helium vdW potential */
	DynamicArray<dVec,1> initialConfig(const Container*, MTRand &, const int); 

    private:
        double sigma;
        double epsilon;
        double Lz;
        double Lzo2;

        double a1x;
        double a1y;
        double a2x;
        double a2y;

        double g1x;
        double g1y;
        double g2x;
        double g2y;

        double b1x;
        double b1y;
        double b2x;
        double b2y;

        double A;

};
#endif

#if NDIM > 2
// ========================================================================  
// GrapheneLUTPotential Class
// ========================================================================  
/** 
 * \brief Returns van der Waals' potential between a helium adatom and a graphene sheet using summation in reciprocal space.
 *
 * Author: Nathan Nichols & Adrian Del Maestro
 * Returns the potential energy resulting from a van der Waals' interaction
 * between a helium adatom and a fixed infinite graphene lattice
 */
class GrapheneLUTPotential: public PotentialBase  {

    public:
        GrapheneLUTPotential(const double, const double, const double, const double, const double, const Container*);
        ~GrapheneLUTPotential();

        /* Return the sum of the van der Waals' interaction energy between the supplied
         * particle and the fixed graphene lattice. */
        inline __attribute__((always_inline)) double V(const dVec &r);

        /* Return the gradient of the sum of the van der Waals' interaction energy between the supplied
         * particle and the fixed graphene lattice. */
        inline __attribute__((always_inline)) dVec gradV(const dVec &r);
        
        /** Initial configuration corresponding to graphene-helium vdW potential */
	DynamicArray<dVec,1> initialConfig(const Container*, MTRand &, const int); 

    private:
        double sigma;
        double epsilon;

        double a1x;
        double a1y;
        double a2x;
        double a2y;

        double g1x;
        double g1y;
        double g2x;
        double g2y;

        double b1x;
        double b1y;
        double b2x;
        double b2y;

        double A;

        double Lzo2;            ///< half the system size in the z-direction
        double Lz;              ///< The size of the system in the z-direction
        double zWall;           ///< The location of the onset of the "hard" wall 
        double invWallWidth;    ///< How fast the wall turns on.
        double V_zmin;          ///< A large potential value used for the cutoff
        
        double q = 2.0; // number of basis vectors
        
        static const int gnum = 3;
        static const int gtot = 16; //pow(gnum + 1,2);

        double dr = 1.0E-5;
        double zmin = 1.5;
        double zmax = 10.0;
        int tableLength;

        /* Array<int,2> karr; */
        /* double garr [gtot]; */

	DynamicArray<int,1> gMagID;    // g-magnitude lookup index
        
        /* The lookup tables */
	DynamicArray<double,2> vg;
	DynamicArray<double,2> gradvg;
};
#endif

#if NDIM > 2
// ========================================================================  
// GrapheneLUT3DPotential Class
// ========================================================================  
/** 
 * \brief Returns van der Waals' potential between a helium adatom and a graphene sheet using summation in reciprocal space.
 *
 * Author: Nathan Nichols
 * Returns the potential energy resulting from a van der Waals' interaction
 * between a helium adatom and a fixed infinite graphene lattice
 */
class GrapheneLUT3DPotential: public PotentialBase  {

    public:
        GrapheneLUT3DPotential(const std::string, const Container*);
        ~GrapheneLUT3DPotential();

        /* Return the sum of the van der Waals' interaction energy between the supplied
         * particle and the fixed graphene lattice. */
        inline __attribute__((always_inline)) double V(const dVec &);

        /* Return the gradient of the sum of the van der Waals' interaction energy between the supplied
         * particle and the fixed graphene lattice. */
        inline __attribute__((always_inline)) dVec gradV(const dVec &);
        inline __attribute__((always_inline)) double grad2V(const dVec &);
        
        /** Initial configuration corresponding to graphene-helium vdW potential */
	DynamicArray<dVec,1> initialConfig(const Container*, MTRand &, const int); 
	DynamicArray<dVec,1> initialConfig1(const Container*, MTRand &, const int); 

        void put_in_uc( dVec &, double, double);
        void cartesian_to_uc( dVec &, double, double, double, double);
        double trilinear_interpolation(DynamicArray<double,3>,dVec,double,double,double);
        double direct_lookup(DynamicArray<double,3>,dVec,double,double,double);

    private:
        double Lzo2;      ///< half the system size in the z-direction
        double zWall;     ///< The location of the onset of the "hard" wall 
        double invWallWidth; ///< How fast the wall turns on.
        
        /* spacing of the lookup tables */
        double dx;
        double dy;
        double dz;

        /* dimensions of the lookup tables */
        double cell_length_a;
        double cell_length_b;
        double zmin;
        double zmax;
        double V_zmin;
        
        /* transfer matrix */
        double A11;
        double A12;
        double A21;
        double A22;

	DynamicArray<double,3> V3d; // Potential lookup table
	DynamicArray<double,3> gradV3d_x; // gradient of potential x direction lookup table
	DynamicArray<double,3> gradV3d_y; // gradient of potential y direction lookup table
	DynamicArray<double,3> gradV3d_z; // gradient of potential z direction lookup table
	DynamicArray<double,3> grad2V3d; // Laplacian of potential
	DynamicArray<double,1> LUTinfo;

};

/****
 * Put inside unit cell
 ***/
inline __attribute__((always_inline)) void GrapheneLUT3DPotential::put_in_uc(dVec &r, double cell_length_a, double cell_length_b) {
    r[0] -= cell_length_a * floor(r[0]/cell_length_a);
    r[1] -= cell_length_b * floor(r[1]/cell_length_b);
}

/****
 * transfer from cartesian to unit cell coordinates
 ***/

inline __attribute__((always_inline)) void GrapheneLUT3DPotential::cartesian_to_uc( dVec &r, double A11, double A12, double A21, double A22) {
    double _x = A11 * r[0] + A12*r[1];
    double _y = A21 * r[0] + A22*r[1];
    r[0] = _x;
    r[1] = _y;
}

inline __attribute__((always_inline)) double GrapheneLUT3DPotential::trilinear_interpolation(DynamicArray<double,3> P,dVec r,double dx,double dy,double dz) {
    double x = r[0];
    double y = r[1];
    double z = r[2];

    double _xidx = floor(x/dx);
    double _yidx = floor(y/dy);
    double _zidx = floor(z/dz);
    int xidx = static_cast<int>(_xidx);
    int yidx = static_cast<int>(_yidx);
    int zidx = static_cast<int>(_zidx);

    double Delta_x = x/dx - _xidx;
    double Delta_y = y/dy - _yidx;
    double Delta_z = z/dz - _zidx;

    double c00 = P(xidx,yidx,zidx)*(1 - Delta_x) + P(xidx + 1,yidx,zidx)*Delta_x;
    double c01 = P(xidx,yidx,zidx + 1)*(1 - Delta_x) + P(xidx + 1,yidx,zidx + 1)*Delta_x;
    double c10 = P(xidx,yidx + 1,zidx)*(1 - Delta_x) + P(xidx + 1,yidx + 1,zidx)*Delta_x;
    double c11 = P(xidx,yidx + 1,zidx + 1)*(1 - Delta_x) + P(xidx + 1,yidx + 1,zidx + 1)*Delta_x;

    double c0 = c00*(1 - Delta_y) + c10*Delta_y;
    double c1 = c01*(1 - Delta_y) + c11*Delta_y;

    double c = c0*(1 - Delta_z) + c1*Delta_z;
    return c;
}

inline __attribute__((always_inline)) double GrapheneLUT3DPotential::direct_lookup(DynamicArray<double,3> P,dVec r,double dx,double dy,double dz) {
    double x = r[0];
    double y = r[1];
    double z = r[2];

    int xidx = static_cast<int>(x/dx);
    int yidx = static_cast<int>(y/dy);
    int zidx = static_cast<int>(z/dz);

    double c0 = P(xidx,yidx,zidx);
    return c0;
}
#endif

#if NDIM > 2
// ========================================================================  
// GrapheneLUT3DPotentialGenerate Class
// ========================================================================  
/** 
 * \brief FIXME Returns van der Waals' potential between a helium adatom and a graphene sheet using summation in reciprocal space.
 *
 * Author: Nathan Nichols
 * Returns the potential energy resulting from a van der Waals' interaction
 * between a helium adatom and a fixed infinite graphene lattice
 */
class GrapheneLUT3DPotentialGenerate: public PotentialBase  {

    public:
        GrapheneLUT3DPotentialGenerate(
                const double, const double, const double, const double,
                const double, const int, const int, const int, const int, const Container*);
        ~GrapheneLUT3DPotentialGenerate();
        
    private:

        double Lzo2;    // half the system size in the z-direction
        
        /* resolution of the lookup tables */
        int xres;
        int yres;
        int zres;

        /* dimensions of the lookup tables */
        /* double zmin; */
        double zmax;
	int k_max;
        /* double V_zmin; */
        
        double Vz_64( double, double, double, int );
        double Vz_64( double, double, double );
        double Vz_64( double, double );
        
        double gradVz_x_64( double, double, double, int );
        double gradVz_x_64( double, double, double );
        double gradVz_x_64( double, double );
        
        double gradVz_y_64( double, double, double, int );
        double gradVz_y_64( double, double, double );
        double gradVz_y_64( double, double );
        
        
        double gradVz_z_64( double, double, double, int );
        double gradVz_z_64( double, double, double );
        double gradVz_z_64( double, double );

        double grad2Vz_64( double, double, double, int );
        double grad2Vz_64( double, double, double );
        double grad2Vz_64( double, double );
        
        double Vg_64( double, double, double, double, double, double, double,
                double, double, double );
        double Vg_64( double, double, double, double, double, double, double,
                double, double );
        
        double gradVg_x_64( double, double, double, double, double, double, double,
                double, double, double );
        double gradVg_x_64( double, double, double, double, double, double, double,
                double, double );

        double gradVg_y_64( double, double, double, double, double, double, double,
                double, double, double );
        double gradVg_y_64( double, double, double, double, double, double, double,
                double, double );

        double gradVg_z_64( double, double, double, double, double, double, double,
                double, double, double );
        double gradVg_z_64( double, double, double, double, double, double, double,
                double, double );

        double grad2Vg_64( double, double, double, double, double, double, double,
                double, double, double );
        double grad2Vg_64( double, double, double, double, double, double, double,
                double, double );

        double V_64( double, double, double, double, double, double,
                std::array<double,2>, std::array<double,2>,
                std::array<double,2>, std::array<double,2>, DynamicArray<int,1>,
                DynamicArray<int,1>, DynamicArray<double,1> );

        double gradV_x_64( double, double, double, double, double, double,
                std::array<double,2>, std::array<double,2>,
                std::array<double,2>, std::array<double,2>, DynamicArray<int,1>,
                DynamicArray<int,1>, DynamicArray<double,1> );

        double gradV_y_64( double, double, double, double, double, double,
                std::array<double,2>, std::array<double,2>,
                std::array<double,2>, std::array<double,2>, DynamicArray<int,1>,
                DynamicArray<int,1>, DynamicArray<double,1> );

        double gradV_z_64( double, double, double, double, double, double,
                std::array<double,2>, std::array<double,2>,
                std::array<double,2>, std::array<double,2>, DynamicArray<int,1>,
                DynamicArray<int,1>, DynamicArray<double,1> );

        double grad2V_64( double, double, double, double, double, double,
                std::array<double,2>, std::array<double,2>,
                std::array<double,2>, std::array<double,2>, DynamicArray<int,1>,
                DynamicArray<int,1>, DynamicArray<double,1> );
        
        std::tuple<
            std::array<double,2>, std::array<double,2>,
            std::array<double,2>, std::array<double,2>,
            std::array<double,2>, std::array<double,2>
                > get_graphene_vectors();
        std::tuple<
            std::array<double,2>, std::array<double,2>,
            std::array<double,2>, std::array<double,2>,
            std::array<double,2>, std::array<double,2>
                > get_graphene_vectors( double );
        std::tuple<
            std::array<double,2>, std::array<double,2>,
            std::array<double,2>, std::array<double,2>,
            std::array<double,2>, std::array<double,2>
                > get_graphene_vectors( double, double, double );
        std::tuple<
            std::array<double,2>, std::array<double,2>,
            std::array<double,2>, std::array<double,2>,
            std::array<double,2>, std::array<double,2>
                > get_graphene_vectors_old( double, double, double );
        std::tuple< DynamicArray<int,1>, DynamicArray<int,1>, DynamicArray<double,1>
            > get_g_magnitudes( std::array<double,2>, std::array<double,2> );

        template <class T> double calculate_magnitude( T vec ) {
            return sqrt(dot(vec,vec));
        }
        
        template <class T1, class T2> double calculate_angle(
                T1 x, T2 y, double magnitude_x, double magnitude_y ) {
            return acos(dot(x,y)/magnitude_x/magnitude_y);
        }
        
        template <class T1, class T2> double calculate_angle( T1 x, T2 y ) {
            return calculate_angle(x,y,calculate_magnitude(x),calculate_magnitude(y));
        }

        void calculate_V3D_64(
                DynamicArray<double,3>, DynamicArray<double,2>, DynamicArray<double,2>,
                DynamicArray<double,1>, double, double,
                double, std::array<double,2>, std::array<double,2>,
                std::array<double,2>, std::array<double,2>,
                DynamicArray<int,1>, DynamicArray<int,1>,
                DynamicArray<double,1> );

        void calculate_gradV3D_x_64(
                DynamicArray<double,3>, DynamicArray<double,2>, DynamicArray<double,2>,
                DynamicArray<double,1>, double, double,
                double, std::array<double,2>, std::array<double,2>,
                std::array<double,2>, std::array<double,2>,
                DynamicArray<int,1>, DynamicArray<int,1>,
                DynamicArray<double,1> );

        void calculate_gradV3D_y_64(
                DynamicArray<double,3>, DynamicArray<double,2>, DynamicArray<double,2>,
                DynamicArray<double,1>, double, double,
                double, std::array<double,2>, std::array<double,2>,
                std::array<double,2>, std::array<double,2>,
                DynamicArray<int,1>, DynamicArray<int,1>,
                DynamicArray<double,1> );

        void calculate_gradV3D_z_64(
                DynamicArray<double,3>, DynamicArray<double,2>, DynamicArray<double,2>,
                DynamicArray<double,1>, double, double,
                double, std::array<double,2>, std::array<double,2>,
                std::array<double,2>, std::array<double,2>,
                DynamicArray<int,1>, DynamicArray<int,1>,
                DynamicArray<double,1> );

        void calculate_grad2V3D_64(
                DynamicArray<double,3>, DynamicArray<double,2>, DynamicArray<double,2>,
                DynamicArray<double,1>, double, double,
                double, std::array<double,2>, std::array<double,2>,
                std::array<double,2>, std::array<double,2>,
                DynamicArray<int,1>, DynamicArray<int,1>,
                DynamicArray<double,1> );

        std::pair<double, double> get_z_min_V_min( 
                double, double, double, double, double,
                std::array<double,2>, std::array<double,2>,
                std::array<double,2>, std::array<double,2>,
                DynamicArray<int,1>, DynamicArray<int,1>, DynamicArray<double,1> );

        std::pair<double, double> get_z_V_to_find( 
                double, double, double,
                std::array<double,2>, std::array<double,2>,
                std::array<double,2>, std::array<double,2>,
                DynamicArray<int,1>, DynamicArray<int,1>, DynamicArray<double,1> );

	DynamicArray<double,3> get_V3D(
                double, double, double, int, int, int, double, double ); 
        std::pair<DynamicArray<double,3> , DynamicArray<double,1>> get_V3D(
                double, double, double, int, int, int, double );
        std::tuple<
            DynamicArray<double,3>, DynamicArray<double,3>, DynamicArray<double,3>,
            DynamicArray<double,3>, DynamicArray<double,3>, DynamicArray<double,2>,
            DynamicArray<double,2>, DynamicArray<double,1>
            > get_V3D_all( double, double, double, int, int, int, double );
};
#endif
#if NDIM > 2
// ========================================================================  
// LeeBenzenePotential Class
// ========================================================================  
/** 
 * @brief Implements an accurate ab-inito calculation for the Benzene-He potential energy
 * @see https://doi.org/10.1063/1.1628217
 *
 * Author: Sutirtha Paul
 * 
 */
class LeeBenzenePotential: public PotentialBase  {

    public:
        LeeBenzenePotential(const Container*);
        ~LeeBenzenePotential();

        /* Return the sum of the Lennard-Jones potential between the supplied
         * particle and the fixed positions found in FILENAME. */
        inline __attribute__((always_inline)) double V(const dVec &r);

    private:
        const Container *boxPtr;
        double Lz;
        double Lx;
        double Ly;
        double Wallcz;
        double Wallcy; 
        double Wallcx;
        double invWallWidth;
        const double r0 = 7.220527;
	const float a = 0.572446;
	const double c3 = 0.058387;
        const double c4 = -6.914851;
        const double c5 = -2.083808;
        const double c6 = 85.429568;
        const double c11 = -16.956997;
        const double c12 = 0.568429;
        const double c22 = 9.073184;
        const double c13 = 1.821920;
        const double c14 = 1.874776;
        const double c23 = -0.098396;
        const double c111 = -3.235122;
        const double c112 = -2.631504;
        const double c122 = 1.423777;
        const double c113 = -2.117152;
        
	DynamicArray<std::array<double,4>,1> fixedParticles;      // The location of the fixed particles
        DynamicArray<std::array<double,2>,1> atomArray;           // The interaction parameters for the fixed particles 
        int numFixedParticles;              // The total number of fixed particles
        int typesofatoms;                   // The various types of atom species

	double w(double);
        double wtilde(double);
        double V2(double);
        double V3(double,double);
        double V4(double,double,double);

};      

inline __attribute__((always_inline)) double LeeBenzenePotential::w(double r) {
        return 1 - exp(-a*(r - r0));
}

inline __attribute__((always_inline)) double LeeBenzenePotential::wtilde(double r) {
        if (r >= r0)
            return w(r);
        else
	    return 0;
}

inline __attribute__((always_inline)) double LeeBenzenePotential::V2(double r) {
    double fval = w(r)*w(r) + c6*(pow(wtilde(r),6));
    fval += c3*(pow(w(r),3)) + c4*(pow(w(r),4)) + c5*(pow(w(r),5));
    return fval;
}

inline __attribute__((always_inline)) double LeeBenzenePotential::V3(double rk, double rl) {
    double fval = c11*w(rk)*w(rl) + c22*(pow(w(rk),2))*(pow(w(rl),2));
    fval += c12*(w(rk)*(pow(w(rl),2)) + w(rl)*(pow(w(rk),2)));
    fval += c13*(w(rk)*(pow(w(rl),3)) + w(rl)*(pow(w(rk),3)));
    fval += c14*(w(rk)*(pow(w(rl),4)) + w(rl)*(pow(w(rk),4)));
    fval += c23*((pow(w(rk),2)*pow(w(rl),3)) + (pow(w(rl),2)*pow(w(rk),3)));
    return fval;
}

inline __attribute__((always_inline)) double LeeBenzenePotential::V4(double rk, double rl, double rm) {
    double fval = c111*w(rk)*w(rl)*w(rm);
    fval += c122*(w(rk)*pow(w(rl),2)*pow(w(rm),2) + w(rl)*pow(w(rm),2)*pow(w(rk),2) + w(rm)*pow(w(rk),2)*pow(w(rl),2));
    fval += c112*(w(rk)*w(rl)*pow(w(rm),2) + w(rk)*pow(w(rl),2)*w(rm) + pow(w(rk),2)*w(rl)*w(rm));
    fval += c113*(w(rk)*w(rl)*pow(w(rm),3) + w(rk)*pow(w(rl),3)*w(rm) + pow(w(rk),3)*w(rl)*w(rm));
    return fval;
}

#endif
#if NDIM > 2
// ========================================================================  
// ShirkovBenzene Potential Class
// ========================================================================  
/** 
 * @brief Implements an accurate ab-inito calculation for the Benzene-He potential energy
 * @see https://doi.org/10.1021/acs.jpca.4c01491
 *
 * Author: Sutirtha Paul
 * 
 */
class ShirkovBenzene: public PotentialBase  {

    public:
        ShirkovBenzene(const Container*);
        ~ShirkovBenzene();

        /* Return the sum of the Lennard-Jones potential between the supplied
         * particle and the fixed positions found in FILENAME. */
        inline __attribute__((always_inline)) double V(const dVec &r);

    private:
        const Container *boxPtr;
        double Lz;
        double Lx;
        double Ly;
        double Wallcz;
        double Wallcy; 
        double Wallcx;
        double invWallWidth;
        
	const double re       = 3.565710000e+00;
	const double ae       = 1.485672859e+00;
	const double c0       = -1.312335657e+02;
	const double v0       = 1.592094647e+03;
	const double c3       = -4.497880474e+00;
	const double c4       = -1.160510603e+00;
	const double c5       = -3.574920218e+00;
	const double c6       = -1.169708775e-01;
	const double c7       = 0.0;
	const double c8       = 0.0;
	const double c12      = -2.894238936e-01;
	const double c112     = 1.465193128e+00;
	const double c1122    = -8.310247790e-01;
	const double c1112    = 2.158143635e+00;
	const double c11112   = 7.135910167e-01;
	const double c11122   = 0.0;
	const double c111112  = 0.0;
	const double c111122  = 0.0;
	const double c111222  = 0.0;
	const double c123     = -1.036553121e+00;
	const double c1123    = -6.340217244e-01;
	const double c11223   = -6.320195633e-02;
	const double c11123   = 0.0;
	const double c112233  = 0.0;

	const double reh     = 3.000000000e+00;
	const double aeh     = 1.920910443e+00;
	const double v0h     = 3.920606202e+02;
	const double d3       = -1.789976234e+00;
	const double d4       = 1.120849120e+00;
	const double d5       = -4.676333910e-01;
	const double d6       = -4.181953639e-02;
	const double d12      = -1.161279955e+00;
	const double d112     = 1.365346870e+00;
	const double d1122    = -1.750118034e+00;
	const double d1112    = 0.0;
	const double dc12   = 0.0;
	const double dc112   = 0.0;

	const double r0       = 5.765366674e+00;
	const double gama     = 1.671477473e+01;
	const double DE       = 8.900000000e+01;

	const double x1       = 2.290626821e+05;
	const double x2       = 4.069217828e+06;
	const double x3       = 9.542607077e+07;
	const double x4       = 0.0;
	const double x5       = -4.775615288e+04;
	const double x6       = -1.263899296e+07;
	const double x7       = -2.959330054e+08;
	const double x8       = 0.0;
	const double x9       = 9.511488779e+05;
	const double x10      = 2.551023972e+08;
	const double x11      = 0.0;
	const double x12      = 1.883887454e+08;

	std::vector<std::vector<double>> generate_benzene_geometry();
	//blitz::Array <blitz::TinyVector<double,4>,1> fixedParticles;      // The location of the fixed particles
        //int numFixedParticles;              // The total number of fixed particles
        //int typesofatoms;                   // The various types of atom species

}; 

inline __attribute__((always_inline)) std::vector<std::vector<double>> ShirkovBenzene::generate_benzene_geometry() {
    double rc = 1.3915;   // Carbon ring radius (Angstrom)
    double rh = 1.080;    // C–H bond length (Angstrom)
    double rh1 = rc + rh; // Hydrogen distance from center

    // Constants for hexagonal arrangement
    std::array<double, 6> ax = {1.0, 0.5, -0.5, -1.0, -0.5, 0.5};
    std::array<double, 6> ay = {
        0.0,
        std::sqrt(3.0) / 2.0,
        std::sqrt(3.0) / 2.0,
        0.0,
        -std::sqrt(3.0) / 2.0,
        -std::sqrt(3.0) / 2.0
    };

    // 3 x 12 coordinate matrix initialized to 0
    std::vector<std::vector<double>> coords(3, std::vector<double>(12, 0.0));

    // Carbon atoms
    for (int i = 0; i < 6; i++) {
        coords[0][i] = ax[i] * rc;  // x
        coords[1][i] = ay[i] * rc;  // y
        // z stays 0
    }

    // Hydrogen atoms
    for (int i = 0; i < 6; i++) {
        coords[0][i + 6] = ax[i] * rh1;  // x
        coords[1][i + 6] = ay[i] * rh1;  // y
        // z stays 0
    }

    return coords;
}
#endif
#if NDIM > 2
// ========================================================================
// GaussianProcess Potential Class
// ========================================================================
/** 
 * @brief Implements a gaussian process trained on accurate ab-inito calculation for the Benzene-He potential energy
 *
 * Author: Sutirtha Paul
 * 
 */
class GPHeBenzenePotential: public PotentialBase  {

    public:
        GPHeBenzenePotential(const Container*, const po::variables_map &);
        ~GPHeBenzenePotential();

        inline __attribute__((always_inline)) double V(const dVec &r);
        inline __attribute__((always_inline)) void V(const dVec*, double*, int);
#ifdef USE_GPU
        void gpuV(const double*, double*, int);
#endif
    private:
        GaussianProcessPotential GPPotential;

        const Container *boxPtr;
        double Lz;
        double Lx;
        double Ly;
        double Wallcz;
        double Wallcy;
        double Wallcx;
        double invWallWidth;
	//Long-range parameters
	const double r0       = 5.765366674e+00;
	const double gama     = 1.671477473e+01;

	const double x1       = 2.290626821e+05;
	const double x2       = 4.069217828e+06;
	const double x3       = 9.542607077e+07;
	const double x4       = 0.0;
	const double x5       = -4.775615288e+04;
	const double x6       = -1.263899296e+07;
	const double x7       = -2.959330054e+08;
	const double x8       = 0.0;
	const double x9       = 9.511488779e+05;
	const double x10      = 2.551023972e+08;
	const double x11      = 0.0;
	const double x12      = 1.883887454e+08;

	double deg2rad(double);
	double rad2deg(double);
	double long_range(double,double,double);
};
inline __attribute__((always_inline)) double GPHeBenzenePotential::deg2rad(double ang) {
    return (ang*M_PI)/180.0;
}
inline __attribute__((always_inline)) double GPHeBenzenePotential::rad2deg(double ang) {
    return (ang*180.0)/M_PI;
}
inline __attribute__((always_inline)) double GPHeBenzenePotential::long_range(double x, double y, double z) {
    double rnorm = std::sqrt(x * x + y * y + z * z);
    double ph0 = std::atan2(y, x);
    double th0 = std::acos(z / rnorm);
    double tt = std::cos(th0);
    double fi = ph0;

    // Angular dependence
    double normm1 = std::sqrt(5.0) / 5.0;
    double normm2 = std::sqrt(13.0) / 13.0;
    double normm3 = 0.1792151994e-4;

    double p1 = -x1 / std::pow(rnorm, 6) - x2 / std::pow(rnorm, 8) - x3 / std::pow(rnorm, 10) - x4 / std::pow(rnorm, 12);
    double p2 = -normm1 * (1.5 * tt * tt - 0.5) *
                (x5 / std::pow(rnorm, 6) + x6 / std::pow(rnorm, 8) + x7 / std::pow(rnorm, 10) + x8 / std::pow(rnorm, 12));
    double p3 = -normm2 * (3.0 / 8.0 + 35.0 / 8.0 * std::pow(tt, 4) - 15.0 / 4.0 * tt * tt);
    double p4 = (x9 / std::pow(rnorm, 8) + x10 / std::pow(rnorm, 10) + x11 / std::pow(rnorm, 12));
    double p5 = -normm3 * 10395.0 * std::pow(1 - tt * tt, 3) * std::cos(6 * fi) * x12 / std::pow(rnorm, 10);
    double vdw = p1 + p2 + p3 * p4 + p5;

    return vdw;
}
#endif

#endif // POTENTIAL_H
