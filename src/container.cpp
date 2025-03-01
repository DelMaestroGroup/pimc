/** 
 * @file container.cpp
 * @author Adrian Del Maestro
 * 
 * @brief Holds all possible container types.
 */

#include "container.h"
#include "constants.h"

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// CONTAINER CLASS -----------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/**************************************************************************//**
 * Initialize all variables.
******************************************************************************/
Container::Container() {
    side     = dVec{};
    sideInv  = dVec{};
    sideInv2 = dVec{};
    pSide    = dVec{};
    periodic.fill(1);
    maxSep   = 0.0;
    volume   = 0.0;
    rcut2    = 0.0;
    name     = "";
    fullyPeriodic = true;

    /* Determine the number of grid boxes in the lookup table */
    numGrid = 1;
    for (int i = 0; i < NDIM; i++)
        numGrid *= NGRIDSEP;
}


/**************************************************************************//**
 * Empty destructor.
******************************************************************************/
Container::~Container() {
}

#if NDIM > 1
/**************************************************************************//**
 *  Given a grid box number, return the associated radius
 *
 *  @param n The grid index
 *  @return The radius squared x^2 + y^2 of the grid box
******************************************************************************/
double Container::gridRadius2(const int n) const {
    iVec _gridIndex;
    for (int i = 0; i < NDIM; i++) {
        int scale = 1;
        for (int j = i+1; j < NDIM; j++) 
            scale *= NGRIDSEP;
        _gridIndex[i] = (n/scale) % NGRIDSEP;
        PIMC_ASSERT(_gridIndex[i]<NGRIDSEP);
    }

    double r2 = 0.0;
    for (int i = 0; i < 2; i++) {
        double ri = -0.5*side[i] + (_gridIndex[i] + 0.5)*gridSize[i];
        r2 += ri*ri;
    }
    return r2;
}
#endif

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// PRISM CLASS ---------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/**************************************************************************//**
 *  Create a NDIM-dimensional hyperprism given density and number of particles.
 * 
 *  We create a NDIM hyperprism with periodic boundary conditions in all
 *  dimensions with a side that is set by the volume, given by the particle
 *  density in atoms/A^3 and the number of particles.
******************************************************************************/
Prism::Prism(double density, int numParticles) {

    /* Here we can only create a cube */
    /* Setup the cube size in each dimension */
    side.fill(pow(1.0*numParticles / density, 1.0/(1.0*NDIM)));
    sideInv = 1.0/side;
    sideInv2 = 2.0*sideInv;

    rcut2 = 0.25*side[NDIM-1]*side[NDIM-1];

    /* The hyper cube has periodic boundary conditions */
    periodic.fill(1); 
    pSide = side;
    fullyPeriodic = true;

    /* Compute the maximum possible separation possible inside the box */
    maxSep = sqrt(dot(side/(periodic + 1u),side/(periodic + 1u)));

    /* Calculate the volume of the cube */
    volume = product(side);

    /* The grid size for the lookup table */
    gridSize = side/NGRIDSEP;

    name = "Prism";
}

/**************************************************************************//**
 *  Create a NDIM-dimensional hyperprism given the edge sizes.
 * 
 *  We create a NDIM hyperprism with periodic boundary conditions in all
 *  dimensions with sides that are set at teh command line.
******************************************************************************/
Prism::Prism(const dVec &_side, const std::array<unsigned int, NDIM> &_periodic) {

    /* Setup the cube size in each dimension */
    side = _side;
    sideInv = 1.0/side;
    sideInv2 = 2.0*sideInv;

    rcut2 = 0.25*side[NDIM-1]*side[NDIM-1];

    /* Setup the periodic boundary conditions */
    periodic = _periodic; 
    pSide = periodic*side;

    /* are there any non-periodic boundary conditions? */
    static constexpr std::array<unsigned int, NDIM> fullPeriodicCondition = make_array<unsigned int, NDIM>(1u);
    fullyPeriodic = periodic == fullPeriodicCondition;

    /* Compute the maximum possible separation possible inside the box */
    maxSep = sqrt(dot(side/(periodic + 1u),side/(periodic + 1u)));

    /* Calculate the volume of the cube */
    volume = product(side);

    /* The grid size for the lookup table */
    gridSize = side/NGRIDSEP;

    name = "Prism";
}

/**************************************************************************//**
 * Empty destructor.
******************************************************************************/
Prism::~Prism() {
}

/**************************************************************************//**
 *  Return a random position inside the cube.
******************************************************************************/
dVec Prism::randPosition(MTRand &random) const {
    dVec randPos;
    for (int i = 0; i < NDIM; i++)
        randPos[i] = side[i]*(-0.5 + random.randExc());
    return randPos;
}

/**************************************************************************//**
 *  Return a random position close to the supplied one.
******************************************************************************/
dVec Prism::randUpdate (MTRand &random, const dVec &pos) const {
    dVec randPos;
    randPos = pos;
    for (int i = 0; i < NDIM; i++) 
        randPos[i] += 4.0*constants()->displaceDelta()*(-0.5 + random.rand());
    putInside(randPos);

    /* Make sure we don't move a particle outside the box */
    if (!fullyPeriodic) {
        for (int i = 0; i < NDIM; i++) {
            if (!periodic[i]) {
                if (randPos[i] >= 0.5*side[i])
                    randPos[i] = 0.5*side[i] - 2*EPS;
                if (randPos[i] < -0.5*side[i]) 
                    randPos[i] = -0.5*side[i] + 2*EPS;
            }
        }
    }
    return randPos;
}

/**************************************************************************//**
 *  Given a particle position, return a single integer which maps to a 
 *  unique grid position.  Used for correlating estimators with individual
 *  particle positions.
 *
 *  @param pos The particle position
 *
 *  @return An integer representing a grid box number
******************************************************************************/
int Prism::gridIndex(const dVec &pos) const {

    int gNumber = 0;
    for (int i = 0; i < NDIM; i++) {  
        int scale = 1;
        for (int j = i+1; j < NDIM; j++) 
            scale *= NGRIDSEP;
        gNumber += scale * 
            static_cast<int>(abs( pos[i] + 0.5*side[i] - EPS ) / (gridSize[i] + EPS));
    }
    PIMC_ASSERT(gNumber<numGrid);
    return gNumber;
}

/**************************************************************************//**
 *  Given a grid index, return the hyper volume of the associated grid box.
 *
 *  @param n The grid index
 *  @return The hyper volume of the grid box
******************************************************************************/
double Prism::gridBoxVolume(const int n) const {
    return std::accumulate(gridSize.begin(), gridSize.end(), 1.0, std::multiplies<>());

}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// CYLINDER CLASS ------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/**************************************************************************//**
 *  Create a cylinder given density, radius and number of particles.
 * 
 *  We create a 3d cylinder (otherwise exit) with periodic boundary conditions
 *  only in the z-direction.  The radius is fixed, and the length is determined
 *  by the density in atoms/A^3.
******************************************************************************/
Cylinder::Cylinder(const double _rho, const double radius, const int numParticles) {

    /* We can only make a cylinder in 3 dimensions */
    if (NDIM != 3) {
        std::cerr << "You can only create a cylinder in 3 dimensions, change NDIM!" 
            << std::endl;
        exit(EXIT_FAILURE);
    }
    else {
        /* We can compute the linear length from the density, number of particles
         * and radius */
        double L = (1.0*numParticles)/(_rho * M_PI * radius * radius);

        /* We check to make sure that our aspect ratio is at least 2:1 */
        if (L < 2.0*radius) 
            std::cerr << "L:r is smaller than 2:1!" << std::endl;

        /* Setup the prism size in each of the three dimensions which the 
         * cylinder will be inscribed insde of.  We make it 2R X 2R X L  */
        side[0] = side[1] = 2.0 * radius;
        side[2] = L;

        rcut2 = 0.25*side[2]*side[2];

        sideInv = 1.0/side;
        sideInv2 = 2.0*sideInv;

        /* setup which dimensions have periodic boundary conditions, 1.0 for yes 
         * and 0.0 for no. */
        periodic[0] = periodic[1] = 0;
        periodic[2] = 1;

        pSide = periodic*side;

        /* Compute the maximum possible separation possible inside the box */
        maxSep = sqrt(dot(side/(periodic + 1),side/(periodic + 1)));

        /* Compute the cylinder volume. We use the radius here instead of the actual
         * side.  This is the 'active' volume */
        volume = M_PI*radius*radius*L;

        name = "Cylinder";

        /* The grid size for the lookup table */
        gridSize = side/NGRIDSEP;

    } // end else
}

/**************************************************************************//**
 *  Create a cylinder given its radius and length.
 * 
 *  We create a 3d cylinder (otherwise exit) with periodic boundary conditions
 *  only in the z-direction.  The radius and length are given in angstroms.
 *  @param L The length of the cylinder
 *  @param radius The radius of the cylinder
******************************************************************************/
Cylinder::Cylinder(const double radius, const double L) {

    /* We can only make a cylinder in 3 dimensions */
    if (NDIM != 3) {
        std::cerr << "You can only create a cylinder in 3 dimensions, change NDIM!" 
            << std::endl;
        exit(EXIT_FAILURE);
    }
    else {
        /* We check to make sure that our aspect ratio is at least 2:1 */
        if (L < 2.0*radius)
            std::cerr << "L:r is smaller than 2:1!" << std::endl;

        /* Setup the prism size in each of the three dimensions which the 
         * cylinder will be inscribed insde of.  We make it 2R X 2R X L 
         * to account for */
        side[0] = side[1] = 2.0 * radius;
        side[2] = L;

        rcut2 = 0.25*side[2]*side[2];

        sideInv = 1.0/side;
        sideInv2 = 2.0*sideInv;

        /* setup which dimensions have periodic boundary conditions, 1 for yes 
         * and 0 for no. */
        periodic[0] = periodic[1] = 0;
        periodic[2] = 1;

        /* Compute the maximum possible separation possible inside the box */
        maxSep = sqrt(dot(side/(periodic + 1.0),side/(periodic + 1.0)));

        pSide = periodic*side;

        /* Compute the cylinder volume. We use the radius here instead of the actual
         * side.  This is the 'active' volume */
        volume = M_PI*radius*radius*L;

        name = "Cylinder";

        /* The grid size for the lookup table */
        gridSize = side/NGRIDSEP;

    } // end else
}

/**************************************************************************//**
 *  Destructor.
******************************************************************************/
Cylinder::~Cylinder() {
}

/**************************************************************************//**
 *  Return a random position inside the cylinder.
 *  @see http://extremelearning.com.au/how-to-generate-uniformly-random-points-on-n-spheres-and-n-balls/
******************************************************************************/
dVec Cylinder::randPosition(MTRand &random) const {

    dVec randPos;
    double r = 0.5*side[0]*sqrt(random.randExc());
    double phi = random.randExc(2.0*M_PI);
    randPos[0] = r * cos(phi);
    randPos[1] = r * sin(phi);
    randPos[2] = side[2]*(-0.5 + random.randExc());

    return randPos;
}

/**************************************************************************//**
 *  Return a random position close to the supplied one.
******************************************************************************/
dVec Cylinder::randUpdate (MTRand &random, const dVec &pos) const {
    dVec randPos;
    if (random.rand() < 0.5) 
        randPos = randUpdateJumpShell(random,pos);
    else 
        randPos = randUpdateSmall(random,pos);
    return randPos;
}

/**************************************************************************//**
 *  Return a random position close to the supplied one. Shell Jump
******************************************************************************/
dVec Cylinder::randUpdateJumpShell (MTRand &random, const dVec &pos) const {
    dVec randPos;
    randPos = pos;
    double theta = atan2(pos[1],pos[0]);
    double oldr = sqrt(pos[0]*pos[0] + pos[1]*pos[1]);
    double newr;

    if (random.rand() > 0.5)
        newr = oldr + (1.00 + 2.50*random.rand());
    else
        newr = oldr - (1.00 + 2.50*random.rand());

    if (newr < 0.0)
        newr *= -1.0;

    double ranTheta = M_PI*(-0.05 + 0.1*random.rand());
    randPos[0] = newr*cos(theta + ranTheta);
    randPos[1] = newr*sin(theta + ranTheta);
    randPos[2] += constants()->displaceDelta()*(-0.5 + random.rand());

    putInside(randPos);
    return randPos;
}

// THIS IS WHAT WAS HERE BEFORE
/**************************************************************************//**
 *  Return a random position close to the supplied one. Same Shell.
******************************************************************************/
dVec Cylinder::randUpdateSmall (MTRand &random, const dVec &pos) const {
    dVec randPos;
    randPos = pos;
    for (int i = 0; i < NDIM; i++) 
        randPos[i] += constants()->displaceDelta()*(-0.5 + random.rand());
    putInside(randPos);
    return randPos;
}
/**************************************************************************//**
 *  Make sure that a supplied vector is put inside the cylinder.
******************************************************************************/
void Cylinder::putInside(dVec &r) const {

    /* For the x and y axis, any outside position is placed near the boundary */
    for (int i = 0; i < NDIM-1; i++) {
        if (r[i] >= 0.5*side[i])
            r[i] = 0.5*side[i] - EPS;
        if (r[i] < -0.5*side[i]) 
            r[i] = -0.5*side[i] + EPS;
    }

    /* Now place the z-coordinate in PBC */
    r[2] -= (r[2] >= 0.5*side[2])*side[2];
    r[2] += (r[2] < -0.5*side[2])*side[2];
}

/**************************************************************************//**
 *  Given a particle position, return a single integer which maps to a 
 *  unique grid position.  Used for correlating estimators with individual
 *  particle positions.
 *
 *  @param pos The particle position
 *
 *  @return An integer representing a grid box number
******************************************************************************/
int Cylinder::gridIndex(const dVec &pos) const {

    // /* Get the r and theta components */
    // double r = sqrt(pos[0]*pos[0]+ pos[1]*pos[1]);
    // double theta = atan2(pos[1],pos[0]);
    // theta += (theta < 0.0)*2.0*M_PI;

    // /* Get the 3d vector index */
    // iVec grid;
    // grid[0] = static_cast<int>(abs(r-EPS)/(gridSize[0]+EPS));
    // grid[1] = static_cast<int>(abs(theta-EPS)/(gridSize[1]+EPS));
    // grid[2] = static_cast<int>(abs(pos[2] + 0.5*side[2] - EPS )/(gridSize[2] + EPS));

    // /* return the flattened index */
    // int gNumber = grid[0]*NGRIDSEP*NGRIDSEP + grid[1]*NGRIDSEP + grid[2];

    // PIMC_ASSERT(gNumber<numGrid);
    // return gNumber;

    int gNumber = 0;
    for (int i = 0; i < NDIM; i++) {  
        int scale = 1;
        for (int j = i+1; j < NDIM; j++) 
            scale *= NGRIDSEP;
        gNumber += scale * 
            static_cast<int>(abs( pos[i] + 0.5*side[i] - EPS ) / (gridSize[i] + EPS));
    }
    PIMC_ASSERT(gNumber<numGrid);

    return gNumber;
}

/**************************************************************************//**
 *  Given a grid index, return the hyper volume of the associated grid box.
 *
 *  @param n The grid index
 *  @return The hyper volume of the grid box
******************************************************************************/
double Cylinder::gridBoxVolume(const int n) const {

    /* Get the first grid box index */
    return std::accumulate(gridSize.begin(), gridSize.end(), 1.0, std::multiplies<>());
}

