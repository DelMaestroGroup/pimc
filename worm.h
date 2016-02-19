/**
 * @file worm.h
 * @author Adrian Del Maestro
 * @date 12.08.2008
 * @brief Worm class definition.
 */

#ifndef WORM_H 
#define WORM_H

#include "common.h"
#include "constants.h"

class Path;
// ========================================================================  
// Worm Class
// ========================================================================  
/** 
 * Contains information on the worm.
 * In a grand-canonical or worm simulation, one world line may lose its
 * periodicity in imaginary time.  This class holds all parameters and
 * methods used to describe worm configurations.
 */

class Worm {

	public:
		Worm(int);
		~Worm();

		beadLocator head;			///< The coordinates of the worm head
		beadLocator tail;			///< The coordinates of the worm tail
		beadLocator special1;		///< Special bead, used in move updates
		beadLocator special2;		///< Special bead, used in move updates
		double maxWormCost;			///< The maximum 'cost' of inserting a worm
		dVec sep;					///< The spatial separation between head and tail
		int length;					///< The length of the worm
		int gap;					///< numTimeSlices - length
		bool isConfigDiagonal;		///< Stores the diagonality of the configuration

		/* Safe get methods for beads */
		inline int beadOn(int,int) const;
		inline int beadOn(const beadLocator&) const;

		/* Get the state of a bead */
		beadState getState (const beadLocator &) const;

		/* The worm-trajectory factor */
		double factor(const beadState, const beadLocator&) const;

		/* Safely add/delete beads */
		inline void delBead(int,int);
		inline void addBead(int,int);
		inline void delBead(const beadLocator&);
		inline void addBead(const beadLocator&);

		/* Reset all worm parameters to null */
		void reset();							

		/* Update all worm parameters */
		void update(Path &, const beadLocator &, const beadLocator &);

		/* Determine if a given beadLocator is on a worm */
		bool foundBead(const Path &, const beadLocator &);

		/** Return true if the worm is too costly */
		inline bool tooCostly() { 
			if (gap == 0) 
				return true;
			else
				return ((dot(sep,sep) * constants()->fourLambdaTauInv() / (1.0*gap)) > maxWormCost);
		}
		/** Return true if the worm is too costly */
		inline bool tooCostly(const dVec& _sep, int _gap) { 
			if (_gap == 0) 
				return true;
			else {
				return ((dot(_sep,_sep) * constants()->fourLambdaTauInv() / (1.0*_gap)) > maxWormCost);
			}
		}

		/** Return the bead list.*/
		const Array <unsigned int, 2> & getBeads() const { return beads; }

		/* Test whether a given bead is on a worm */
		/** Return the number of active beads. */
		int getNumBeadsOn() const {return numBeadsOn;}
		/** Increment the number of active beads. */
		void incNumBeadsOn() {++numBeadsOn;}
		/** Decrement the number of active beads. */
		void decNumBeadsOn() {--numBeadsOn;}
		/** Reset the number of active beads. */
		void resetNumBeadsOn() {numBeadsOn = sum(beads);}

		friend class Path;						// Path needs access to beads
		friend class PathIntegralMonteCarlo;	// Friends for I/O

	private:
		Array <unsigned int,2> beads;			// Is a bead present?
		int numBeadsOn;							// How many beads are present
};

// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// INLINE FUNCTION DEFINITIONS
// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
/** Get the state of the supplied bead? */
inline beadState Worm::getState(const beadLocator &beadIndex) const {
	if (all(beadIndex==head) || all(beadIndex==tail))
		return HEADTAIL;
	else if (all(beadIndex==special1) || all(beadIndex==special2))
		return SPECIAL;
	else
		return NONE;
}

// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

/** Safely delete a bead (int indexed) */
inline void Worm::delBead(int slice, int ptcl) {
	PIMC_ASSERT((slice >= 0) && (slice < constants()->numTimeSlices()));
	beads(slice,ptcl) = 0;
}

/** Safely delete a bead (beadLocator indexed) */
inline void Worm::delBead(const beadLocator &beadIndex) {
	PIMC_ASSERT((beadIndex[0] >= 0) && (beadIndex[0] < constants()->numTimeSlices()));
	beads(beadIndex) = 0;
}

// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
/** Safely add a bead (int indexed) */
inline void Worm::addBead(int slice, int ptcl) {
	PIMC_ASSERT((slice >= 0) && (slice < constants()->numTimeSlices()));
	beads(slice,ptcl) = 1;
}

/** Safely add a bead (beadLocator indexed) */
inline void Worm::addBead(const beadLocator &beadIndex) {
	PIMC_ASSERT((beadIndex[0] >= 0) && (beadIndex[0] < constants()->numTimeSlices()));
	beads(beadIndex) = 1;
}

// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
/** Safely get a bead (int indexed) */
inline int Worm::beadOn(int slice, int ptcl) const {
	PIMC_ASSERT((slice >=0) && (slice < constants()->numTimeSlices()));
	return beads(slice,ptcl);
}

/** Safely get a bead (beadLocator indexed) */
inline int Worm::beadOn(const beadLocator &beadIndex) const {
	return beads(beadIndex);
}

#endif
