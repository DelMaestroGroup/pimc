/**
 * @file move.h
 * @author Adrian Del Maestro
 * @date 10.21.2008
 *
 * @brief Move class definitions.
 */

#include "common.h"

#ifndef MOVE_H 
#define MOVE_H

//#define DEBUG_WORM
//#define DEBUG_MOVE

class Path;
class ActionBase;

// ========================================================================  
// MoveBase Class
// ========================================================================  
/** 
 * The base class that all moves will be derived from.
 * 
 * There will be a bunch of different types of moves, those which
 * move individual particles, the center of mass of an entire worldine
 * loop etc.  They will all share the basic functionality defined here.
 */
class MoveBase {

	public:
		MoveBase (Path &, ActionBase *, MTRand &, string _name="", 
                ensemble _operateOnConfig=ANY);
		virtual ~MoveBase();

		string name;	                ///< The name of the estimator
        ensemble operateOnConfig;       ///< What configurations do we operate on?

		/** Get the acceptance ratio. */
		double getAcceptanceRatio() {
			return (numAttempted == 0 ? 0.0 : 1.0*numAccepted/(1.0*numAttempted));
		}

		/** Get the total acceptance ratio. */
		double getTotAcceptanceRatio() {
			return (totAttempted == 0 ? 0.0 : 1.0*totAccepted/(1.0*totAttempted));
		}

		/** Get the acceptance ratio by level. */
		double getAcceptanceRatioLevel(int n) { 
			return (numAttemptedLevel(n) == 0 ? 0.0 : 
					1.0*numAcceptedLevel(n)/(1.0*numAttemptedLevel(n)));
		}

		/** Get the number of moves attempted */
		int getNumAttempted() { return numAttempted; }
		/** Get the number of moves accepted */
		int getNumAccepted() { return numAccepted; }
		/** Get the number of moves attempted by level */
		int getNumAttemptedLevel(int n) { return numAttemptedLevel(n); }
		/** Get the number of moves accepted by level */
		int getNumAcceptedLevel(int n) { return numAcceptedLevel(n); }

		/** Attempt the move (will be overloaded). */
		virtual bool attemptMove() = 0;	

		/** Reset the total accepted counter */
		void resetTotAccept() { totAccepted = totAttempted = 0; }
		/** Reset the number accepted counter */
		void resetAccept() { numAccepted = numAttempted = 0; }

	protected:
		friend class PathIntegralMonteCarlo;	// Friends for I/O

		Path &path;						///< A reference to the paths
		ActionBase *actionPtr;			///< A base pointer to the action
		MTRand &random;					///< A reference to the RNG

		bool success;					///< Did we sucessfully perform a move?

		uint32 numAccepted;				///< The number of accepted moves
		uint32 numAttempted;			///< The number of attempted moves
		int numToMove;					///< The number of particles moved

		static uint32 totAccepted;		///< The total number of  moves accepted
		static uint32 totAttempted;		///< The total number of  moves attempted

		Array <uint32,1> numAcceptedLevel;	///< The number of moves accepted at each level
		Array <uint32,1> numAttemptedLevel;	///< The number of moves attempted at each level

		Array <dVec,1> originalPos;		///< The original particle positions

		double oldAction;				///< The original potential action
		double newAction;				///< The new potential action

		double sqrt2LambdaTau;			///< sqrt(2 * Lambda * tau)

		dVec neighborPos; 				///< Staging neighbor position
		dVec newRanPos; 				///< Staing random position

		/** Keep the move */
		virtual void keepMove();		
		/** undo the move */
		virtual void undoMove() = 0;	

		/* Returns a new bead position based on the staging algorithm */
		dVec newStagingPosition(const beadLocator &, const beadLocator &, const int, const int);

		/* Return a new bead position which samples the free particle density matrix */
		dVec newFreeParticlePosition(const beadLocator &);

		///@cond DEBUG
		double newK,oldK;				// The old and new action pieces
		double newV,oldV;

		/* Debugging methods */
		void printMoveState(string);
		void checkMove(int,double);
		///@endcond DEBUG
};

// ========================================================================  
// CenterOfMass Move Class 
// ========================================================================  
/** 
 * A derived class which performs a simple displacement of the center of
 * mass of the entire wordline for a particle.
 */
class CenterOfMassMove: public MoveBase {

	public:
		CenterOfMassMove(Path &, ActionBase *, MTRand &, 
                string _name="center of mass", ensemble _operateOnConfig=ANY);
		~CenterOfMassMove();

		bool attemptMove();

	private:
		beadLocator startBead,endBead;	// The start and end beads
		void undoMove();				// revert everything back
};


// ========================================================================  
// Staging Move Class 
// ========================================================================  
/** 
 * A derived class which performs a staging move, which exactly samples
 * the kinetic action.
 */
class StagingMove: public MoveBase {

	public:
		StagingMove(Path &, ActionBase *, MTRand &, string _name="staging", 
                ensemble _operateOnConfig=ANY);
		~StagingMove();

		bool attemptMove();

	private:
		beadLocator startBead,endBead;		// The start and end of the stage
		bool foundTail;						// Have we found the worm tail?

		void undoMove(); 					// Undo the move
};

// ========================================================================  
// Open Move Class 
// ========================================================================  
/** 
 * A derived class which performs an open move, creating a worm with a well
 * defined head and tail.
 */
class OpenMove: public MoveBase {

	public:
		OpenMove(Path &, ActionBase *, MTRand &, string _name="open",
                ensemble _operateOnConfig=DIAGONAL);
		~OpenMove();

		bool attemptMove();
		bool attemptMove1();

	private:
		beadLocator headBead, tailBead;	// The temporary head and tail locatores
		int gapLength;					// The proposed WL length to remove 
		int numLevels;					// The 2^numLevels = num slices moved

		void undoMove();				// Undo a move
		void keepMove();				// keep the move

};

// ========================================================================  
// Close Move Class 
// ========================================================================  
/** 
 * A derived class which performs a close move, creating a diagonal world
 * line configuration.
 */
class CloseMove: public MoveBase {

	public:
		CloseMove(Path &, ActionBase *, MTRand &, string _name="close",
                ensemble _operateOnConfig=OFFDIAGONAL);
		~CloseMove();

		bool attemptMove();
		bool attemptMove1();

	private:
		beadLocator headBead,tailBead;	// The temporary head and tail slices
		int healLength;					// The proposed WL length to close
		int numLevels;						// The 2^numLevels = num slices moved

		Array <int,1> oldBeadOn;		// The old and new bead states

		void undoMove();				// Undo a move
		void keepMove();				// keep the move
};

// ========================================================================  
// Insert Move Class 
// ========================================================================  
/** 
 * A derived class which performs an insert move, creating an off-diagonal 
 * world line configuration with a single worm.
 */
class InsertMove: public MoveBase {

	public:
		InsertMove(Path &, ActionBase *, MTRand &, string _name="insert",
                ensemble _operateOnConfig=DIAGONAL);
		~InsertMove();

		bool attemptMove();
		bool attemptMove1();

	private:
		beadLocator headBead,tailBead;	// The temporary head and tail beads

		int wormLength;					// The proposed WL length to close
		int numLevels;					// The 2^numLevels = num slices moved

		void undoMove();				// Undo a move
		void keepMove();				// keep the move
};

// ========================================================================  
// Remove Move Class 
// ========================================================================  
/** 
 * A derived class which performs a remove move, creating a diagonal 
 * world line configuration by destroying a single worm.
 */
class RemoveMove: public MoveBase {

	public:
		RemoveMove(Path &, ActionBase *, MTRand &, string _name="remove",
                ensemble _operateOnConfig=OFFDIAGONAL);
		~RemoveMove();

		
		bool attemptMove();
		bool attemptMove1();

	private:
		int numLevels;					// The 2^numLevels = num slices moved

		void undoMove();				// Undo a move
		void keepMove();				// keep the move
};

// ========================================================================  
// Advance Head Move Class 
// ========================================================================  
/** 
 * A derived class which performs an advance head move, causing the head of 
 * a worm in a off-diagonal configuration to advance in imaginary time.
 */
class AdvanceHeadMove: public MoveBase {

	public:
		AdvanceHeadMove(Path &, ActionBase *, MTRand &, 
                string _name="advance head", 
                ensemble _operateOnConfig=OFFDIAGONAL);
		~AdvanceHeadMove();
		
		bool attemptMove();
		bool attemptMove1();

	private:
		beadLocator headBead;			// The temporary new head

		int advanceLength;				// The proposed WL length to advance
		int numLevels;					// The 2^numLevels = num slices moved

		void undoMove();				// Undo a move
		void keepMove();				// keep the move

		int healLength;
		beadLocator startBead;			
		Array <dVec,1> newPos;					// The modified particle positions
		Array <unsigned int,1> oldBeadOn;		// The old and new bead states

};

// ========================================================================  
// Advance Tail Move Class 
// ========================================================================  
/** 
 * A derived class which performs an advance tail move, causing the tail of 
 * a worm in a off-diagonal configuration to advance in imaginary time, resulting
 * in a shorter worm.
 */
class AdvanceTailMove: public MoveBase {

	public:
        AdvanceTailMove(Path &, ActionBase *, MTRand &, 
                string _name="advance tail",
                ensemble _operateOnConfig=OFFDIAGONAL);
		~AdvanceTailMove();

		bool attemptMove();
		bool attemptMove1();

	private:
		beadLocator tailBead;			// The temporary new tail

		int advanceLength;				// The proposed WL length to advance
		int numLevels;					// The 2^numLevels = num slices moved

		void undoMove();				// Undo a move
		void keepMove();				// keep the move
};

// ========================================================================  
// Recede Head Move Class 
// ========================================================================  
/** 
 * A derived class which performs a recede move on the head, causing a worm 
 * head to propagate backwards in imaginary time by removing beads and links.
 */
class RecedeHeadMove: public MoveBase {

	public:
		RecedeHeadMove(Path &, ActionBase *, MTRand &,
                string _name="recede head",
                ensemble _operateOnConfig=OFFDIAGONAL);
		~RecedeHeadMove();
		
		bool attemptMove();
		bool attemptMove1();

	private:
		beadLocator headBead;			// The proposed new head position
		int recedeLength;				// The number of slices to recede by
		int numLevels;					// The 2^numLevels = num slices moved

		void undoMove();				// Undo a move
		void keepMove();				// keep the move
};

// ========================================================================  
// Recede Tail Move Class 
// ========================================================================  
/** 
 * A derived class which performs a recede move on the tail, causing a worm 
 * tail to propagate backwards in imaginary time by adding beads and links.
 */
class RecedeTailMove: public MoveBase {

	public:
		RecedeTailMove(Path &, ActionBase *, MTRand &,
                string _name="recede tail",
                ensemble _operateOnConfig=OFFDIAGONAL);
		~RecedeTailMove();
		
		bool attemptMove();
		bool attemptMove1();

	private:
		beadLocator tailBead;			// The proposed new head position
		int recedeLength;				// The number of slices to recede by
		int numLevels;					// The 2^numLevels = num slices moved

		void undoMove();				// Undo a move
		void keepMove();				// keep the move
};

// ========================================================================  
// Swap Move Base Class 
// ========================================================================  
/** 
 * A derived class which forms the base of a swap head and swap tail move
 * class.
 */
class SwapMoveBase: public MoveBase {

	public:
		SwapMoveBase(Path &, ActionBase *, MTRand &, string _name="swap",
                ensemble _operateOnConfig=OFFDIAGONAL);
		~SwapMoveBase();

	protected:
		int swapLength;						///< The length of worldLine to be moved
		int numLevels;						///< The number of bisection levels

		vector <double> cumulant;			///< The cumulant array used in selecting a pivot

		beadLocator pivot;					///< The pivot bead
		beadLocator swap;					///< The swap bead

		double SigmaSwap;					///< Probability normalization factor

		/* Returns the normalization factor for the probability dist. */
		double getNorm(const beadLocator&);
		
		/* Gets the bead where the swap will pivot. */
		beadLocator selectPivotBead();
};

// ========================================================================  
// Swap Head Move Class 
// ========================================================================  
/** 
 * A derived class which performs a swap head move, which mixes up worldlines
 * by reconnecting the worm head and is essential for systems with permutation 
 * symmetry (such as bosons).
 */
class SwapHeadMove: public SwapMoveBase {

	public:
		SwapHeadMove(Path &, ActionBase *, MTRand &,
                string _name="swap head",
                ensemble _operateOnConfig=OFFDIAGONAL);
		~SwapHeadMove();

		bool attemptMove();

	private:

		double SigmaHead;		// The probability normalization factor

		beadLocator nextSwap;	// Used for re-linking

		void undoMove();		// Undo a move
		void keepMove();		// keep the move
};

// ========================================================================  
// Swap Tail Move Class 
// ========================================================================  
/** 
 * A derived class which performs a swap tail move, which mixes up worldlines
 * by reconnecting the worm tail and is essential for systems with permutation 
 * symmetry (such as our bosons).
 */
class SwapTailMove: public SwapMoveBase {

	public:
		SwapTailMove(Path &, ActionBase *, MTRand &,
                string _name="swap tail",
                ensemble _operateOnConfig=OFFDIAGONAL);
		~SwapTailMove();
		
		bool attemptMove();

	private:

		double SigmaTail;				// The probability normalization factor

		beadLocator prevSwap;			// Used for re-linking

		void undoMove();				// Undo a move
		void keepMove();				// keep the move
};

#endif
