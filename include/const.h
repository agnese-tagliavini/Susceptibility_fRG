
/******************************************************************************************//** @file
 *  		
 * 	file: 		const.h
 * 	contents:  	Holds relevant constants for the grids and the numerical calculations
 * 
 ****************************************************************************************************/


#pragma once

#include <complex>
#include <string>

const std::complex<double> I( 0.0, 1.0 );			///< Imaginary unit
const double PI = 3.14159265358979323846; 			///< PI
const double LN_10 = 2.30258509299;				///< Natural log of 10
const double CHOP_ERR = 1E-16; 					///< chopping sets cosine, sine and weight to zero under this value


const int COUNT = 2; 

// ----- SE dimensions

const int POS_FFREQ_COUNT_SIG = 4 * COUNT;			///< Amount of positive frequencies in self-energy grid
const int FFREQ_COUNT_SIG = 2 * POS_FFREQ_COUNT_SIG;		///< Amount of frequencies in self-energy grid

// ----- phi dimensions

const int POS_FFREQ_COUNT_PHI = COUNT;				///< Amount of positive fermionic frequencies in phi grid 
const int FFREQ_COUNT_PHI = 2 * POS_FFREQ_COUNT_PHI;		///< Amount of fermionic frequencies in phi grid

const int POS_BFREQ_COUNT_PHI = 2 * COUNT; 			///< Amount of positive bosonic frequencies in phi grid 
const int BFREQ_COUNT_PHI = 2 * POS_BFREQ_COUNT_PHI + 1;	///< Amount of bosonic frequencies in phi grid

// ----- P dimensions

const int POS_FFREQ_COUNT_P = 4 * COUNT;			///< Amount of positive fermionic frequencies in P grid
const int FFREQ_COUNT_P = 2 * POS_FFREQ_COUNT_P;		///< Amount of fermionic frequencies in P grid
                                                                                                                          
const int POS_BFREQ_COUNT_P = 3 * POS_FFREQ_COUNT_P / 2;		///< Amount of positive bosonic frequencies in P grid 
const int BFREQ_COUNT_P = 2 * POS_BFREQ_COUNT_P + 1;		///< Amount of bosonic frequencies in P grid

// ----- internal integration range and green function grid

const int POS_INT_RANGE = 4 * FFREQ_COUNT_SIG;			///< Positive range for internal integrations
const int TAIL_LENGTH = POS_INT_RANGE / 5; 			///< Length of tail used for fitting matsubara sum
const int FIT_ORDER = 4; 					///< Fit tail function has exponents one lower than this constant

// ----- chi dimensions

const int POS_BFREQ_COUNT_CHI = 2 * POS_INT_RANGE; // !!! DECREASED FOR ffcount=7      2 * POS_INT_RANGE; 		///< Amount of positive bosonic frequencies in chi grid
const int BFREQ_COUNT_CHI = 2 * POS_BFREQ_COUNT_CHI + 1;	///< Amount of bosonic frequencies in chi grid

const int POS_BFREQ_COUNT_SUSCEPT = 8 * COUNT; 
const int BFREQ_COUNT_SUSCEPT = 2 * POS_BFREQ_COUNT_SUSCEPT + 1; 


// ----- Green function and single scale propagator dimension

const int POS_1P_RANGE = POS_BFREQ_COUNT_SUSCEPT + POS_INT_RANGE; 	///< Positive range for Green functions and Single scale propagator

// ----- Output ranges

const int POS_PLOT_RANGE_PHI = 1.2 * POS_FFREQ_COUNT_PHI; 	///< Amount of positive frequencies in phi output grid
const int POS_PLOT_RANGE_VERT = POS_PLOT_RANGE_PHI; 		///< Amount of positive frequencies in vertex output grid


#ifdef NO_MOMENTA
const int MAX_KPOS = 1;
const int PATCH_COUNT = 1;					///< Amount of k-patches

const int NN_COUNT = 0;
const int FFACTOR_COUNT = 1; 

#else
#ifdef UNIFORM_GRID
const int MAX_KPOS = 4;						///< Points on positive kx,ky,kz directions 
const int PATCH_COUNT = 4*MAX_KPOS*MAX_KPOS;		///< Amount of k-patches
#elif defined SIMPLE_PATCH
const int PATCH_COUNT = 8;
#endif
const int NN_COUNT = 1;						///< Next-neighbour considered
const int FFACTOR_COUNT = 5; ///< Amount of form factors->changed accordingly with NN_COUNT in grid.h 					
const int REAL_GRID_FF_SHELL= FFACTOR_COUNT;

const int FFREAL_DIM = 3;                                       ///< Adjust accordingly with the NN_COUNT (e.g NN_COUNT=0 -> FFREAL_DIM = 1, NN_COUNT=1,2,3 -> FFREAL_DIM = 3) 


const int REAL_GRID = 13; 
const int MAX_PROJ_R_GRID = 4;
const int FFT_DIM = 16; 					///< Dimension of the FFT grid in every directions

#endif

const int QN_COUNT = 1;						///< Amount of possible tuples of the discrete quantum numbers


