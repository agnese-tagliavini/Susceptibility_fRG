
/*******************************************************************************************//** @file
 *  		
 * 	file: 		frg.h
 * 	contents:  	Definition of Green function ( i.e. system ) and fRG specific functions 
 * 
 ****************************************************************************************************/


#pragma once

#include <def.h>
#include <params.h>

extern gf_1p_t gf_Sig_read; 

extern gf_phi_t gf_phi_pp_read;
extern gf_phi_t gf_phi_ph_read;
extern gf_phi_t gf_phi_xph_read;

extern gf_P_t gf_P_pp_read;
extern gf_P_t gf_P_ph_read;
extern gf_P_t gf_P_xph_read;

extern gf_chi_t gf_chi_pp_read;
extern gf_chi_t gf_chi_ph_read;
extern gf_chi_t gf_chi_xph_read;

extern const int POS_FERM_PHI_IN, POS_BOS_PHI_IN, POS_SIG_IN; 
extern const int POS_FERM_P_IN, POS_BOS_P_IN, POS_BOS_CHI_IN; 

void read_exact(); 
#ifdef NO_MOMENTA
MatQN Gam( double w );				///< Hybridization function according to the compile flags in makefile
#endif

MatQN G0inv( double w, double kx, double ky);					///< Inverse non-interacting Greens function ( scale-independent! )

MatQN G( double w, double kx, double ky,  const MatQN& selfEn );	///< Scale-dependent Greens function, introduce regulator here!
dcomplex vert_bare( int s1_in, int s2_in, int s1_out, int s2_out );  	///< Initial vertex values in Orbital basis
dcomplex vert_bare( const idx_2p_t& idx );				///< Return initial vertex value for given index object

dcomplex Sig_init( const idx_1p_t& idx );					///< Return initial self-energy value for given index object
dcomplex phi_init_pp( const idx_phi_t& idx ); 
dcomplex phi_init_ph( const idx_phi_t& idx ); 
dcomplex phi_init_xph( const idx_phi_t& idx ); 

dcomplex P_init_pp( const idx_P_t& idx ); 
dcomplex P_init_ph( const idx_P_t& idx ); 
dcomplex P_init_xph( const idx_P_t& idx ); 

dcomplex chi_init_pp( const idx_chi_t& idx );				///< Return initial vertex value for chiasch functions
dcomplex chi_init_ph( const idx_chi_t& idx );				///< Return initial vertex value for chiasch functions
dcomplex chi_init_xph( const idx_chi_t& idx );				///< Return initial vertex value for chiasch functions

dcomplex suscept_init_sc( const idx_suscept_t& idx );				///< Return initial vertex value for chiasch functions
dcomplex suscept_init_d( const idx_suscept_t& idx );				///< Return initial vertex value for chiasch functions
dcomplex suscept_init_m( const idx_suscept_t& idx );				///< Return initial vertex value for chiasch functions

dcomplex tri_init( const idx_tri_t& idx );				///< Return initial vertex value for chiasch functions
dcomplex tri_bare( int n_in, int n_out, int s1_in, int s2_in, int s1_out, int s2_out);				///< Return initial vertex value for chiasch functions

dcomplex asympt_GG_pp( int W );	///< Give estimate for the integral 1/2/PI * G(W/2-w-1,Lam) * G(W/2+w,Lam) for the outside region N \ [ -POS_INT_RANGE, POS_INT_RANGE )
dcomplex asympt_GG_ph( int W );	///< Give estimate for the integral 1/2/PI * G(w-W/2,Lam) * G(w+W/2,Lam) for the outside region N \ [ -POS_INT_RANGE, POS_INT_RANGE )

inline double w_val( int w_idx )					///< Return value of fermionic matsubara frequency for a given index
{
   return PI / BETA * ( 2 * w_idx + 1 ); 
}

inline double W_val( int W_idx )					///< Return value of bosonic matsubara frequency for a given index
{
   return 2.0 * PI / BETA * W_idx; 
}
