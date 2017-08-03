
/*******************************************************************************************//** @file
 *  		
 * 	file: 		symmetries.h
 * 	contents:  	List of symmetry functions
 * 
 ****************************************************************************************************/

#pragma once

#include <def.h>
#include <grid.h>
#include <symmetry_group.h>
#include <translate.h>

// Symmetries - Two-particle ( Copied from Vertex code )
template<typename notation> operation exch_in( idx_2p_t& idx );			///< Exchange ingoing lines
template<typename notation> operation exch_out( idx_2p_t& idx );		///< Exchange outgoing lines
template<typename notation> operation compl_conj( idx_2p_t& idx );		///< Complex conjugation 
template<typename notation> operation time_rev( idx_2p_t& idx );		///< Time reversal symmetry
template<typename notation> operation particle_hole( idx_2p_t& idx );		///< Particle hole symmetry in nambu notation
template<typename notation> operation spin_symm( idx_2p_t& idx );		///< Spin symmetry in nambu notation

template<typename notation> operation rot_k( idx_2p_t& idx );			///< Rotate all momenta by 90 degrees
template<typename notation> operation mirror_vert( idx_2p_t& idx );		///< Mirror all momenta vertically
template<typename notation> operation mirror_diag( idx_2p_t& idx );		///< Mirror all momenta diagonally

//// Translation of Two-particle symmetries to P ( problematic )
//using symm_func_t = operation (*)( idx_2p_t& idx ); 
//template<symm_func_t symm_func> operation symm_P( idx_P_t& idx )
//{
//   idx_phi_t idx_phi = iP_to_iphi( idx ); 
//   operation op = symm_func( idx_phi ); 
//   idx = iphi_to_iP( idx_phi ); 
//}

// symmetries from symmetries.cpp , lattice/system independent
// Two-particle - Mixed Notation ( analytic derivation )
operation hmirror_phi_pp( idx_phi_t& idx );	///< horizontal mirroring - exchange of both incoming/outcoming legs
operation hmirror_phi_ph( idx_phi_t& idx );	///< horizontal mirroring - exchange of both incoming/outcoming legs
operation hmirror_phi_xph( idx_phi_t& idx );	///< horizontal mirroring - exchange of both incoming/outcoming legs
operation cconj_phi_pp( idx_phi_t& idx );	///< complex conjugate -exchange incoming with outgoing
operation cconj_phi_ph( idx_phi_t& idx );	///< complex conjugate -exchange incoming with outgoing
operation cconj_phi_xph( idx_phi_t& idx );	///< complex conjugate -exchange incoming with outgoing
operation timerev_phi_pp( idx_phi_t& idx);	///< time reversal
operation timerev_phi_ph( idx_phi_t& idx);	///< time reversal
operation timerev_phi_xph( idx_phi_t& idx);	///< time reversal

// Specific symmetries of P functions - Observed
operation hmirror_P_pp( idx_P_t& idx );	///< horizontal mirroring - exchange of both incoming/outcoming legs
operation timerev_P_ph( idx_P_t& idx );	///< horizontal mirroring - exchange of both incoming/outcoming legs
operation cconj_timerev_P_pp( idx_P_t& idx );	///< complex conjugate -exchange incoming with outgoing
operation cconj_P_ph( idx_P_t& idx );	///< complex conjugate -exchange incoming with outgoing
operation hmirror_cconj_P_xph( idx_P_t& idx );	///< complex conjugate -exchange incoming with outgoing
operation hmirror_timerev_P_xph( idx_P_t& idx);	///< time reversal

// Specific symmetries of chi functions - Observed


operation hmirror_chi_pp( idx_chi_t& idx );	///< horizontal mirroring - exchange of both incoming/outcoming legs
operation hmirror_chi_ph( idx_chi_t& idx );	///< horizontal mirroring - exchange of both incoming/outcoming legs
operation hmirror_chi_xph( idx_chi_t& idx );	///< horizontal mirroring - exchange of both incoming/outcoming legs
operation cconj_chi_pp( idx_chi_t& idx );	///< complex conjugate -exchange incoming with outgoing
operation cconj_chi_ph( idx_chi_t& idx );	///< complex conjugate -exchange incoming with outgoing
operation cconj_chi_xph( idx_chi_t& idx );	///< complex conjugate -exchange incoming with outgoing
operation timerev_chi_pp( idx_chi_t& idx);	///< time reversal
operation timerev_chi_ph( idx_chi_t& idx);	///< time reversal
operation timerev_chi_xph( idx_chi_t& idx);	///< time reversal

operation hmirror_suscept_pp( idx_suscept_t& idx );	///< horizontal mirroring - exchange of both incoming/outcoming legs
operation hmirror_suscept_ph( idx_suscept_t& idx );	///< horizontal mirroring - exchange of both incoming/outcoming legs
operation hmirror_suscept_xph( idx_suscept_t& idx );	///< horizontal mirroring - exchange of both incoming/outcoming legs
operation cconj_suscept_pp( idx_suscept_t& idx );	///< complex conjugate -exchange incoming with outgoing
operation cconj_suscept_ph( idx_suscept_t& idx );	///< complex conjugate -exchange incoming with outgoing
operation cconj_suscept_xph( idx_suscept_t& idx );	///< complex conjugate -exchange incoming with outgoing
operation timerev_suscept_pp( idx_suscept_t& idx);	///< time reversal
operation timerev_suscept_ph( idx_suscept_t& idx);	///< time reversal
operation timerev_suscept_xph( idx_suscept_t& idx);	///< time reversal
// Specific symmetries of Self energy
operation cconj_sig( idx_1p_t& idx );	///< Change sign of bosonic frequency + compl conj

// Specific symmetries of the bubble

operation swap_bubble(  idx_bubble_mat_t& idx );	///< Swap of the form factor dependence
operation cconj_bubble(  idx_bubble_mat_t& idx );	///< complex conjugate for the bubble (same for ph and pp)
operation hmirror_bubble_pp(  idx_bubble_mat_t& idx );	///< Horizontal mirroring operation (valid just for the pp channel)
operation hmirror_timerev_bubble_ph(  idx_bubble_mat_t& idx );	///< Horizontal mirroring operation (valid just for the pp channel)
// symmetries from symmetries_lattice.cpp, adjusted for 3D cubic lattice
// --- phi --

operation rot_pi2_x_phi( idx_phi_t& idx );
operation rot_pi2_y_phi( idx_phi_t& idx );
operation rot_pi2_z_phi( idx_phi_t& idx ); 
operation rot_2pi3_111_phi( idx_phi_t& idx ); 
operation rot_2pi3_1m11_phi( idx_phi_t& idx ); 
operation rot_2pi3_m111_phi( idx_phi_t& idx ); 
operation rot_2pi3_m1m11_phi( idx_phi_t& idx );
operation rot_pi_110_phi( idx_phi_t& idx );
operation rot_pi_1m10_phi( idx_phi_t& idx );
operation rot_pi_011_phi( idx_phi_t& idx ); 
operation rot_pi_01m1_phi( idx_phi_t& idx ); 
operation rot_pi_101_phi( idx_phi_t& idx );
operation rot_pi_10m1_phi( idx_phi_t& idx ); 
operation mirror_y_phi( idx_phi_t& idx ); 

// --- P --

operation rot_pi2_x_P( idx_P_t& idx );
operation rot_pi2_y_P( idx_P_t& idx );
operation rot_pi2_z_P( idx_P_t& idx ); 
operation rot_2pi3_111_P( idx_P_t& idx ); 
operation rot_2pi3_1m11_P( idx_P_t& idx ); 
operation rot_2pi3_m111_P( idx_P_t& idx ); 
operation rot_2pi3_m1m11_P( idx_P_t& idx );
operation rot_pi_110_P( idx_P_t& idx );
operation rot_pi_1m10_P( idx_P_t& idx );
operation rot_pi_011_P( idx_P_t& idx ); 
operation rot_pi_01m1_P( idx_P_t& idx ); 
operation rot_pi_101_P( idx_P_t& idx );
operation rot_pi_10m1_P( idx_P_t& idx ); 
operation mirror_y_P( idx_P_t& idx ); 

// --- chi --

operation rot_pi2_x_chi( idx_chi_t& idx );
operation rot_pi2_y_chi( idx_chi_t& idx );
operation rot_pi2_z_chi( idx_chi_t& idx ); 
operation rot_2pi3_111_chi( idx_chi_t& idx ); 
operation rot_2pi3_1m11_chi( idx_chi_t& idx ); 
operation rot_2pi3_m111_chi( idx_chi_t& idx ); 
operation rot_2pi3_m1m11_chi( idx_chi_t& idx );
operation rot_pi_110_chi( idx_chi_t& idx );
operation rot_pi_1m10_chi( idx_chi_t& idx );
operation rot_pi_011_chi( idx_chi_t& idx ); 
operation rot_pi_01m1_chi( idx_chi_t& idx ); 
operation rot_pi_101_chi( idx_chi_t& idx );
operation rot_pi_10m1_chi( idx_chi_t& idx ); 
operation mirror_y_chi( idx_chi_t& idx ); 

// --- suscept --

operation rot_pi2_x_suscept( idx_suscept_t& idx );
operation rot_pi2_y_suscept( idx_suscept_t& idx );
operation rot_pi2_z_suscept( idx_suscept_t& idx ); 
operation rot_2pi3_111_suscept( idx_suscept_t& idx ); 
operation rot_2pi3_1m11_suscept( idx_suscept_t& idx ); 
operation rot_2pi3_m111_suscept( idx_suscept_t& idx ); 
operation rot_2pi3_m1m11_suscept( idx_suscept_t& idx );
operation rot_pi_110_suscept( idx_suscept_t& idx );
operation rot_pi_1m10_suscept( idx_suscept_t& idx );
operation rot_pi_011_suscept( idx_suscept_t& idx ); 
operation rot_pi_01m1_suscept( idx_suscept_t& idx ); 
operation rot_pi_101_suscept( idx_suscept_t& idx );
operation rot_pi_10m1_suscept( idx_suscept_t& idx ); 
operation mirror_y_suscept( idx_suscept_t& idx ); 

// --- Sig --

operation rot_pi2_x_sig( idx_1p_t& idx );
operation rot_pi2_y_sig( idx_1p_t& idx );
operation rot_pi2_z_sig( idx_1p_t& idx ); 
operation rot_2pi3_111_sig( idx_1p_t& idx ); 
operation rot_2pi3_1m11_sig( idx_1p_t& idx ); 
operation rot_2pi3_m111_sig( idx_1p_t& idx ); 
operation rot_2pi3_m1m11_sig( idx_1p_t& idx );
operation rot_pi_110_sig( idx_1p_t& idx );
operation rot_pi_1m10_sig( idx_1p_t& idx );
operation rot_pi_011_sig( idx_1p_t& idx ); 
operation rot_pi_01m1_sig( idx_1p_t& idx ); 
operation rot_pi_101_sig( idx_1p_t& idx );
operation rot_pi_10m1_sig( idx_1p_t& idx ); 
operation mirror_y_sig( idx_1p_t& idx ); 

// --- Vert BARE FF --

operation rot_pi2_x_vert_ff(idx_vert_bare_ff_t& idx );
operation rot_pi2_y_vert_ff(idx_vert_bare_ff_t& idx );
operation rot_pi2_z_vert_ff(idx_vert_bare_ff_t& idx ); 
operation rot_2pi3_111_vert_ff(idx_vert_bare_ff_t& idx ); 
operation rot_2pi3_1m11_vert_ff(idx_vert_bare_ff_t& idx ); 
operation rot_2pi3_m111_vert_ff(idx_vert_bare_ff_t& idx ); 
operation rot_2pi3_m1m11_vert_ff(idx_vert_bare_ff_t& idx );
operation rot_pi_110_vert_ff(idx_vert_bare_ff_t& idx );
operation rot_pi_1m10_vert_ff(idx_vert_bare_ff_t& idx );
operation rot_pi_011_vert_ff(idx_vert_bare_ff_t& idx ); 
operation rot_pi_01m1_vert_ff(idx_vert_bare_ff_t& idx ); 
operation rot_pi_101_vert_ff(idx_vert_bare_ff_t& idx );
operation rot_pi_10m1_vert_ff(idx_vert_bare_ff_t& idx ); 

// --- PROJ MATRICES --

operation rot_pi2_x_projmat(idx_proj_matrix_t& idx );
operation rot_pi2_y_projmat(idx_proj_matrix_t& idx );
operation rot_pi2_z_projmat(idx_proj_matrix_t& idx ); 
operation rot_2pi3_111_projmat(idx_proj_matrix_t& idx ); 
operation rot_2pi3_1m11_projmat(idx_proj_matrix_t& idx ); 
operation rot_2pi3_m111_projmat(idx_proj_matrix_t& idx ); 
operation rot_2pi3_m1m11_projmat(idx_proj_matrix_t& idx );
operation rot_pi_110_projmat(idx_proj_matrix_t& idx );
operation rot_pi_1m10_projmat(idx_proj_matrix_t& idx );
operation rot_pi_011_projmat(idx_proj_matrix_t& idx ); 
operation rot_pi_01m1_projmat(idx_proj_matrix_t& idx ); 
operation rot_pi_101_projmat(idx_proj_matrix_t& idx );
operation rot_pi_10m1_projmat(idx_proj_matrix_t& idx ); 
operation mirror_y_projmat(idx_proj_matrix_t& idx ); 

// One-particle .
//operation compl_conj( idx_1p_t& ind );	///< Complex conjugation 
//operation time_rev( idx_1p_t& ind );		///< Time reversal symmetry
//operation particle_hole( idx_1p_t& ind );	///< Particle hole symmetry
//operation spin_symm( idx_1p_t& ind );		///< Spin symmetry in nambu notation

//operation rot_k( idx_1p_t& ind );		///< Rotate all momenta by 90 degrees
//operation mirror_vert( idx_1p_t& ind );	///< Mirror all momenta vertically

// -- Tools

/* rotation functions (inline calling the matrix)   -> symmetries.h
 * in symmetries.h: #include <grid.h>  already done
 * */
inline int rot_pi2_x( int& idx, int bf[3])
{
   bf[0]= K_Grid::rot_pi2_x_bfx_matrix[idx]; 
   bf[1]= K_Grid::rot_pi2_x_bfy_matrix[idx]; 
   bf[2]= K_Grid::rot_pi2_x_bfz_matrix[idx]; 
   idx = K_Grid::rot_pi2_x_matrix[idx];
}

inline void rot_pi2_y( int& idx, int bf[3] )
{
   bf[0]= K_Grid::rot_pi2_y_bfx_matrix[idx]; 
   bf[1]= K_Grid::rot_pi2_y_bfy_matrix[idx]; 
   bf[2]= K_Grid::rot_pi2_y_bfz_matrix[idx]; 
   idx = K_Grid::rot_pi2_y_matrix[idx];
}

inline void rot_pi2_z( int& idx, int bf[3] )
{
   bf[0]= K_Grid::rot_pi2_z_bfx_matrix[idx]; 
   bf[1]= K_Grid::rot_pi2_z_bfy_matrix[idx]; 
   bf[2]= K_Grid::rot_pi2_z_bfz_matrix[idx]; 
   idx = K_Grid::rot_pi2_z_matrix[idx]; 
}

// Rotation around the diagonals

inline void rot_2pi3_111( int& idx, int bf[3] )
{
   bf[0]= K_Grid::rot_2pi3_111_bfx_matrix[idx]; 
   bf[1]= K_Grid::rot_2pi3_111_bfy_matrix[idx]; 
   bf[2]= K_Grid::rot_2pi3_111_bfz_matrix[idx]; 
   idx = K_Grid::rot_2pi3_111_matrix[idx];
}

inline void rot_2pi3_1m11( int& idx, int bf[3] )
{
   bf[0]= K_Grid::rot_2pi3_1m11_bfx_matrix[idx]; 
   bf[1]= K_Grid::rot_2pi3_1m11_bfy_matrix[idx]; 
   bf[2]= K_Grid::rot_2pi3_1m11_bfz_matrix[idx]; 
   idx = K_Grid::rot_2pi3_1m11_matrix[idx];
}

inline void rot_2pi3_m111( int& idx, int bf[3] )
{
   bf[0]= K_Grid::rot_2pi3_m111_bfx_matrix[idx]; 
   bf[1]= K_Grid::rot_2pi3_m111_bfy_matrix[idx]; 
   bf[2]= K_Grid::rot_2pi3_m111_bfz_matrix[idx]; 
   idx = K_Grid::rot_2pi3_m111_matrix[idx];
}

inline void rot_2pi3_m1m11( int& idx, int bf[3] )
{
   bf[0]= K_Grid::rot_2pi3_m1m11_bfx_matrix[idx]; 
   bf[1]= K_Grid::rot_2pi3_m1m11_bfy_matrix[idx]; 
   bf[2]= K_Grid::rot_2pi3_m1m11_bfz_matrix[idx]; 
   idx = K_Grid::rot_2pi3_m1m11_matrix[idx];
}

// Rotations around the plane diagonals

inline void rot_pi_110( int& idx, int bf[3] )
{
   bf[0]= K_Grid::rot_pi_110_bfx_matrix[idx]; 
   bf[1]= K_Grid::rot_pi_110_bfy_matrix[idx]; 
   bf[2]= K_Grid::rot_pi_110_bfz_matrix[idx]; 
   idx = K_Grid::rot_pi_110_matrix[idx];
}

inline void rot_pi_1m10( int& idx, int bf[3] )
{
   bf[0]= K_Grid::rot_pi_1m10_bfx_matrix[idx]; 
   bf[1]= K_Grid::rot_pi_1m10_bfy_matrix[idx]; 
   bf[2]= K_Grid::rot_pi_1m10_bfz_matrix[idx]; 
   idx = K_Grid::rot_pi_1m10_matrix[idx];
}

inline void rot_pi_011( int& idx, int bf[3] )
{
   bf[0]= K_Grid::rot_pi_011_bfx_matrix[idx]; 
   bf[1]= K_Grid::rot_pi_011_bfy_matrix[idx]; 
   bf[2]= K_Grid::rot_pi_011_bfz_matrix[idx]; 
   idx = K_Grid::rot_pi_011_matrix[idx];
}

inline void rot_pi_01m1( int& idx, int bf[3] )
{
   bf[0]= K_Grid::rot_pi_01m1_bfx_matrix[idx]; 
   bf[1]= K_Grid::rot_pi_01m1_bfy_matrix[idx]; 
   bf[2]= K_Grid::rot_pi_01m1_bfz_matrix[idx]; 
   idx = K_Grid::rot_pi_01m1_matrix[idx];
}

inline void rot_pi_101( int& idx, int bf[3] )
{
   bf[0]= K_Grid::rot_pi_101_bfx_matrix[idx]; 
   bf[1]= K_Grid::rot_pi_101_bfy_matrix[idx]; 
   bf[2]= K_Grid::rot_pi_101_bfz_matrix[idx]; 
   idx = K_Grid::rot_pi_101_matrix[idx];
}

inline void rot_pi_10m1( int& idx, int bf[3] )
{
   bf[0]= K_Grid::rot_pi_10m1_bfx_matrix[idx]; 
   bf[1]= K_Grid::rot_pi_10m1_bfy_matrix[idx]; 
   bf[2]= K_Grid::rot_pi_10m1_bfz_matrix[idx]; 
   idx = K_Grid::rot_pi_10m1_matrix[idx];
}

inline void mirror_y( int& idx, int bf[3] )
{
   bf[0]= K_Grid::mirror_y_bfx_matrix[idx]; 
   bf[1]= K_Grid::mirror_y_bfy_matrix[idx]; 
   bf[2]= K_Grid::mirror_y_bfz_matrix[idx]; 
   idx = K_Grid::mirror_y_matrix[idx];
}

inline void mirror_mom_pipipi( int& idx, int bf[3] )
{ 
   int k_pipipi[3]= {MAX_KPOS-1,MAX_KPOS-1,MAX_KPOS-1};
   unsigned int patch_pipipi = K_Grid::get_patch(k_pipipi);
   bf[0]= K_Grid::dif_k_bfx_matrix[patch_pipipi][idx]; 
   bf[1]= K_Grid::dif_k_bfy_matrix[patch_pipipi][idx]; 
   bf[2]= K_Grid::dif_k_bfz_matrix[patch_pipipi][idx]; 
   idx = K_Grid::dif_k_matrix[patch_pipipi][idx];
}

//#include <symmetries_impl.h> 	// contains implementations of template functions TODO: check whether it is really used!
