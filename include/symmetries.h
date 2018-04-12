
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

// Two-particle "LEFT DIAGRAMS" used for the MULTILOOP IMPLEMENTATION - Mixed Notation ( analytic derivation )
operation hmirror_phi_pp_L( idx_phi_t& idx );	///< horizontal mirroring - exchange of both incoming/outcoming legs
operation timerev_phi_ph_L( idx_phi_t& idx );	///< horizontal mirroring - exchange of both incoming/outcoming legs
operation cconj_timerev_phi_pp_L( idx_phi_t& idx );	///< complex conjugate -exchange incoming with outgoing
operation cconj_phi_ph_L( idx_phi_t& idx );	///< complex conjugate -exchange incoming with outgoing
operation hmirror_cconj_phi_xph_L( idx_phi_t& idx );	///< complex conjugate -exchange incoming with outgoing
operation hmirror_timerev_phi_xph_L( idx_phi_t& idx);	///< time reversal

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

// Specific symmetries of Self energy
operation cconj_sig( idx_1p_t& idx );	///< Change sign of bosonic frequency + compl conj

// Specific symmetries of the bubble

operation swap_bubble(  idx_bubble_mat_t& idx );	///< Swap of the form factor dependence
operation cconj_bubble(  idx_bubble_mat_t& idx );	///< complex conjugate for the bubble (same for ph and pp)
operation hmirror_bubble_pp(  idx_bubble_mat_t& idx );	///< Horizontal mirroring operation (valid just for the pp channel)
operation hmirror_timerev_ph( idx_bubble_mat_t& idx);   ///< combined horizontal mirroring and timerev (valid only for the ph channel)

// Specific symmetries of Trileg functions - Observed
operation hmirror_tri_sc( idx_tri_t& idx );	///< horizontal mirroring - exchange of both incoming/outcoming legs
operation timerev_tri_d( idx_tri_t& idx );	///< horizontal mirroring - exchange of both incoming/outcoming legs
operation cconj_timerev_tri_sc( idx_tri_t& idx );	///< complex conjugate -exchange incoming with outgoing
operation cconj_tri_d( idx_tri_t& idx );	///< complex conjugate -exchange incoming with outgoing
operation hmirror_cconj_tri_m( idx_tri_t& idx );	///< complex conjugate -exchange incoming with outgoing
operation hmirror_timerev_tri_m( idx_tri_t& idx);	///< time reversal

// Specific symmetries of Trileg functions - Observed

operation hmirror_asytri_sc( idx_asytri_t& idx );	///< horizontal mirroring - exchange of both incoming/outcoming legs
operation timerev_asytri_d( idx_asytri_t& idx );	///< horizontal mirroring - exchange of both incoming/outcoming legs
operation cconj_timerev_asytri_sc( idx_asytri_t& idx );	///< complex conjugate -exchange incoming with outgoing
operation cconj_asytri_d( idx_asytri_t& idx );	///< complex conjugate -exchange incoming with outgoing
operation hmirror_cconj_asytri_m( idx_asytri_t& idx );	///< complex conjugate -exchange incoming with outgoing
operation hmirror_timerev_asytri_m( idx_asytri_t& idx);	///< time reversal

// Specific symmetries of Susceptibility functions - Observed

operation hmirror_susc_sc( idx_susc_t& idx );	///< horizontal mirroring - exchange of both incoming/outcoming legs
operation hmirror_susc_d( idx_susc_t& idx );	///< horizontal mirroring - exchange of both incoming/outcoming legs
operation hmirror_susc_m( idx_susc_t& idx );	///< horizontal mirroring - exchange of both incoming/outcoming legs
operation cconj_susc_sc( idx_susc_t& idx );	///< complex conjugate -exchange incoming with outgoing
operation cconj_susc_d( idx_susc_t& idx );	///< complex conjugate -exchange incoming with outgoing
operation cconj_susc_m( idx_susc_t& idx );	///< complex conjugate -exchange incoming with outgoing
operation timerev_susc_sc( idx_susc_t& idx);	///< time reversal
operation timerev_susc_d( idx_susc_t& idx);	///< time reversal
operation timerev_susc_m( idx_susc_t& idx);	///< time reversal


// symmetries from symmetries_lattice.cpp, adjusted for 3D cubic lattice
// --- phi --

operation rot_pi2_z_phi( idx_phi_t& idx ); 
operation mirror_y_phi( idx_phi_t& idx ); 

// --- P --

operation rot_pi2_z_P( idx_P_t& idx ); 
operation mirror_y_P( idx_P_t& idx ); 

// --- chi --

operation rot_pi2_z_chi( idx_chi_t& idx ); 
operation mirror_y_chi( idx_chi_t& idx ); 

// --- Sig --

operation rot_pi2_z_sig( idx_1p_t& idx ); 
operation mirror_y_sig( idx_1p_t& idx ); 

// --- Vert BARE FF --

operation rot_pi2_z_vert_ff(idx_vert_bare_ff_t& idx ); 
operation mirror_y_vert_ff(idx_vert_bare_ff_t& idx ); 

// --- PROJ MATRICES --

operation rot_pi2_z_projmat(idx_proj_matrix_t& idx ); 
operation mirror_y_projmat(idx_proj_matrix_t& idx ); 
operation mirror_diagonal_projmat(idx_proj_matrix_t& idx ); 

// --- Trileg --

operation rot_pi2_z_tri( idx_tri_t& idx ); 
operation mirror_y_tri( idx_tri_t& idx ); 

// --- Asymptotics Trileg --

operation rot_pi2_z_asytri( idx_asytri_t& idx ); 
operation mirror_y_asytri( idx_asytri_t& idx ); 
// --- Susc --

operation rot_pi2_z_susc( idx_susc_t& idx ); 
operation mirror_y_susc( idx_susc_t& idx ); 

// -- Tools

/* rotation functions (inline calling the matrix)   -> symmetries.h
 * in symmetries.h: #include <grid.h>  already done
 * */

// Rotation of pi/2 of the XY plane around the z axis
inline void rot_pi2_z( int& idx, int bf[2] )
{
   bf[0]= K_Grid::rot_pi2_z_bfx_matrix[idx]; 
   bf[1]= K_Grid::rot_pi2_z_bfy_matrix[idx]; 
   idx  = K_Grid::rot_pi2_z_matrix[idx]; 
}

// Reflection of the XY plane with respect to the yz plane

inline void mirror_y( int& idx, int bf[2] )
{
   bf[0]= K_Grid::mirror_y_bfx_matrix[idx]; 
   bf[1]= K_Grid::mirror_y_bfy_matrix[idx]; 
   idx = K_Grid::mirror_y_matrix[idx];
}

// reflection with respect to the diagonal (x<->y)
inline void mirror_diagonal( int& idx)
{
   idx = K_Grid::mirror_diagonal_matrix[idx];
}

#ifdef UNIFORM_GRID
inline void mirror_mom_pipi( int& idx, int bf[2] )
{ 
   int k_pipi[2]= {MAX_KPOS-1,MAX_KPOS-1};
   unsigned int patch_pipi = K_Grid::get_patch(k_pipi);
   bf[0]=K_Grid::dif_k_bfx_matrix[patch_pipi][idx];
   bf[1]=K_Grid::dif_k_bfy_matrix[patch_pipi][idx];
   idx = K_Grid::dif_k_matrix[patch_pipi][idx];
}

#elif defined SIMPLE_PATCH
inline void mirror_mom_pipi( int& idx )
{ 
   idx = K_Grid::mirror_mom_pipi_matrix[idx];
}
#endif
//#include <symmetries_impl.h> 	// contains implementations of template functions TODO: check whether it is really used!
