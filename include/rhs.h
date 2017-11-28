
/*******************************************************************************************//** @file
 *  		
 * 	file: 		rhs.h
 * 	contents:  	Functor specifying the calculation of the rhs of the flow equation
 * 
 ****************************************************************************************************/


#pragma once

#include <state.h>
#include <grid.h>
#include <symmetry_group.h>
#include <symmetries.h>
#include <fftw3.h>

class rhs_t 		///< Functor to specify the rhs calculation for both the self-energy and the vertex
{
   public:
      void operator() ( const state_t& state_vec, state_t &dfdl ); 	///< Overload call operator to calculate full rhs

      /********************* RHS evaluation ********************/
      static dcomplex eval_diag_suscept_trip( const idx_suscept_t& idx, const state_t& state_vec, const gf_bubble_mat_t& bubble_pp );		///< Calculate particle-particle diagram for the flow
      static dcomplex eval_diag_suscept_s( const idx_suscept_t& idx, const state_t& state_vec, const gf_bubble_mat_t& bubble_pp );		///< Calculate particle-hole diagram for the flow
      static dcomplex eval_diag_suscept_d( const idx_suscept_t& idx, const state_t& state_vec, const gf_bubble_mat_t& bubble_ph );		///< Calculate particle-hole diagram for the flow
      static dcomplex eval_diag_suscept_m( const idx_suscept_t& idx, const state_t& state_vec, const gf_bubble_mat_t& bubble_ph );		///< Calculate particle-hole diagram for the flow
      
      static dcomplex eval_diag_tri_sc( const    idx_tri_t& idx, const state_t& state_vec, const gf_bubble_mat_t& bubble_pp );		///< Calculate particle-hole diagram for the flow
      static dcomplex eval_diag_tri_d( const    idx_tri_t& idx, const state_t& state_vec, const gf_bubble_mat_t& bubble_ph );		///< Calculate particle-hole diagram for the flow
      static dcomplex eval_diag_tri_m( const    idx_tri_t& idx, const state_t& state_vec, const gf_bubble_mat_t& bubble_ph );		///< Calculate particle-hole diagram for the flow

      static MatPatch eval_diag_bubble_pp( const idx_bubble_mat_t& idx, const gf_1p_mat_real_t& Gvec_real, fftw_plan p_b );		///< Calculate particle-hole diagram for the flow
      static MatPatch eval_diag_bubble_ph( const idx_bubble_mat_t& idx, const gf_1p_mat_real_t& Gvec_real, fftw_plan p_b );		///< Calculate particle-hole diagram for the flow
      static MatReal eval_Gvec_real( const idx_1p_mat_real_t& idx, const gf_1p_mat_t& Gvec, fftw_plan p_g );		///< Calculate particle-hole diagram for the flow

   private:
      
      static symmetry_grp_t<dcomplex,8> symm_grp_chi_tripl; 
      static symmetry_grp_t<dcomplex,8> symm_grp_chi_singl; 
      static symmetry_grp_t<dcomplex,8> symm_grp_chi_dens; 
      static symmetry_grp_t<dcomplex,8> symm_grp_chi_mag; 
      
      static symmetry_grp_t<dcomplex,6> symm_grp_chi_pp; 
      static symmetry_grp_t<dcomplex,6> symm_grp_chi_ph; 
      static symmetry_grp_t<dcomplex,6> symm_grp_chi_xph; 
      
      static symmetry_grp_t<dcomplex,9> symm_grp_tri_sc; 
      static symmetry_grp_t<dcomplex,9> symm_grp_tri_dens; 
      static symmetry_grp_t<dcomplex,9> symm_grp_tri_mag; 
      
      
      static symmetry_grp_t<MatPatch,8> symm_grp_bubble_pp;
      static symmetry_grp_t<MatPatch,8> symm_grp_bubble_ph;
      
      static symmetry_grp_t<MatReal,3> symm_grp_Gvec_real;
      
      static gf_weight_vec_t weight_vec; 
      static gf<double, 2> weight_vec_2d; 
};
