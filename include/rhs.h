
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
      static dcomplex eval_rhs_Sig_old( const idx_1p_t& se_idx, const state_t& state_vec, const gf_1p_mat_t& Gvec ); 					///< Calculate rhs of self-energy flow for specific index
      
      static dcomplex eval_rhs_Sig( const idx_1p_t& se_idx, const state_t& state_vec, const gf_1p_mat_t& Gvec, gf_bubble_mat_t& bubble_ph , gf_bubble_mat_t& bubble_pp ,F_factors ffactor_mom ); 					///< Calculate rhs of self-energy flow for specific index
      
      static dcomplex eval_rhs_Sig_ph( const idx_1p_t& se_idx, const state_t& state_vec, const gf_1p_mat_t& Gvec, gf_bubble_mat_t& bubble_ph , F_factors ffactor_mom ); 					///< Calculate rhs of self-energy flow for specific index
      
      static dcomplex eval_rhs_Sig_pp( const idx_1p_t& se_idx, const state_t& state_vec, const gf_1p_mat_t& Gvec, gf_bubble_mat_t& bubble_pp , F_factors ffactor_mom ); 					///< Calculate rhs of self-energy flow for specific index
      
      static dcomplex eval_rhs_Sig_xph( const idx_1p_t& se_idx, const state_t& state_vec, const gf_1p_mat_t& Gvec, gf_bubble_mat_t& bubble_ph , F_factors ffactor_mom ); 					///< Calculate rhs of self-energy flow for specific index
      static dcomplex eval_diag_susc_sc( const idx_susc_t& idx, const state_t& state_vec, const gf_bubble_mat_t& bubble_pp );		///< Calculate particle-particle diagram for the flow
      static dcomplex eval_diag_susc_d( const idx_susc_t& idx, const state_t& state_vec, const gf_bubble_mat_t& bubble_ph );		///< Calculate particle-hole diagram for the flow
      static dcomplex eval_diag_susc_m( const idx_susc_t& idx, const state_t& state_vec, const gf_bubble_mat_t& bubble_ph );		///< Calculate particle-hole diagram for the flow
      
      static dcomplex eval_diag_tri_sc( const idx_tri_t& idx, const state_t& state_vec, const gf_bubble_mat_t& bubble_pp );			///< Calculate particle-hole diagram for the flow
      static dcomplex eval_diag_tri_d( const idx_tri_t& idx, const state_t& state_vec, const gf_bubble_mat_t& bubble_ph );			///< Calculate particle-hole diagram for the flow
      static dcomplex eval_diag_tri_m( const idx_tri_t& idx, const state_t& state_vec, const gf_bubble_mat_t& bubble_ph );			///< Calculate particle-hole diagram for the flow

      static MatPatch eval_diag_bubble_pp( const idx_bubble_mat_t& idx, const gf_1p_mat_real_t& Gvec_real, fftw_plan p_b );			///< Calculate particle-hole diagram for the flow
      static MatPatch eval_diag_bubble_ph( const idx_bubble_mat_t& idx, const gf_1p_mat_real_t& Gvec_real, fftw_plan p_b );			///< Calculate particle-hole diagram for the flow
      static MatReal eval_Gvec_real( const idx_1p_mat_real_t& idx, const gf_1p_mat_t& Gvec, fftw_plan p_g );					///< Calculate particle-hole diagram for the flow

   private:
      
      static symmetry_grp_t<dcomplex,4> symm_grp_sig; 
      
      static symmetry_grp_t<dcomplex,8> symm_grp_susc_sc; 
      static symmetry_grp_t<dcomplex,8> symm_grp_susc_d; 
      static symmetry_grp_t<dcomplex,8> symm_grp_susc_m; 
      
      static symmetry_grp_t<dcomplex,7> symm_grp_asytri_sc; 
      static symmetry_grp_t<dcomplex,7> symm_grp_asytri_d; 
      static symmetry_grp_t<dcomplex,7> symm_grp_asytri_m; 
      
      static symmetry_grp_t<dcomplex,9> symm_grp_tri_sc; 
      static symmetry_grp_t<dcomplex,9> symm_grp_tri_d; 
      static symmetry_grp_t<dcomplex,9> symm_grp_tri_m; 
      
      static symmetry_grp_t<MatPatch,8> symm_grp_bubble_pp;
      static symmetry_grp_t<MatPatch,8> symm_grp_bubble_ph;
      
      static symmetry_grp_t<MatReal,3> symm_grp_Gvec_real;
      
      static gf_weight_vec_t weight_vec; 
      static gf<double, 2> weight_vec_2d; 
};
