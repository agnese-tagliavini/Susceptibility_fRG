
/******************************************************************************************//** @file
 *  		
 * 	file: 		state.h
 * 	contents:  	Definition of ODE state_t class including its norm
 * 
 ****************************************************************************************************/


#pragma once 

#include <arithmetic_tuple.h>
#include <grid.h>
#include <boost/numeric/odeint.hpp>
#include <def.h>
#include <projection.h>  
#include <symmetry_group.h>
#include <symmetries.h>

class state_t: public ReaK::arithmetic_tuple< gf_1p_t, gf_phi_t, gf_phi_t, gf_phi_t, gf_P_t, gf_P_t, gf_P_t, gf_chi_t, gf_chi_t, gf_chi_t, gf_suscept_t, gf_suscept_t, gf_suscept_t,gf_tri_t, gf_tri_t, gf_tri_t, gf_suscept_t, gf_suscept_t, gf_suscept_t > ///< Type of the state vector of the ODE solver
{
   public:
      using base_t = arithmetic_tuple< gf_1p_t, gf_phi_t, gf_phi_t, gf_phi_t, gf_P_t, gf_P_t, gf_P_t, gf_chi_t, gf_chi_t, gf_chi_t, gf_suscept_t, gf_suscept_t, gf_suscept_t, gf_tri_t, gf_tri_t, gf_tri_t, gf_suscept_t, gf_suscept_t, gf_suscept_t >; 
      using base_base_t = base_t::base_t; 
      using Sig_t = gf_1p_t; 
      using phi_t = gf_phi_t; 
      using P_t = gf_P_t; 
      using chi_t = gf_chi_t; 
      using suscept_t = gf_suscept_t; 
      using tri_t = gf_tri_t; 

      inline Sig_t& gf_Sig()			        { return std::get<0>( *this ); }
      inline const Sig_t& gf_Sig() const	        { return std::get<0>( *this ); }

      inline phi_t& gf_phi_pp()			        { return std::get<1>( *this ); }
      inline const phi_t& gf_phi_pp() const	        { return std::get<1>( *this ); }

      inline phi_t& gf_phi_ph()			        { return std::get<2>( *this ); }
      inline const phi_t& gf_phi_ph() const	        { return std::get<2>( *this ); }

      inline phi_t& gf_phi_xph()		        { return std::get<3>( *this ); }
      inline const phi_t& gf_phi_xph() const	        { return std::get<3>( *this ); }

      inline P_t& gf_P_pp()			        { return std::get<4>( *this ); }
      inline const P_t& gf_P_pp() const		        { return std::get<4>( *this ); }

      inline P_t& gf_P_ph()			        { return std::get<5>( *this ); }
      inline const P_t& gf_P_ph() const		        { return std::get<5>( *this ); }

      inline P_t& gf_P_xph()			        { return std::get<6>( *this ); }
      inline const P_t& gf_P_xph() const	        { return std::get<6>( *this ); }

      inline chi_t& gf_chi_pp()			        { return std::get<7>( *this ); }
      inline const chi_t& gf_chi_pp() const	        { return std::get<7>( *this ); }

      inline chi_t& gf_chi_ph()			        { return std::get<8>( *this ); }
      inline const chi_t& gf_chi_ph() const	        { return std::get<8>( *this ); }

      inline chi_t& gf_chi_xph()		        { return std::get<9>( *this ); }
      inline const chi_t& gf_chi_xph() const	        { return std::get<9>( *this ); }

      inline suscept_t& gf_suscept_sc()			{ return std::get<10>( *this ); }
      inline const suscept_t& gf_suscept_sc() const	{ return std::get<10>( *this ); }
      
      //inline suscept_t& gf_suscept_s()			{ return std::get<11>( *this ); }
      //inline const suscept_t& gf_suscept_s() const	{ return std::get<11>( *this ); }
      
      inline suscept_t& gf_suscept_d()			{ return std::get<11>( *this ); }
      inline const suscept_t& gf_suscept_d() const	{ return std::get<11>( *this ); }
      
      inline suscept_t& gf_suscept_m()			{ return std::get<12>( *this ); }
      inline const suscept_t& gf_suscept_m() const	{ return std::get<12>( *this ); }
      
      inline tri_t& gf_tri_sc()				{ return std::get<13>( *this ); }
      inline const tri_t& gf_tri_sc() const		{ return std::get<13>( *this ); }
      
      inline tri_t& gf_tri_d()				{ return std::get<14>( *this ); }
      inline const tri_t& gf_tri_d() const		{ return std::get<14>( *this ); }
      
      inline tri_t& gf_tri_m()				{ return std::get<15>( *this ); }
      inline const tri_t& gf_tri_m() const		{ return std::get<15>( *this ); }
      
      inline suscept_t& gf_asytri_sc()			{ return std::get<16>( *this ); }
      inline const suscept_t& gf_asytri_sc() const	{ return std::get<16>( *this ); }
      
      inline suscept_t& gf_asytri_d()			{ return std::get<17>( *this ); }
      inline const suscept_t& gf_asytri_d() const	{ return std::get<17>( *this ); }
      
      inline suscept_t& gf_asytri_m()			{ return std::get<18>( *this ); }
      inline const suscept_t& gf_asytri_m() const	{ return std::get<18>( *this ); }
      
      
      static gf_vert_bare_ff_t proj_vert_bare;

      static gf_proj_matrix_t proj_matrix_ph_to_pp;
      static gf_proj_matrix_t proj_matrix_pp_to_ph;
      static gf_proj_matrix_t proj_matrix_ph_to_xph;
      
      static symmetry_grp_t<dcomplex,6> symmetry_grp_vert_ff; 
      static symmetry_grp_t<dcomplex,6> symmetry_grp_proj_matrix_ph_to_pp; 
      static symmetry_grp_t<dcomplex,6> symmetry_grp_proj_matrix_pp_to_ph; 
      static symmetry_grp_t<dcomplex,6> symmetry_grp_proj_matrix_ph_to_xph; 
      
      
      F_factors ffactor_mom;

      state_t():
	 base_t(),
	 ffactor_mom()
   {}; 

      state_t( const state_t& state_vec ):
	 base_t( state_vec ), 
	 ffactor_mom()
   {}; 

      state_t( state_t&& state_vec ):
	 base_t( std::move(state_vec) ), 
	 ffactor_mom()
   {}; 

      state_t& operator=( const state_t& state_vec )
      {
	 base_t::operator=( state_vec );
      }

      state_t& operator=( state_t&& state_vec )
      {
	 base_t::operator=( std::move( state_vec ) );
      }

      /********************* Interfacing gf containers  ********************/

      dcomplex Sig( int w, int k, int s_in, int s_out ) const; 		///< Return self-energy value for a specific set of index
      dcomplex Sig_out( int w, int k, int s_in, int s_out ) const; 		///< Return self-energy value for a specific set of index

      MatQN SigMat( int w, int k ) const; 				///< Return self-energy quantum number matrix for specific momentum and frequency given the current state vector
      
      dcomplex vertx( int w1_in, int w2_in, int w1_out, int k1_in, int k2_in, int k1_out, int s1_in, int s2_in, int s1_out, int s2_out ) const;	///< Return vertex-tensor element for up, down, up, down
      dcomplex vertx_pp( int W, int w_in, int w_out, int K, int n_in, int n_out, int s1_in, int s2_in, int s1_out, int s2_out ) const;	///< Return vertex-tensor element
      dcomplex vertx_ph( int W, int w_in, int w_out, int K, int n_in, int n_out, int s1_in, int s2_in, int s1_out, int s2_out ) const;	///< Return vertex-tensor element
      dcomplex vertx_xph( int W, int w_in, int w_out, int K, int n_in, int n_out, int s1_in, int s2_in, int s1_out, int s2_out ) const;	///< Return vertex-tensor element
      
      dcomplex vertx_t( int W, int w_in, int w_out, int K, int n_in, int n_out, int s1_in, int s2_in, int s1_out, int s2_out ) const;	///< Return vertex-tensor element
      dcomplex vertx_s( int W, int w_in, int w_out, int K, int n_in, int n_out, int s1_in, int s2_in, int s1_out, int s2_out ) const;	///< Return vertex-tensor element
      dcomplex vertx_d( int W, int w_in, int w_out, int K, int n_in, int n_out, int s1_in, int s2_in, int s1_out, int s2_out ) const;	///< Return vertex-tensor element
      dcomplex vertx_m( int W, int w_in, int w_out, int K, int n_in, int n_out, int s1_in, int s2_in, int s1_out, int s2_out ) const;	///< Return vertex-tensor element
      dcomplex mom_phi_pp( int w1_in, int w2_in, int w1_out, int k1_in, int k2_in, int k1_out,int s1_in, int s2_in, int s1_out, int s2_out ) const;
      dcomplex mom_phi_ph( int w1_in, int w2_in, int w1_out, int k1_in, int k2_in, int k1_out,int s1_in, int s2_in, int s1_out, int s2_out ) const;
      dcomplex mom_phi_xph( int w1_in, int w2_in, int w1_out, int k1_in, int k2_in, int k1_out,int s1_in, int s2_in, int s1_out, int s2_out ) const;
      
      dcomplex phi_pp( int W, int w_in, int w_out, int K, int n_in, int n_out, int s1_in, int s2_in, int s1_out, int s2_out ) const;	///< Return Phi function for PP-channel
      dcomplex phi_ph( int W, int w_in, int w_out, int K, int n_in, int n_out, int s1_in, int s2_in, int s1_out, int s2_out ) const;	///< Return Phi function for PH-channel
      dcomplex phi_xph( int W, int w_in, int w_out, int K, int n_in, int n_out, int s1_in, int s2_in, int s1_out, int s2_out ) const;	///< Return Phi function for PH-channel

      dcomplex phi_pp_outside( int W, int w_in, int w_out, int K, int n_in, int n_out, int s1_in, int s2_in, int s1_out, int s2_out ) const;	///< Asymptotic part of Phi function for PP-channel
      dcomplex phi_ph_outside( int W, int w_in, int w_out, int K, int n_in, int n_out, int s1_in, int s2_in, int s1_out, int s2_out ) const;	///< Asymptotic part of Phi function for PH-channel
      dcomplex phi_xph_outside( int W, int w_in, int w_out, int K, int n_in, int n_out, int s1_in, int s2_in, int s1_out, int s2_out ) const;	///< Asymptotic part of Phi function for PH-channel

      dcomplex chi_pp( int W, int K, int s1_in, int s2_in, int s1_out, int s2_out ) const;	///< Return chi function for PP-channel
      dcomplex chi_ph( int W, int K, int s1_in, int s2_in, int s1_out, int s2_out ) const;	///< Return chi function for PH-channel
      dcomplex chi_xph( int W, int K, int s1_in, int s2_in, int s1_out, int s2_out ) const;	///< Return chi function for PH-channel

      dcomplex P_pp( int W, int w, int K, int n, int s1_in, int s2_in, int s1_out, int s2_out ) const;	///< Return P function for PP-channel
      dcomplex P_ph( int W, int w, int K, int n, int s1_in, int s2_in, int s1_out, int s2_out ) const;	///< Return P function for PH-channel
      dcomplex P_xph( int W, int w, int K, int n, int s1_in, int s2_in, int s1_out, int s2_out ) const;	///< Return P function for PH-channel

      dcomplex R_pp( int W, int w_in, int w_out, int K, int n_in, int n_out, int s1_in, int s2_in, int s1_out, int s2_out ) const;	///< Return R function for PP-channel
      dcomplex R_ph( int W, int w_in, int w_out, int K, int n_in, int n_out, int s1_in, int s2_in, int s1_out, int s2_out ) const;	///< Return R function for PH-channel
      dcomplex R_xph( int W, int w_in, int w_out, int K, int n_in, int n_out, int s1_in, int s2_in, int s1_out, int s2_out ) const;	///< Return R function for PH-channel
      
      dcomplex suscept_sc( int W, int K, int n_in, int n_out, int s1_in, int s2_in, int s1_out, int s2_out ) const;	///< Return chi function for PH-channel
      dcomplex suscept_d( int W, int K, int n_in, int n_out, int s1_in, int s2_in, int s1_out, int s2_out ) const;	///< Return chi function for PH-channel
      dcomplex suscept_m( int W, int K, int n_in, int n_out, int s1_in, int s2_in, int s1_out, int s2_out ) const;	///< Return chi function for PH-channel
      
      dcomplex tri_sc( int W, int w, int K, int n_in, int n_out, int s1_in, int s2_in, int s1_out, int s2_out ) const;	///< Return chi function for PH-channel
      dcomplex tri_d( int W, int w, int K, int n_in, int n_out, int s1_in, int s2_in, int s1_out, int s2_out ) const;	///< Return chi function for PH-channel
      dcomplex tri_m( int W, int w, int K, int n_in, int n_out, int s1_in, int s2_in, int s1_out, int s2_out ) const;	///< Return chi function for PH-channel
   
      dcomplex asytri_sc( int W, int K, int n_in, int n_out, int s1_in, int s2_in, int s1_out, int s2_out ) const;	///< Return chi function for PH-channel
      dcomplex asytri_d( int W, int K, int n_in, int n_out, int s1_in, int s2_in, int s1_out, int s2_out ) const;	///< Return chi function for PH-channel
      dcomplex asytri_m( int W, int K, int n_in, int n_out, int s1_in, int s2_in, int s1_out, int s2_out ) const;	///< Return chi function for PH-channel
      
      dcomplex proj_pp_phi_ph( int W, int w_in, int w_out, int K, int n_in, int n_out, int s1_in, int s2_in, int s1_out, int s2_out ) const;	///< Return phi_ph projected in pp
      dcomplex proj_pp_phi_xph( int W, int w_in, int w_out, int K, int n_in, int n_out, int s1_in, int s2_in, int s1_out, int s2_out ) const;	///< Return phi_xph projected in pp
      dcomplex proj_ph_phi_pp( int W, int w_in, int w_out, int K, int n_in, int n_out, int s1_in, int s2_in, int s1_out, int s2_out ) const;	///< Return phi_pp projected in ph
      dcomplex proj_ph_phi_xph( int W, int w_in, int w_out, int K, int n_in, int n_out, int s1_in, int s2_in, int s1_out, int s2_out ) const;	///< Return phi_xph projected in ph
      dcomplex proj_xph_phi_pp( int W, int w_in, int w_out, int K, int n_in, int n_out, int s1_in, int s2_in, int s1_out, int s2_out ) const;	///< Return phi_pp projected in xph
      dcomplex proj_xph_phi_ph( int W, int w_in, int w_out, int K, int n_in, int n_out, int s1_in, int s2_in, int s1_out, int s2_out ) const;	///< Return phi_ph projected in xph
      

      dcomplex Gval( const idx_1p_t& Gvec ) const; 

      MatQN GMat( const idx_1p_mat_t& Gvec) const; 

      class precalc
      {
      public:
	 precalc()
	 {
	    symmetry_grp_vert_ff.init(proj_vert_bare,[](const idx_vert_bare_ff_t& idx){return generate_proj_vert_bare(idx); } );
	    symmetry_grp_proj_matrix_ph_to_pp.init(proj_matrix_ph_to_pp,[](const idx_proj_matrix_t& idx){return generate_proj_matrix_ph_to_pp(idx); } );
	    symmetry_grp_proj_matrix_pp_to_ph.init(proj_matrix_pp_to_ph,[](const idx_proj_matrix_t& idx){return generate_proj_matrix_pp_to_ph(idx); } );
	    symmetry_grp_proj_matrix_ph_to_xph.init(proj_matrix_ph_to_xph,[](const idx_proj_matrix_t& idx){return generate_proj_matrix_ph_to_xph(idx); } );
	 }
      };

      static precalc precalculation;
      
      inline dcomplex Sig_out( const idx_1p_t& idx ) const
      {
         return Sig_out( idx(0), idx(1), idx(2), idx(3) ); 
      }
      
      inline dcomplex vertx_pp( const idx_phi_t& idx ) const
      {
         return vertx_pp( idx(0), idx(1), idx(2), idx(3), idx(4), idx(5), idx(6), idx(7), idx(8), idx(9) ); 
      }

      inline dcomplex vertx_ph( const idx_phi_t& idx ) const
      {
         return vertx_ph( idx(0), idx(1), idx(2), idx(3), idx(4), idx(5), idx(6), idx(7), idx(8), idx(9) ); 
      }

      inline dcomplex vertx_xph( const idx_phi_t& idx ) const
      {
         return vertx_xph( idx(0), idx(1), idx(2), idx(3), idx(4), idx(5), idx(6), idx(7), idx(8), idx(9) ); 
      }
      
      inline dcomplex phi_pp( const idx_phi_t& idx ) const
      {
         return phi_pp( idx(0), idx(1), idx(2), idx(3), idx(4), idx(5), idx(6), idx(7), idx(8), idx(9) ); 
      }

      inline dcomplex phi_ph( const idx_phi_t& idx ) const
      {
         return phi_ph( idx(0), idx(1), idx(2), idx(3), idx(4), idx(5), idx(6), idx(7), idx(8), idx(9) ); 
      }

      inline dcomplex phi_xph( const idx_phi_t& idx ) const
      {
         return phi_xph( idx(0), idx(1), idx(2), idx(3), idx(4), idx(5), idx(6), idx(7), idx(8), idx(9) ); 
      }
      
      inline dcomplex R_pp( const idx_phi_t& idx ) const
      {
         return R_pp( idx(0), idx(1), idx(2), idx(3), idx(4), idx(5), idx(6), idx(7), idx(8), idx(9) ); 
      }

      inline dcomplex R_ph( const idx_phi_t& idx ) const
      {
         return R_ph( idx(0), idx(1), idx(2), idx(3), idx(4), idx(5), idx(6), idx(7), idx(8), idx(9) ); 
      }

      inline dcomplex R_xph( const idx_phi_t& idx ) const
      {
         return R_xph( idx(0), idx(1), idx(2), idx(3), idx(4), idx(5), idx(6), idx(7), idx(8), idx(9) ); 
      }

      
      inline dcomplex suscept_sc( const idx_suscept_t& idx ) const
      {
	 return suscept_sc( idx(0), idx(1), idx(2), idx(3), idx(4), idx(5), idx(6), idx(7) ); 
      }
      
      inline dcomplex suscept_d( const idx_suscept_t& idx ) const
      {
	 return suscept_d( idx(0), idx(1), idx(2), idx(3), idx(4), idx(5), idx(6), idx(7) ); 
      }
      
      inline dcomplex suscept_m( const idx_suscept_t& idx ) const
      {
	 return suscept_m( idx(0), idx(1), idx(2), idx(3), idx(4), idx(5), idx(6), idx(7) ); 
      }

}; 


namespace boost {
   namespace numeric {
      namespace odeint {
	 template<>
	    struct vector_space_norm_inf< state_t >
	    {
	       typedef double result_type;
	       double operator()( const state_t &state_vec ) const
	       {
		  using namespace std; 
		  return norm( state_vec ); 
	       }
	    };
      }
   }
}
