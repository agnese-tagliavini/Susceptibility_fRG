

/******************************************************************************************//** @file
 *  		
 * 	file: 		state.cpp
 * 	contents:  	See state.h
 * 
 ****************************************************************************************************/


#include <state.h>
#include <mymath.h>
#include <frg.h>
#include <projection.h>
//#include <cubature.h>

using namespace std; //can be deleted after debugging
/********************* Interfacing gf containers  ********************/


symmetry_grp_t<dcomplex,6> state_t::symmetry_grp_vert_ff(gf_vert_bare_ff_t(),{ });

//------ with symmetries 
symmetry_grp_t<dcomplex,6> state_t::symmetry_grp_proj_matrix_ph_to_pp(gf_proj_matrix_t(),{
      rot_pi2_z_projmat,
      mirror_y_projmat,
      mirror_diagonal_projmat 
      });
symmetry_grp_t<dcomplex,6> state_t::symmetry_grp_proj_matrix_pp_to_ph(gf_proj_matrix_t(),{ 
      rot_pi2_z_projmat,
      mirror_y_projmat,
      mirror_diagonal_projmat 
      });
symmetry_grp_t<dcomplex,6> state_t::symmetry_grp_proj_matrix_ph_to_xph(gf_proj_matrix_t(),{ 
      rot_pi2_z_projmat,
      mirror_y_projmat,
      mirror_diagonal_projmat 
      });
//-------no symmetries
//symmetry_grp_t<dcomplex,6> state_t::symmetry_grp_proj_matrix_ph_to_pp(gf_proj_matrix_t(),{ });
//symmetry_grp_t<dcomplex,6> state_t::symmetry_grp_proj_matrix_pp_to_ph(gf_proj_matrix_t(),{ });
//symmetry_grp_t<dcomplex,6> state_t::symmetry_grp_proj_matrix_ph_to_xph(gf_proj_matrix_t(),{ });

gf_vert_bare_ff_t state_t::proj_vert_bare;

gf_proj_matrix_t state_t::proj_matrix_ph_to_pp;
gf_proj_matrix_t state_t::proj_matrix_pp_to_ph;
gf_proj_matrix_t state_t::proj_matrix_ph_to_xph;

state_t::precalc state_t::precalculation;



dcomplex state_t::Sig( int w, int k, int s_in, int s_out ) const
{
   if ( w < -POS_FFREQ_COUNT_SIG || w > POS_FFREQ_COUNT_SIG - 1 ) 
      return 0.0; 

   return gf_Sig()[w][k][s_in][s_out]; 
}


MatQN state_t::SigMat( int w, int k ) const
{
   if ( w < -POS_FFREQ_COUNT_SIG || w > POS_FFREQ_COUNT_SIG - 1 ) 
      return MatQN::Zero(); 
   return Eigen::Map<const MatQN>( &(gf_Sig()[w][k][0][0]) );  
}


dcomplex state_t::vertx_pp( int W, int w_in, int w_out, int K, int n_in, int n_out, int s1_in, int s2_in, int s1_out, int s2_out ) const
{
      return phi_pp( W, w_in, w_out, K, n_in, n_out, s1_in, s2_in, s1_out, s2_out ) + 
      proj_pp_phi_ph( W, w_in, w_out, K, n_in, n_out, s1_in, s2_in, s1_out, s2_out ) +
      proj_pp_phi_xph( W, w_in, w_out, K, n_in, n_out, s1_in, s2_in, s1_out, s2_out ) +
      proj_vert_bare[n_in][n_out][s1_in][s2_in][s1_out][s2_out];
}

dcomplex state_t::vertx_ph( int W, int w_in, int w_out, int K, int n_in, int n_out, int s1_in, int s2_in, int s1_out, int s2_out ) const
{

   return phi_ph( W, w_in, w_out, K, n_in, n_out, s1_in, s2_in, s1_out, s2_out ) + 
      proj_ph_phi_pp( W, w_in, w_out, K, n_in, n_out, s1_in, s2_in, s1_out, s2_out ) +
      proj_ph_phi_xph( W, w_in, w_out, K, n_in, n_out, s1_in, s2_in, s1_out, s2_out ) +
      proj_vert_bare[n_in][n_out][s1_in][s2_in][s1_out][s2_out];
}

dcomplex state_t::vertx_xph( int W, int w_in, int w_out, int K, int n_in, int n_out, int s1_in, int s2_in, int s1_out, int s2_out ) const
{
   
   return phi_xph( W, w_in, w_out, K, n_in, n_out, s1_in, s2_in, s1_out, s2_out ) + 
      proj_xph_phi_pp( W, w_in, w_out, K, n_in, n_out, s1_in, s2_in, s1_out, s2_out ) +
      proj_xph_phi_ph( W, w_in, w_out, K, n_in, n_out, s1_in, s2_in, s1_out, s2_out ) +
      proj_vert_bare[n_in][n_out][s1_in][s2_in][s1_out][s2_out];

}
/*************************************** VEERTEX IN PHYSICAL BASIS D/M/S/T ****************************************************************/

dcomplex state_t::vertx_t(int W, int w_in, int w_out, int K, int n_in, int n_out, int s1_in, int s2_in, int s1_out, int s2_out ) const
{
   return vertx_pp(W, w_in, w_out, K, n_in, n_out, s1_in, s2_in, s1_out, s2_out)-(double)Parity(n_out) * vertx_pp(W, w_in, -w_out-1, K, n_in, n_out, s1_in, s2_in, s1_out, s2_out);

}

dcomplex state_t::vertx_s(int W, int w_in, int w_out, int K, int n_in, int n_out, int s1_in, int s2_in, int s1_out, int s2_out ) const
{
   return vertx_pp(W, w_in, w_out, K, n_in, n_out, s1_in, s2_in, s1_out, s2_out)+(double)Parity(n_out) * vertx_pp(W, w_in, -w_out-1, K, n_in, n_out, s1_in, s2_in, s1_out, s2_out);

}

dcomplex state_t::vertx_d(int W, int w_in, int w_out, int K, int n_in, int n_out, int s1_in, int s2_in, int s1_out, int s2_out ) const
{
   return 2.0 * vertx_ph(W, w_in, w_out, K, n_in, n_out, s1_in, s2_in, s1_out, s2_out)-vertx_xph(W, w_in, w_out, K, n_in, n_out, s1_in, s2_in, s1_out, s2_out);

}

dcomplex state_t::vertx_m(int W, int w_in, int w_out, int K, int n_in, int n_out, int s1_in, int s2_in, int s1_out, int s2_out ) const
{
   return -vertx_xph(W, w_in, w_out, K, n_in, n_out, s1_in, s2_in, s1_out, s2_out);

}


/************************************* PHIs IN ALL CHANNELS *****************************************************************************/

dcomplex state_t::phi_pp( int W, int w_in, int w_out, int K, int n_in, int n_out, int s1_in, int s2_in, int s1_out, int s2_out ) const
{
   if ( W < -POS_BFREQ_COUNT_PHI || W > POS_BFREQ_COUNT_PHI ||
	 w_in < -POS_FFREQ_COUNT_PHI || w_in > POS_FFREQ_COUNT_PHI - 1 ||
	 w_out < -POS_FFREQ_COUNT_PHI || w_out > POS_FFREQ_COUNT_PHI - 1 )
      return phi_pp_outside( W, w_in, w_out, K, n_in, n_out, s1_in, s2_in, s1_out, s2_out );

   return gf_phi_pp()[W][w_in][w_out][K][n_in][n_out][s1_in][s2_in][s1_out][s2_out]; 
}

dcomplex state_t::phi_ph( int W, int w_in, int w_out, int K, int n_in, int n_out, int s1_in, int s2_in, int s1_out, int s2_out ) const
{
   if ( W < -POS_BFREQ_COUNT_PHI || W > POS_BFREQ_COUNT_PHI ||
	 w_in < -POS_FFREQ_COUNT_PHI || w_in > POS_FFREQ_COUNT_PHI - 1 ||
	 w_out < -POS_FFREQ_COUNT_PHI || w_out > POS_FFREQ_COUNT_PHI - 1 )
      return phi_ph_outside( W, w_in, w_out, K, n_in, n_out, s1_in, s2_in, s1_out, s2_out );

   return gf_phi_ph()[W][w_in][w_out][K][n_in][n_out][s1_in][s2_in][s1_out][s2_out]; 
}

dcomplex state_t::phi_xph( int W, int w_in, int w_out, int K, int n_in, int n_out, int s1_in, int s2_in, int s1_out, int s2_out ) const
{
   if ( W < -POS_BFREQ_COUNT_PHI || W > POS_BFREQ_COUNT_PHI ||
	 w_in < -POS_FFREQ_COUNT_PHI || w_in > POS_FFREQ_COUNT_PHI - 1 ||
	 w_out < -POS_FFREQ_COUNT_PHI || w_out > POS_FFREQ_COUNT_PHI - 1 )
      return phi_xph_outside( W, w_in, w_out, K, n_in, n_out, s1_in, s2_in, s1_out, s2_out );

   return gf_phi_xph()[W][w_in][w_out][K][n_in][n_out][s1_in][s2_in][s1_out][s2_out];
}

dcomplex state_t::phi_pp_outside( int W, int w_in, int w_out, int K, int n_in, int n_out, int s1_in, int s2_in, int s1_out, int s2_out ) const
{
   return chi_pp( W, K, s1_in, s2_in, s1_out, s2_out ) * (double)((n_in==0) * (n_out==0))  + 
      P_pp( W, w_in, K, n_in, s1_in, s2_in, s1_out, s2_out ) * double(n_out==0) + 
      P_pp( W, w_out, K, n_out, s1_out, s2_out, s1_in, s2_in ) * double(n_in==0);  // time reversal symmetry used
}

dcomplex state_t::phi_ph_outside( int W, int w_in, int w_out, int K, int n_in, int n_out, int s1_in, int s2_in, int s1_out, int s2_out ) const
{
   int bf[2] = {0,0};
   int minusK = neg_k(K,bf);
   int bf_sign_x = (1-(1-translate_2pi_x(n_in)*translate_2pi_x(n_out))*bf[0]);
   int bf_sign_y = (1-(1-translate_2pi_y(n_in)*translate_2pi_y(n_out))*bf[1]);

   return chi_ph( W, K, s1_in, s2_in, s1_out, s2_out ) * double((n_in==0) * (n_out==0)) +
      P_ph( W, w_in, K, n_in, s1_in, s2_in, s1_out, s2_out ) * double(n_out==0) + 
      P_ph( -W, w_out, minusK, n_out, s2_in, s1_in, s2_out, s1_out ) * double(bf_sign_x * bf_sign_y * (n_in==0));  // flip diagram horizontally
}

dcomplex state_t::phi_xph_outside( int W, int w_in, int w_out, int K, int n_in, int n_out, int s1_in, int s2_in, int s1_out, int s2_out ) const
{
   int bf[2] = {0,0};
   int minusK = neg_k(K,bf);
   int bf_sign_x = (1-(1-translate_2pi_x(n_in)*translate_2pi_x(n_out))*bf[0]);
   int bf_sign_y = (1-(1-translate_2pi_y(n_in)*translate_2pi_y(n_out))*bf[1]);
   return chi_xph( W, K, s1_in, s2_in, s1_out, s2_out ) * double((n_in==0) * (n_out==0)) + 
      P_xph( W, w_in, K,n_in, s1_in, s2_in, s1_out, s2_out ) * double(n_out==0) + 
      P_xph( -W, w_out, minusK, n_out, s2_in, s1_in, s2_out, s1_out ) * double(bf_sign_x * bf_sign_y * (n_in==0));  // flip diagram horizontally
}

dcomplex state_t::chi_pp( int W, int K, int s1_in, int s2_in, int s1_out, int s2_out ) const
{
   if ( W < -POS_BFREQ_COUNT_CHI || W > POS_BFREQ_COUNT_CHI ) 
      return 0.0; 

   return gf_chi_pp()[W][K][s1_in][s2_in][s1_out][s2_out]; 
}

dcomplex state_t::chi_ph( int W, int K, int s1_in, int s2_in, int s1_out, int s2_out ) const
{
   if ( W < -POS_BFREQ_COUNT_CHI || W > POS_BFREQ_COUNT_CHI ) 
      return 0.0; 

   return gf_chi_ph()[W][K][s1_in][s2_in][s1_out][s2_out]; 
}

dcomplex state_t::chi_xph( int W, int K, int s1_in, int s2_in, int s1_out, int s2_out ) const
{   
   if ( W < -POS_BFREQ_COUNT_CHI || W > POS_BFREQ_COUNT_CHI ) 
      return 0.0; 

   return gf_chi_xph()[W][K][s1_in][s2_in][s1_out][s2_out]; 
}

dcomplex state_t::P_pp( int W, int w, int K, int n, int s1_in, int s2_in, int s1_out, int s2_out ) const
{
   if ( W < -POS_BFREQ_COUNT_P || W > POS_BFREQ_COUNT_P ||
	 w < -POS_FFREQ_COUNT_P || w > POS_FFREQ_COUNT_P - 1 )
      return 0.0; 

   return gf_P_pp()[W][w][K][n][s1_in][s2_in][s1_out][s2_out]; 
}

dcomplex state_t::P_ph( int W, int w, int K, int n, int s1_in, int s2_in, int s1_out, int s2_out ) const
{
   if ( W < -POS_BFREQ_COUNT_P || W > POS_BFREQ_COUNT_P ||
	 w < -POS_FFREQ_COUNT_P || w > POS_FFREQ_COUNT_P - 1 )
      return 0.0;

   return gf_P_ph()[W][w][K][n][s1_in][s2_in][s1_out][s2_out]; 
}

dcomplex state_t::P_xph( int W, int w, int K, int n, int s1_in, int s2_in, int s1_out, int s2_out ) const
{
   if ( W < -POS_BFREQ_COUNT_P || W > POS_BFREQ_COUNT_P ||
	 w < -POS_FFREQ_COUNT_P || w > POS_FFREQ_COUNT_P - 1 )
      return 0.0;

   return gf_P_xph()[W][w][K][n][s1_in][s2_in][s1_out][s2_out]; 
}

dcomplex state_t::R_pp( int W, int w_in, int w_out, int K, int n_in, int n_out, int s1_in, int s2_in, int s1_out, int s2_out ) const
{
   if ( W < -POS_BFREQ_COUNT_PHI || W > POS_BFREQ_COUNT_PHI ||
	 w_in < -POS_FFREQ_COUNT_PHI || w_in > POS_FFREQ_COUNT_PHI - 1 ||
	 w_out < -POS_FFREQ_COUNT_PHI || w_out > POS_FFREQ_COUNT_PHI - 1 )
      return 0.0;

   return phi_pp( W, w_in, w_out, K, n_in, n_out, s1_in, s2_in, s1_out, s2_out ) -
      phi_pp_outside( W, w_in, w_out, K, n_in, n_out, s1_in, s2_in, s1_out, s2_out );  
}

dcomplex state_t::R_ph( int W, int w_in, int w_out, int K, int n_in, int n_out, int s1_in, int s2_in, int s1_out, int s2_out ) const
{
   if ( W < -POS_BFREQ_COUNT_PHI || W > POS_BFREQ_COUNT_PHI ||
	 w_in < -POS_FFREQ_COUNT_PHI || w_in > POS_FFREQ_COUNT_PHI - 1 ||
	 w_out < -POS_FFREQ_COUNT_PHI || w_out > POS_FFREQ_COUNT_PHI - 1 )
      return 0.0;

   return phi_ph( W, w_in, w_out, K, n_in, n_out, s1_in, s2_in, s1_out, s2_out ) -
      phi_ph_outside( W, w_in, w_out, K, n_in, n_out, s1_in, s2_in, s1_out, s2_out ); 
}

dcomplex state_t::R_xph( int W, int w_in, int w_out, int K, int n_in, int n_out, int s1_in, int s2_in, int s1_out, int s2_out ) const
{
   if ( W < -POS_BFREQ_COUNT_PHI || W > POS_BFREQ_COUNT_PHI ||
	 w_in < -POS_FFREQ_COUNT_PHI || w_in > POS_FFREQ_COUNT_PHI - 1 ||
	 w_out < -POS_FFREQ_COUNT_PHI || w_out > POS_FFREQ_COUNT_PHI - 1 )
      return 0.0;

   return phi_xph( W, w_in, w_out, K, n_in, n_out, s1_in, s2_in, s1_out, s2_out ) -
      phi_xph_outside( W, w_in, w_out, K, n_in, n_out, s1_in, s2_in, s1_out, s2_out );  
}


dcomplex state_t::Gval( const idx_1p_t& idx) const
{
   return G( w_val( idx(0) ), get_kx(idx(1)), get_ky(idx(1)), SigMat( idx(0), idx(1) ) )( idx(2), idx(3) ); 
}

MatQN state_t::GMat( const idx_1p_mat_t& idx) const
{
   int idx_k = p_to_k_patch(idx(1));
   return G( w_val( idx(0) ), get_px(idx(1)), get_py(idx(1)), SigMat( idx(0), idx_k ) ); 
}

/******************************************** SUSCEPTIBILITIES *************************************************************************************/

dcomplex state_t::suscept_trip( int W, int K, int n_in, int n_out, int s1_in, int s2_in, int s1_out, int s2_out ) const
{
   if ( W < -POS_BFREQ_COUNT_SUSCEPT || W > POS_BFREQ_COUNT_SUSCEPT ) 
      return 0.0; 

   return gf_suscept_trip()[W][K][n_in][n_out][s1_in][s2_in][s1_out][s2_out]; 
}

dcomplex state_t::suscept_s( int W, int K, int n_in, int n_out, int s1_in, int s2_in, int s1_out, int s2_out ) const
{
   if ( W < -POS_BFREQ_COUNT_SUSCEPT || W > POS_BFREQ_COUNT_SUSCEPT ) 
      return 0.0; 

   return gf_suscept_s()[W][K][n_in][n_out][s1_in][s2_in][s1_out][s2_out]; 
}

dcomplex state_t::suscept_d( int W, int K,int n_in, int n_out, int s1_in, int s2_in, int s1_out, int s2_out ) const
{
   if ( W < -POS_BFREQ_COUNT_SUSCEPT || W > POS_BFREQ_COUNT_SUSCEPT ) 
      return 0.0; 

   return gf_suscept_d()[W][K][n_in][n_out][s1_in][s2_in][s1_out][s2_out]; 
}

dcomplex state_t::suscept_m( int W, int K, int n_in, int n_out, int s1_in, int s2_in, int s1_out, int s2_out ) const
{
   if ( W < -POS_BFREQ_COUNT_SUSCEPT || W > POS_BFREQ_COUNT_SUSCEPT ) 
      return 0.0; 

   return gf_suscept_m()[W][K][n_in][n_out][s1_in][s2_in][s1_out][s2_out]; 
}

//PROJECTION FUNCTIONS
//TO PP
dcomplex state_t::proj_pp_phi_ph( int W, int w_in, int w_out, int K, int n_in, int n_out, int s1_in, int s2_in, int s1_out, int s2_out ) const
{
   dcomplex val( 0.0, 0.0 );

   int W_ph = w_out - w_in;
   for(int m = 0; m<FFACTOR_COUNT; ++m)
     for(int mp = 0; mp<FFACTOR_COUNT; ++mp)
      for(int Q = 0; Q<PATCH_COUNT; ++Q)
	val+= proj_matrix_ph_to_pp[K][Q][n_in][n_out][m][mp] * phi_ph(W_ph, w_in + div2_ceil( W ) + div2_floor( W_ph ), div2_floor( W ) - w_out + div2_floor( W_ph ) - 1, Q, m, mp, s1_in, s2_in, s1_out, s2_out ); 
   
   val *= 1.0/PATCH_COUNT; //TODO: FFACTOR_COUNT?
   return val;
}

dcomplex state_t::proj_pp_phi_xph( int W, int w_in, int w_out, int K, int n_in, int n_out, int s1_in, int s2_in, int s1_out, int s2_out ) const
{
   dcomplex val( 0.0, 0.0 );

   int W_xph = - ( W + 100000 ) % 2 - w_out - w_in - 1;
   for(int m = 0; m<FFACTOR_COUNT; ++m)
     for(int mp = 0; mp<FFACTOR_COUNT; ++mp)
      for(int Q = 0; Q<PATCH_COUNT; ++Q)
	val+= (double)Parity(n_out) * proj_matrix_ph_to_pp[K][Q][n_in][n_out][m][mp] * phi_xph(W_xph, w_in + div2_ceil( W ) + div2_floor( W_xph ), w_out + div2_ceil( W ) + div2_floor( W_xph ), Q, m, mp, s1_in, s2_in, s1_out, s2_out ); // symmetry of projection matrix according to paper TUfRG
   
   val *= 1.0/PATCH_COUNT; //TODO: FFACTOR_COUNT?
   return val;
}

//TO PH
dcomplex state_t::proj_ph_phi_pp( int W, int w_in, int w_out, int K, int n_in, int n_out, int s1_in, int s2_in, int s1_out, int s2_out ) const
{
   dcomplex val( 0.0, 0.0 );

   int W_pp = w_in + w_out + ( W + 100000 ) % 2 + 1;
   for(int m = 0; m<FFACTOR_COUNT; ++m)
     for(int mp = 0; mp<FFACTOR_COUNT; ++mp)
      for(int Q = 0; Q<PATCH_COUNT; ++Q)
	val+= proj_matrix_pp_to_ph[K][Q][n_in][n_out][m][mp] * phi_pp( W_pp, w_in - div2_floor( W ) - div2_ceil( W_pp ), w_in + div2_ceil( W ) - div2_ceil( W_pp ), Q, m, mp, s1_in, s2_in, s1_out, s2_out ); 
   
   val *= 1.0/PATCH_COUNT; //TODO: FFACTOR_COUNT?
   return val;
}

dcomplex state_t::proj_ph_phi_xph( int W, int w_in, int w_out, int K, int n_in, int n_out, int s1_in, int s2_in, int s1_out, int s2_out ) const
{
   dcomplex val( 0.0, 0.0 );

   int W_xph = w_out - w_in;
   for(int m = 0; m<FFACTOR_COUNT; ++m)
     for(int mp = 0; mp<FFACTOR_COUNT; ++mp)
      for(int Q = 0; Q<PATCH_COUNT; ++Q)
	val+= proj_matrix_ph_to_xph[K][Q][n_in][n_out][m][mp] * phi_xph( W_xph, w_in - div2_floor( W ) + div2_floor( W_xph ), w_out + div2_ceil( W ) - div2_ceil( W_xph ), Q, m, mp, s1_in, s2_in, s1_out, s2_out ); // symmetry of projection matrix according to paper TUfRG
   
   val *= 1.0/PATCH_COUNT; //TODO: FFACTOR_COUNT?
   return val;
}

//TO XPH
dcomplex state_t::proj_xph_phi_pp( int W, int w_in, int w_out, int K, int n_in, int n_out, int s1_in, int s2_in, int s1_out, int s2_out ) const
{
   dcomplex val( 0.0, 0.0 );

   int W_pp = w_in + w_out + ( W + 100000 ) % 2 + 1;
   for(int m = 0; m<FFACTOR_COUNT; ++m)
     for(int mp = 0; mp<FFACTOR_COUNT; ++mp)
      for(int Q = 0; Q<PATCH_COUNT; ++Q)
	val+=  (double)Parity(mp) * proj_matrix_pp_to_ph[K][Q][n_in][n_out][m][mp] * phi_pp( W_pp, w_in - div2_floor( W ) - div2_ceil( W_pp ), w_out - div2_floor( W ) - div2_ceil( W_pp ), Q, m, mp, s1_in, s2_in, s1_out, s2_out );// symmetry of projection matrix according to paper TUfRG 
   
   val *= 1.0/PATCH_COUNT; //TODO: FFACTOR_COUNT?
   return val;
}

dcomplex state_t::proj_xph_phi_ph( int W, int w_in, int w_out, int K, int n_in, int n_out, int s1_in, int s2_in, int s1_out, int s2_out ) const
{
   dcomplex val( 0.0, 0.0 );

   int W_ph = w_out - w_in;
   for(int m = 0; m<FFACTOR_COUNT; ++m)
     for(int mp = 0; mp<FFACTOR_COUNT; ++mp)
      for(int Q = 0; Q<PATCH_COUNT; ++Q)
	val+= proj_matrix_ph_to_xph[K][Q][n_in][n_out][m][mp] * phi_ph(W_ph, w_in - div2_floor( W ) + div2_floor( W_ph ), w_out + div2_ceil( W ) - div2_ceil( W_ph ), Q, m, mp, s1_in, s2_in, s1_out, s2_out ); 
   
   val *= 1.0/PATCH_COUNT; //TODO: FFACTOR_COUNT?
   return val;
}
