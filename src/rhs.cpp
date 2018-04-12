
/************************************************************************************************//**
 *  		
 * 	file: 		rhs.cpp
 * 	contents:   	See rhs.h
 * 
 ****************************************************************************************************/


#include <rhs.h>
#include <frg.h>
#include <params.h>
#include <cmath>
#include <translate.h>
#include <mymath.h>
#include <projection.h>
#include <complex>

using namespace std; 
using boost::bind;

// --- Create static objects: symmetry groups, weight vector

// ALL SYMM

symmetry_grp_t<dcomplex,4> rhs_t::symm_grp_sig( gf_1p_t(),{
       cconj_sig
      , rot_pi2_z_sig
      , mirror_y_sig
      } );

symmetry_grp_t<dcomplex,8> rhs_t::symm_grp_susc_sc( gf_susc_t(), { 
      hmirror_susc_sc
      , cconj_susc_sc
      , timerev_susc_sc
      , rot_pi2_z_susc
      , mirror_y_susc
      } ); 

symmetry_grp_t<dcomplex,8> rhs_t::symm_grp_susc_d( gf_susc_t(), { 
      hmirror_susc_d
      , cconj_susc_d
      , timerev_susc_d
      , rot_pi2_z_susc
      , mirror_y_susc
      } ); 

symmetry_grp_t<dcomplex,8> rhs_t::symm_grp_susc_m( gf_susc_t(), { 
      hmirror_susc_m
      , cconj_susc_m
      , timerev_susc_m
      , rot_pi2_z_susc
      , mirror_y_susc
      } ); 

symmetry_grp_t<dcomplex,9> rhs_t::symm_grp_tri_sc( gf_tri_t(),{
      hmirror_tri_sc
      , cconj_timerev_tri_sc
      , rot_pi2_z_tri
      , mirror_y_tri
      } ); 

symmetry_grp_t<dcomplex,9> rhs_t::symm_grp_tri_d( gf_tri_t(),{ 
      timerev_tri_d
      , cconj_tri_d
      , rot_pi2_z_tri
      , mirror_y_tri
      } ); 

symmetry_grp_t<dcomplex,9> rhs_t::symm_grp_tri_m( gf_tri_t(),{ 
       hmirror_cconj_tri_m 
      , hmirror_timerev_tri_m
      , rot_pi2_z_tri
      , mirror_y_tri
      } ); 

symmetry_grp_t<dcomplex,7> rhs_t::symm_grp_asytri_sc( gf_asytri_t(),{
      hmirror_asytri_sc
      , cconj_timerev_asytri_sc
      , rot_pi2_z_asytri
      , mirror_y_asytri
      } ); 

symmetry_grp_t<dcomplex,7> rhs_t::symm_grp_asytri_d( gf_asytri_t(),{ 
      timerev_asytri_d
      , cconj_asytri_d
      , rot_pi2_z_asytri
      , mirror_y_asytri
      } ); 

symmetry_grp_t<dcomplex,7> rhs_t::symm_grp_asytri_m( gf_asytri_t(),{ 
       hmirror_cconj_asytri_m 
      , hmirror_timerev_asytri_m
      , rot_pi2_z_asytri
      , mirror_y_asytri
      } ); 


symmetry_grp_t<MatReal,3> rhs_t::symm_grp_Gvec_real(gf_1p_mat_real_t(),{});

symmetry_grp_t<MatPatch,8> rhs_t::symm_grp_bubble_pp(gf_bubble_mat_t(),{swap_bubble, cconj_bubble, hmirror_bubble_pp});
symmetry_grp_t<MatPatch,8> rhs_t::symm_grp_bubble_ph(gf_bubble_mat_t(),{swap_bubble, cconj_bubble});


gf_weight_vec_t rhs_t::weight_vec( generate_weights( POS_INT_RANGE-TAIL_LENGTH, TAIL_LENGTH, FIT_ORDER ) ); 
gf<double, 2> rhs_t::weight_vec_2d( generate_2d_weights( POS_INT_RANGE-TAIL_LENGTH, TAIL_LENGTH, FIT_ORDER ) ); 

void rhs_t::operator() ( const state_t& state_vec, state_t& dfdl )
{

   // Precalculate Green function and single-scale propagator on frequency grid

   gf_1p_mat_t Gvec( POS_1P_RANGE, FFT_DIM*FFT_DIM ); 
   
   Gvec.init( bind( &state_t::GMat, boost::cref(state_vec), _1 ) ); // Initialize big Green function vector 
   
   F_factors ffactor_mom;
   
   //gf_1p_mat_t Gvec_big( POS_1P_RANGE_SE, FFT_DIM*FFT_DIM ); 
   
   //Gvec_big.init( bind( &state_t::GMat, boost::cref(state_vec), _1 ) ); // Initialize big Green function vector 

   // Precalculate Green's functions in real space by FFT
   gf_1p_mat_real_t Gvec_real;
   gf_bubble_mat_t bubble_pp;
   gf_bubble_mat_t bubble_ph;
   
   cout << " ... bubble PT1 " << endl;
   fftw_complex *input, *output;	// input output functions for the FFTW
   fftw_plan p, p_b; 			// plan for FFTW
   
   input = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * FFT_DIM * FFT_DIM);
   output = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * FFT_DIM * FFT_DIM);

   p = fftw_plan_dft_2d(FFT_DIM, FFT_DIM, input, output, FFTW_BACKWARD, FFTW_ESTIMATE); //plan for FFTW, same for Gvec and Svec (BACKWARD), different for bubble (FORWARD)
   
   
   symm_grp_Gvec_real.init( Gvec_real, [&Gvec, p]( const idx_1p_mat_real_t& idx ){ return eval_Gvec_real( idx, Gvec, p); } );
   
   
   cout << " ... bubble PT2 " << endl;
   
   p_b = fftw_plan_dft_2d(FFT_DIM, FFT_DIM, input, output, FFTW_FORWARD, FFTW_ESTIMATE); //plan for FFTW, same for Gvec and Svec (BACKWARD), different for bubble (FORWARD)
   
   symm_grp_bubble_pp.init( bubble_pp, [&Gvec_real, p_b]( const idx_bubble_mat_t& idx ){ return eval_diag_bubble_pp( idx, Gvec_real, p_b ); } ); 
   symm_grp_bubble_ph.init( bubble_ph, [&Gvec_real, p_b]( const idx_bubble_mat_t& idx ){ return eval_diag_bubble_ph( idx, Gvec_real, p_b ); } ); 

   fftw_destroy_plan(p);
   fftw_destroy_plan(p_b);
   fftw_free(input); fftw_free(output);

   cout << " Rhs calculation " << endl;

   /*********************  Open MP parallelized RHS calculation  ********************/

   // --- Conventional flow
   //cout << " ... Self Energy SDE " << endl;
   //symm_grp_sig.init( dfdl.gf_Sig_old(), [&state_vec, &Gvec]( const idx_1p_t& idx ){ return eval_rhs_Sig_old( idx, state_vec, Gvec ); } );
   //cout << " ... Self Energy SDE " << endl;
   //symm_grp_sig.init( dfdl.gf_Sig(), [&state_vec, &Gvec, &bubble_ph, &bubble_pp, &ffactor_mom]( const idx_1p_t& idx ){ return eval_rhs_Sig( idx, state_vec, Gvec, bubble_ph, bubble_pp, ffactor_mom); } );
   //cout << " ... Self Energy PH SDE " << endl;
   //symm_grp_sig.init( dfdl.gf_Sig_ph(), [&state_vec, &Gvec, &bubble_ph, &ffactor_mom]( const idx_1p_t& idx ){ return eval_rhs_Sig_ph( idx, state_vec, Gvec, bubble_ph, ffactor_mom ); } );
   //cout << " ... Self Energy PP SDE " << endl;
   //symm_grp_sig.init( dfdl.gf_Sig_pp(), [&state_vec, &Gvec, &bubble_pp, &ffactor_mom]( const idx_1p_t& idx ){ return eval_rhs_Sig_pp( idx, state_vec, Gvec, bubble_pp, ffactor_mom ); } );
   //cout << " ... Self Energy XPH SDE " << endl;
   //symm_grp_sig.init( dfdl.gf_Sig_xph(), [&state_vec, &Gvec, &bubble_ph, &ffactor_mom]( const idx_1p_t& idx ){ return eval_rhs_Sig_xph( idx, state_vec, Gvec, bubble_ph, ffactor_mom ); } );
   
   cout << " ... Susceptibilities " << endl;
   
   cout << " ... SC " << endl;
   symm_grp_susc_sc.init( dfdl.gf_susc_sc(), [&state_vec, &bubble_pp]( const idx_susc_t& idx ){ return eval_diag_susc_sc( idx, state_vec, bubble_pp ); } ); 
   cout << " ... Density " << endl;
   symm_grp_susc_d.init( dfdl.gf_susc_d(),     [&state_vec, &bubble_ph]( const idx_susc_t& idx ){ return eval_diag_susc_d( idx, state_vec, bubble_ph ); } ); 
   cout << " ... Magnetic " << endl;
   symm_grp_susc_m.init( dfdl.gf_susc_m(),      [&state_vec, &bubble_ph]( const idx_susc_t& idx ){ return eval_diag_susc_m( idx, state_vec, bubble_ph ); } ); 
   
   cout << " ... trileg asymptotics " << endl;
   cout << " ... Sc " << endl;
   symm_grp_asytri_sc.init( dfdl.gf_asytri_sc(),    [&state_vec, &bubble_pp]( const idx_asytri_t& idx ){ return eval_diag_tri_sc( iasy_to_itri(idx), state_vec, bubble_pp ); } ); 
   cout << " ... Density " << endl;
   symm_grp_asytri_d.init( dfdl.gf_asytri_d(),     [&state_vec, &bubble_ph]( const idx_asytri_t& idx ){ return eval_diag_tri_d( iasy_to_itri(idx), state_vec, bubble_ph ); } ); 
   cout << " ... Magnetic " << endl;
   symm_grp_asytri_m.init( dfdl.gf_asytri_m(),      [&state_vec, &bubble_ph]( const idx_asytri_t& idx ){ return eval_diag_tri_m( iasy_to_itri(idx), state_vec, bubble_ph ); } ); 
   
   cout << " ... trileg " << endl;
   cout << " ... Sc " << endl;
   symm_grp_tri_sc.init( dfdl.gf_tri_sc(),    [&state_vec, &bubble_pp]( const idx_tri_t& idx ){ return eval_diag_tri_sc( idx, state_vec, bubble_pp ); } ); 
   cout << " ... Density " << endl;
   symm_grp_tri_d.init( dfdl.gf_tri_d(),     [&state_vec, &bubble_ph]( const idx_tri_t& idx ){ return eval_diag_tri_d( idx, state_vec, bubble_ph ); } ); 
   cout << " ... Magnetic " << endl;
   symm_grp_tri_m.init( dfdl.gf_tri_m(),      [&state_vec, &bubble_ph]( const idx_tri_t& idx ){ return eval_diag_tri_m( idx, state_vec, bubble_ph ); } ); 

}


// --- RHS DIAGRAMS ---

dcomplex rhs_t::eval_rhs_Sig_old( const idx_1p_t& idx, const state_t& state_vec, const gf_1p_mat_t& Gvec )
{
   dcomplex val( 0.0, 0.0 );
   //int count =0;
   for( int s1 = 0; s1 < QN_COUNT; ++s1 )
      for( int s1p = 0; s1p < QN_COUNT; ++s1p )
         for( int s2 = 0; s2 < QN_COUNT; ++s2 )
            for( int s2p = 0; s2p < QN_COUNT; ++s2p )
               for( int s3 = 0; s3 < QN_COUNT; ++s3 )
         	  for( int s3p = 0; s3p < QN_COUNT; ++s3p )
   	 	     for( int k = 0; k < PATCH_COUNT; ++k )
   	                for( int kp = 0; kp < PATCH_COUNT; ++kp )
   	                   for( int w = -POS_INT_RANGE; w < POS_INT_RANGE; ++w )
   	                      for( int wp = -POS_INT_RANGE; wp < POS_INT_RANGE; ++wp ){
				int p1 = k_to_p_patch(k);
   	             		int p3 = k_to_p_patch(kp);
                     		int bf[2];// if this is even(0) -> backfolding with +sign
   		     		int k_plus_kp = add_k(k, kp, bf); 
   		     		int k_plus_kp_minus_kin = dif_k(k_plus_kp,idx(I1P::k), bf);
   		     		int p2 = k_to_p_patch(k_plus_kp_minus_kin);  
				val += vert_bare(idx( I1P::s_in ), s3, s2p, s1p) * 
   		                       Gvec[w][p1](s1,s1p) * Gvec[wp][p3](s2,s2p) *
   			               Gvec[w+wp-idx(I1P::w)][p2](s3,s3p) * 
   			               weight_vec_2d[w][wp] *
   		                       (state_vec.vertx( w, wp, w+wp-idx(I1P::w), k, kp, k_plus_kp_minus_kin, s1, s2,idx( I1P::s_out ), s3p ));  
	    }

   val *= 1.0/BETA/BETA/PATCH_COUNT/PATCH_COUNT;
   //cout << "Sig val SDE : "  <<  -val  << endl;
   //cout << "Sig val: "  <<  gf_Sig_read[idx( I1P::w )][idx( I1P::k )][idx( I1P::s_in )][idx( I1P::s_out )] << endl; 
   return -val;
}


// --- New way of calculating the SE taking into account all channels as three contributions in their native channel

dcomplex rhs_t::eval_rhs_Sig( const idx_1p_t& idx, const state_t& state_vec, const gf_1p_mat_t& Gvec, gf_bubble_mat_t& bubble_ph, gf_bubble_mat_t& bubble_pp, F_factors ffactor_mom )
{
   int k  = idx(I1P::k);
   int w_in  = idx(I1P::w);
   dcomplex val( 0.0, 0.0 );
   for( int s1 = 0; s1 < QN_COUNT; ++s1 )
      for( int s1p = 0; s1p < QN_COUNT; ++s1p )
         for( int s2 = 0; s2 < QN_COUNT; ++s2 )
            for( int s2p = 0; s2p < QN_COUNT; ++s2p )
               for( int s3 = 0; s3 < QN_COUNT; ++s3 )
         	  for( int s3p = 0; s3p < QN_COUNT; ++s3p )
   	 	     for( int q = 0; q < PATCH_COUNT; ++q )
   	                for( int m = 0; m < FFACTOR_COUNT ; ++m )
   	                   for( int n = 0; n < FFACTOR_COUNT ; ++n )
   	                      for( int w = -POS_INT_RANGE; w < POS_INT_RANGE; ++w )
   	                         for( int wp = -POS_INT_RANGE; wp < POS_INT_RANGE; ++wp ){
				   
                     		  int bf[2];// if this is even(0) -> backfolding with +sign
   		     		  int k_plus_q = add_k(k, q, bf);
				  int k_minus_q = dif_k(k, q, bf);
				  int p_xph = k_to_p_patch(k_minus_q); 
				  int p_ph = k_to_p_patch(k_plus_q);
   		     		  int q_minus_k = dif_k(q, k, bf); 
				  int p_pp = k_to_p_patch(q_minus_k);
				  int qx = get_kx(q); 
				  int qy = get_ky(q);
				  int kx = get_kx(k); 
				  int ky = get_ky(k); 
				   
				   val -= vert_bare(idx( I1P::s_in ), s3, s1p, s2p) * 
				       Gvec[w+idx(I1P::w)][p_ph](s1p,s1) * ffactor_mom[m](kx+0.5*(qx),ky+0.5*(qy)) *
				       weight_vec_2d[w][wp] *
	       		    	       bubble_ph[w][wp][0][n][s3][s3p][s2][s2p](q) * //TODO: to check the order of s
				       (state_vec.phi_ph( w, wp, w_in+div2_ceil(w), q, n, m, s1, s2,idx( I1P::s_out ), s3p ) +
					state_vec.phi_xph( w, wp, w_in+div2_ceil(w), q, n, m, s1, s2,idx( I1P::s_out ), s3p ) +
					state_vec.proj_vert_bare[0][m][s1][s2][idx(I1P::s_out)][s3p])+
				       vert_bare(idx( I1P::s_in ), s3, s1p, s2p) * 
				       Gvec[w-idx(I1P::w)][p_pp](s1p,s1) * ffactor_mom[m](0.5*(qx)-kx,0.5*(qy)-ky) *
				       weight_vec_2d[w][wp] *
	       		    	       bubble_pp[w][wp][0][n][s3][s3p][s2][s2p](q) * //TODO: to check the order of s
				       state_vec.phi_pp( w, wp, div2_floor(w)-w_in, q, n, m, s1, s2,idx( I1P::s_out ), s3p );  
				      
	    }

   val *= 1.0/BETA/BETA/PATCH_COUNT/2.0/PI;
   return val;
}

dcomplex rhs_t::eval_rhs_Sig_ph( const idx_1p_t& idx, const state_t& state_vec, const gf_1p_mat_t& Gvec, gf_bubble_mat_t& bubble_ph, F_factors ffactor_mom )
{
   int k  = idx(I1P::k);
   int w_in  = idx(I1P::w);
   dcomplex val( 0.0, 0.0 );
   //int count =0;
   for( int s1 = 0; s1 < QN_COUNT; ++s1 )
      for( int s1p = 0; s1p < QN_COUNT; ++s1p )
         for( int s2 = 0; s2 < QN_COUNT; ++s2 )
            for( int s2p = 0; s2p < QN_COUNT; ++s2p )
               for( int s3 = 0; s3 < QN_COUNT; ++s3 )
         	  for( int s3p = 0; s3p < QN_COUNT; ++s3p )
   	 	     for( int q = 0; q < PATCH_COUNT; ++q )
   	                for( int m = 0; m < FFACTOR_COUNT ; ++m )
   	                   for( int n = 0; n < FFACTOR_COUNT ; ++n )
   	                      for( int w = -POS_INT_RANGE; w < POS_INT_RANGE; ++w )
   	                         for( int wp = -POS_INT_RANGE; wp < POS_INT_RANGE; ++wp ){
				   
                     		  int bf[2];// if this is even(0) -> backfolding with +sign
   		     		  int k_plus_q = add_k(k, q, bf); 
				  int p = k_to_p_patch(k_plus_q);
				  int qx = get_kx(q); 
				  int qy = get_ky(q);
				  int kx = get_kx(k); 
				  int ky = get_ky(k); 
				   
				   val -= vert_bare(idx( I1P::s_in ), s3, s1p, s2p) * 
				       Gvec[w+idx(I1P::w)][p](s1p,s1) * ffactor_mom[m](kx+0.5*(qx),ky+0.5*(qy)) *
				       weight_vec_2d[w][wp] *
	       		    	       bubble_ph[w][wp][0][n][s3][s3p][s2][s2p](q) * //TODO: to check the order of s
				       state_vec.vertx_ph( w, wp, w_in+div2_ceil(w), q, n, m, s1, s2,idx( I1P::s_out ), s3p );  
	    }

   val *= 1.0/BETA/BETA/PATCH_COUNT/2.0/PI;
   cout << "Sig val SDE : "  <<  val  << endl;
   cout << "Sig val: "  <<  gf_Sig_read[idx( I1P::w )][idx( I1P::k )][idx( I1P::s_in )][idx( I1P::s_out )] << endl; 
   //return -val + gf_Sig_read[idx( I1P::w )][idx( I1P::k )][idx( I1P::s_in )][idx( I1P::s_out )];
   return val;
}

dcomplex rhs_t::eval_rhs_Sig_pp( const idx_1p_t& idx, const state_t& state_vec, const gf_1p_mat_t& Gvec, gf_bubble_mat_t& bubble_pp , F_factors ffactor_mom )
{
   int k  = idx(I1P::k);
   int w_in  = idx(I1P::w);
   dcomplex val( 0.0, 0.0 );
   //int count =0;
   for( int s1 = 0; s1 < QN_COUNT; ++s1 )
      for( int s1p = 0; s1p < QN_COUNT; ++s1p )
         for( int s2 = 0; s2 < QN_COUNT; ++s2 )
            for( int s2p = 0; s2p < QN_COUNT; ++s2p )
               for( int s3 = 0; s3 < QN_COUNT; ++s3 )
         	  for( int s3p = 0; s3p < QN_COUNT; ++s3p )
   	 	     for( int q = 0; q < PATCH_COUNT; ++q )
   	                for( int m = 0; m < FFACTOR_COUNT ; ++m )
   	                   for( int n = 0; n < FFACTOR_COUNT ; ++n )
   	                      for( int w = -POS_INT_RANGE; w < POS_INT_RANGE; ++w )
   	                         for( int wp = -POS_INT_RANGE; wp < POS_INT_RANGE; ++wp ){
				   
                     		  int bf[2];// if this is even(0) -> backfolding with +sign
   		     		  int q_minus_k = dif_k(q, k, bf); 
				  int p = k_to_p_patch(q_minus_k);
				  int qx = get_kx(q); 
				  int qy = get_ky(q);
				  int kx = get_kx(k); 
				  int ky = get_ky(k); 
				   
				   val -= vert_bare(idx( I1P::s_in ), s3, s1p, s2p) * 
				       Gvec[w-idx(I1P::w)][p](s1p,s1) * ffactor_mom[m](0.5*(qx)-kx,0.5*(qy)-ky) *
				       weight_vec_2d[w][wp] *
	       		    	       bubble_pp[w][wp][0][n][s3][s3p][s2][s2p](q) * //TODO: to check the order of s
				       state_vec.vertx_pp( w, wp, div2_floor(w)-w_in, q, n, m, s1, s2,idx( I1P::s_out ), s3p );  
	    }

   val *= 1.0/BETA/BETA/PATCH_COUNT/2.0/PI;
   cout << "Sig val SDE : "  <<  val  << endl;
   cout << "Sig val: "  <<  gf_Sig_read[idx( I1P::w )][idx( I1P::k )][idx( I1P::s_in )][idx( I1P::s_out )] << endl; 
   //return -val + gf_Sig_read[idx( I1P::w )][idx( I1P::k )][idx( I1P::s_in )][idx( I1P::s_out )];
   return val;
}

dcomplex rhs_t::eval_rhs_Sig_xph( const idx_1p_t& idx, const state_t& state_vec, const gf_1p_mat_t& Gvec, gf_bubble_mat_t& bubble_ph , F_factors ffactor_mom )
{
   int k  = idx(I1P::k);
   int w_in  = idx(I1P::w);
   dcomplex val( 0.0, 0.0 );
   //int count =0;
   for( int s1 = 0; s1 < QN_COUNT; ++s1 )
      for( int s1p = 0; s1p < QN_COUNT; ++s1p )
         for( int s2 = 0; s2 < QN_COUNT; ++s2 )
            for( int s2p = 0; s2p < QN_COUNT; ++s2p )
               for( int s3 = 0; s3 < QN_COUNT; ++s3 )
         	  for( int s3p = 0; s3p < QN_COUNT; ++s3p )
   	 	     for( int q = 0; q < PATCH_COUNT; ++q )
   	                for( int m = 0; m < FFACTOR_COUNT ; ++m )
   	                   for( int n = 0; n < FFACTOR_COUNT ; ++n )
   	                      for( int w = -POS_INT_RANGE; w < POS_INT_RANGE; ++w )
   	                         for( int wp = -POS_INT_RANGE; wp < POS_INT_RANGE; ++wp ){
				   
                     		  int bf[2];// if this is even(0) -> backfolding with +sign
   		     		  int k_plus_q = add_k(k, q, bf); 
				  int p = k_to_p_patch(k_plus_q);
				  int qx = get_kx(q); 
				  int qy = get_ky(q);
				  int kx = get_kx(k); 
				  int ky = get_ky(k); 
				   
				   val -= vert_bare(idx( I1P::s_in ), s3, s1p, s2p) * 
				       Gvec[w+idx(I1P::w)][p](s1p,s1) * ffactor_mom[m](kx+0.5*(qx),ky+0.5*(qy)) *
				       weight_vec_2d[w][wp] *
	       		    	       bubble_ph[w][wp][0][n][s3][s3p][s2][s2p](q) * //TODO: to check the order of s
				       state_vec.vertx_xph( w, wp, w_in+div2_ceil(w), q, n, m, s1, s2,idx( I1P::s_out ), s3p );  
	    }

   val *= 1.0/BETA/BETA/PATCH_COUNT/2.0/PI;
   cout << "Sig val SDE : "  <<  val  << endl;
   cout << "Sig val: "  <<  gf_Sig_read[idx( I1P::w )][idx( I1P::k )][idx( I1P::s_in )][idx( I1P::s_out )] << endl; 
   return val;
}

MatReal rhs_t::eval_Gvec_real( const idx_1p_mat_real_t& idx, const gf_1p_mat_t& Gvec, fftw_plan p_g)
{
   int w = idx(0);
   int s1 = idx(1);
   int s2 = idx(2);  

   MatReal val_g;
   fftw_complex *input_g, *output_g;
   input_g = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * FFT_DIM * FFT_DIM );
   output_g = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * FFT_DIM * FFT_DIM );
   
   for(int p=0;p<FFT_DIM*FFT_DIM; ++p){
      input_g[p][0] = Gvec[w][p](s1,s2).real(); 
      input_g[p][1] = Gvec[w][p](s1,s2).imag(); 
   }
   
   // plan is created outside and passed to function because only fftw_execute_dft is thread save
   fftw_execute_dft(p_g,input_g,output_g); 
   
   for(int R=0;R<FFT_DIM*FFT_DIM; ++R){
      val_g(R) =  1.0/FFT_DIM/FFT_DIM * (output_g[Matidxshift(get_px_idx(R),get_py_idx(R))][0] + I * output_g[Matidxshift(get_px_idx(R),get_py_idx(R))][1]);
   }

   fftw_free(input_g); fftw_free(output_g);
   return val_g; 
}


MatPatch rhs_t::eval_diag_bubble_pp( const idx_bubble_mat_t& idx, const gf_1p_mat_real_t& Gvec_real, fftw_plan p_b )
{
   int W = idx(0);
   int w = idx(1);
   int m = idx(2);
   int n = idx(3);
   int s4 = idx(4);
   int s4p = idx(5);
   int s3 = idx(6);
   int s3p = idx(7);
   
   MatReal val = MatReal::Zero();
   MatPatch val_patch = MatPatch::Zero();
   fftw_complex *input_b, *output_b;
   input_b = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * FFT_DIM * FFT_DIM);
   output_b = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * FFT_DIM * FFT_DIM);
   
   for(int r=0; r< REAL_GRID; ++r){
    if(abs(F_factors_real::fft_fxf_weight[r][m][n]) > CHOP_ERR){// 1E-16 is error on the weights exponents, cosines,...
       int rx = std::get<0>(R_Grid()[r]);
       int ry = std::get<1>(R_Grid()[r]);
       for(int Rx=0;Rx<FFT_DIM; ++Rx)
       for(int Ry=0;Ry<FFT_DIM; ++Ry){	 
	       dcomplex temp = Gvec_real[ w + div2_ceil( W ) ][s4][s4p](Matidx(Rx,Ry)) * Gvec_real[ div2_floor( W ) - w - 1 ][s3][s3p](Matidx(Rx+rx,Ry+ry));
		  input_b[Matidxshift(Rx,Ry)][0] = temp.real();
	          input_b[Matidxshift(Rx,Ry)][1] = temp.imag();
	    }
      // plan is created outside and passed to function because only fftw_execute_dft is thread save
	fftw_execute_dft(p_b, input_b, output_b);
	// calculate weight and multiply it
       for (int K=0; K < PATCH_COUNT; ++K){
	  int Kx_idx = std::get<0>(K_Grid::get_indices(K));  
	  int Ky_idx = std::get<1>(K_Grid::get_indices(K));  
	  int Kx_shift = (Kx_idx*FFT_DIM/2/MAX_KPOS + 10*FFT_DIM) % FFT_DIM;
	  int Ky_shift = (Ky_idx*FFT_DIM/2/MAX_KPOS + 10*FFT_DIM) % FFT_DIM;
	  dcomplex weight = F_factors_real::fft_fxf_weight[r][m][n];
	  weight *= (cos_k(Kx_idx*rx + Ky_idx*ry) - I * sin_k(Kx_idx*rx + Ky_idx*ry));
          val_patch(K) += (output_b[Matidx(Kx_shift,Ky_shift)][0]+ I * output_b[Matidx(Kx_shift,Ky_shift)][1]) * weight; 
       }  
    }
   }
  
  fftw_free(input_b); fftw_free(output_b);
  return val_patch;

} 

MatPatch rhs_t::eval_diag_bubble_ph( const idx_bubble_mat_t& idx, const gf_1p_mat_real_t& Gvec_real, fftw_plan p_b)
{
   int W = idx(0);
   int w = idx(1);
   int m = idx(2);
   int n = idx(3);
   int s4 = idx(4);
   int s4p = idx(5);
   int s3 = idx(6);
   int s3p = idx(7);
   
   MatReal val = MatReal::Zero();
   MatPatch val_patch = MatPatch::Zero();
   fftw_complex *input_b, *output_b;
   input_b = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * FFT_DIM * FFT_DIM);
   output_b = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * FFT_DIM * FFT_DIM);
   
   for(int r=0; r< REAL_GRID; ++r){
     if(abs(F_factors_real::fft_fxf_weight[r][m][n])>  CHOP_ERR){// 1E-16 is error on the weights exponents, cosines,...
	int rx = std::get<0>(R_Grid()[r]);
        int ry = std::get<1>(R_Grid()[r]);
        for(int Rx=0;Rx<FFT_DIM; ++Rx)
        for(int Ry=0;Ry<FFT_DIM; ++Ry){ 
	    dcomplex temp = Gvec_real[ w - div2_floor( W ) ][s4][s4p](Matidx(FFT_DIM-Rx,FFT_DIM-Ry)) * Gvec_real[ div2_ceil( W ) + w ][s3][s3p](Matidx(Rx-rx,Ry-ry));
	    input_b[Matidxshift(Rx,Ry)][0] = temp.real();
	    input_b[Matidxshift(Rx,Ry)][1] = temp.imag();
	}
	// plan is created outside and passed to function because only fftw_execute_dft is thread save
	fftw_execute_dft(p_b,input_b,output_b);
	// calculate weight and multiply it
	for (int K=0; K < PATCH_COUNT; ++K){
	int Kx_idx = std::get<0>(K_Grid::get_indices(K));  
	int Ky_idx = std::get<1>(K_Grid::get_indices(K));  
	int Kx_shift = (Kx_idx*FFT_DIM/2/MAX_KPOS + 10*FFT_DIM) % FFT_DIM;
	int Ky_shift = (Ky_idx*FFT_DIM/2/MAX_KPOS + 10*FFT_DIM) % FFT_DIM;
	dcomplex weight = F_factors_real::fft_fxf_weight[r][m][n];
	weight *= (cos_k(Kx_idx*rx + Ky_idx*ry) + I * sin_k(Kx_idx*rx + Ky_idx*ry));
	val_patch(K) += (output_b[Matidx(Kx_shift,Ky_shift)][0]+ I * output_b[Matidx(Kx_shift,Ky_shift)][1]) * weight;
	}
     }
   }

  fftw_free(input_b); fftw_free(output_b);
  return val_patch;
}


dcomplex rhs_t::eval_diag_susc_sc( const idx_susc_t& idx, const state_t& state_vec, const gf_bubble_mat_t& bubble_pp)
{
   dcomplex val( 0.0, 0.0 );

   // Introduce help variables
   int W = idx( ISUSC::W );

   int K = idx( ISUSC::K );

   int n_in = idx( ISUSC::n_in);
   int n_out = idx( ISUSC::n_out);

   int s1_in = idx( ISUSC::s1_in );
   int s2_in = idx( ISUSC::s2_in );
   int s1_out = idx( ISUSC::s1_out );
   int s2_out = idx( ISUSC::s2_out );
   
   for( int s3 = 0; s3 < QN_COUNT; ++s3 )
      for( int s3p = 0; s3p < QN_COUNT; ++s3p )
	 for( int s4 = 0; s4 < QN_COUNT; ++s4 )
	    for( int s4p = 0; s4p < QN_COUNT; ++s4p )
	       for( int m = 0; m < FFACTOR_COUNT; ++m )
		  for( int mp = 0; mp < FFACTOR_COUNT; ++mp )
		     for( int w = -POS_INT_RANGE; w < POS_INT_RANGE; ++w )
		       for( int wp = -POS_INT_RANGE; wp < POS_INT_RANGE; ++wp )
	       	       {
	       		val += 1.0/BETA/BETA/4.0/PI/PI	* // Two pairs of equivalent lines -> Factor 1/4
			    bubble_pp[W][w][n_in][mp][s4][s2_in][s3][s1_in](K) * 
			    state_vec.vertx_pp( W, w, wp, K, mp, m, s3, s4, s3p, s4p ) *  
	       		    bubble_pp[W][wp][m][n_out][s2_out][s4p][s1_out][s3p](K) *
			    weight_vec_2d[w][wp];
	       	       }
   for( int w = -POS_INT_RANGE; w < POS_INT_RANGE; ++w )
   {
      val += 1.0/BETA *	//Two equivalent lines -> factor 1/2
	 bubble_pp[W][w][n_in][n_out][s2_out][s2_in][s1_out][s1_in](K) *
	 weight_vec[w];
   }

   return val; 
}

dcomplex rhs_t::eval_diag_susc_d( const idx_susc_t& idx, const state_t& state_vec, const gf_bubble_mat_t& bubble_ph)
{
   dcomplex val( 0.0, 0.0 );

   // Introduce help variables
   int W = idx( ISUSC::W );

   int K = idx( ISUSC::K );
   
   int n_in = idx( ISUSC::n_in);
   int n_out = idx( ISUSC::n_out);

   int s1_in = idx( ISUSC::s1_in );
   int s2_in = idx( ISUSC::s2_in );
   int s1_out = idx( ISUSC::s1_out );
   int s2_out = idx( ISUSC::s2_out );
   
   for( int s3 = 0; s3 < QN_COUNT; ++s3 )
      for( int s3p = 0; s3p < QN_COUNT; ++s3p )
	 for( int s4 = 0; s4 < QN_COUNT; ++s4 )
	    for( int s4p = 0; s4p < QN_COUNT; ++s4p )
	       for( int m = 0; m < FFACTOR_COUNT; ++m )
		  for( int mp = 0; mp < FFACTOR_COUNT; ++mp )
		     for( int w = -POS_INT_RANGE; w < POS_INT_RANGE; ++w )
		       for( int wp = -POS_INT_RANGE; wp < POS_INT_RANGE; ++wp )
	       	       {
			//cout << "sum" << endl;  
	       		val += 1.0/BETA/BETA/4.0/PI/PI	 * //Two internal loops
			    bubble_ph[W][w][n_in][mp][s4][s2_in][s3][s1_in](K) * 
			    state_vec.vertx_d( W, w, wp, K, mp, m, s3, s4, s3p, s4p ) *  
	       		    bubble_ph[W][wp][m][n_out][s2_out][s4p][s1_out][s3p](K) *
			    weight_vec_2d[w][wp];
	       	       }
   for( int w = -POS_INT_RANGE; w < POS_INT_RANGE; ++w )
   {
      val -= 1/BETA * // One internal loop
	 bubble_ph[W][w][n_in][n_out][s2_out][s2_in][s1_out][s1_in](K) *
	 weight_vec[w];
   }

   return val; 
}

dcomplex rhs_t::eval_diag_susc_m( const idx_susc_t& idx, const state_t& state_vec, const gf_bubble_mat_t& bubble_ph)
{
   dcomplex val( 0.0, 0.0 );

   // Introduce help variables
   int W = idx( ISUSC::W );

   int K = idx( ISUSC::K );
   
   int n_in = idx( ISUSC::n_in);
   int n_out = idx( ISUSC::n_out);

   int s1_in = idx( ISUSC::s1_in );
   int s2_in = idx( ISUSC::s2_in );
   int s1_out = idx( ISUSC::s1_out );
   int s2_out = idx( ISUSC::s2_out );
   
   for( int s3 = 0; s3 < QN_COUNT; ++s3 )
      for( int s3p = 0; s3p < QN_COUNT; ++s3p )
	 for( int s4 = 0; s4 < QN_COUNT; ++s4 )
	    for( int s4p = 0; s4p < QN_COUNT; ++s4p )
	       for( int m = 0; m < FFACTOR_COUNT; ++m )
		  for( int mp = 0; mp < FFACTOR_COUNT; ++mp )
		     for( int w = -POS_INT_RANGE; w < POS_INT_RANGE; ++w )
		       for( int wp = -POS_INT_RANGE; wp < POS_INT_RANGE; ++wp )
	       	       {
	       		val += 1.0/BETA/BETA/4.0/PI/PI *
			    bubble_ph[W][w][n_in][mp][s4][s2_in][s3][s1_in](K) * 
			    state_vec.vertx_m( W, w, wp, K, mp, m, s3, s4, s3p, s4p ) *  
	       		    bubble_ph[W][wp][m][n_out][s2_out][s4p][s1_out][s3p](K) *
			    weight_vec_2d[w][wp];
	       	       }
   for( int w = -POS_INT_RANGE; w < POS_INT_RANGE; ++w )
   {
      val -= 1/BETA * // One internal loop
	 bubble_ph[W][w][n_in][n_out][s2_out][s2_in][s1_out][s1_in](K) *
	 weight_vec[w];
   }

   return val; 
}


dcomplex rhs_t::eval_diag_tri_sc( const idx_tri_t& idx, const state_t& state_vec, const gf_bubble_mat_t& bubble_pp)
{
   dcomplex val( 0.0, 0.0 );
   
   // Introduce help variables
   int W    = idx( ITRI::W );
   int w_in = idx( ITRI::w );

   int K = idx( ITRI::K );

   int n_in  = idx( ITRI::n_in);
   int n_out = idx( ITRI::n_out);

   int s1_in = idx(  ITRI::s1_in );
   int s2_in = idx(  ITRI::s2_in );
   int s1_out = idx( ITRI::s1_out );
   int s2_out = idx( ITRI::s2_out );
   
   for( int s3 = 0; s3 < QN_COUNT; ++s3 )
      for( int s3p = 0; s3p < QN_COUNT; ++s3p )
	 for( int s4 = 0; s4 < QN_COUNT; ++s4 )
	    for( int s4p = 0; s4p < QN_COUNT; ++s4p )
	       for( int m = 0; m < FFACTOR_COUNT; ++m )
		  for( int mp = 0; mp < FFACTOR_COUNT; ++mp ){
		     for( int w = -POS_INT_RANGE- abs(W/2) - (W+100000) % 2; w < POS_INT_RANGE + abs(W/2); ++w )
	       	       {
	       		val += 1.0/BETA/4.0/PI/PI * // Two pairs of equivalent lines -> Factor 1/4
			    tri_bare(n_in, m, s1_in, s2_in, s3p, s4p) *
			    bubble_pp[W][w][m][mp][s4][s4p][s3][s3p](K) * 
			    state_vec.vertx_pp( W, w, w_in, K, mp, n_out, s3, s4, s1_out, s2_out );
		       } 
			int w = POS_INT_RANGE + abs(W/2);
	      		val += 1.0/4.0/PI/PI * double(s3 == s3p) * double(s4 == s4p) * double(m == mp) * asympt_GG_pp( W ) *
 				tri_bare(n_in, m, s1_in, s2_in, s3p, s4p) *
				state_vec.vertx_pp( W, w, w_in, K, mp, n_out, s3, s4, s1_out, s2_out );		
			    //weight_vec[w];
	       	       }
   
   return val + double(n_in == n_out); 
}

dcomplex rhs_t::eval_diag_tri_d( const idx_tri_t& idx, const state_t& state_vec, const gf_bubble_mat_t& bubble_ph)
{
   dcomplex val( 0.0, 0.0 );
   
   // Introduce help variables
   int W    = idx( ITRI::W );
   int w_in = idx( ITRI::w );

   int K = idx( ITRI::K );

   int n_in  = idx( ITRI::n_in);
   int n_out = idx( ITRI::n_out);

   int s1_in = idx(  ITRI::s1_in );
   int s2_in = idx(  ITRI::s2_in );
   int s1_out = idx( ITRI::s1_out );
   int s2_out = idx( ITRI::s2_out );
   
   for( int s3 = 0; s3 < QN_COUNT; ++s3 )
      for( int s3p = 0; s3p < QN_COUNT; ++s3p )
	 for( int s4 = 0; s4 < QN_COUNT; ++s4 )
	    for( int s4p = 0; s4p < QN_COUNT; ++s4p )
	       for( int m = 0; m < FFACTOR_COUNT; ++m )
		  for( int mp = 0; mp < FFACTOR_COUNT; ++mp ){
		     for( int w = -POS_INT_RANGE - abs(W/2) - (W + 100000) % 2; w < POS_INT_RANGE+ abs(W/2); ++w )
	       	       {
	       		val -= 1.0/BETA/4.0/PI/PI * // Two pairs of equivalent lines -> Factor 1/4
			    tri_bare(n_in, m, s1_in, s2_in, s3p, s4p) *
			    bubble_ph[W][w][m][mp][s4][s4p][s3][s3p](K) * 
			    state_vec.vertx_d( W, w, w_in, K, mp, n_out, s3, s4, s1_out, s2_out ); 
		       }
			int w = POS_INT_RANGE + abs(W/2);
	      		val -= 1.0/4.0/PI/PI * double(s3 == s3p) * double(s4 == s4p) * double(m == mp) * asympt_GG_ph( W ) *
 				tri_bare(n_in, m, s1_in, s2_in, s3p, s4p) *
				state_vec.vertx_d( W, w, w_in, K, mp, n_out, s3, s4, s1_out, s2_out );		
			    //weight_vec[w];
	       	       }


   return val + double(n_in == n_out); 
}

dcomplex rhs_t::eval_diag_tri_m( const idx_tri_t& idx, const state_t& state_vec, const gf_bubble_mat_t& bubble_ph)
{
   dcomplex val( 0.0, 0.0 );
   
   // Introduce help variables
   int W    = idx( ITRI::W );
   int w_in = idx( ITRI::w );

   int K = idx( ITRI::K );

   int n_in  = idx( ITRI::n_in);
   int n_out = idx( ITRI::n_out);

   int s1_in = idx(  ITRI::s1_in );
   int s2_in = idx(  ITRI::s2_in );
   int s1_out = idx( ITRI::s1_out );
   int s2_out = idx( ITRI::s2_out );
   
   for( int s3 = 0; s3 < QN_COUNT; ++s3 )
      for( int s3p = 0; s3p < QN_COUNT; ++s3p )
	 for( int s4 = 0; s4 < QN_COUNT; ++s4 )
	    for( int s4p = 0; s4p < QN_COUNT; ++s4p )
	       for( int m = 0; m < FFACTOR_COUNT; ++m )
		  for( int mp = 0; mp < FFACTOR_COUNT; ++mp ){
		     for( int w = -POS_INT_RANGE-abs(W/2)-(W + 100000) % 2; w < POS_INT_RANGE + abs(W/2); ++w )
	       	       {
	       		val -= 1.0/BETA/4.0/PI/PI * // Two pairs of equivalent lines -> Factor 1/4
			    tri_bare(n_in, m, s1_in, s2_in, s3p, s4p) *
			    bubble_ph[W][w][m][mp][s4][s4p][s3][s3p](K) * 
			    state_vec.vertx_m( W, w, w_in, K, mp, n_out, s3, s4, s1_out, s2_out );
		       }
		       
		        int w = POS_INT_RANGE + abs(W/2);
	      		val -= 1.0/4.0/PI/PI * double(s3 == s3p) * double(s4 == s4p) * double(m == mp) * asympt_GG_ph( W ) *
 				tri_bare(n_in, m, s1_in, s2_in, s3p, s4p) *
				state_vec.vertx_m( W, w, w_in, K, mp, n_out, s3, s4, s1_out, s2_out );		
			    //weight_vec[w];
	       	       }


   return val + double(n_in == n_out); 
}
