
/************************************************************************************************//**
 *  		
 * 	file: 		rhs.cpp
 * 	contents:   	See rhs.h
 * 
 ****************************************************************************************************/


#include <rhs.h>
#include <frg.h>
#include <params.h>
//#include <observables.h>
#include <cmath>
#include <translate.h>
#include <mymath.h>
#include <projection.h>
#include <complex>

using namespace std; 
using boost::bind;

// --- Create static objects: symmetry groups, weight vector
// ALL SYMM

//TODO: include mirro_y_symmetry

symmetry_grp_t<dcomplex,8> rhs_t::symm_grp_chi_tripl( gf_suscept_t(), { hmirror_suscept_pp, cconj_suscept_pp, timerev_suscept_pp,rot_pi2_x_suscept,rot_pi2_y_suscept,rot_pi2_z_suscept,rot_2pi3_111_suscept,rot_2pi3_1m11_suscept,rot_2pi3_m111_suscept,rot_2pi3_m1m11_suscept,rot_pi_110_suscept,rot_pi_1m10_suscept,rot_pi_011_suscept,rot_pi_01m1_suscept,rot_pi_101_suscept,rot_pi_10m1_suscept,mirror_y_suscept  });
symmetry_grp_t<dcomplex,8> rhs_t::symm_grp_chi_singl( gf_suscept_t(), { hmirror_suscept_pp, cconj_suscept_pp, timerev_suscept_pp,rot_pi2_x_suscept,rot_pi2_y_suscept,rot_pi2_z_suscept,rot_2pi3_111_suscept,rot_2pi3_1m11_suscept,rot_2pi3_m111_suscept,rot_2pi3_m1m11_suscept,rot_pi_110_suscept,rot_pi_1m10_suscept,rot_pi_011_suscept,rot_pi_01m1_suscept,rot_pi_101_suscept,rot_pi_10m1_suscept,mirror_y_suscept });
symmetry_grp_t<dcomplex,8> rhs_t::symm_grp_chi_dens(  gf_suscept_t(),{ hmirror_suscept_ph, cconj_suscept_ph, timerev_suscept_ph,rot_pi2_x_suscept,rot_pi2_y_suscept,rot_pi2_z_suscept,rot_2pi3_111_suscept,rot_2pi3_1m11_suscept,rot_2pi3_m111_suscept,rot_2pi3_m1m11_suscept,rot_pi_110_suscept,rot_pi_1m10_suscept,rot_pi_011_suscept,rot_pi_01m1_suscept,rot_pi_101_suscept,rot_pi_10m1_suscept,mirror_y_suscept }); 
symmetry_grp_t<dcomplex,8> rhs_t::symm_grp_chi_mag(   gf_suscept_t(),{ hmirror_suscept_xph, cconj_suscept_xph, timerev_suscept_xph,rot_pi2_x_suscept,rot_pi2_y_suscept,rot_pi2_z_suscept,rot_2pi3_111_suscept,rot_2pi3_1m11_suscept,rot_2pi3_m111_suscept,rot_2pi3_m1m11_suscept,rot_pi_110_suscept,rot_pi_1m10_suscept,rot_pi_011_suscept,rot_pi_01m1_suscept,rot_pi_101_suscept,rot_pi_10m1_suscept,mirror_y_suscept }); 

symmetry_grp_t<MatReal,3> rhs_t::symm_grp_Gvec_real(gf_1p_mat_real_t(),{});

symmetry_grp_t<MatPatch,8> rhs_t::symm_grp_bubble_pp(gf_bubble_mat_t(),{swap_bubble, cconj_bubble, hmirror_bubble_pp});
symmetry_grp_t<MatPatch,8> rhs_t::symm_grp_bubble_ph(gf_bubble_mat_t(),{swap_bubble, cconj_bubble, hmirror_timerev_bubble_ph});



//-----------------NO SYMMETRIES ----------------------------------

//symmetry_grp_t<dcomplex,4> rhs_t::symm_grp_sig( gf_1p_t(), {});
//symmetry_grp_t<dcomplex,10> rhs_t::symm_grp_phi_pp( gf_phi_t(), {}); 
//symmetry_grp_t<dcomplex,10> rhs_t::symm_grp_phi_ph( gf_phi_t(), {}); 
//symmetry_grp_t<dcomplex,10> rhs_t::symm_grp_phi_xph( gf_phi_t(),{});
//
//symmetry_grp_t<dcomplex,8> rhs_t::symm_grp_P_pp( gf_P_t(), {});
//symmetry_grp_t<dcomplex,8> rhs_t::symm_grp_P_ph( gf_P_t(), {});   
//symmetry_grp_t<dcomplex,8> rhs_t::symm_grp_P_xph( gf_P_t(),{}); 
//
//symmetry_grp_t<dcomplex,6> rhs_t::symm_grp_chi_pp( gf_chi_t(), {});
//symmetry_grp_t<dcomplex,6> rhs_t::symm_grp_chi_ph( gf_chi_t(), {});
//symmetry_grp_t<dcomplex,6> rhs_t::symm_grp_chi_xph( gf_chi_t(),{}); 
//
//symmetry_grp_t<MatReal,3> rhs_t::symm_grp_Gvec_real(gf_1p_mat_real_t(),{});
//
//symmetry_grp_t<MatPatch,7> rhs_t::symm_grp_GG_pp(gf_GG_mat_t(),{});
//symmetry_grp_t<MatPatch,7> rhs_t::symm_grp_GG_ph(gf_GG_mat_t(),{});
//
//symmetry_grp_t<MatPatch,8> rhs_t::symm_grp_bubble_pp(gf_bubble_mat_t(),{});
//symmetry_grp_t<MatPatch,8> rhs_t::symm_grp_bubble_ph(gf_bubble_mat_t(),{});


gf_weight_vec_t rhs_t::weight_vec( generate_weights( POS_INT_RANGE-TAIL_LENGTH, TAIL_LENGTH, FIT_ORDER ) ); 
gf<double, 2> rhs_t::weight_vec_2d( generate_2d_weights( POS_INT_RANGE-TAIL_LENGTH, TAIL_LENGTH, FIT_ORDER ) ); 

void rhs_t::operator() ( const state_t& state_vec, state_t& dfdl)
{

   // Precalculate Green function and single-scale propagator on frequency grid

   gf_1p_mat_t Gvec( POS_1P_RANGE, FFT_DIM*FFT_DIM*FFT_DIM ); 
   
   Gvec.init( bind( &state_t::GMat, boost::cref(state_vec), _1) ); // Initialize big Green function vector 

   // Precalculate Green's functions in real space by FFT
   gf_1p_mat_real_t Gvec_real;
   gf_bubble_mat_t bubble_pp;
   gf_bubble_mat_t bubble_ph;
   
   cout << " ... bubble PT1 " << endl;
   fftw_complex *input, *output;	// input output functions for the FFTW
   fftw_plan p, p_b; 			// plan for FFTW
   
   input = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * FFT_DIM * FFT_DIM * FFT_DIM);
   output = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * FFT_DIM * FFT_DIM * FFT_DIM);

   p = fftw_plan_dft_3d(FFT_DIM, FFT_DIM, FFT_DIM, input, output, FFTW_BACKWARD, FFTW_ESTIMATE); //plan for FFTW, same for Gvec and Svec (BACKWARD), different for bubble (FORWARD)
   
   
   symm_grp_Gvec_real.init( Gvec_real, [&Gvec, p]( const idx_1p_mat_real_t& idx ){ return eval_Gvec_real( idx, Gvec, p); } );
   
   
   cout << " ... bubble PT2 " << endl;
   
   p_b = fftw_plan_dft_3d(FFT_DIM, FFT_DIM, FFT_DIM, input, output, FFTW_FORWARD, FFTW_ESTIMATE); //plan for FFTW, same for Gvec and Svec (BACKWARD), different for bubble (FORWARD)
   
   symm_grp_bubble_pp.init( bubble_pp, [&Gvec_real, p_b]( const idx_bubble_mat_t& idx ){ return eval_diag_bubble_pp( idx, Gvec_real, p_b ); } ); 
   symm_grp_bubble_ph.init( bubble_ph, [&Gvec_real, p_b]( const idx_bubble_mat_t& idx ){ return eval_diag_bubble_ph( idx, Gvec_real, p_b ); } ); 
   // Calculate right hand side for vertex and self-energy

   fftw_destroy_plan(p);
   fftw_destroy_plan(p_b);
   fftw_free(input); fftw_free(output);
   cout << " Rhs calculation " << endl;

   /*********************  Open MP parallelized RHS calculation  ********************/

   // --- Conventional flow
   cout << " ... susceptibilities " << endl;
   cout << " ... Triplet " << endl;
   symm_grp_chi_tripl.init( dfdl.gf_suscept_trip(), [&state_vec, &bubble_pp]( const idx_suscept_t& idx ){ return eval_diag_suscept_trip( idx, state_vec, bubble_pp ); } ); 
   cout << " ... Singlet " << endl;
   symm_grp_chi_singl.init( dfdl.gf_suscept_s(),    [&state_vec, &bubble_pp]( const idx_suscept_t& idx ){ return eval_diag_suscept_s( idx, state_vec, bubble_pp ); } ); 
   cout << " ... Density " << endl;
   symm_grp_chi_dens.init( dfdl.gf_suscept_d(),     [&state_vec, &bubble_ph]( const idx_suscept_t& idx ){ return eval_diag_suscept_d( idx, state_vec, bubble_ph ); } ); 
   cout << " ... Magnetic " << endl;
   symm_grp_chi_mag.init( dfdl.gf_suscept_m(),      [&state_vec, &bubble_ph]( const idx_suscept_t& idx ){ return eval_diag_suscept_m( idx, state_vec, bubble_ph ); } ); 

}


MatReal rhs_t::eval_Gvec_real( const idx_1p_mat_real_t& idx, const gf_1p_mat_t& Gvec, fftw_plan p_g )
{
   int w = idx(0);
   int s1 = idx(1);
   int s2 = idx(2);  

   MatReal val_g = MatReal::Zero();
   fftw_complex *input_g, *output_g;
   input_g = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * FFT_DIM * FFT_DIM * FFT_DIM);
   output_g = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * FFT_DIM * FFT_DIM * FFT_DIM);
   
   for(int p=0;p<FFT_DIM*FFT_DIM*FFT_DIM; ++p){
      input_g[p][0] = Gvec[w][p](s1,s2).real(); 
      input_g[p][1] = Gvec[w][p](s1,s2).imag(); 
   }
   
   // plan is created outside and passed to function because only fftw_execute_dft is thread save
   fftw_execute_dft(p_g,input_g,output_g); 
   
   for(int R=0;R<FFT_DIM*FFT_DIM*FFT_DIM; ++R){
      val_g(R) =  1.0/FFT_DIM/FFT_DIM/FFT_DIM * (output_g[Matidxshift(get_px_idx(R),get_py_idx(R),get_pz_idx(R))][0] + I * output_g[Matidxshift(get_px_idx(R),get_py_idx(R),get_pz_idx(R))][1]);
   }

   fftw_free(input_g); fftw_free(output_g);
   return val_g; 
}




MatPatch rhs_t::eval_diag_bubble_pp( const idx_bubble_mat_t& idx, const gf_1p_mat_real_t& Gvec_real, fftw_plan p_b)
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
   input_b = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * FFT_DIM * FFT_DIM * FFT_DIM);
   output_b = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * FFT_DIM * FFT_DIM * FFT_DIM);
   
   for(int r=0; r< REAL_GRID; ++r){
      if(abs(F_factors_real::fft_fxf_weight[r][m][n])> CHOP_ERR){// 1E-16 is error on the weights exponents, cosines,...
      int rx = std::get<0>(R_Grid()[r]);
      int ry = std::get<1>(R_Grid()[r]);
      int rz = std::get<2>(R_Grid()[r]);
      for(int Rx=0;Rx<FFT_DIM; ++Rx)
         for(int Ry=0;Ry<FFT_DIM; ++Ry)
            for(int Rz=0;Rz<FFT_DIM; ++Rz){
	       dcomplex temp = Gvec_real[ w + div2_ceil( W ) ][s4][s4p](Matidx(Rx,Ry,Rz)) * Gvec_real[ div2_floor( W ) - w - 1 ][s3][s3p](Matidx(Rx+rx,Ry+ry,Rz+rz));
		  input_b[Matidxshift(Rx,Ry,Rz)][0] = temp.real();
	          input_b[Matidxshift(Rx,Ry,Rz)][1] = temp.imag();
	    }
      
      // plan is created outside and passed to function because only fftw_execute_dft is thread save
      fftw_execute_dft(p_b, input_b, output_b);

	 for (int K=0; K < PATCH_COUNT; ++K){
	    int Kx_idx = std::get<0>(K_Grid::get_indices(K));  
	    int Ky_idx = std::get<1>(K_Grid::get_indices(K));
	    int Kz_idx = std::get<2>(K_Grid::get_indices(K));  
	    int Kx_shift = (Kx_idx*FFT_DIM/2/MAX_KPOS + 10*FFT_DIM) % FFT_DIM;
            int Ky_shift = (Ky_idx*FFT_DIM/2/MAX_KPOS + 10*FFT_DIM) % FFT_DIM;
            int Kz_shift = (Kz_idx*FFT_DIM/2/MAX_KPOS + 10*FFT_DIM) % FFT_DIM;
	    dcomplex weight = F_factors_real::fft_fxf_weight[r][m][n];
	    weight *= (cos_k(Kx_idx*rx + Ky_idx*ry+Kz_idx*rz) - I * sin_k(Kx_idx*rx + Ky_idx*ry+Kz_idx*rz));
		 val_patch(K) += (output_b[Matidx(Kx_shift,Ky_shift,Kz_shift)][0]+ I * output_b[Matidx(Kx_shift,Ky_shift,Kz_shift)][1]) * weight; 
	  }
      }
   }
  
  
  fftw_free(input_b); fftw_free(output_b);
  return val_patch;

}

MatPatch rhs_t::eval_diag_bubble_ph( const idx_bubble_mat_t& idx, const gf_1p_mat_real_t& Gvec_real, fftw_plan p_b )
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
   input_b = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * FFT_DIM * FFT_DIM * FFT_DIM);
   output_b = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * FFT_DIM * FFT_DIM * FFT_DIM);
   
   for(int r=0; r< REAL_GRID; ++r){
      if(abs(F_factors_real::fft_fxf_weight[r][m][n])>  CHOP_ERR){// 1E-16 is error on the weights exponents, cosines,...
      int rx = std::get<0>(R_Grid()[r]);
      int ry = std::get<1>(R_Grid()[r]);
      int rz = std::get<2>(R_Grid()[r]);
      for(int Rx=0;Rx<FFT_DIM; ++Rx)
         for(int Ry=0;Ry<FFT_DIM; ++Ry)
            for(int Rz=0;Rz<FFT_DIM; ++Rz){
	       dcomplex temp = Gvec_real[ w - div2_floor( W ) ][s4][s4p](Matidx(FFT_DIM-Rx,FFT_DIM-Ry,FFT_DIM-Rz)) * Gvec_real[ div2_ceil( W ) + w ][s3][s3p](Matidx(Rx-rx,Ry-ry,Rz-rz));
	       input_b[Matidxshift(Rx,Ry,Rz)][0] = temp.real();
	       input_b[Matidxshift(Rx,Ry,Rz)][1] = temp.imag();
	    }
      
      // plan is created outside and passed to function because only fftw_execute_dft is thread save
      fftw_execute_dft(p_b,input_b,output_b);
      
      for (int K = 0; K < PATCH_COUNT; ++K){
         int Kx_idx = std::get<0>(K_Grid::get_indices(K));  
         int Ky_idx = std::get<1>(K_Grid::get_indices(K)); 
         int Kz_idx = std::get<2>(K_Grid::get_indices(K));  
         int Kx_shift = (Kx_idx*FFT_DIM/2/MAX_KPOS + 10*FFT_DIM) % FFT_DIM;
         int Ky_shift = (Ky_idx*FFT_DIM/2/MAX_KPOS + 10*FFT_DIM) % FFT_DIM;
         int Kz_shift = (Kz_idx*FFT_DIM/2/MAX_KPOS + 10*FFT_DIM) % FFT_DIM;
	    dcomplex weight = F_factors_real::fft_fxf_weight[r][m][n];
	    weight *= (cos_k(Kx_idx*rx + Ky_idx*ry+Kz_idx*rz) + I * sin_k(Kx_idx*rx + Ky_idx*ry+Kz_idx*rz));
	    val_patch(K) += (output_b[Matidx(Kx_shift,Ky_shift,Kz_shift)][0]+ I * output_b[Matidx(Kx_shift,Ky_shift,Kz_shift)][1]) * weight;
	  }
      }
   }

  fftw_free(input_b); fftw_free(output_b);
  return val_patch;
}

dcomplex rhs_t::eval_diag_suscept_trip( const idx_suscept_t& idx, const state_t& state_vec, const gf_bubble_mat_t& bubble_pp)
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
	       		val += 0.25/BETA/BETA * // Two pairs of equivalent lines -> Factor 1/4
			    bubble_pp[W][w][n_in][mp][s4][s2_in][s3][s1_in](K) * 
			    state_vec.vertx_t( W, w, wp, K, mp, m, s3, s4, s3p, s4p ) *  
	       		    bubble_pp[W][wp][m][n_out][s2_out][s4p][s1_out][s3p](K) *
			    weight_vec_2d[w][wp];
	       	       }
   for( int w = -POS_INT_RANGE; w < POS_INT_RANGE; ++w )
   {
      val += 0.5/BETA *	//Two equivalent lines -> factor 1/2
	 bubble_pp[W][w][n_in][n_out][s2_out][s2_in][s1_out][s1_in](K) *
	 weight_vec[w];
   }

   return val; 
}

dcomplex rhs_t::eval_diag_suscept_s( const idx_suscept_t& idx, const state_t& state_vec, const gf_bubble_mat_t& bubble_pp)
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
	       		val += 0.25/BETA/BETA * // Two pairs of equivalent lines -> Factor 1/4
			    bubble_pp[W][w][n_in][mp][s4][s2_in][s3][s1_in](K) * 
			    state_vec.vertx_s( W, w, wp, K, mp, m, s3, s4, s3p, s4p ) *  
	       		    bubble_pp[W][wp][m][n_out][s2_out][s4p][s1_out][s3p](K) *
			    weight_vec_2d[w][wp];
	       	       }
   for( int w = -POS_INT_RANGE; w < POS_INT_RANGE; ++w )
   {
      val += 0.5/BETA * // Two equivalent lines -> factor 1/2
	 bubble_pp[W][w][n_in][n_out][s2_out][s2_in][s1_out][s1_in](K) *
	 weight_vec[w];
   }

   return val; 
}

dcomplex rhs_t::eval_diag_suscept_d( const idx_suscept_t& idx, const state_t& state_vec, const gf_bubble_mat_t& bubble_ph)
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
	       		val += 1.0/BETA/BETA * //Two internal loops
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

dcomplex rhs_t::eval_diag_suscept_m( const idx_suscept_t& idx, const state_t& state_vec, const gf_bubble_mat_t& bubble_ph)
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
	       		val += 1.0/BETA/BETA *
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
