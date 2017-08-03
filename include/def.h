
/******************************************************************************************//** @file
 *  		
 * 	file: 		def.h
 * 	contents:  	Definition of used correlation function containers ( wrapper around gf container )
 * 
 ****************************************************************************************************/

#pragma once

#define EIGEN_STACK_ALLOCATION_LIMIT 0
#include <complex>
#include <Eigen/Core>
#include <vector>
#include <boost/multi_array.hpp>

#include <const.h>
#include <gf.h>

using dcomplex = std::complex<double>;						///< Complex double type
using MatQN = Eigen::Matrix<dcomplex, QN_COUNT, QN_COUNT, Eigen::RowMajor>;	///< Complex matrix representing the discrete quantum number structure
using MatQNQN = Eigen::Matrix<dcomplex, QN_COUNT*QN_COUNT, QN_COUNT*QN_COUNT, Eigen::RowMajor>;	///< Complex matrix representing the discrete quantum number structure of two-particle function

using MatReal = Eigen::Matrix<dcomplex, FFT_DIM * FFT_DIM * FFT_DIM, 1>;	///< Complex matrix representing the real grid structure
using MatPatch = Eigen::Matrix<dcomplex, PATCH_COUNT, 1>;	///< Complex matrix representing the real grid structure

#define INSERT_COPY_AND_ASSIGN(X) 					\
X( const X & gf_obj ):    						\
   base_t( gf_obj )							\
{}       								\
X( X && gf_obj ):							\
   base_t( std::move(gf_obj) )						\
{}      								\
X & operator=( const X & gf_obj )					\
{									\
   base_t::operator=( gf_obj ); 					\
   return *this; 							\
} 									\
X & operator=( X && gf_obj )						\
{									\
   base_t::operator=( std::move( gf_obj) ); 				\
   return *this; 							\
} 


// Container and index types
class gf_1p_mat_t : public gf< MatQN, 2 > 			///< Matrix-valued container type for one-particle correlation functions
{
   public:
      using base_t = gf< MatQN, 2 >; 

      gf_1p_mat_t( int pos_freq_count_= POS_FFREQ_COUNT_SIG, int patch_count_ = PATCH_COUNT ):
	 base_t( boost::extents[ffreq(pos_freq_count_)][patch_count_] ) 
   {}
      INSERT_COPY_AND_ASSIGN(gf_1p_mat_t)
}; 
using idx_1p_mat_t = gf_1p_mat_t::idx_t;   

class gf_1p_mat_real_t : public gf< MatReal, 3 > 			///< Matrix-valued container type for one-particle correlation functions
{
   public:
      using base_t = gf< MatReal, 3 >; 

      gf_1p_mat_real_t( int pos_freq_count_= POS_1P_RANGE ):
	 base_t( boost::extents[ffreq(pos_freq_count_)][QN_COUNT][QN_COUNT] ) 
   {}
      INSERT_COPY_AND_ASSIGN(gf_1p_mat_real_t)
}; 
using idx_1p_mat_real_t = gf_1p_mat_real_t::idx_t;   

class gf_1p_real_t : public gf< dcomplex, 4 > 			///< Matrix-valued container type for one-particle correlation functions
{
   public:
      using base_t = gf< dcomplex, 4 >; 

      gf_1p_real_t( int pos_freq_count_= POS_1P_RANGE ):
	 base_t( boost::extents[ffreq(pos_freq_count_)][FFT_DIM*FFT_DIM*FFT_DIM][QN_COUNT][QN_COUNT] ) 
   {}
      INSERT_COPY_AND_ASSIGN(gf_1p_real_t)
}; 
using idx_1p_real_t = gf_1p_real_t::idx_t;   

class gf_GG_mat_t : public gf< MatPatch, 7 > 			///< Matrix-valued container type for one-particle correlation functions
{
   public:
      using base_t = gf< MatPatch, 7 >; 

      gf_GG_mat_t( int pos_bfreq_count_= POS_BFREQ_COUNT_CHI, int pos_freq_count_= POS_INT_RANGE+POS_BFREQ_COUNT_CHI/2):
	 base_t( boost::extents[bfreq(pos_bfreq_count_)][ffreq(pos_freq_count_)][REAL_GRID][QN_COUNT][QN_COUNT][QN_COUNT][QN_COUNT] ) 
   {}
      INSERT_COPY_AND_ASSIGN(gf_GG_mat_t)
}; 
using idx_GG_mat_t = gf_GG_mat_t::idx_t;   

class gf_GG_t : public gf< dcomplex, 8 > 			///< Matrix-valued container type for one-particle correlation functions
{
   public:
      using base_t = gf< dcomplex, 8 >; 

      gf_GG_t( int pos_bfreq_count_= POS_BFREQ_COUNT_CHI, int pos_freq_count_= POS_INT_RANGE+POS_BFREQ_COUNT_CHI/2):
	 base_t( boost::extents[bfreq(pos_bfreq_count_)][ffreq(pos_freq_count_)][PATCH_COUNT][REAL_GRID][QN_COUNT][QN_COUNT][QN_COUNT][QN_COUNT] ) 
   {}
      INSERT_COPY_AND_ASSIGN(gf_GG_t)
}; 
using idx_GG_t = gf_GG_t::idx_t;   

class gf_bubble_mat_t : public gf< MatPatch, 8 > 			///< Matrix-valued container type for one-particle correlation functions
{
   public:
      using base_t = gf< MatPatch, 8 >; 

      gf_bubble_mat_t( int pos_bfreq_count_= POS_BFREQ_COUNT_SUSCEPT, int pos_freq_count_= POS_INT_RANGE+POS_BFREQ_COUNT_SUSCEPT/2):
	 base_t( boost::extents[bfreq(pos_bfreq_count_)][ffreq(pos_freq_count_)][FFACTOR_COUNT][FFACTOR_COUNT][QN_COUNT][QN_COUNT][QN_COUNT][QN_COUNT] ) 
   {}
      INSERT_COPY_AND_ASSIGN(gf_bubble_mat_t)
}; 
using idx_bubble_mat_t = gf_bubble_mat_t::idx_t;   


enum class I1P{ w, k, s_in, s_out }; 
class gf_1p_t : public  gf< dcomplex, 4 > 		///< Container type for one-particle correlation functions, holds dcomplex
{
   public:
      using base_t = gf< dcomplex, 4 >; 

      gf_1p_t( int pos_freq_count_ = POS_FFREQ_COUNT_SIG ):
	 base_t( boost::extents[ffreq(pos_freq_count_)][PATCH_COUNT][QN_COUNT][QN_COUNT] )
   {}
      INSERT_COPY_AND_ASSIGN(gf_1p_t)
}; 
using idx_1p_t = gf_1p_t::idx_t; 

enum class I2P{ w1_in, w2_in, w1_out, k1_in, k2_in, k1_out, s1_in, s2_in, s1_out, s2_out }; 
class gf_2p_t : public gf< dcomplex, 10 > 		///< Container type for two-particle correlation functions, holds dcomplex
{
   public:
      using base_t = gf< dcomplex, 10 >; 

      gf_2p_t( int pos_freq_count_ = POS_FFREQ_COUNT_PHI ):
	 base_t( boost::extents[ffreq(pos_freq_count_)][ffreq(pos_freq_count_)][ffreq(pos_freq_count_)]
	       [PATCH_COUNT][PATCH_COUNT][PATCH_COUNT]
	       [QN_COUNT][QN_COUNT][QN_COUNT][QN_COUNT] )
   {}
      INSERT_COPY_AND_ASSIGN(gf_2p_t)
}; 
using idx_2p_t = gf_2p_t::idx_t; 

enum class IPHI{ W, w_in, w_out, K, n_in, n_out, s1_in, s2_in, s1_out, s2_out }; 
class gf_phi_t : public gf< dcomplex, 10 > 		///< Container type for two-particle correlation functions, holds ind_cpl_t
{
   public:
      using base_t = gf< dcomplex, 10 >;

      gf_phi_t( int pos_bfreq_count_ = POS_BFREQ_COUNT_PHI, int pos_freq_count_ = POS_FFREQ_COUNT_PHI ):
	 base_t( boost::extents[bfreq(pos_bfreq_count_)][ffreq(pos_freq_count_)][ffreq(pos_freq_count_)]
	       [PATCH_COUNT][FFACTOR_COUNT][FFACTOR_COUNT]
	       [QN_COUNT][QN_COUNT][QN_COUNT][QN_COUNT] )
   {}
      INSERT_COPY_AND_ASSIGN(gf_phi_t)
}; 
using idx_phi_t = gf_phi_t::idx_t; 

enum class IP{ W, w, K, n, s1_in, s2_in, s1_out, s2_out }; 
class gf_P_t : public gf< dcomplex, 8 > 		///< Container type for two-particle correlation functions, holds dcomplex
{
   public:
      using base_t = gf< dcomplex, 8 >; 

      gf_P_t( int pos_bfreq_count_ = POS_BFREQ_COUNT_P,int pos_freq_count_ = POS_FFREQ_COUNT_P ):
	 base_t( boost::extents[bfreq(pos_bfreq_count_)][ffreq(pos_freq_count_)]
	       [PATCH_COUNT][FFACTOR_COUNT]
	       [QN_COUNT][QN_COUNT][QN_COUNT][QN_COUNT] )
   {}
      INSERT_COPY_AND_ASSIGN(gf_P_t)
}; 
using idx_P_t = gf_P_t::idx_t; 

enum class ICHI{ W, K, s1_in, s2_in, s1_out, s2_out }; 
class gf_chi_t : public gf< dcomplex, 6 > 		///< Container type for two-particle correlation functions
{
   public:
      using base_t = gf< dcomplex, 6 >; 

      gf_chi_t( int pos_bfreq_count_ = POS_BFREQ_COUNT_CHI ):
	 base_t( boost::extents[bfreq(pos_bfreq_count_)]
	       [PATCH_COUNT][QN_COUNT][QN_COUNT][QN_COUNT][QN_COUNT] )
   {}
      INSERT_COPY_AND_ASSIGN(gf_chi_t)
}; 
using idx_chi_t = gf_chi_t::idx_t; 


enum class ISUSC{ W, K, n_in, n_out, s1_in, s2_in, s1_out, s2_out }; 
class gf_suscept_t : public gf< dcomplex, 8 > 		///< Container type for two-particle correlation functions
{
   public:
      using base_t = gf< dcomplex, 8 >; 

      gf_suscept_t( int pos_bfreq_count = POS_BFREQ_COUNT_SUSCEPT ):
	 base_t( boost::extents[bfreq(pos_bfreq_count)]
	       [PATCH_COUNT][FFACTOR_COUNT][FFACTOR_COUNT]
	       [QN_COUNT][QN_COUNT][QN_COUNT][QN_COUNT] )
   {}
      INSERT_COPY_AND_ASSIGN(gf_suscept_t)
}; 
using idx_suscept_t = gf_suscept_t::idx_t; 

// Class for the bare vertex in the form factors + band basis -> only needed for input to the frg flow 
enum class VERT_BARE_FF{n_in, n_out, s1_in, s2_in, s1_out, s2_out }; 
class gf_vert_bare_ff_t : public gf< dcomplex, 6 > 		///< Container type for two-particle correlation functions
{
   public:
      using base_t = gf< dcomplex, 6 >; 

      gf_vert_bare_ff_t():
	 base_t( boost::extents[FFACTOR_COUNT][FFACTOR_COUNT]
	       [QN_COUNT][QN_COUNT][QN_COUNT][QN_COUNT] )
   {}
      INSERT_COPY_AND_ASSIGN(gf_vert_bare_ff_t)
}; 
using idx_vert_bare_ff_t = gf_vert_bare_ff_t::idx_t; 

// Class for the precalculation of the projection matrices -> only needed for input to the frg flow 
enum class PROJ_MATRIX{ K_in, K_out, m, n, mp, np}; 
class gf_proj_matrix_t : public gf< dcomplex, 6 > 		///< Container type for two-particle correlation functions
{
   public:
      using base_t = gf< dcomplex, 6>; 

      gf_proj_matrix_t():
	 base_t( boost::extents[PATCH_COUNT][PATCH_COUNT]
	       [FFACTOR_COUNT][FFACTOR_COUNT][FFACTOR_COUNT][FFACTOR_COUNT])
   {}
      INSERT_COPY_AND_ASSIGN(gf_proj_matrix_t)
}; 
using idx_proj_matrix_t = gf_proj_matrix_t::idx_t; 

class gf_weight_vec_t : public gf< double, 1 > 			///< Matrix-valued container type for one-particle correlation functions
{
   public:
      using base_t = gf< double, 1 >; 

      gf_weight_vec_t( int pos_freq_count_= POS_INT_RANGE ):
	 base_t( boost::extents[ffreq(pos_freq_count_)] ) 
   {}
      INSERT_COPY_AND_ASSIGN(gf_weight_vec_t)
}; 

enum class K_SUM_MATRIX{ K_in, K_out}; 
class gf_k_sum_matrix_t : public gf< unsigned int, 2 > 		///< Container type for two-particle correlation functions
{
   public:
      using base_t = gf< unsigned int, 2>; 

      gf_k_sum_matrix_t( ):
	 base_t( boost::extents[PATCH_COUNT][PATCH_COUNT])
   {}
      INSERT_COPY_AND_ASSIGN(gf_k_sum_matrix_t)
}; 
using idx_k_sum_matrix_t = gf_k_sum_matrix_t::idx_t; 

enum class K_ROT_MATRIX{ K }; 
class gf_k_rot_matrix_t : public gf< unsigned int, 1 > 		///< Container type for two-particle correlation functions
{
   public:
      using base_t = gf< unsigned int, 1>; 

      gf_k_rot_matrix_t( ):
	 base_t( boost::extents[PATCH_COUNT])
   {}
      INSERT_COPY_AND_ASSIGN(gf_k_rot_matrix_t)
}; 
using idx_k_rot_matrix_t = gf_k_rot_matrix_t::idx_t; 

enum class K_4MAXKPOS_MATRIX{ K }; 
class gf_k_4MAXKPOS_matrix_t: public gf< double, 1 > 		///< Container type for two-particle correlation functions
{
   public:
      using base_t = gf< double, 1>; 

      gf_k_4MAXKPOS_matrix_t( ):
	 base_t( boost::extents[4*MAX_KPOS])
   {}
      INSERT_COPY_AND_ASSIGN(gf_k_4MAXKPOS_matrix_t)
}; 
using idx_k_4MAXKPOS_matrix_t = gf_k_4MAXKPOS_matrix_t::idx_t; 

enum class FFT_FxF_WEIGHT{ R, m, n }; 
class gf_fft_fxf_weight_t : public gf< dcomplex, 3 > 		///< Container type for two-particle correlation functions
{
   public:
      using base_t = gf< dcomplex, 3>; 

      gf_fft_fxf_weight_t( ):
	 base_t( boost::extents[REAL_GRID][FFACTOR_COUNT][FFACTOR_COUNT])
   {}
      INSERT_COPY_AND_ASSIGN(gf_fft_fxf_weight_t)
}; 
using idx__fft_fxf_weight_t = gf_fft_fxf_weight_t::idx_t; 
