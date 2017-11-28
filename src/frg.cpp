
/************************************************************************************************//**
 *  		
 * 	file: 		frg.cpp
 * 	contents:  	for further documentation see frg.h
 * 
 ****************************************************************************************************/


#include <frg.h>
#include <Eigen/Eigenvalues>
#include <cmath>
#include <H5Tools.h>
#include <params.h>
#include <mymath.h>
#include <iostream>

using namespace Eigen;
using namespace std;

#ifdef READIN

const int POS_FERM_PHI_IN = 12;
const int POS_BOS_PHI_IN = 8;

const int POS_FERM_P_IN = 16;
const int POS_BOS_P_IN = 24;

const int POS_BOS_CHI_IN = 256;

const int POS_SIG_IN = 16; 
// Read the 2P reducible vertex on the form factor basis

gf_phi_t gf_phi_pp_read( POS_BOS_PHI_IN, POS_FERM_PHI_IN); 
gf_phi_t gf_phi_ph_read( POS_BOS_PHI_IN, POS_FERM_PHI_IN);
gf_phi_t gf_phi_xph_read(POS_BOS_PHI_IN, POS_FERM_PHI_IN);

// Read the K2 on the form factor basis
gf_P_t gf_P_pp_read( POS_BOS_P_IN,POS_FERM_P_IN); 
gf_P_t gf_P_ph_read( POS_BOS_P_IN,POS_FERM_P_IN);
gf_P_t gf_P_xph_read(POS_BOS_P_IN,POS_FERM_P_IN);

// Read the K1 on the form factor basis
gf_chi_t gf_chi_pp_read( POS_BOS_CHI_IN); 
gf_chi_t gf_chi_ph_read( POS_BOS_CHI_IN);
gf_chi_t gf_chi_xph_read(POS_BOS_CHI_IN);

//Read the SE
gf_1p_t gf_Sig_read(POS_SIG_IN); 


void read_exact()
{
   using namespace H5;
   H5File input_file( "dat/dat_U2_Beta5_PFCB256_HUB_OMFL_SU2_2D_4PISYMM_mkp2_ONLY_K_SYMM_WITHSE.h5", H5F_ACC_RDONLY );

   cout << "Got file" << endl;

   Group phi_group =  input_file.openGroup("/phi_func");
   
   read( gf_phi_pp_read, phi_group,"_PP");
   cout << "Got phi pp data" << endl; 
   read( gf_phi_ph_read, phi_group,"_PH");
   cout << "Got phi ph data" << endl; 
   read( gf_phi_xph_read, phi_group,"_XPH");
   cout << "Got phi xph data" << endl; 
   
   Group P_group =  input_file.openGroup("/P_func");
   
   read( gf_P_pp_read, P_group,"_PP");
   cout << "Got P pp data" << endl; 
   read( gf_P_ph_read, P_group,"_PH");
   cout << "Got P ph data" << endl; 
   read( gf_P_xph_read, P_group,"_XPH");
   cout << "Got P xph data" << endl; 
   
   Group chi_group =  input_file.openGroup("/chi_func");
   
   read( gf_chi_pp_read, chi_group,"_PP");
   cout << "Got chi pp data" << endl; 
   read( gf_chi_ph_read, chi_group,"_PH");
   cout << "Got chi ph data" << endl; 
   read( gf_chi_xph_read, chi_group,"_XPH");
   cout << "Got chi xph data" << endl; 
   
   Group Sig_group = input_file.openGroup("/Sig") ; 
   read( gf_Sig_read, Sig_group );
   cout << "SE" << endl;
   
}

#endif


/********************* Hybridization function ********************/

#ifdef ED_BATH

MatQN Gam( double w )
{
   dcomplex val = 0.0; 
   for( int ed_idx = 0; ed_idx < energies.size(); ++ed_idx )
      val += hybridizations[ed_idx]*hybridizations[ed_idx] /( I*w - energies[ed_idx]) ; 

   MatQN Gam; 
   Gam <<
      val;
   return Gam; 
}

#elif QMC_BATH

MatQN Gam( double w )
{
   dcomplex val = 0.0;
   val = I / 2.0 * ( w - sgn( w ) * sqrt( 4.0 + w*w ) ); 

   MatQN Gam; 
   Gam <<
      val; 
   return Gam; 
}

#else // Wide band 

MatQN Gam( double w )
{
   double sq = sqrt( w*w + DEL*DEL );
   double dt = DEL / sq * DD;

   MatQN Gam; 
   Gam << 
      -I*w/sq; 
   return Gam; 
}

#endif


/*********************  INTERACTION FlOW  ********************/


MatQN G( double w, double kx, double ky, const MatQN& selfEn )
{
   return ( G0inv( w, kx, ky) - selfEn ).inverse();
}

dcomplex asympt_GG_pp( int W_int )
{
   double PIR = POS_INT_RANGE + abs(W_int/2); 
   double W = W_int; 

   if( W_int == 0 )
      return BETA/(2.0*pow(PI,2.0)*PIR); 
   
   return (BETA*atanh(W/(2.0*PIR)))/(pow(PI,2.0)*W); 
}


/************************ Common to all implemented cutoff schemes ********************************/

#ifdef NO_MOMENTA

// SIAM / SQDJJ, G0 scale independent for interaction cutoff, Matrix in Nambu basis
MatQN G0inv( double w, double kx, double ky )
{
   MatQN ginv; 
   ginv <<		
      I*w - EPS - B;

   return ginv - Gam( w ); 
}

#else

// Hubbard model, G0 scale independent for interaction cutoff, Matrix in Nambu basis

/*****************MULTIORBITAL MATRIX EXPRESSION****************************
MatQN G0inv( double w, double kx, double ky, double kz )
{
   double cos_kx  = cos( kx );
   double cos_ky  = cos( ky );
   double cos_kz  = cos( kz );

   MatQN G0inv; 
   
   G0inv <<		
	I*w+1.5*(cos_kx + cos_ky)-MU, 			sqrt(3.)*0.5*(cos_ky - cos_kx),
        sqrt(3.)*0.5*(cos_ky - cos_kx),  	I*w+0.5*(cos_kx + cos_ky)+2.*cos_kz-MU;
   
   return G0inv; 
}

****************************************************************************/

// Hubbard model, G0 scale independent for interaction cutoff, Matrix in Nambu basis
MatQN G0inv( double w, double kx, double ky )
{
   double cos_kx  = cos( kx );
   double cos_ky  = cos( ky );

   MatQN G0inv; 
   G0inv <<		
      I*w + 2.0 * ( cos_kx + cos_ky ) + 4.0 * (cos_kx * cos_ky ) * T_PRIME - B - MU;

   return G0inv; 
}

#endif

dcomplex asympt_GG_ph( int W )
{
   return -asympt_GG_pp( W ); 
}

// ---- Initial values 

/*********************** MULTIORBITAL INTERACTION MATRIX *********************************

dcomplex vert_bare( int s1_in, int s2_in, int s1_out, int s2_out ) // bare vertex values in Nambu basis
{
   if ((s1_in == s2_in) && (s1_in == s1_out) && (s1_in == s2_out))
      return -UINT;

   else if ((s1_in == s2_out) && (s2_in == s1_out))
      return -J;

   else if ((s1_in == s1_out) && (s2_in == s2_out))
      return -UINTP;

   else if ((s1_in == s2_in) && (s1_out == s2_out))
      return -JP;
   
   else
   return dcomplex(0.0,0.0);
   
}

*****************************************************************************************/

dcomplex vert_bare( const idx_2p_t& idx ) // initial vertex values in Nambu basis
{
   return vert_bare( idx( I2P::s1_in ), idx( I2P::s2_in ), idx( I2P::s1_out ), idx( I2P::s2_out ) ); // return bare interaction 
}
// Unitary transformation from orbital to band basis

dcomplex vert_bare( int s1_in, int s2_in, int s1_out, int s2_out ) // bare vertex values in Nambu basis
{
   return -UINT;
}

dcomplex Sig_init( const idx_1p_t& idx ) // initial self-energy values in Nambu basis   
{
#ifdef READIN
   return gf_Sig_read[idx(0)][idx(1)][idx(2)][idx(3)]; 
#endif
   return 0.0; // vanishes for weak coupling flows
}

dcomplex phi_init_pp( const idx_phi_t& idx ) // initial self-energy values in Nambu basis
{
#ifdef READIN
   return gf_phi_pp_read[idx(IPHI::W)][idx(IPHI::w_in)][idx(IPHI::w_out)][idx(IPHI::K)][idx(IPHI::n_in)][idx(IPHI::n_out)][idx(IPHI::s1_in)][idx(IPHI::s2_in)][idx(IPHI::s1_out)][idx(IPHI::s2_out)] ; 
#endif
   return 0.0; // vanishes for weak coupling flows
}

dcomplex phi_init_ph( const idx_phi_t& idx ) // initial self-energy values in Nambu basis
{
#ifdef READIN
   return gf_phi_ph_read[idx(IPHI::W)][idx(IPHI::w_in)][idx(IPHI::w_out)][idx(IPHI::K)][idx(IPHI::n_in)][idx(IPHI::n_out)][idx(IPHI::s1_in)][idx(IPHI::s2_in)][idx(IPHI::s1_out)][idx(IPHI::s2_out)]; 
#endif
   return 0.0; // vanishes for weak coupling flows
}

dcomplex phi_init_xph( const idx_phi_t& idx ) // initial self-energy values in Nambu basis
{
#ifdef READIN
   return gf_phi_xph_read[idx(IPHI::W)][idx(IPHI::w_in)][idx(IPHI::w_out)][idx(IPHI::K)][idx(IPHI::n_in)][idx(IPHI::n_out)][idx(IPHI::s1_in)][idx(IPHI::s2_in)][idx(IPHI::s1_out)][idx(IPHI::s2_out)]; 
#endif
   return 0.0; // vanishes for weak coupling flows
}

dcomplex P_init_pp( const idx_P_t& idx ) // initial self-energy values in Nambu basis
{
#ifdef READIN
   return gf_P_pp_read[idx(IP::W)][idx(IP::w)][idx(IP::K)][idx(IP::n)][idx(IP::s1_in)][idx(IP::s2_in)][idx(IP::s1_out)][idx(IP::s2_out)]; 
#endif
   return 0.0; // vanishes for weak coupling flows
}

dcomplex P_init_ph( const idx_P_t& idx ) // initial self-energy values in Nambu basis
{
#ifdef READIN
   return gf_P_ph_read[idx(IP::W)][idx(IP::w)][idx(IP::K)][idx(IP::n)][idx(IP::s1_in)][idx(IP::s2_in)][idx(IP::s1_out)][idx(IP::s2_out)]; 
#endif
   return 0.0; // vanishes for weak coupling flows
}

dcomplex P_init_xph( const idx_P_t& idx ) // initial self-energy values in Nambu basis
{
#ifdef READIN
   return gf_P_xph_read[idx(IP::W)][idx(IP::w)][idx(IP::K)][idx(IP::n)][idx(IP::s1_in)][idx(IP::s2_in)][idx(IP::s1_out)][idx(IP::s2_out)]; 
#endif
   return 0.0; // vanishes for weak coupling flows
}

dcomplex chi_init_pp( const idx_chi_t& idx ) // inital Karrasch functions in Nambu basis
{
#ifdef READIN
   return gf_chi_pp_read[idx(ICHI::W)][idx(ICHI::K)][idx(ICHI::s1_in)][idx(ICHI::s2_in)][idx(ICHI::s1_out)][idx(ICHI::s2_out)]; 
#endif
   return 0.0; // vanishes for weak coupling flows
}

dcomplex chi_init_ph( const idx_chi_t& idx ) // inital Karrasch functions in Nambu basis
{
#ifdef READIN
   return gf_chi_ph_read[idx(ICHI::W)][idx(ICHI::K)][idx(ICHI::s1_in)][idx(ICHI::s2_in)][idx(ICHI::s1_out)][idx(ICHI::s2_out)]; 
#endif
   return 0.0; // vanishes for weak coupling flows
}

dcomplex chi_init_xph( const idx_chi_t& idx ) // inital Karrasch functions in Nambu basis
{
#ifdef READIN
   return gf_chi_xph_read[idx(ICHI::W)][idx(ICHI::K)][idx(ICHI::s1_in)][idx(ICHI::s2_in)][idx(ICHI::s1_out)][idx(ICHI::s2_out)]; 
#endif
   return 0.0; // vanishes for weak coupling flows
}
dcomplex suscept_init_t( const idx_suscept_t& idx ) // inital Karrasch functions in Nambu basis
{
   return 0.0; // vanishes for weak coupling flows
}

dcomplex suscept_init_s( const idx_suscept_t& idx ) // inital Karrasch functions in Nambu basis
{
   return 0.0; // vanishes for weak coupling flows
}

dcomplex suscept_init_d( const idx_suscept_t& idx ) // inital Karrasch functions in Nambu basis
{
   return 0.0; // vanishes for weak coupling flows
}

dcomplex suscept_init_m( const idx_suscept_t& idx ) // inital Karrasch functions in Nambu basis
{
   return 0.0; // vanishes for weak coupling flows
}

dcomplex tri_init( const idx_tri_t& idx ) // inital Karrasch functions in Nambu basis
{
  return tri_bare(idx(ITRI::n_in),idx(ITRI::n_out), idx(ITRI::s1_in), idx(ITRI::s2_in), idx(ITRI::s1_out), idx(ITRI::s2_out)); // vanishes for weak coupling flows
  //return 0.0;
}

dcomplex tri_bare(int n_in, int n_out, int s1_in, int s2_in, int s1_out, int s2_out )
{
  return double(n_in == n_out);
}
