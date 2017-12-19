
/************************************************************************************************//**
 *  		
 * 	file: 		output.cpp
 * 	contents:   	See output.h
 * 
 ****************************************************************************************************/


#include <output.h>
#include <params.h>
#include <const.h>
#include <frg.h>
#include <grid.h>
#include <ctime>
#include <H5Tools.h>

//using namespace boost;
using namespace H5; 
using namespace std; 


void init( H5File& file )
{
   time_t curr_time = time( nullptr );
   char* time_str = asctime( localtime(&curr_time )); 
      
   DataSpace dsp = DataSpace( H5S_SCALAR );

   StrType strdatatype( PredType::C_S1, H5T_VARIABLE ); // String type with variable length

   hid_t time_att_id = H5Acreate2( file.getId(), "DATE_TIME", strdatatype.getId(), dsp.getId(), H5P_DEFAULT, H5P_DEFAULT ); // Adds attribute to root group of file - find c++ equivalent! 

   Attribute time_att( time_att_id );
   time_att.write( strdatatype, &time_str );
}

void write_config( H5File& file )
{
   Group group( file.createGroup("/Config"));

   vector<pair<string, double>> config_scalar_list = 
   { 
      { "POS_FFREQ_COUNT_SIG", POS_FFREQ_COUNT_SIG }, 
      { "POS_FFREQ_COUNT_PHI", POS_FFREQ_COUNT_PHI }, 
      { "POS_BFREQ_COUNT_PHI", POS_BFREQ_COUNT_PHI }, 
      { "POS_FFREQ_COUNT_P", POS_FFREQ_COUNT_P }, 
      { "POS_BFREQ_COUNT_P", POS_BFREQ_COUNT_P }, 
      { "POS_BFREQ_COUNT_CHI", POS_BFREQ_COUNT_CHI }, 
      { "POS_INT_RANGE", POS_INT_RANGE }, 
      { "PATCH_COUNT", PATCH_COUNT }, 
      { "QN_COUNT", QN_COUNT }, 
      { "FFACTOR_COUNT", FFACTOR_COUNT},
   };

   for( auto conf : config_scalar_list )
      write( conf.second, group, conf.first); 
   
   DataSpace config_dsp = DataSpace ( H5S_SCALAR );

   int katanin_int =
#ifdef KATANIN
      1;
#else
      0; 
#endif
   write( katanin_int, group, "KATANIN" ); 
   
   int forced_zeros_int =
#ifdef FORCED_ZEROS
      1;
#else
      0; 
#endif
   write( forced_zeros_int, group, "FORCED_ZEROS" ); 
}

void write_params( H5File& file )
{

   Group group( file.createGroup("/Params"));

vector<pair<string, double>> par_lst = 
{ 
   { "UINT", UINT },
/************** MULTIORBITAL CASE********** 
   { "UINTP", UNITP },
   { "J", J },
   { "JP", JP }, 
*****************************************/
   { "BETA", BETA },
   { "B", B }, 
   { "MU", MU } 
};

   for( auto par : par_lst )
      write( par.second, group, par.first); 

}

void write_suscept_func( H5File& file, const state_t& state_vec )
{
   Group group( file.createGroup("/suscept_func") );

   write( state_vec.gf_suscept_sc(), group, "_SC" ); 
   write( state_vec.gf_suscept_d(), group, "_DENSITY" ); 
   write( state_vec.gf_suscept_m(), group, "_MAGNETIC" ); 

   write( Bos_Grid( POS_BFREQ_COUNT_SUSCEPT, 2.0*PI / BETA ), group );
}

void write_tri_func( H5File& file, const state_t& state_vec )
{
   Group group( file.createGroup("/tri_func") );

   write( state_vec.gf_tri_sc(), group, "_SC" ); 
   write( state_vec.gf_tri_d(), group, "_DENSITY" ); 
   write( state_vec.gf_tri_m(), group, "_MAGNETIC" ); 

   write( Bos_Grid( POS_BFREQ_COUNT_TRI, 2.0*PI / BETA ), group );
   write( F_Grid( POS_FFREQ_COUNT_TRI, 2.0*PI / BETA ), group );
}

void write_asytri_func( H5File& file, const state_t& state_vec )
{
   Group group( file.createGroup("/asytri_func") );

   write( state_vec.gf_asytri_sc(), group, "_SC" ); 
   write( state_vec.gf_asytri_d(), group, "_DENSITY" ); 
   write( state_vec.gf_asytri_m(), group, "_MAGNETIC" ); 

   write( Bos_Grid( POS_BFREQ_COUNT_SUSCEPT, 2.0*PI / BETA ), group );
}

void write_Sig_tensor( H5::H5File& file, const state_t& state_vec )
{
   Group group( file.createGroup("/Sig") );
   write( state_vec.gf_Sig(), group ); 
   write( F_Grid( POS_FFREQ_COUNT_SIG, 2.0*PI / BETA ), group );
}


void write_vert_func( H5File& file, const state_t& state_vec )
{
   Group group( file.createGroup("/vert_func") );

   gf_phi_t gf_vert_pp_plot( POS_BFREQ_COUNT_PHI, POS_PLOT_RANGE_PHI ); 
   gf_phi_t gf_vert_ph_plot( POS_BFREQ_COUNT_PHI, POS_PLOT_RANGE_PHI ); 
   gf_phi_t gf_vert_xph_plot( POS_BFREQ_COUNT_PHI, POS_PLOT_RANGE_PHI ); 

   gf_vert_pp_plot.init( bind( &state_t::vertx_pp, boost::cref(state_vec), _1 ) ); 
   gf_vert_ph_plot.init( bind( &state_t::vertx_ph, boost::cref(state_vec), _1 ) ); 
   gf_vert_xph_plot.init( bind( &state_t::vertx_xph, boost::cref(state_vec), _1 ) ); 

   write( gf_vert_pp_plot, group, "_PP" ); 
   write( gf_vert_ph_plot, group, "_PH" ); 
   write( gf_vert_xph_plot, group, "_XPH" ); 

   write( Bos_Grid( POS_BFREQ_COUNT_PHI, 2.0*PI / BETA ), group );
   write( F_Grid( POS_PLOT_RANGE_PHI, 2.0*PI / BETA ), group );

}

void write_phi_func( H5File& file, const state_t& state_vec )
{
   Group group( file.createGroup("/phi_func") );

   gf_phi_t gf_phi_pp_plot( POS_BFREQ_COUNT_PHI, POS_PLOT_RANGE_PHI ); 
   gf_phi_t gf_phi_ph_plot( POS_BFREQ_COUNT_PHI, POS_PLOT_RANGE_PHI ); 
   gf_phi_t gf_phi_xph_plot( POS_BFREQ_COUNT_PHI, POS_PLOT_RANGE_PHI ); 

   gf_phi_pp_plot.init( bind( &state_t::phi_pp, boost::cref(state_vec), _1 ) ); 
   gf_phi_ph_plot.init( bind( &state_t::phi_ph, boost::cref(state_vec), _1 ) ); 
   gf_phi_xph_plot.init( bind( &state_t::phi_xph, boost::cref(state_vec), _1 ) ); 

   write( gf_phi_pp_plot, group, "_PP" ); 
   write( gf_phi_ph_plot, group, "_PH" ); 
   write( gf_phi_xph_plot, group, "_XPH" ); 

   write( Bos_Grid( POS_BFREQ_COUNT_PHI, 2.0*PI / BETA ), group );
   write( F_Grid( POS_PLOT_RANGE_PHI, 2.0*PI / BETA ), group );

}

void write_chi_func( H5File& file, const state_t& state_vec )
{
   Group group( file.createGroup("/chi_func") );

   write( state_vec.gf_chi_pp(), group, "_PP" ); 
   write( state_vec.gf_chi_ph(), group, "_PH" ); 
   write( state_vec.gf_chi_xph(), group, "_XPH" ); 

   write( Bos_Grid( POS_BFREQ_COUNT_CHI, 2.0*PI / BETA ), group );
}

void write_P_func( H5File& file, const state_t& state_vec )
{
   Group group( file.createGroup("/P_func") );

   write( state_vec.gf_P_pp(), group, "_PP" ); 
   write( state_vec.gf_P_ph(), group, "_PH" ); 
   write( state_vec.gf_P_xph(), group, "_XPH" ); 

   write( Bos_Grid( POS_BFREQ_COUNT_P, 2.0*PI / BETA ), group );
   write( F_Grid( POS_FFREQ_COUNT_P, 2.0*PI / BETA ), group );
}

