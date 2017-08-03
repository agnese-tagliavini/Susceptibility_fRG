
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

   write( state_vec.gf_suscept_trip(), group, "_TRIPLET" ); 
   write( state_vec.gf_suscept_s(), group, "_SINGLET" ); 
   write( state_vec.gf_suscept_d(), group, "_DENSITY" ); 
   write( state_vec.gf_suscept_m(), group, "_MAGNETIC" ); 

   write( Bos_Grid( POS_BFREQ_COUNT_SUSCEPT, 2.0*PI / BETA ), group );
}

