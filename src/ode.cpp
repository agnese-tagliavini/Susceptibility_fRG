
/************************************************************************************************//**
 *  		
 * 	file: 		ode.cpp
 * 	contents:   	main routine, sets up and runs ODE solver, then outputs results
 * 
 ****************************************************************************************************/


#include <ode.h>
#include <boost/numeric/odeint.hpp>
#include <rhs.h>
#include <params.h>
#include <frg.h>
#include <output.h>
#include <state.h>
#include <projection.h>
#include <grid.h>

using namespace boost::numeric::odeint;
using namespace H5; 
using namespace std; 

int main ( int argc, char * argv[])
{
   // ------ Set output format

   cout << scientific << setprecision(4); 

   // ------ Read in parameters from file

   if( argc == readIn_lst.size() + 1 )
   {
      for( int i = 0; i < readIn_lst.size(); i++ )
	 *readIn_lst[i] = atof( argv[i+1]);
      update_dep_params(); 
   }
   else
   {
      cout << "Wrong amount of arguments, using defaults." << endl << endl; 
   }

   // ------ Initialize rhs object and observer list

   cout << " Initializing RHS... " << endl << endl;
   cout << endl;
   cout << " WARNING: THE INPUT FILE SHOULD HAVE THE SAME PARAMETERS SPECIFIED IN CONST.H!!!" << endl;
   rhs_t rhs;

#ifdef READIN
   // ------ Reading fRG data as specified in fRG.cpp
   cout << "Reading from file:" << endl;
   read_exact(); 
#endif

   // ------ Write initial values

   cout << " Writing initial conditions... " << endl << endl; 
   state_t state_vec; 

   // Write initial values for Sig, phi's, P's and chi's
   state_vec.gf_Sig().init( Sig_init ); 

   state_vec.gf_phi_pp().init( phi_init_pp ); 
   state_vec.gf_phi_ph().init( phi_init_ph ); 
   state_vec.gf_phi_xph().init( phi_init_xph );

   state_vec.gf_P_pp().init( P_init_pp ); 
   state_vec.gf_P_ph().init( P_init_ph ); 
   state_vec.gf_P_xph().init( P_init_xph ); 

   state_vec.gf_chi_pp().init( chi_init_pp ); 
   state_vec.gf_chi_ph().init( chi_init_ph ); 
   state_vec.gf_chi_xph().init( chi_init_xph );

   state_vec.gf_suscept_sc().init( suscept_init_sc ); 
   state_vec.gf_suscept_d().init( suscept_init_d );
   state_vec.gf_suscept_m().init( suscept_init_m );
   
   state_vec.gf_tri_sc().init( tri_init ); 
   state_vec.gf_tri_d().init( tri_init );
   state_vec.gf_tri_m().init( tri_init );


/************************** CALC SUSCEPTIBILITIES **************************************/
   
   state_t state_vec_old;
   
   state_vec_old = state_vec;

   rhs( state_vec_old, state_vec );
   
   // ------ Write output file

   H5std_string FILE_NAME("dat/suscept");
   vector<pair<string, double>> fname_params = {{ "U", UINT }, { "Beta", BETA }, { "PFCB", POS_BFREQ_COUNT_CHI }}; //{ "PFC", POS_FFREQ_COUNT_P }}; //, { "Eps", EPS }};
   for( auto parStr: fname_params )
   {
      FILE_NAME.append("_" + parStr.first );
      string valStr = to_string( ( long double ) parStr.second );
      valStr.erase( valStr.find_last_not_of('0') + 1, string::npos ); // delete trailing zeros
      valStr.erase( valStr.find_last_not_of('.') + 1, string::npos ); // delete trailing dot
      replace( valStr.begin(), valStr.end(), '.', 'p');	// replace dot with p
      FILE_NAME.append( valStr ); 
   }

#ifdef NO_MOMENTA

#ifdef QMC_BATH
   FILE_NAME.append("_QMCB"); 
#elif ED_BATH
   FILE_NAME.append("_EDB"); 
#endif

#else
   FILE_NAME.append("_HUB"); 
#endif

#ifdef KATANIN
   FILE_NAME.append("_KAT"); 
#endif

#ifdef TWOLOOP
   FILE_NAME.append("_2LOOP"); 
#endif

   FILE_NAME.append("_SU2"); 
   FILE_NAME.append("_2D");
   FILE_NAME.append("_8LOOP");

   FILE_NAME.append("_OMFL.h5"); 

   //const H5std_string	FILE_NAME("dat/dat.h5");
   H5File file( FILE_NAME, H5F_ACC_TRUNC );

   init( file );		
   write_config( file );		
   write_params( file );
   write_suscept_func( file, state_vec );
   write_tri_func( file, state_vec );
   write_asytri_func( file, state_vec );
   //write_Sig_tensor( file, state_vec );
   //write_vert_func( file, state_vec ); 
   //write_phi_func( file, state_vec ); 
   //write_chi_func( file, state_vec ); 
   //write_P_func( file, state_vec ); 
   
   return 0;

}
