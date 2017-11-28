
/*******************************************************************************************//** @file
 *  		
 * 	file: 		translate.cpp
 * 	contents:  	Functions that allow for translations between different notations
 * 
 ****************************************************************************************************/


#include <translate.h>
#include <grid.h>
#include <cmath>
#include <mymath.h>


idx_phi_t iP_to_iphi( const idx_P_t& idx )
{
   return  idx_phi_t( { idx( IP::W ), idx( IP::w ), 5*BFREQ_COUNT_CHI, idx( IP::K ), idx( IP::n ), 0, idx( IP::s1_in ), idx( IP::s2_in ), idx( IP::s1_out ), idx( IP::s2_out )} ); 
}

idx_phi_t ichi_to_iphi( const idx_chi_t& idx )
{
   return  idx_phi_t( { idx( ICHI::W ), 5*BFREQ_COUNT_CHI, 5*BFREQ_COUNT_CHI, idx( ICHI::K ), 0, 0, idx( ICHI::s1_in ), idx( ICHI::s2_in ), idx( ICHI::s1_out ), idx( ICHI::s2_out )} );
}

idx_chi_t iP_to_ichi( const idx_P_t& idx )
{
   return  idx_chi_t( { idx( IP::W ), idx( IP::K ), idx( IP::s1_in ), idx( IP::s2_in ), idx( IP::s1_out ), idx( IP::s2_out )} ); 
}

idx_tri_t iasy_to_itri( const idx_chi_t& idx )
{
   return  idx_tri_t( { idx( ICHI::W ), 5*BFREQ_COUNT_CHI, idx( ICHI::K ), 0, 0, idx( ICHI::s1_in ), idx( ICHI::s2_in ), idx( ICHI::s1_out ), idx( ICHI::s2_out )} ); 
}


