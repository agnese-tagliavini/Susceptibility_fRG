

/*******************************************************************************************//** @file
 *  		
 * 	file: 		projection.h
 * 	contents:  	Definition of projection-related quantities 
 * 
****************************************************************************************************/


#pragma once

#include <def.h>
#include <grid.h>
#include <params.h>
//#include <state.h>

using proj_matrix_t = gf_proj_matrix_t;
using vert_bare_ff_t = gf_vert_bare_ff_t;

dcomplex generate_proj_matrix_ph_to_pp( const idx_proj_matrix_t& idx );
dcomplex generate_proj_matrix_pp_to_ph( const idx_proj_matrix_t& idx );
dcomplex generate_proj_matrix_ph_to_xph( const idx_proj_matrix_t& idx  );

proj_matrix_t generate_proj_matrix_ph_to_pp( );
proj_matrix_t generate_proj_matrix_pp_to_ph( );
proj_matrix_t generate_proj_matrix_ph_to_xph( );

dcomplex generate_proj_vert_bare(const idx_vert_bare_ff_t& idx);

