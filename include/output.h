
/*******************************************************************************************//** @file
 *  		
 * 	file: 		output.h
 * 	contents:  	Functions to write output to hdf5 file
 * 
 ****************************************************************************************************/


#pragma once

#include <H5Cpp.h>
#include <rhs.h>

void init( H5::H5File& file );			///< Initialize with time and date
void write_config( H5::H5File& file );		///< Write configuration to output file
void write_params( H5::H5File& file );		///< Write parameters to output file
void write_susc_func( H5::H5File& file, const state_t& state_vec );  	///<	Write 1PI two-particle vertex to output file
void write_tri_func( H5::H5File& file, const state_t& state_vec );  	///<	Write 1PI two-particle vertex to output file
void write_asytri_func( H5::H5File& file, const state_t& state_vec );  	///<	Write 1PI two-particle vertex to output file
void write_Sig_tensor( H5::H5File& file, const state_t& state_vec );		///<	Write self-energy to output file
void write_Sig_old_tensor( H5::H5File& file, const state_t& state_vec );		///<	Write self-energy to output file
void write_Sig_ph_tensor( H5::H5File& file, const state_t& state_vec );		///<	Write self-energy to output file
void write_Sig_pp_tensor( H5::H5File& file, const state_t& state_vec );		///<	Write self-energy to output file
void write_Sig_xph_tensor( H5::H5File& file, const state_t& state_vec );		///<	Write self-energy to output file
void write_vert_func( H5::H5File& file, const state_t& state_vec ); 		///< 	Write phi functions to output file
void write_phi_func( H5::H5File& file, const state_t& state_vec ); 		///< 	Write phi functions to output file
void write_chi_func( H5::H5File& file, const state_t& state_vec );  	///<	Write Karrasch functions to output file
void write_P_func( H5::H5File& file, const state_t& state_vec ); 		///<	Write P functions to output file
void write_Giw_tensor( H5::H5File& file, const state_t& state_vec );		///<	Write self-energy to output file
