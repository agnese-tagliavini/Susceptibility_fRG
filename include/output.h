
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
void write_suscept_func( H5::H5File& file, const state_t& state_vec );  	///<	Write 1PI two-particle vertex to output file
