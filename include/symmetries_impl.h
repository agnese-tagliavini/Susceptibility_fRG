
#include <grid.h>

using namespace std;

// --- Two particle symmetries ( Copied from Code working with vertex )

template<typename notation>
operation exch_in( idx_2p_t& idx_not )
{
   idx_2p_t idx = to_ferm<notation>( idx_not ); 
   swap( idx( I2P::w1_in ), idx( I2P::w2_in ) );
   swap( idx( I2P::k1_in ), idx( I2P::k2_in ) );
   swap( idx( I2P::s1_in ), idx( I2P::s2_in ) );

   idx_not = from_ferm<notation>( idx ); 
   return operation( true, false );
}

template<typename notation>
operation exch_out( idx_2p_t& idx_not ) 
{
   idx_2p_t idx = to_ferm<notation>( idx_not ); 
   int w2_out = idx( I2P::w1_in ) + idx( I2P::w2_in ) - idx( I2P::w1_out ); // calculate w2_out by means of frequency conservation
   idx( I2P::w1_out ) = w2_out;
   idx( I2P::k1_out ) = dif_k( add_k( idx( I2P::k1_in ), idx( I2P::k2_in ) ), idx( I2P::k1_out ) ) ; // calculate k2_out and assign to k1_out
   swap( idx( I2P::s1_out ), idx( I2P::s2_out ) );

   idx_not = from_ferm<notation>( idx ); 
   return operation( true, false );
}

template<typename notation>
operation compl_conj( idx_2p_t& idx_not )
{
   idx_2p_t idx = to_ferm<notation>( idx_not ); 
   int w2_out = idx( I2P::w1_in ) + idx( I2P::w2_in ) - idx( I2P::w1_out ); // calculate w2_out by means of frequency conservation
   swap( idx( I2P::w1_in ), idx( I2P::w1_out ) );
   idx( I2P::w2_in ) = w2_out;

   idx( I2P::w1_in ) = -idx( I2P::w1_in ) - 1;
   idx( I2P::w2_in ) = -idx( I2P::w2_in ) - 1;
   idx( I2P::w1_out ) = -idx( I2P::w1_out ) - 1;

   // Changing all momenta signs is unnecessary since equal to rotating twice by 90 degrees

   int k2_out = dif_k( add_k( idx( I2P::k1_in ), idx( I2P::k2_in ) ), idx( I2P::k1_out ) ) ; // calculate k2_out 
   swap( idx( I2P::k1_in ), idx( I2P::k1_out ) ); 	// swap k1_in and k1_out
   idx( I2P::k2_in ) = k2_out; 		// swap k2_in and k2_out

   swap( idx( I2P::s1_in ), idx( I2P::s1_out ) );
   swap( idx( I2P::s2_in ), idx( I2P::s2_out ) );

   idx_not = from_ferm<notation>( idx ); 
   return operation( false, true );
}

template<typename notation>
operation time_rev( idx_2p_t& idx_not )
{
   idx_2p_t idx = to_ferm<notation>( idx_not ); 
   int w2_out = idx( I2P::w1_in ) + idx( I2P::w2_in ) - idx( I2P::w1_out ); // calculate w2_out by means of frequency conservation
   swap( idx( I2P::w1_in ), idx( I2P::w1_out ) );
   idx( I2P::w2_in ) = w2_out;

   // Changing all momenta signs is unnecessary since equal to rotating twice by 90 degrees

   int k2_out = dif_k( add_k( idx( I2P::k1_in ), idx( I2P::k2_in ) ), idx( I2P::k1_out ) ) ; // calculate k2_out 
   swap( idx( I2P::k1_in ), idx( I2P::k1_out ) ); 	// swap k1_in and k1_out
   idx( I2P::k2_in ) = k2_out; 		// swap k2_in and k2_out

   swap( idx( I2P::s1_in ), idx( I2P::s1_out ) );
   swap( idx( I2P::s2_in ), idx( I2P::s2_out ) );

   idx_not = from_ferm<notation>( idx ); 
   return operation( false, false );
}

template<typename notation>
operation particle_hole( idx_2p_t& idx_not )
{
   idx_2p_t idx = to_ferm<notation>( idx_not ); 
#ifndef NO_MOMENTA
   mirror_mom_pipipi( idx( I2P::k1_in ) );
   mirror_mom_pipipi( idx( I2P::k2_in ) );
   mirror_mom_pipipi( idx( I2P::k1_out ) );
#endif

   idx_not = from_ferm<notation>( idx ); 
   return operation( false, true );
}

template<typename notation>
operation spin_symm( idx_2p_t& idx_not )
{
   idx_2p_t idx = to_ferm<notation>( idx_not ); 

   // ----- NAMBU VERSION ( first part not necessary as long as complex conjugation is present as well )
   int w2_out = idx( I2P::w1_in ) + idx( I2P::w2_in ) - idx( I2P::w1_out ); // calculate w2_out by means of frequency conservation
   swap( idx( I2P::w1_in ), idx( I2P::w1_out ) );
   idx( I2P::w2_in ) = w2_out;

   idx( I2P::w1_in ) = -idx( I2P::w1_in ) - 1;
   idx( I2P::w2_in ) = -idx( I2P::w2_in ) - 1;
   idx( I2P::w1_out ) = -idx( I2P::w1_out ) - 1;

   // Changing all momenta signs is unnecessary since equal to rotating twice by 90 degrees

   int k2_out = dif_k( add_k( idx( I2P::k1_in ), idx( I2P::k2_in ) ), idx( I2P::k1_out ) ); // calculate k2_out 
   swap( idx( I2P::k1_in ), idx( I2P::k1_out ) ); 	// swap k1_in and k1_out
   idx( I2P::k2_in ) = k2_out; 				// swap k2_in and k2_out

   swap( idx( I2P::s1_in ), idx( I2P::s1_out ) );
   swap( idx( I2P::s2_in ), idx( I2P::s2_out ) );

   flip_spin( idx( I2P::s1_in ) );
   flip_spin( idx( I2P::s2_in ) );
   flip_spin( idx( I2P::s1_out ) );
   flip_spin( idx( I2P::s2_out ) );

   idx_not = from_ferm<notation>( idx ); 
   if ( ( idx( I2P::s1_in ) + idx( I2P::s2_in ) + !( idx( I2P::s1_out ) ) + !( idx( I2P::s2_out ) ) ) % 2 != 0 ) // Minus sign for uneven amount of creation operators
      return operation( true, false );

   return operation( false, false );
}


template<typename notation>
operation rot_pi2_x( idx_2p_t& idx_not )
{
   idx_2p_t idx = to_ferm<notation>( idx_not ); 
   rot_pi2_x( idx( I2P::k1_in ) );
   rot_pi2_x( idx( I2P::k2_in ) );
   rot_pi2_x( idx( I2P::k1_out ) );

   idx_not = from_ferm<notation>( idx ); 
   return operation( false, false );
}

   
template<typename notation>
operation rot_pi2_y( idx_2p_t& idx_not )
{
   idx_2p_t idx = to_ferm<notation>( idx_not ); 
   rot_pi2_y( idx( I2P::k1_in ) );
   rot_pi2_y( idx( I2P::k2_in ) );
   rot_pi2_y( idx( I2P::k1_out ) );

   idx_not = from_ferm<notation>( idx ); 
   return operation( false, false );
}


template<typename notation>
operation rot_pi2_z( idx_2p_t& idx_not )
{
   idx_2p_t idx = to_ferm<notation>( idx_not ); 
   rot_pi2_z( idx( I2P::k1_in ) );
   rot_pi2_z( idx( I2P::k2_in ) );
   rot_pi2_z( idx( I2P::k1_out ) );

   idx_not = from_ferm<notation>( idx ); 
   return operation( false, false );
}


template<typename notation>
operation rot_2pi3_111( idx_2p_t& idx_not )
{
   idx_2p_t idx = to_ferm<notation>( idx_not ); 
   rot_2pi3_111( idx( I2P::k1_in ) );
   rot_2pi3_111( idx( I2P::k2_in ) );
   rot_2pi3_111( idx( I2P::k1_out ) );

   idx_not = from_ferm<notation>( idx ); 
   return operation( false, false );
}

   
template<typename notation>
operation rot_2pi3_1m11( idx_2p_t& idx_not )
{
   idx_2p_t idx = to_ferm<notation>( idx_not ); 
   rot_2pi3_1m11( idx( I2P::k1_in ) );
   rot_2pi3_1m11( idx( I2P::k2_in ) );
   rot_2pi3_1m11( idx( I2P::k1_out ) );

   idx_not = from_ferm<notation>( idx ); 
   return operation( false, false );
}


template<typename notation>
operation rot_2pi3_m111( idx_2p_t& idx_not )
{
   idx_2p_t idx = to_ferm<notation>( idx_not ); 
   rot_2pi3_m111( idx( I2P::k1_in ) );
   rot_2pi3_m111( idx( I2P::k2_in ) );
   rot_2pi3_m111( idx( I2P::k1_out ) );

   idx_not = from_ferm<notation>( idx ); 
   return operation( false, false );
}


template<typename notation>
operation rot_2pi3_m1m11( idx_2p_t& idx_not )
{
   idx_2p_t idx = to_ferm<notation>( idx_not ); 
   rot_2pi3_m1m11( idx( I2P::k1_in ) );
   rot_2pi3_m1m11( idx( I2P::k2_in ) );
   rot_2pi3_m1m11( idx( I2P::k1_out ) );

   idx_not = from_ferm<notation>( idx ); 
   return operation( false, false );
}


template<typename notation>
operation rot_pi_110( idx_2p_t& idx_not )
{
   idx_2p_t idx = to_ferm<notation>( idx_not ); 
   rot_pi_110( idx( I2P::k1_in ) );
   rot_pi_110( idx( I2P::k2_in ) );
   rot_pi_110( idx( I2P::k1_out ) );

   idx_not = from_ferm<notation>( idx ); 
   return operation( false, false );
}

template<typename notation>
operation rot_pi_1m10( idx_2p_t& idx_not )
{
   idx_2p_t idx = to_ferm<notation>( idx_not ); 
   rot_pi_1m10( idx( I2P::k1_in ) );
   rot_pi_1m10( idx( I2P::k2_in ) );
   rot_pi_1m10( idx( I2P::k1_out ) );

   idx_not = from_ferm<notation>( idx ); 
   return operation( false, false );
}

template<typename notation>
operation rot_pi_011( idx_2p_t& idx_not )
{
   idx_2p_t idx = to_ferm<notation>( idx_not ); 
   rot_pi_011( idx( I2P::k1_in ) );
   rot_pi_011( idx( I2P::k2_in ) );
   rot_pi_011( idx( I2P::k1_out ) );

   idx_not = from_ferm<notation>( idx ); 
   return operation( false, false );
}


template<typename notation>
operation rot_pi_01m1( idx_2p_t& idx_not )
{
   idx_2p_t idx = to_ferm<notation>( idx_not ); 
   rot_pi_01m1( idx( I2P::k1_in ) );
   rot_pi_01m1( idx( I2P::k2_in ) );
   rot_pi_01m1( idx( I2P::k1_out ) );

   idx_not = from_ferm<notation>( idx ); 
   return operation( false, false );
}

template<typename notation>
operation rot_pi_101( idx_2p_t& idx_not )
{
   idx_2p_t idx = to_ferm<notation>( idx_not ); 
   rot_pi_101( idx( I2P::k1_in ) );
   rot_pi_101( idx( I2P::k2_in ) );
   rot_pi_101( idx( I2P::k1_out ) );

   idx_not = from_ferm<notation>( idx ); 
   return operation( false, false );
}

template<typename notation>
operation rot_pi_10m1( idx_2p_t& idx_not )
{
   idx_2p_t idx = to_ferm<notation>( idx_not ); 
   rot_pi_10m1( idx( I2P::k1_in ) );
   rot_pi_10m1( idx( I2P::k2_in ) );
   rot_pi_10m1( idx( I2P::k1_out ) );

   idx_not = from_ferm<notation>( idx ); 
   return operation( false, false );
}

