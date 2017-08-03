
/******************************************************************************************//** @file
 *  		
 * 	file: 		mymath.h
 * 	contents:  	Contains useful mathematical functions / routines
 * 
 ****************************************************************************************************/

#pragma once

#include <def.h>
#include <gsl/gsl_errno.h>
#include <complex>
#include <cmath>

//---------------- MatReal get indices 1d array to 3d array -----------------------------------

inline int Matidx(int rx, int ry){
   return (ry+10*FFT_DIM) % FFT_DIM+FFT_DIM*((rx+10*FFT_DIM) % FFT_DIM);
}

inline int Matidxshift(int rx, int ry){
   return (ry+FFT_DIM/2 + 10*FFT_DIM) % FFT_DIM+FFT_DIM*((rx+FFT_DIM/2 + 10*FFT_DIM) % FFT_DIM );

}


// ---- Useful functions

inline int div2_ceil( int W )
{
   return (W - 1000000)/2 + 500000; 
}

inline int div2_floor( int W )
{
   return (W + 1000000)/2 - 500000; 
}

template <typename T> double sgn(T val) 
{
   return (T(0) < val) - (val < T(0));
}

template <typename T> double Theta(T val)
{
   return 0.5 * ( sgn<T>( val ) + 1 ); 
}

template <typename T1, typename T2> double Theta( T1 val1, T2 val2 )
{
   return Theta( val1 ) * Theta( val2 ); 
}

double kdel( double a, double b );	///< Kronecker delta for two double values


// --- Self Consistency ---

typedef void(*sc_Func )( const std::complex<double>[], std::complex<double>[], void * q ); ///< General function for self_consistency
double dist( const dcomplex a[], const dcomplex b[], int len );	///< Arbitrary distance measure on array of length len
/**
 * 	Perform self consistency loop using general function type defined above, till the distance measure dist(..) between two
 * 	steps fullfills the errorbound err	
 */
int self_consistency( std::complex<double> init[], std::complex<double> fin[], int len, sc_Func fw, double err, void * q );

void rotation_matrix_2d( int n, double rot[2][2]);
void mirror_matrix_y(  double rot[2][2]);
void mirror_matrix_diagonal(  double rot[2][2]);
void rotation_matrix_3d( double axis[3], int n, double rot[3][3]);

// --- ARRAY INTEGRATION -------

/**
 *	Trapezoidal integration of function y( x ) with array of x and array of y values as input ( both of length len )
 *	Nonequidistant grid for x possible
 */
std::complex<double> arrInt( double x[], std::complex<double> y[], int len );
/**
 *	Trapezoidal integration of function y( x ) with array of y values as input
 *	Equidistant grid with spacing h for x grid assumed
 */
std::complex<double> arrInt( double h, std::complex<double> y[], int len );
/**
 *	Extended Simpson integration of function y( x ) with array of x and array of y values as input ( both of length len )
 *	Nonequidistant grid for x possible
 */
std::complex<double> arrIntSimpson( double h, std::complex<double> y[], int len );
/**
 *	Extended Simpson integration of function y( x ) with array of y values as input
 *	Equidistant grid with spacing h for x grid assumed
 */
std::complex<double> arrIntSimpson( double h, std::complex<double> y[], int lenx, int leny );




// --- Integrations using ODE solver ---------

typedef int(*intFuncP )( double, const double[], double[], void*);
typedef std::complex<double>(*intFuncPCPLX )( double, void*);

double integrate( intFuncP fw, void* q, double err );
double integrate( intFuncP fw, double err );

std::complex<double> integrateCPLX( intFuncPCPLX fw, void* q, double err );
int realFunc( double w, const double y[], double g[], void* q );
int imagFunc( double w, const double y[], double g[], void* q );

struct intCPLXList{
    intFuncPCPLX fw;
    void* q;
    intCPLXList( intFuncPCPLX fw_val, void* q_val ){fw = fw_val; q = q_val;}
};


// --- INTEGRATION QAGS ---------

typedef double(*intFuncP_QAGS )( double, void*);

double integrate_QAGS( intFuncP_QAGS fw, void* q, double err );
double integrate_QAGS( intFuncP_QAGS fw, double err );

std::complex<double> integrateCPLX( intFuncPCPLX fw, void* q, double err );
double realFunc_QAGS( double w, void* q );
double imagFunc_QAGS( double w, void* q );

// -- Oscilatory

double integrate_QAGS_COS( intFuncP_QAGS fw, void* q, double err );
double integrate_QAGS_SIN( intFuncP_QAGS fw, void* q, double err );

// generate_sum_function provides the fitting of a finite matsubara sum in order to get the extension of the sum up to infinity
// The Least-Square fitting uses the following fitting function : sum_{i=0}^{fit_order} a_i (x)^{-i}
// This requires the minimization of the chisquare which has been done in a semi-analytic way (see notes by Georg Roehringer)
// generate_sum_func returns sumfit of a given single-variable function 
std::vector<double> generate_tail_weights( int iMin, int tail_length, int fit_order ); 
gf_weight_vec_t generate_weights( int iMin, int tail_length, int fit_order ); 
gf<double, 2> generate_2d_weights( int iMin, int tail_length, int fit_order ); 

template< typename value_t >
boost::function< value_t ( boost::function< value_t ( int i ) > ) > generate_sum_func( int iMin, int tail_length, int fit_order )
{
   std::vector<double> tail_weights = generate_tail_weights( iMin, tail_length, fit_order ); 

   return [tail_weights,iMin,tail_length]( boost::function< value_t ( int i ) > func )
   {
      value_t val; 
      for( int i = -iMin; i <= iMin; ++i )
	 val += func(i); 

      for( int i = 1; i < tail_length; ++i )
	 val += ( func(iMin + i) + func(-iMin -i) ) * tail_weights[i]; 

      return val; 
   }; 
}
