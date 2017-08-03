
/************************************************************************************************//**
 *  		
 * 	file: 		mymath.cpp
 * 	contents:   	See mymath.h
 * 
 ****************************************************************************************************/


#include "mymath.h"
#include <math.h>
#include <iostream>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
//#include <gsl/gsl_odeiv2.h>
//#include <gsl/gsl_integration.h>

using std::cout; using std::endl; 

// ROTATION MATRIX //

void rotation_matrix( double axis[3], int n, double rot[3][3])
{
       const double PI = 3.14159265358979323846; 
       double theta = 2.0*PI/n;
       double vx = axis[0];
       double vy = axis[1];
       double vz = axis[2]; 
       double cos2  = cos(theta/2.0)*cos(theta/2.0);
       double sin2  = sin(theta/2.0)*sin(theta/2.0);
       
       rot[0][0] = cos2 + sin2*(2.0*vx*vx-1.0);
       rot[0][1] = 2.0*vx*vy*sin2+vz*sin(theta);
       rot[0][2] = 2.0*vx*vz*sin2-vy*sin(theta);

       rot[1][0] = 2.0*vy*vx*sin2-vz*sin(theta);	
       rot[1][1] = cos2 + sin2*(2.0*vy*vy-1.0);
       rot[1][2] = 2.0*vy*vz*sin2+vx*sin(theta);

       rot[2][0] = 2.0*vz*vx*sin2+vy*sin(theta);
       rot[2][1] = 2.0*vz*vy*sin2-vx*sin(theta);
       rot[2][2] = cos2 + sin2*(2.0*vz*vz-1.0); 
}        

//REFLECTION OPERATION WITH RESPECT TO THE XZ PLANE (Y -> -Y)

void mirror_matrix_y( double rot[3][3])
{
       rot[0][0] = 1.;
       rot[0][1] = 0.;
       rot[0][2] = 0.;

       rot[1][0] = 0.;
       rot[1][1] = -1.;
       rot[1][2] = 0.;
       
       rot[2][0] = 0.;	
       rot[2][1] = 0.;
       rot[2][2] = 1.;

}        
// CONST //

size_t wspace_size = 10000;

double kdel( double a, double b )
{
    if ( a == b ){ return 1.0;}
    return 0.0;
}

// -------------- Self Consistency ----------------
double dist( const dcomplex a[], const dcomplex b[], int len ) // maximumsnorm
{
    double err = 0.0;
    for( int i = 0; i < len; i++)
    {
        if( abs( a[i] - b[i]) > err )
        {
            err = abs( a[i] - b[i]);
        }
    }
    return err;
}

int self_consistency( dcomplex init[], dcomplex fin[], int len, sc_Func fw, double err, void * q )
{
    fw( init, fin, q );
    int cnt = 1;
    while( dist( init, fin, len ) > err && cnt < 100 )
    {
        for( int i = 0; i < len; i++)
        {
            init[i] = fin[i];
        }
        fw( init, fin, q );
        cnt++;
    }

    if ( cnt == 99 ){return -1;}
    return cnt;
}


// -------------- ARRAY INTEGRATION ---------------

dcomplex arrInt( double x[], dcomplex y[], int len )
{
    dcomplex val = 0.0;

    if ( len < 2 ){ return val; }

    val += ( y[0] * ( x[1] - x[0]) + y[len - 1] * ( x[len-1] - x[len-2])); // /2.0; !!!!!!!!!!!!!!! OPEN ENDED INTEGRATION!!

    for ( int i = 1; i < len-1; i++)
    {
        val += y[i] * ( x[i+1] - x[i-1]) / 2.0;
    }

    return val;
}

dcomplex arrInt( double h, dcomplex y[], int len )
{
    dcomplex val = 0.0;

    if ( len < 2 ){ return val; }

    val += ( y[0] + y[len - 1] )/2.0;

    for ( int i = 1; i < len-1; i++)
    {
        val += y[i] ;
    }

    return h*val;
}

dcomplex arrIntSimpson( double h, dcomplex y[], int len ) // len > 2 und Len ungerade!!!
{
    dcomplex val = 0.0;

    if ( len % 2 == 0 || len < 3 ){ return 0.0; cout << " LENGTH PROBLEM IntSimpson " << endl; }

    val += ( y[0] - y[len-1] )/ 3.0; // Korrigiert folgenden LOOP

    for ( int i = 1; i < len-1; i += 2 )
    {
        val += ( 4.0 * y[i] + 2.0 * y[i+1] ) / 3.0;
    }

    return h*val;
}

dcomplex arrIntSimpson( double h, dcomplex y[], int lenx, int leny )
{
    dcomplex val = 0.0;

    if ( lenx % 2 == 0 || lenx < 3 ){ return 0.0; cout << " LENGTH PROBLEM IntSimpson " << endl; }

    val += ( arrIntSimpson( h, y, leny ) -  arrIntSimpson( h, y + ( lenx-1 )*leny, leny ) )/ 3.0; // Korrigiert folgenden LOOP

    for ( int i = 1; i < lenx-1; i += 2 )
    {
        val += ( 4.0 * arrIntSimpson( h, y + i*leny, leny ) + 2.0 * arrIntSimpson( h, y + ( i+1 )*leny, leny ) ) / 3.0;
    }

    return h*val;
}

// -------------- INTEGRATION ODE -----------------

//double integrate( intFuncP fw, void* q, double err ) // Integrate function from 0 to infinity
//{
    //gsl_odeiv2_system sysJC = {fw, 0, 1, q};
    //gsl_odeiv2_driver * drvJC = gsl_odeiv2_driver_alloc_y_new (&sysJC, gsl_odeiv2_step_rkf45, 1e-5, err, err*1e+2 );

////    double wi1 = -DELTA_START, wf1 = 0.0, wi2 = 0.0, wf2 = DELTA_START; //1e+1*err,
    //double wi = 0.0, wf = 1e+10; //1e+1*err,
    //double z[1] = { 0.0 };

    //int status = gsl_odeiv2_driver_apply ( drvJC, &wi, wf, z );

    //if ( status != GSL_SUCCESS )
    //{
        //printf ("error, return value=%d\n", status );
    //}
    //gsl_odeiv2_driver_free ( drvJC );

    //return z[0];
//}

//double integrate( intFuncP fw, double err ) // Integrate function from 0 to infinity
//{
    //return integrate( fw, NULL, err );
//}

//dcomplex integrateCPLX( intFuncPCPLX fw, void* q, double err )
//{
    //intCPLXList intL( fw, q );
    //return dcomplex( integrate( realFunc, &intL, err ), integrate( imagFunc, &intL, err ));
    ////return complex<double>( integrate_QAGS( realFunc_QAGS, &intL, err ), integrate_QAGS( imagFunc_QAGS, &intL, err ));
    ////return complex<double>( integrate_QAGS_COS( realFunc_QAGS, &intL, err ), integrate_QAGS_SIN( imagFunc_QAGS, &intL, err )); //CAUTION, Integrates f * sin or f * cos
//}

//int realFunc( double w, const double y[], double g[], void* q )
//{
    //intCPLXList* intL = ( intCPLXList*) q;
    //g[0] = real( intL->fw( w, intL->q ) + intL->fw(-w, intL->q ));
    //return GSL_SUCCESS;
//}

//int imagFunc( double w, const double y[], double g[], void* q )
//{
    //intCPLXList* intL = ( intCPLXList*) q;
    //g[0] = imag( intL->fw( w, intL->q ) + intL->fw(-w, intL->q ));
    //return GSL_SUCCESS;
//}

// -------------- INTEGRATION QAGS -----------------

//double integrate_QAGS( intFuncP_QAGS fw, void* q, double err ) // Integrate function from 0 to infinity
//{
    //gsl_integration_workspace * w  = gsl_integration_workspace_alloc ( wspace_size );

    //double result, error;

    //gsl_function F;
    //F.function = fw;
    //F.params = q;

    //gsl_integration_qagiu(&F, 0.0, err, err*1e+2, wspace_size, w, &result, &error );
    ////gsl_integration_qags (&F, 0.0, 1e+10, err, err*1e+2, wspace_size, w, &result, &error );

    ////cout << " Error " << error << " relative " << error/result << endl;getchar();

    //gsl_integration_workspace_free ( w );

    //return result;
//}

//double integrate_QAGS( intFuncP fw, double err ) // Integrate function from 0 to infinity
//{
    //return integrate( fw, NULL, err );
//}

double realFunc_QAGS( double w, void* q )
{
    intCPLXList* intL = ( intCPLXList*) q;
    return real( intL->fw( w, intL->q ) + intL->fw(-w, intL->q ));
}

double imagFunc_QAGS( double w, void* q )
{
    intCPLXList* intL = ( intCPLXList*) q;
    return imag( intL->fw( w, intL->q ) + intL->fw(-w, intL->q ));
}

// ---------------- OSCILATORY INTEGRATION ---------------

//double integrate_QAGS_COS( intFuncP_QAGS fw, void* q, double err ) // Integrate function from 0 to infinity
//{
    //gsl_integration_workspace * w  = gsl_integration_workspace_alloc ( wspace_size );

    //double result, error;

    //gsl_function F;
    //F.function = fw;
    //F.params = q;

     ////Integrate oscilatory functions
    //gsl_integration_qawo_table * wf = gsl_integration_qawo_table_alloc ( 1e-10, 1e+10, GSL_INTEG_COSINE, ( int ) 10.0*log( wspace_size )/log( 2 ) );
    //gsl_integration_qawo (&F, 0.0, err, err*1e+2, wspace_size, w, wf, &result, &error );
    //gsl_integration_qawo_table_free ( wf );

    ////cout << " Error " << error << " relative " << error/result << endl;getchar();

    //gsl_integration_workspace_free ( w );

    //return result;
//}

//double integrate_QAGS_SIN( intFuncP_QAGS fw, void* q, double err ) // Integrate function from 0 to infinity
//{
    //gsl_integration_workspace * w  = gsl_integration_workspace_alloc ( wspace_size );

    //double result, error;

    //gsl_function F;
    //F.function = fw;
    //F.params = q;

     ////Integrate oscilatory functions
    //gsl_integration_qawo_table * wf = gsl_integration_qawo_table_alloc ( 1e-10, 1e+10, GSL_INTEG_SINE, ( int ) 10*log( wspace_size )/log( 2 ) );
    //gsl_integration_qawo (&F, 0.0, err, err*1e+2, wspace_size, w, wf, &result, &error );
    //gsl_integration_qawo_table_free ( wf );

    ////cout << " Error " << error << " relative " << error/result << endl;getchar();

    //gsl_integration_workspace_free ( w );

    //return result;
//}

std::vector<double> generate_tail_weights( int iMin, int tail_length, int fit_order )
{
   Eigen::MatrixXd MatFit(fit_order,fit_order); 
   MatFit.setZero(); 

   int iMax = tail_length + iMin; 

   for( int i = 0; i < fit_order; ++i )
      for( int j = 0; j < fit_order; ++j )
	 for( int k = iMin; k < iMax; ++k )
	    MatFit(i,j) += 1.0 / pow(k, i + j); 

   Eigen::MatrixXd MatInv = MatFit.inverse(); 

   std::vector<double> w_array(tail_length,0.0); 
   for( int i = 0; i < tail_length; ++i )
      for( int l = 0; l < fit_order; ++l )
	 w_array[i] += MatInv(0,l) / pow(i + iMin, l); 

   std::vector<double> tail_weights(tail_length,0.0);
   for( int i = 0; i < tail_length; ++i )
      for( int l = i; l < tail_length; ++l )
	 tail_weights[i] += w_array[l]; 

   return tail_weights; 
}

gf_weight_vec_t generate_weights( int iMin, int tail_length, int fit_order )
{
   std::vector<double> tail_weights = generate_tail_weights( iMin, tail_length, fit_order ); 

   gf_weight_vec_t weight_vec; 
   weight_vec.init( []( const gf_weight_vec_t::idx_t& idx ){ return 1.0; } ); 

   for( int i = 0; i < tail_length; ++i )
   {
      weight_vec[iMin + i] = tail_weights[i]; 
      weight_vec[-iMin - i - 1] = tail_weights[i]; 
   }

   return weight_vec; 
}

gf<double, 2> generate_2d_weights( int iMin, int tail_length, int fit_order )
{
   std::vector<double> tail_weights = generate_tail_weights( iMin, tail_length, fit_order ); 

   gf<double, 2> weight_vec( boost::extents[ffreq(POS_INT_RANGE)][ffreq(POS_INT_RANGE)] ); 
   weight_vec.init( []( const gf<double, 2>::idx_t& idx ){ return 1.0; } ); 

   for( int i = 0; i < tail_length; ++i )
   {
      for( int j = -iMin - i - 1; j <= iMin + i; ++j )
      {
	 weight_vec[iMin + i][j] = tail_weights[i]; 
	 weight_vec[-iMin - i - 1][j] = tail_weights[i]; 
	 weight_vec[j][iMin + i] = tail_weights[i]; 
	 weight_vec[j][-iMin - i - 1] = tail_weights[i]; 
      }
   }

   return weight_vec; 
}
