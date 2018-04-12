
/****************************************************************************************************
 *  		
 * 	file: 		grid.cpp
 * 	contents:   	See grid.h
 * 
 ****************************************************************************************************/


#include <grid.h>
#include <fstream>
#include <string>
#include <const.h>
#include <tuple>

using namespace std; 

/********************* F_Grid Class  ********************/

F_Grid::F_Grid( unsigned int pos_freq_count_, double step_size_ ) :
   pos_freq_count( pos_freq_count_ ), step_size( step_size_ )
{
   for( int i = -pos_freq_count; i < pos_freq_count ; ++i )
      grid_points.push_back( step_size * i + step_size / 2.0 );
}

void F_Grid::print( string fname )
{
   fstream data;
   data.open("dat/" + fname, ios::out | ios::trunc );
   for( double freq : grid_points )
      data << freq << endl;
   data.flush();
   data.close();
}

unsigned int F_Grid::get_pos_freq_count() const
{
   return pos_freq_count;
}

double F_Grid::get_step_size() const
{
   return step_size;
}

/********************* Bos_Grid Class  ********************/

Bos_Grid::Bos_Grid( unsigned int pos_freq_count_, double step_size_ ) :
   pos_freq_count( pos_freq_count_ ), step_size( step_size_ )
{
   for( int i = -pos_freq_count; i < pos_freq_count + 1 ; ++i )
      grid_points.push_back( step_size * i );
}

void Bos_Grid::print( string fname )
{
   fstream data;
   data.open("dat/" + fname, ios::out | ios::trunc );
   for( double freq : grid_points )
      data << freq << endl;
   data.flush();
   data.close();
}

unsigned int Bos_Grid::get_pos_freq_count() const
{
   return pos_freq_count;
}

double Bos_Grid::get_step_size() const
{
   return step_size;
}

R_Grid::R_Grid()
{
//   cout<< "beginning K_Grid()"<<endl;
   if(NN_COUNT == 0)
      push_back( tuple<int, int>(0,0));
   else if(NN_COUNT == 1){
      push_back( tuple<int, int>(0,0)); 
      push_back( tuple<int, int>(-1,0));
      push_back( tuple<int, int>(1,0));
      push_back( tuple<int, int>(0,-1));
      push_back( tuple<int, int>(0,1));
      //Until here FFACTOR_COUNT
      push_back( tuple<int, int>(1,1)); 
      push_back( tuple<int, int>(-1,1));
      push_back( tuple<int, int>(1,-1));
      push_back( tuple<int, int>(-1,-1));
      //Next cube
      push_back( tuple<int, int>(-2,0));
      push_back( tuple<int, int>(2,0));
      push_back( tuple<int, int>(0,-2));
      push_back( tuple<int, int>(0,2));
   }
}

#ifdef UNIFORM_GRID

K_Grid::K_Grid()
{
//   cout<< "beginning K_Grid()"<<endl;
   for( int kx =-MAX_KPOS+1; kx < MAX_KPOS+1 ; ++kx ){
	for( int ky =-MAX_KPOS+1; ky < MAX_KPOS+1; ++ky ){
		push_back( tuple<double, double>( PI*kx/MAX_KPOS,PI*ky/MAX_KPOS));
	}
   }
}

int K_Grid::k_origin[2] = {0,0};
unsigned int K_Grid::patch_origin = K_Grid::get_patch(k_origin);

gf_k_sum_matrix_t K_Grid::add_k_matrix(K_Grid().generate_add_k_matrix());
gf_k_sum_matrix_t K_Grid::add_k_bfx_matrix(K_Grid().generate_add_k_bfx_matrix());
gf_k_sum_matrix_t K_Grid::add_k_bfy_matrix(K_Grid().generate_add_k_bfy_matrix());
gf_k_sum_matrix_t K_Grid::dif_k_matrix(K_Grid().generate_dif_k_matrix());
gf_k_sum_matrix_t K_Grid::dif_k_bfx_matrix(K_Grid().generate_dif_k_bfx_matrix());
gf_k_sum_matrix_t K_Grid::dif_k_bfy_matrix(K_Grid().generate_dif_k_bfy_matrix());

gf_k_rot_matrix_t K_Grid::rot_pi2_z_matrix(K_Grid().generate_rot_pi2_z_matrix());
gf_k_rot_matrix_t K_Grid::rot_pi2_z_bfx_matrix(K_Grid().generate_rot_pi2_z_bfx_matrix());
gf_k_rot_matrix_t K_Grid::rot_pi2_z_bfy_matrix(K_Grid().generate_rot_pi2_z_bfy_matrix());

gf_k_rot_matrix_t K_Grid::mirror_y_matrix(K_Grid().generate_mirror_y_matrix());
gf_k_rot_matrix_t K_Grid::mirror_y_bfx_matrix(K_Grid().generate_mirror_y_bfx_matrix());
gf_k_rot_matrix_t K_Grid::mirror_y_bfy_matrix(K_Grid().generate_mirror_y_bfy_matrix());

gf_k_rot_matrix_t K_Grid::mirror_diagonal_matrix(K_Grid().generate_mirror_diagonal_matrix());  //backfolding never needed!

gf_k_4MAXKPOS_matrix_t K_Grid::COS(K_Grid().generate_cos_matrix());
gf_k_4MAXKPOS_matrix_t K_Grid::SIN(K_Grid().generate_sin_matrix());

unsigned int K_Grid::get_patch(int ind[2] )
{
   int x_patch, y_patch; 
   x_patch = (ind[0]+MAX_KPOS-1+1000*MAX_KPOS)%(2*MAX_KPOS);
   y_patch = (ind[1]+MAX_KPOS-1+1000*MAX_KPOS)%(2*MAX_KPOS);

   return x_patch*2*MAX_KPOS+y_patch;      
}

unsigned int K_Grid::get_bf(int ind )
{
   return abs((2*ind-1) / (2*MAX_KPOS))%2;      
   //return (ind+MAX_KPOS-1) / (2*MAX_KPOS);      
}

std::tuple<int, int> K_Grid::get_indices(double kx,double ky)
{
   return tuple<int, int>(std::round(kx*MAX_KPOS/PI),std::round(ky*MAX_KPOS/PI));
}

std::tuple<int, int> K_Grid::get_indices(unsigned int patch)
{
   int x = (patch / (2*MAX_KPOS))  -MAX_KPOS +1;
   int y = (patch % (2*MAX_KPOS))  -MAX_KPOS +1;
   return tuple<int, int>(x,y);
}

std::tuple<int, int> K_Grid::get_indices_enum(unsigned int patch)
{
   int x = (patch / (2*MAX_KPOS));
   int y = (patch % (2*MAX_KPOS));
   return tuple<int, int>(x,y);
}

gf_k_4MAXKPOS_matrix_t K_Grid::generate_cos_matrix(){
   cout <<"in generate cos matrix"<< endl;
   gf_k_4MAXKPOS_matrix_t matrix;
   for( int i= 0; i< 4*MAX_KPOS; ++i){
      matrix[i] = cos(0.5 * double(i * PI/MAX_KPOS));
      if (abs(matrix[i])<CHOP_ERR) matrix[i]=0;
   } cout << "End of generate cos k:"<< endl;    
   return matrix;
}
   
gf_k_4MAXKPOS_matrix_t K_Grid::generate_sin_matrix(){
     cout <<"in generate sin matrix"<< endl;
     gf_k_4MAXKPOS_matrix_t matrix;
     for( int i= 0; i< 4*MAX_KPOS; ++i){
    	   matrix[i] = sin(0.5 * double(i * PI/MAX_KPOS));
    	   if (abs(matrix[i])<CHOP_ERR) matrix[i]=0;
       } cout << "End of generate sin k:"<< endl;    
     return matrix;
  }
  
//const double K_Grid::COS(K_Grid().get_cosine());  // cosine( 1/2 * Q_x * PI/MAX_KPOS ) starting from 0 up to 4*MAX_KPOS-1
//const double K_Grid::SIN(K_Grid().get_sine());  // sine( 1/2 * Q_x * PI/MAX_KPOS )     !pay attention if input becomes minus!

int K_Grid::get_add_k(unsigned int i, unsigned int j)
{
      int xi = std::get<0>(get_indices(i));
      int yi = std::get<1>(get_indices(i));
      int xj = std::get<0>(get_indices(j));
      int yj = std::get<1>(get_indices(j));
      int vec_add[2] = {xi+xj, yi+yj};
      return  get_patch(vec_add);
}  
int K_Grid::get_add_k_bfx(unsigned int i, unsigned int j)
{
      int xi = std::get<0>(get_indices(i));
      int xj = std::get<0>(get_indices(j));
      return get_bf(xi+xj);
}
int K_Grid::get_add_k_bfy(unsigned int i, unsigned int j)
{
      int yi = std::get<1>(get_indices(i));
      int yj = std::get<1>(get_indices(j));
      return get_bf(yi+yj);
}

int K_Grid::get_dif_k(unsigned int i, unsigned int j)
{
      int xi = std::get<0>(get_indices(i));
      int yi = std::get<1>(get_indices(i));
      int xj = std::get<0>(get_indices(j));
      int yj = std::get<1>(get_indices(j));
      int vec_dif[2] = {xi-xj, yi-yj};

      return get_patch(vec_dif);
}
int K_Grid::get_dif_k_bfx(unsigned int i, unsigned int j)
{
      int xi = std::get<0>(get_indices(i));
      int xj = std::get<0>(get_indices(j));
      return get_bf(xi-xj);
}
int K_Grid::get_dif_k_bfy(unsigned int i, unsigned int j)
{
      int yi = std::get<1>(get_indices(i));
      int yj = std::get<1>(get_indices(j));
      return get_bf(yi-yj);
}

unsigned int K_Grid::get_rot_k(int n, unsigned int p)
{
   double rotation[2][2];
   int  vec_rot[2];
   std::tuple<int, int> indices;
   
   rotation_matrix_2d( n, rotation);
   
   indices = get_indices(p);
   for(int j = 0; j<2; ++j){
      vec_rot[j] = std::round(rotation[j][0]*(std::get<0>(indices))
			     +rotation[j][1]*(std::get<1>(indices)));
   }
   return get_patch(vec_rot);
}
unsigned int K_Grid::get_rot_k_bfx(int n, unsigned int p)
{
   double rotation[2][2];
   int  vec_rot[2];
   std::tuple<int, int> indices;
   
   rotation_matrix_2d( n, rotation);
   
   indices = get_indices(p);
   for(int j = 0; j<2; ++j){
      vec_rot[j] = std::round(rotation[j][0]*(std::get<0>(indices))
			     +rotation[j][1]*(std::get<1>(indices)));
   }
   return get_bf(vec_rot[0]);
}
unsigned int K_Grid::get_rot_k_bfy(int n, unsigned int p)
{
   double rotation[2][2];
   int  vec_rot[2];
   std::tuple<int, int> indices;
   
   rotation_matrix_2d( n, rotation);
   
   indices = get_indices(p);
   for(int j = 0; j<2; ++j){
      vec_rot[j] = std::round(rotation[j][0]*(std::get<0>(indices))
			     +rotation[j][1]*(std::get<1>(indices)));
   }
   return get_bf(vec_rot[1]);
}

unsigned int K_Grid::get_mirror_k(unsigned int p)
{
   double rotation[2][2];
   int  vec_rot[2];
   std::tuple<int, int> indices;
   
   mirror_matrix_y(rotation);
   
   indices = get_indices(p);
   for(int j = 0; j<2; ++j){
      vec_rot[j] = std::round(rotation[j][0]*(std::get<0>(indices))
			     +rotation[j][1]*(std::get<1>(indices)));
   }
   return get_patch(vec_rot);
}
unsigned int K_Grid::get_mirror_k_bfx(unsigned int p)
{
   double rotation[2][2];
   int  vec_rot[2];
   std::tuple<int, int> indices;
   
   mirror_matrix_y(rotation);
   
   indices = get_indices(p);
   for(int j = 0; j<2; ++j){
      vec_rot[j] = std::round(rotation[j][0]*(std::get<0>(indices))
			     +rotation[j][1]*(std::get<1>(indices)));
   }
   return get_bf(vec_rot[0]);
}
unsigned int K_Grid::get_mirror_k_bfy(unsigned int p)
{
   double rotation[2][2];
   int  vec_rot[2];
   std::tuple<int, int> indices;
   
   mirror_matrix_y(rotation);
   
   indices = get_indices(p);
   for(int j = 0; j<2; ++j){
      vec_rot[j] = std::round(rotation[j][0]*(std::get<0>(indices))
			     +rotation[j][1]*(std::get<1>(indices)));
   }
   return get_bf(vec_rot[1]);
}

unsigned int K_Grid::get_mirror_diagonal_k(unsigned int p)
{
   double rotation[2][2];
   int  vec_rot[2];
   std::tuple<int, int> indices;
   
   mirror_matrix_diagonal(rotation);
   
   indices = get_indices(p);
   for(int j = 0; j<2; ++j){
      vec_rot[j] = std::round(rotation[j][0]*(std::get<0>(indices))
			     +rotation[j][1]*(std::get<1>(indices)));
   }
   return get_patch(vec_rot);
}
/**********************GENERATION FUNCTIONS FOR ADD_K_MATRIX AND DIF_K_MATRIX ********************************/

gf_k_sum_matrix_t K_Grid::generate_add_k_matrix()
{  
   cout <<"in generate add k matrix"<< endl;
   gf_k_sum_matrix_t add_k_matrix;
   for( int i= 0; i< PATCH_COUNT; ++i){
      for( int k= 0; k< PATCH_COUNT; ++k){
	   add_k_matrix[i][k] = get_add_k(i,k);
      }
   } cout << "End of generate add_k:"<< endl;    
     return add_k_matrix;
}
gf_k_sum_matrix_t K_Grid::generate_add_k_bfx_matrix()
{
   cout <<"in generate add k bfx matrix"<< endl;
   gf_k_sum_matrix_t add_k_bfx_matrix;
   for( int i= 0; i< PATCH_COUNT; ++i){
      for( int k= 0; k< PATCH_COUNT; ++k){
	   add_k_bfx_matrix[i][k] = get_add_k_bfx(i,k);
      }
   } cout << "End of generate add_k_bfx:"<< endl; 
     return add_k_bfx_matrix;
}
gf_k_sum_matrix_t K_Grid::generate_add_k_bfy_matrix()
{
   cout <<"in generate add k bfy matrix"<< endl;
   gf_k_sum_matrix_t add_k_bfy_matrix;
   for( int i= 0; i< PATCH_COUNT; ++i){
      for( int k= 0; k< PATCH_COUNT; ++k){
	   add_k_bfy_matrix[i][k] = get_add_k_bfy(i,k);
      }
   } cout << "End of generate add_k_bfy:"<< endl; 
     return add_k_bfy_matrix;
}

gf_k_sum_matrix_t K_Grid::generate_dif_k_matrix()
{
   cout <<"in generate dif k matrix"<< endl;
   gf_k_sum_matrix_t dif_k_matrix;
   for( int i= 0; i< PATCH_COUNT; ++i){
      for( int k= 0; k< PATCH_COUNT; ++k){
	   dif_k_matrix[i][k] = get_dif_k(i,k);
      }
   } cout << "End of generate dif_k:"<< endl; 
     return dif_k_matrix;
}
gf_k_sum_matrix_t K_Grid::generate_dif_k_bfx_matrix()
{
   cout <<"in generate dif k bfx matrix"<< endl;
   gf_k_sum_matrix_t dif_k_bfx_matrix;
   for( int i= 0; i< PATCH_COUNT; ++i){
      for( int k= 0; k< PATCH_COUNT; ++k){
	   dif_k_bfx_matrix[i][k] = get_dif_k_bfx(i,k);
      }
   } cout << "End of generate dif_k_bfx:"<< endl; 
     return dif_k_bfx_matrix;
}
gf_k_sum_matrix_t K_Grid::generate_dif_k_bfy_matrix()
{
   cout <<"in generate dif k bfy matrix"<< endl;
   gf_k_sum_matrix_t dif_k_bfy_matrix;
   for( int i= 0; i< PATCH_COUNT; ++i){
      for( int k= 0; k< PATCH_COUNT; ++k){
	   dif_k_bfy_matrix[i][k] = get_dif_k_bfy(i,k);
      }
   } cout << "End of generate dif_k_bfy:"<< endl; 
     return dif_k_bfy_matrix;
}

//Rotation functions around the main axes
//clockwise rotation

//Rotation of the plane XY of pi/2 around the z axis
gf_k_rot_matrix_t K_Grid::generate_rot_pi2_z_matrix()
{
   cout <<"in generate rot pi2 z matrix"<< endl;
   gf_k_rot_matrix_t rot_pi2_z_matrix;
   for(unsigned int i= 0; i< PATCH_COUNT; ++i){
      rot_pi2_z_matrix[i] = get_rot_k(4,i);
   }  
   return rot_pi2_z_matrix;
}
gf_k_rot_matrix_t K_Grid::generate_rot_pi2_z_bfx_matrix()
{
   cout <<"in generate rot pi2 z bfx matrix"<< endl;
   gf_k_rot_matrix_t rot_pi2_z_bfx_matrix;
   for(unsigned int i= 0; i< PATCH_COUNT; ++i){
      rot_pi2_z_bfx_matrix[i] = get_rot_k_bfx(4,i);
   }  
   return rot_pi2_z_bfx_matrix;
}
gf_k_rot_matrix_t K_Grid::generate_rot_pi2_z_bfy_matrix()
{
   cout <<"in generate rot pi2 z bfy matrix"<< endl;
   gf_k_rot_matrix_t rot_pi2_z_bfy_matrix;
   for(unsigned int i= 0; i< PATCH_COUNT; ++i){
      rot_pi2_z_bfy_matrix[i] = get_rot_k_bfy(4,i);
   }  
   return rot_pi2_z_bfy_matrix;
}


// Reflection of the XY plane around the yz plane
gf_k_rot_matrix_t K_Grid::generate_mirror_y_matrix()
{
   cout <<"in generate mirror y matrix"<< endl;
   gf_k_rot_matrix_t mirror_y_matrix;
   for(unsigned int i= 0; i< PATCH_COUNT; ++i){
      mirror_y_matrix[i] = get_mirror_k(i);
   }  
   return mirror_y_matrix;
}
gf_k_rot_matrix_t K_Grid::generate_mirror_y_bfx_matrix()
{
   cout <<"in generate mirror y bfx matrix"<< endl;
   gf_k_rot_matrix_t mirror_y_bfx_matrix;
   for(unsigned int i= 0; i< PATCH_COUNT; ++i){
      mirror_y_bfx_matrix[i] = get_mirror_k_bfx(i);
   }  
   return mirror_y_bfx_matrix;
}
gf_k_rot_matrix_t K_Grid::generate_mirror_y_bfy_matrix()
{
   cout <<"in generate mirror y bfy matrix"<< endl;
   gf_k_rot_matrix_t mirror_y_bfy_matrix;
   for(unsigned int i= 0; i< PATCH_COUNT; ++i){
      mirror_y_bfy_matrix[i] = get_mirror_k_bfy(i);
   }  
   return mirror_y_bfy_matrix;
}

// Reflection of the XY plane at the diagonal (x<->y)
gf_k_rot_matrix_t K_Grid::generate_mirror_diagonal_matrix()
{
   cout <<"in generate mirror diagonal matrix"<< endl;
   gf_k_rot_matrix_t mirror_y_matrix;
   for(unsigned int i= 0; i< PATCH_COUNT; ++i){
      mirror_y_matrix[i] = get_mirror_diagonal_k(i);
   }  
   return mirror_y_matrix;
}

#endif
#ifdef SIMPLE_PATCH
K_Grid::K_Grid()
{
      push_back( pair<double, double>( 0.0, 	0.0 ));		// PATCH 0
      push_back( pair<double, double>( PI/2.0, 	PI/2.0 ));	// PATCH 1
      push_back( pair<double, double>( PI/2.0, 	-PI/2.0 ));	// PATCH 2
      push_back( pair<double, double>(-PI/2.0, 	-PI/2.0 ));	// PATCH 3
      push_back( pair<double, double>(-PI/2.0, 	PI/2.0 ));	// PATCH 4
      push_back( pair<double, double>( 0.0, 	PI ));		// PATCH 5
      push_back( pair<double, double>( PI, 	PI ));		// PATCH 6
      push_back( pair<double, double>( PI, 	0.0 ));		// PATCH 7
      patch_count = size();
}

#endif

unsigned int K_Grid::get_patch_count()
{
   return patch_count;
}

void K_Grid::print( string fname )
{
      fstream data;
         data.open("dat/" + fname, ios::out | ios::trunc );
	    data << "# kx\t\tky" << endl;         
	       for( tuple<double, double> mom : *this )
		        data << std::get<0>(mom) << "\t\t" << std::get<1>(mom) << endl;
	          data.flush();
		     data.close();
}


/******************* FORM FACTORS***********************/
//Defined in grid.h class F_factors -> TODO: to be updated in the implementation of next-nearest neighbors



F_factors::F_factors()
    {
       //LOCAL FORM FACTORS
        if(NN_COUNT == 0){
//	   cout << "Just s-wave" << endl;
	   push_back([](double kx, double ky){return 1./sqrt(4*PI*PI);}); //s-wave
	   //parity_arr[0] = 1.0;
	}
       //FIRST NEAREST NEIGHBOR
	if(NN_COUNT == 1){
	   //cout << " s-wave and p-wave" << endl;
	   push_back([](double kx, double ky){return 1./sqrt(4*PI*PI);}); //s-wave}      
           push_back([](double kx, double ky){return 1./sqrt(2*PI*PI)*cos(kx);});
           push_back([](double kx, double ky){return 1./sqrt(2*PI*PI)*cos(ky);});
           //push_back([](double kx, double ky){return 1./sqrt(2*PI*PI)*sin(kx);});
           //push_back([](double kx, double ky){return 1./sqrt(2*PI*PI)*sin(ky);});
    }
        ffactor_count = size();
    }


unsigned int F_factors::get_ffactor_count()
{
   return ffactor_count;
}

#ifdef LOCAL
int F_factors::parity_arr[FFACTOR_COUNT] = {1};
int F_factors::translate_2pi_arr[FFACTOR_COUNT][2]= { 1, 1} ;
//#elif defined FIRST_NN
//int F_factors::parity_arr[FFACTOR_COUNT] = {1,1,1,-1,-1};
//int F_factors::translate_2pi_arr[FFACTOR_COUNT][2] = { { 1, 1}, {-1, 1 },  { 1,-1}, {-1, 1}, { 1,-1} }; 
#elif defined FIRST_NN // USED FOR THE s+d form-factors
int F_factors::parity_arr[FFACTOR_COUNT] = {1,1,1};
int F_factors::translate_2pi_arr[FFACTOR_COUNT][2] = { { 1, 1}, {-1, 1 },  { 1,-1} };

#endif

constexpr int F_factors::FF_rot_pi2_z_arr[FFACTOR_COUNT]; 
constexpr int F_factors::FF_mirror_y_arr[FFACTOR_COUNT]; 
constexpr int F_factors::FF_mirror_diagonal_arr[FFACTOR_COUNT]; 
constexpr int F_factors::FF_rot_sign_pi2_z_arr[FFACTOR_COUNT]; 
constexpr int F_factors::FF_mirror_sign_y_arr[FFACTOR_COUNT]; 

F_factors_real::F_factors_real()
    {
#ifdef LOCAL
       //LOCAL FORM FACTORS
        if(NN_COUNT == 0){
//	   cout << "Just s-wave" << endl;
	   push_back(std::array<std::array<dcomplex,FFREAL_DIM>,FFREAL_DIM>{{{1}}}); //s-wave
	}
       //FIRST NEAREST NEIGHBOR
#elif defined FIRST_NN
	if(NN_COUNT == 1){
//	   cout << " s-wave and p-wave" << endl;
	   push_back(array<array<dcomplex,FFREAL_DIM>,FFREAL_DIM>{{{0,0,0},{0,1,0},{0,0,0}}}); // s-wave
	   push_back(array<array<dcomplex,FFREAL_DIM>,FFREAL_DIM>{{{0,1./sqrt(2.),0},{0,0,0},{0,1./sqrt(2.),0}}}); // cos x
	   push_back(array<array<dcomplex,FFREAL_DIM>,FFREAL_DIM>{{{0,0,0},{1./sqrt(2.),0,1./sqrt(2.)},{0,0,0}}}); // cos y
	   //push_back(array<array<dcomplex,FFREAL_DIM>,FFREAL_DIM>{{{0,-I/sqrt(2.),0},{0,0,0},{0,I/sqrt(2.),0}}}); // sin x
	   //push_back(array<array<dcomplex,FFREAL_DIM>,FFREAL_DIM>{{{0,0,0},{-I/sqrt(2.),0,I/sqrt(2.)},{0,0,0}}}); // sin y
    }
#endif
        ffactor_real_count = size();
    }


unsigned int F_factors_real::get_ffactor_real_count()
{
   return ffactor_real_count;
}

gf_fft_fxf_weight_t F_factors_real::fft_fxf_weight(F_factors_real().generate_fft_fxf_weight());


gf_fft_fxf_weight_t F_factors_real::generate_fft_fxf_weight()
{
 gf_fft_fxf_weight_t fft_fxf_weight;
 for(int R=0; R<REAL_GRID;++R)
 {
      int rx = std::get<0>(R_Grid()[R]);
      int ry = std::get<1>(R_Grid()[R]);
    for(int m=0; m< FFACTOR_COUNT; ++m)
       for(int n=0; n< FFACTOR_COUNT; ++n){
	  dcomplex weight(0.0,0.0);
	  for(int rp=0; rp< REAL_GRID_FF_SHELL; ++rp){
	     int rpx = std::get<0>(R_Grid()[rp]);
	     int rpy = std::get<1>(R_Grid()[rp]);
	     int diffx = rpx-rx;
   	     int diffy = rpy-ry;
	     if(abs(diffx)<= (FFREAL_DIM-1)/2 && abs(diffy)<= (FFREAL_DIM-1)/2){
		weight += conj(ffactor_real(m, diffx, diffy))*ffactor_real(n, rpx, rpy);
	     }
	  }

	  fft_fxf_weight[R][m][n] = weight;
       }
 }
 return fft_fxf_weight;
}
