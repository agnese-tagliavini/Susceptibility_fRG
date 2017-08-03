
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
      push_back( tuple<int, int, int>(0,0,0));
   else if(NN_COUNT == 1){
      push_back( tuple<int, int, int>(0,0,0)); 
      push_back( tuple<int, int, int>(-1,0,0));
      push_back( tuple<int, int, int>(1,0,0));
      push_back( tuple<int, int, int>(0,-1,0));
      push_back( tuple<int, int, int>(0,1,0));
      push_back( tuple<int, int, int>(0,0,-1));
      push_back( tuple<int, int, int>(0,0,1));      
      //Until here REAL_GRID_FF_SHELL(=FFACTOR_COUNT normally)
      push_back( tuple<int, int, int>(0,1,1)); 
      push_back( tuple<int, int, int>(0,-1,1));
      push_back( tuple<int, int, int>(0,1,-1));
      push_back( tuple<int, int, int>(0,-1,-1));
      push_back( tuple<int, int, int>(1,1,0));
      push_back( tuple<int, int, int>(-1,1,0));
      push_back( tuple<int, int, int>(1,-1,0));
      push_back( tuple<int, int, int>(-1,-1,0));      
      push_back( tuple<int, int, int>(1,0,1));
      push_back( tuple<int, int, int>(-1,0,1));
      push_back( tuple<int, int, int>(1,0,-1));
      push_back( tuple<int, int, int>(-1,0,-1));
      //Next cube
      push_back( tuple<int, int, int>(-2,0,0));
      push_back( tuple<int, int, int>(2,0,0));
      push_back( tuple<int, int, int>(0,-2,0));
      push_back( tuple<int, int, int>(0,2,0));
      push_back( tuple<int, int, int>(0,0,-2));
      push_back( tuple<int, int, int>(0,0,2));      
   }
}

#ifdef UNIFORM_GRID

K_Grid::K_Grid()
{
//   cout<< "beginning K_Grid()"<<endl;
   for( int kx =-MAX_KPOS+1; kx < MAX_KPOS+1 ; ++kx ){
	for( int ky =-MAX_KPOS+1; ky < MAX_KPOS+1; ++ky ){
		for( int kz =-MAX_KPOS+1; kz < MAX_KPOS+1; ++kz ){ 
		   	push_back( tuple<double, double, double>( PI*kx/MAX_KPOS,PI*ky/MAX_KPOS, PI*kz/MAX_KPOS));
		}
	}
   }
}

int K_Grid::k_origin[3] = {0,0,0};
unsigned int K_Grid::patch_origin = K_Grid::get_patch(k_origin);

gf_k_sum_matrix_t K_Grid::add_k_matrix(K_Grid().generate_add_k_matrix());
gf_k_sum_matrix_t K_Grid::add_k_bfx_matrix(K_Grid().generate_add_k_bfx_matrix());
gf_k_sum_matrix_t K_Grid::add_k_bfy_matrix(K_Grid().generate_add_k_bfy_matrix());
gf_k_sum_matrix_t K_Grid::add_k_bfz_matrix(K_Grid().generate_add_k_bfz_matrix());

gf_k_sum_matrix_t K_Grid::dif_k_matrix(K_Grid().generate_dif_k_matrix());
gf_k_sum_matrix_t K_Grid::dif_k_bfx_matrix(K_Grid().generate_dif_k_bfx_matrix());
gf_k_sum_matrix_t K_Grid::dif_k_bfy_matrix(K_Grid().generate_dif_k_bfy_matrix());
gf_k_sum_matrix_t K_Grid::dif_k_bfz_matrix(K_Grid().generate_dif_k_bfz_matrix());

gf_k_rot_matrix_t K_Grid::rot_pi2_x_matrix(K_Grid().generate_rot_pi2_x_matrix());
gf_k_rot_matrix_t K_Grid::rot_pi2_x_bfx_matrix(K_Grid().generate_rot_pi2_x_bfx_matrix());
gf_k_rot_matrix_t K_Grid::rot_pi2_x_bfy_matrix(K_Grid().generate_rot_pi2_x_bfy_matrix());
gf_k_rot_matrix_t K_Grid::rot_pi2_x_bfz_matrix(K_Grid().generate_rot_pi2_x_bfz_matrix());
gf_k_rot_matrix_t K_Grid::rot_pi2_y_matrix(K_Grid().generate_rot_pi2_y_matrix());
gf_k_rot_matrix_t K_Grid::rot_pi2_y_bfx_matrix(K_Grid().generate_rot_pi2_y_bfx_matrix());
gf_k_rot_matrix_t K_Grid::rot_pi2_y_bfy_matrix(K_Grid().generate_rot_pi2_y_bfy_matrix());
gf_k_rot_matrix_t K_Grid::rot_pi2_y_bfz_matrix(K_Grid().generate_rot_pi2_y_bfz_matrix());
gf_k_rot_matrix_t K_Grid::rot_pi2_z_matrix(K_Grid().generate_rot_pi2_z_matrix());
gf_k_rot_matrix_t K_Grid::rot_pi2_z_bfx_matrix(K_Grid().generate_rot_pi2_z_bfx_matrix());
gf_k_rot_matrix_t K_Grid::rot_pi2_z_bfy_matrix(K_Grid().generate_rot_pi2_z_bfy_matrix());
gf_k_rot_matrix_t K_Grid::rot_pi2_z_bfz_matrix(K_Grid().generate_rot_pi2_z_bfz_matrix());

gf_k_rot_matrix_t K_Grid::rot_2pi3_111_matrix(K_Grid().generate_rot_2pi3_111_matrix());
gf_k_rot_matrix_t K_Grid::rot_2pi3_111_bfx_matrix(K_Grid().generate_rot_2pi3_111_bfx_matrix());
gf_k_rot_matrix_t K_Grid::rot_2pi3_111_bfy_matrix(K_Grid().generate_rot_2pi3_111_bfy_matrix());
gf_k_rot_matrix_t K_Grid::rot_2pi3_111_bfz_matrix(K_Grid().generate_rot_2pi3_111_bfz_matrix());
gf_k_rot_matrix_t K_Grid::rot_2pi3_1m11_matrix(K_Grid().generate_rot_2pi3_1m11_matrix());
gf_k_rot_matrix_t K_Grid::rot_2pi3_1m11_bfx_matrix(K_Grid().generate_rot_2pi3_1m11_bfx_matrix());
gf_k_rot_matrix_t K_Grid::rot_2pi3_1m11_bfy_matrix(K_Grid().generate_rot_2pi3_1m11_bfy_matrix());
gf_k_rot_matrix_t K_Grid::rot_2pi3_1m11_bfz_matrix(K_Grid().generate_rot_2pi3_1m11_bfz_matrix());
gf_k_rot_matrix_t K_Grid::rot_2pi3_m111_matrix(K_Grid().generate_rot_2pi3_m111_matrix());
gf_k_rot_matrix_t K_Grid::rot_2pi3_m111_bfx_matrix(K_Grid().generate_rot_2pi3_m111_bfx_matrix());
gf_k_rot_matrix_t K_Grid::rot_2pi3_m111_bfy_matrix(K_Grid().generate_rot_2pi3_m111_bfy_matrix());
gf_k_rot_matrix_t K_Grid::rot_2pi3_m111_bfz_matrix(K_Grid().generate_rot_2pi3_m111_bfz_matrix());
gf_k_rot_matrix_t K_Grid::rot_2pi3_m1m11_matrix(K_Grid().generate_rot_2pi3_m1m11_matrix());
gf_k_rot_matrix_t K_Grid::rot_2pi3_m1m11_bfx_matrix(K_Grid().generate_rot_2pi3_m1m11_bfx_matrix());
gf_k_rot_matrix_t K_Grid::rot_2pi3_m1m11_bfy_matrix(K_Grid().generate_rot_2pi3_m1m11_bfy_matrix());
gf_k_rot_matrix_t K_Grid::rot_2pi3_m1m11_bfz_matrix(K_Grid().generate_rot_2pi3_m1m11_bfz_matrix());


gf_k_rot_matrix_t K_Grid::rot_pi_110_matrix(K_Grid().generate_rot_pi_110_matrix());
gf_k_rot_matrix_t K_Grid::rot_pi_110_bfx_matrix(K_Grid().generate_rot_pi_110_bfx_matrix());
gf_k_rot_matrix_t K_Grid::rot_pi_110_bfy_matrix(K_Grid().generate_rot_pi_110_bfy_matrix());
gf_k_rot_matrix_t K_Grid::rot_pi_110_bfz_matrix(K_Grid().generate_rot_pi_110_bfz_matrix());
gf_k_rot_matrix_t K_Grid::rot_pi_1m10_matrix(K_Grid().generate_rot_pi_1m10_matrix());
gf_k_rot_matrix_t K_Grid::rot_pi_1m10_bfx_matrix(K_Grid().generate_rot_pi_1m10_bfx_matrix());
gf_k_rot_matrix_t K_Grid::rot_pi_1m10_bfy_matrix(K_Grid().generate_rot_pi_1m10_bfy_matrix());
gf_k_rot_matrix_t K_Grid::rot_pi_1m10_bfz_matrix(K_Grid().generate_rot_pi_1m10_bfz_matrix());
gf_k_rot_matrix_t K_Grid::rot_pi_011_matrix(K_Grid().generate_rot_pi_011_matrix());
gf_k_rot_matrix_t K_Grid::rot_pi_011_bfx_matrix(K_Grid().generate_rot_pi_011_bfx_matrix());
gf_k_rot_matrix_t K_Grid::rot_pi_011_bfy_matrix(K_Grid().generate_rot_pi_011_bfy_matrix());
gf_k_rot_matrix_t K_Grid::rot_pi_011_bfz_matrix(K_Grid().generate_rot_pi_011_bfz_matrix());
gf_k_rot_matrix_t K_Grid::rot_pi_01m1_matrix(K_Grid().generate_rot_pi_01m1_matrix());
gf_k_rot_matrix_t K_Grid::rot_pi_01m1_bfx_matrix(K_Grid().generate_rot_pi_01m1_bfx_matrix());
gf_k_rot_matrix_t K_Grid::rot_pi_01m1_bfy_matrix(K_Grid().generate_rot_pi_01m1_bfy_matrix());
gf_k_rot_matrix_t K_Grid::rot_pi_01m1_bfz_matrix(K_Grid().generate_rot_pi_01m1_bfz_matrix());
gf_k_rot_matrix_t K_Grid::rot_pi_101_matrix(K_Grid().generate_rot_pi_101_matrix());
gf_k_rot_matrix_t K_Grid::rot_pi_101_bfx_matrix(K_Grid().generate_rot_pi_101_bfx_matrix());
gf_k_rot_matrix_t K_Grid::rot_pi_101_bfy_matrix(K_Grid().generate_rot_pi_101_bfy_matrix());
gf_k_rot_matrix_t K_Grid::rot_pi_101_bfz_matrix(K_Grid().generate_rot_pi_101_bfz_matrix());
gf_k_rot_matrix_t K_Grid::rot_pi_10m1_matrix(K_Grid().generate_rot_pi_10m1_matrix());
gf_k_rot_matrix_t K_Grid::rot_pi_10m1_bfx_matrix(K_Grid().generate_rot_pi_10m1_bfx_matrix());
gf_k_rot_matrix_t K_Grid::rot_pi_10m1_bfy_matrix(K_Grid().generate_rot_pi_10m1_bfy_matrix());
gf_k_rot_matrix_t K_Grid::rot_pi_10m1_bfz_matrix(K_Grid().generate_rot_pi_10m1_bfz_matrix());

gf_k_rot_matrix_t K_Grid::mirror_y_matrix(K_Grid().generate_mirror_y_matrix());
gf_k_rot_matrix_t K_Grid::mirror_y_bfx_matrix(K_Grid().generate_mirror_y_bfx_matrix());
gf_k_rot_matrix_t K_Grid::mirror_y_bfy_matrix(K_Grid().generate_mirror_y_bfy_matrix());
gf_k_rot_matrix_t K_Grid::mirror_y_bfz_matrix(K_Grid().generate_mirror_y_bfz_matrix());


gf_k_4MAXKPOS_matrix_t K_Grid::COS(K_Grid().generate_cos_matrix());
gf_k_4MAXKPOS_matrix_t K_Grid::SIN(K_Grid().generate_sin_matrix());

unsigned int K_Grid::get_patch(int ind[3] )
{
   int x_patch, y_patch, z_patch; 
   x_patch = (ind[0]+MAX_KPOS-1+1000*MAX_KPOS)%(2*MAX_KPOS);
   y_patch = (ind[1]+MAX_KPOS-1+1000*MAX_KPOS)%(2*MAX_KPOS);
   z_patch = (ind[2]+MAX_KPOS-1+1000*MAX_KPOS)%(2*MAX_KPOS);

   return x_patch*4*MAX_KPOS*MAX_KPOS+y_patch*2*MAX_KPOS+z_patch;      
}

std::tuple<int, int, int> K_Grid::get_indices(double kx,double ky,double kz)
{
   return tuple<int, int, int>(std::round(kx*MAX_KPOS/PI),std::round(ky*MAX_KPOS/PI),std::round(kz*MAX_KPOS/PI));
}

std::tuple<int, int, int> K_Grid::get_indices(unsigned int patch)
  {
     int x = (patch) /(4* MAX_KPOS*MAX_KPOS) -MAX_KPOS +1;
     int y = (patch/(2*MAX_KPOS)) % (2*MAX_KPOS) -MAX_KPOS +1;
     int z = ((patch % (4*MAX_KPOS*MAX_KPOS))) -MAX_KPOS +1;
     return tuple<int, int, int>(x,y,z);
  }

std::tuple<int, int, int> K_Grid::get_indices_enum(unsigned int patch)
  {
     int x = (patch / (4*MAX_KPOS*MAX_KPOS));
     int y = (patch / (2*MAX_KPOS));
     int z = (patch % (4*MAX_KPOS*MAX_KPOS));
     return tuple<int, int, int>(x,y,z);
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

unsigned int K_Grid::get_bf(int ind )
{
   return abs((2*ind-1) / (2*MAX_KPOS))%2;      
}

int K_Grid::get_add_k(unsigned int i, unsigned int j)
{
      int xi = std::get<0>(get_indices(i));
      int yi = std::get<1>(get_indices(i));
      int zi = std::get<2>(get_indices(i));
      int xj = std::get<0>(get_indices(j));
      int yj = std::get<1>(get_indices(j));
      int zj = std::get<2>(get_indices(j));
      int vec_add[3] = {xi+xj, yi+yj, zi+zj};
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
int K_Grid::get_add_k_bfz(unsigned int i, unsigned int j)
{
      int zi = std::get<2>(get_indices(i));
      int zj = std::get<2>(get_indices(j));
      return get_bf(zi+zj);
}

int K_Grid::get_dif_k(unsigned int i, unsigned int j)
{
      int xi = std::get<0>(get_indices(i));
      int yi = std::get<1>(get_indices(i));
      int zi = std::get<2>(get_indices(i));
      	    
      int xj = std::get<0>(get_indices(j));
      int yj = std::get<1>(get_indices(j));
      int zj = std::get<2>(get_indices(j));
      int vec_dif[3] = {xi-xj, yi-yj, zi-zj};
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
int K_Grid::get_dif_k_bfz(unsigned int i, unsigned int j)
{
      int zi = std::get<2>(get_indices(i));
      int zj = std::get<2>(get_indices(j));
      return get_bf(zi-zj);
}

unsigned int K_Grid::get_rot_k(double axis[3], int n, unsigned int p)
{
   double rotation[3][3];
   int  vec_rot[3];
   std::tuple<int, int, int> indices;
   
   rotation_matrix(axis, n, rotation);
   
   indices = get_indices(p);
   for(int j = 0; j<3; ++j){
      vec_rot[j] = std::round(rotation[j][0]*(std::get<0>(indices))
	                     +rotation[j][1]*(std::get<1>(indices))
			     +rotation[j][2]*(std::get<2>(indices)));
   }
   return get_patch(vec_rot);
}
unsigned int K_Grid::get_rot_k_bfx(double axis[3], int n, unsigned int p)
{
   double rotation[3][3];
   int  vec_rot[3];
   std::tuple<int, int, int> indices;
   
   rotation_matrix(axis, n, rotation);
   
   indices = get_indices(p);
   for(int j = 0; j<3; ++j){
      vec_rot[j] = std::round(rotation[j][0]*(std::get<0>(indices))
			     +rotation[j][1]*(std::get<1>(indices))
			     +rotation[j][2]*(std::get<2>(indices)));
   }
   return get_bf(vec_rot[0]);
}
unsigned int K_Grid::get_rot_k_bfy(double axis[3], int n, unsigned int p)
{
   double rotation[3][3];
   int  vec_rot[3];
   std::tuple<int, int, int> indices;
   
   rotation_matrix(axis, n, rotation);
   
   indices = get_indices(p);
   for(int j = 0; j<3; ++j){
      vec_rot[j] = std::round(rotation[j][0]*(std::get<0>(indices))
			     +rotation[j][1]*(std::get<1>(indices))
			     +rotation[j][2]*(std::get<2>(indices)));
   }
   return get_bf(vec_rot[1]);
}
unsigned int K_Grid::get_rot_k_bfz(double axis[3], int n, unsigned int p)
{
   double rotation[3][3];
   int  vec_rot[3];
   std::tuple<int, int, int> indices;
   
   rotation_matrix(axis, n, rotation);
   
   indices = get_indices(p);
   for(int j = 0; j<3; ++j){
      vec_rot[j] = std::round(rotation[j][0]*(std::get<0>(indices))
			     +rotation[j][1]*(std::get<1>(indices))
			     +rotation[j][2]*(std::get<2>(indices)));
   }
   return get_bf(vec_rot[2]);
}

unsigned int K_Grid::get_mirror_y(unsigned int p)
{
   double rotation[3][3];
   int vec_rot[3];
   std::tuple<int, int, int> indices;
   
   mirror_matrix_y(rotation);
  
  cout << "in get_mirror_y p="<<p<<endl; 
   indices = get_indices(p);
  cout << "                p="<< (std::get<0>(indices))<<","<<(std::get<1>(indices))<<","<<(std::get<2>(indices))<<endl; 
   for(int j = 0; j<3; ++j){
      vec_rot[j] = std::round(rotation[j][0]*(std::get<0>(indices))
			     +rotation[j][1]*(std::get<1>(indices))
			     +rotation[j][2]*(std::get<2>(indices)));
   }
  cout << "           mirr p="<<vec_rot[0]<<","<<vec_rot[1]<<","<<vec_rot[2]<<endl; 
  cout << "           mirr p="<<get_patch(vec_rot)<<endl; 
   return get_patch(vec_rot);
}
unsigned int K_Grid::get_mirror_y_bfx(unsigned int p)
{
   double rotation[3][3];
   int  vec_rot[3];
   std::tuple<int, int, int> indices;
   
   mirror_matrix_y(rotation);
   
   indices = get_indices(p);
   for(int j = 0; j<3; ++j){
      vec_rot[j] = std::round(rotation[j][0]*(std::get<0>(indices))
			     +rotation[j][1]*(std::get<1>(indices))
			     +rotation[j][2]*(std::get<2>(indices)));
   }
   return get_bf(vec_rot[0]);
}
unsigned int K_Grid::get_mirror_y_bfy(unsigned int p)
{
   double rotation[3][3];
   int  vec_rot[3];
   std::tuple<int, int, int> indices;
   
   mirror_matrix_y(rotation);
   
   indices = get_indices(p);
   for(int j = 0; j<3; ++j){
      vec_rot[j] = std::round(rotation[j][0]*(std::get<0>(indices))
			     +rotation[j][1]*(std::get<1>(indices))
			     +rotation[j][2]*(std::get<2>(indices)));
   }
   return get_bf(vec_rot[1]);
}
unsigned int K_Grid::get_mirror_y_bfz(unsigned int p)
{
   double rotation[3][3];
   int  vec_rot[3];
   std::tuple<int, int, int> indices;
   
   mirror_matrix_y(rotation);
   
   indices = get_indices(p);
   for(int j = 0; j<3; ++j){
      vec_rot[j] = std::round(rotation[j][0]*(std::get<0>(indices))
			     +rotation[j][1]*(std::get<1>(indices))
			     +rotation[j][2]*(std::get<2>(indices)));
   }
   return get_bf(vec_rot[2]);
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
gf_k_sum_matrix_t K_Grid::generate_add_k_bfz_matrix()
{
   cout <<"in generate add k bfz matrix"<< endl;
   gf_k_sum_matrix_t add_k_bfz_matrix;
   for( int i= 0; i< PATCH_COUNT; ++i){
      for( int k= 0; k< PATCH_COUNT; ++k){
	   add_k_bfz_matrix[i][k] = get_add_k_bfz(i,k);
      }
   } cout << "End of generate add_k_bfz:"<< endl; 
     return add_k_bfz_matrix;
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
   for( int i= 0; i< PATCH_COUNT; ++i)
   {
      for( int k= 0; k< PATCH_COUNT; ++k){
	   dif_k_bfy_matrix[i][k] = get_dif_k_bfy(i,k);
      }
   } cout << "End of generate dif_k_bfy:"<< endl; 
     return dif_k_bfy_matrix;
}
gf_k_sum_matrix_t K_Grid::generate_dif_k_bfz_matrix()
{
   cout <<"in generate dif k bfz matrix"<< endl;
   gf_k_sum_matrix_t dif_k_bfz_matrix;
   for( int i= 0; i< PATCH_COUNT; ++i){
      for( int k= 0; k< PATCH_COUNT; ++k){
	   dif_k_bfz_matrix[i][k] = get_dif_k_bfz(i,k);
      }
   } cout << "End of generate dif_k_bfz:"<< endl; 
     return dif_k_bfz_matrix;
}

//Rotation functions around the main axes
//clockwise rotation
gf_k_rot_matrix_t K_Grid::generate_rot_pi2_x_matrix()
{
   cout <<"in generate rot pi2 x matrix"<< endl;
   gf_k_rot_matrix_t rot_pi2_x_matrix;
   double axis[3]={1.0,0.,0.};
   for(unsigned int i= 0; i< PATCH_COUNT; ++i){
      rot_pi2_x_matrix[i] = get_rot_k(axis,4,i);
   }  
   return rot_pi2_x_matrix;
}
gf_k_rot_matrix_t K_Grid::generate_rot_pi2_x_bfx_matrix()
{
   cout <<"in generate rot pi2 x bfx matrix"<< endl;
   gf_k_rot_matrix_t rot_pi2_x_bfx_matrix;
   double axis[3]={1.0,0.,0.};
   for(unsigned int i= 0; i< PATCH_COUNT; ++i){
      rot_pi2_x_bfx_matrix[i] = get_rot_k_bfx(axis,4,i);
   }  
   return rot_pi2_x_bfx_matrix;
}
gf_k_rot_matrix_t K_Grid::generate_rot_pi2_x_bfy_matrix()
{
   cout <<"in generate rot pi2 x bfy matrix"<< endl;
   gf_k_rot_matrix_t rot_pi2_x_bfy_matrix;
   double axis[3]={1.0,0.,0.};
   for(unsigned int i= 0; i< PATCH_COUNT; ++i){
      rot_pi2_x_bfy_matrix[i] = get_rot_k_bfy(axis,4,i);
   }  
   return rot_pi2_x_bfy_matrix;
}
gf_k_rot_matrix_t K_Grid::generate_rot_pi2_x_bfz_matrix()
{
   cout <<"in generate rot pi2 x bfz matrix"<< endl;
   gf_k_rot_matrix_t rot_pi2_x_bfz_matrix;
   double axis[3]={1.0,0.,0.};
   for(unsigned int i= 0; i< PATCH_COUNT; ++i){
      rot_pi2_x_bfz_matrix[i] = get_rot_k_bfz(axis,4,i);
   }  
   return rot_pi2_x_bfz_matrix;
}

gf_k_rot_matrix_t K_Grid::generate_rot_pi2_y_matrix()
{
   cout <<"in generate rot pi2 y matrix"<< endl;
   gf_k_rot_matrix_t rot_pi2_y_matrix;
   double axis[3]={0.,1.0,0.};
   for(unsigned int i= 0; i< PATCH_COUNT; ++i){
      rot_pi2_y_matrix[i] = get_rot_k(axis,4,i);
   }  
   return rot_pi2_y_matrix;
}
gf_k_rot_matrix_t K_Grid::generate_rot_pi2_y_bfx_matrix()
{
   cout <<"in generate rot pi2 y bfx matrix"<< endl;
   gf_k_rot_matrix_t rot_pi2_y_bfx_matrix;
   double axis[3]={0.,1.0,0.};
   for(unsigned int i= 0; i< PATCH_COUNT; ++i){
      rot_pi2_y_bfx_matrix[i] = get_rot_k_bfx(axis,4,i);
   }  
   return rot_pi2_y_bfx_matrix;
}
gf_k_rot_matrix_t K_Grid::generate_rot_pi2_y_bfy_matrix()
{
   cout <<"in generate rot pi2 y bfy matrix"<< endl;
   gf_k_rot_matrix_t rot_pi2_y_bfy_matrix;
   double axis[3]={0.,1.0,0.};
   for(unsigned int i= 0; i< PATCH_COUNT; ++i){
      rot_pi2_y_bfy_matrix[i] = get_rot_k_bfy(axis,4,i);
   }  
   return rot_pi2_y_bfy_matrix;
}
gf_k_rot_matrix_t K_Grid::generate_rot_pi2_y_bfz_matrix()
{
   cout <<"in generate rot pi2 y bfz matrix"<< endl;
   gf_k_rot_matrix_t rot_pi2_y_bfz_matrix;
   double axis[3]={0.,1.0,0.};
   for(unsigned int i= 0; i< PATCH_COUNT; ++i){
      rot_pi2_y_bfz_matrix[i] = get_rot_k_bfz(axis,4,i);
   }  
   return rot_pi2_y_bfz_matrix;
}

gf_k_rot_matrix_t K_Grid::generate_rot_pi2_z_matrix()
{
   cout <<"in generate rot pi2 z matrix"<< endl;
   gf_k_rot_matrix_t rot_pi2_z_matrix;
   double axis[3]={0.,0.,1.};
   for(unsigned int i= 0; i< PATCH_COUNT; ++i){
      rot_pi2_z_matrix[i] = get_rot_k(axis,4,i);
   }  
   return rot_pi2_z_matrix;
}
gf_k_rot_matrix_t K_Grid::generate_rot_pi2_z_bfx_matrix()
{
   cout <<"in generate rot pi2 y bfx matrix"<< endl;
   gf_k_rot_matrix_t rot_pi2_z_bfx_matrix;
   double axis[3]={0.,0.,1.};
   for(unsigned int i= 0; i< PATCH_COUNT; ++i){
      rot_pi2_z_bfx_matrix[i] = get_rot_k_bfx(axis,4,i);
   }  
   return rot_pi2_z_bfx_matrix;
}
gf_k_rot_matrix_t K_Grid::generate_rot_pi2_z_bfy_matrix()
{
   cout <<"in generate rot pi2 y bfy matrix"<< endl;
   gf_k_rot_matrix_t rot_pi2_z_bfy_matrix;
   double axis[3]={0.,0.,1.};
   for(unsigned int i= 0; i< PATCH_COUNT; ++i){
      rot_pi2_z_bfy_matrix[i] = get_rot_k_bfy(axis,4,i);
   }  
   return rot_pi2_z_bfy_matrix;
}
gf_k_rot_matrix_t K_Grid::generate_rot_pi2_z_bfz_matrix()
{
   cout <<"in generate rot pi2 y bfz matrix"<< endl;
   gf_k_rot_matrix_t rot_pi2_z_bfz_matrix;
   double axis[3]={0.,0.,1.};
   for(unsigned int i= 0; i< PATCH_COUNT; ++i){
      rot_pi2_z_bfz_matrix[i] = get_rot_k_bfz(axis,4,i);
   }  
   return rot_pi2_z_bfz_matrix;
}

// Rotation around the diagonals
gf_k_rot_matrix_t K_Grid::generate_rot_2pi3_111_matrix()
{
   cout <<"in generate rot 2pi3 111 matrix"<< endl;
   gf_k_rot_matrix_t rot_2pi3_111_matrix;
   double axis[3]={1./sqrt(3),1./sqrt(3),1./sqrt(3)};
   for(unsigned int i= 0; i< PATCH_COUNT; ++i){
      rot_2pi3_111_matrix[i] = get_rot_k(axis,3,i);
   }  
   return rot_2pi3_111_matrix;
}
gf_k_rot_matrix_t K_Grid::generate_rot_2pi3_111_bfx_matrix()
{
   cout <<"in generate rot 2pi3 111 bfx matrix"<< endl;
   gf_k_rot_matrix_t rot_2pi3_111_bfx_matrix;
   double axis[3]={1./sqrt(3),1./sqrt(3),1./sqrt(3)};
   for(unsigned int i= 0; i< PATCH_COUNT; ++i){
      rot_2pi3_111_bfx_matrix[i] = get_rot_k_bfx(axis,3,i);
   }  
   return rot_2pi3_111_bfx_matrix;
}
gf_k_rot_matrix_t K_Grid::generate_rot_2pi3_111_bfy_matrix()
{
   cout <<"in generate rot 2pi3 111 bfy matrix"<< endl;
   gf_k_rot_matrix_t rot_2pi3_111_bfy_matrix;
   double axis[3]={1./sqrt(3),1./sqrt(3),1./sqrt(3)};
   for(unsigned int i= 0; i< PATCH_COUNT; ++i){
      rot_2pi3_111_bfy_matrix[i] = get_rot_k_bfy(axis,3,i);
   }  
   return rot_2pi3_111_bfy_matrix;
}
gf_k_rot_matrix_t K_Grid::generate_rot_2pi3_111_bfz_matrix()
{
   cout <<"in generate rot 2pi3 111 bfz matrix"<< endl;
   gf_k_rot_matrix_t rot_2pi3_111_bfz_matrix;
   double axis[3]={1./sqrt(3),1./sqrt(3),1./sqrt(3)};
   for(unsigned int i= 0; i< PATCH_COUNT; ++i){
      rot_2pi3_111_bfz_matrix[i] = get_rot_k_bfz(axis,3,i);
   }  
   return rot_2pi3_111_bfz_matrix;
}

gf_k_rot_matrix_t K_Grid::generate_rot_2pi3_1m11_matrix()
{
   cout <<"in generate rot pi 1m11 matrix"<< endl;
   gf_k_rot_matrix_t rot_2pi3_1m11_matrix;
   double axis[3]={1./sqrt(3),-1./sqrt(3),1./sqrt(3)};
   for(unsigned int i= 0; i< PATCH_COUNT; ++i){
      rot_2pi3_1m11_matrix[i] = get_rot_k(axis,3,i);
   }  
   return rot_2pi3_1m11_matrix;
}
gf_k_rot_matrix_t K_Grid::generate_rot_2pi3_1m11_bfx_matrix()
{
   cout <<"in generate rot pi 1m11 bfx matrix"<< endl;
   gf_k_rot_matrix_t rot_2pi3_1m11_bfx_matrix;
   double axis[3]={1./sqrt(3),-1./sqrt(3),1./sqrt(3)};
   for(unsigned int i= 0; i< PATCH_COUNT; ++i){
      rot_2pi3_1m11_bfx_matrix[i] = get_rot_k_bfx(axis,3,i);
   }  
   return rot_2pi3_1m11_bfx_matrix;
}
gf_k_rot_matrix_t K_Grid::generate_rot_2pi3_1m11_bfy_matrix()
{
   cout <<"in generate rot pi 1m11 bfy matrix"<< endl;
   gf_k_rot_matrix_t rot_2pi3_1m11_bfy_matrix;
   double axis[3]={1./sqrt(3),-1./sqrt(3),1./sqrt(3)};
   for(unsigned int i= 0; i< PATCH_COUNT; ++i){
      rot_2pi3_1m11_bfy_matrix[i] = get_rot_k_bfy(axis,3,i);
   }  
   return rot_2pi3_1m11_bfy_matrix;
}
gf_k_rot_matrix_t K_Grid::generate_rot_2pi3_1m11_bfz_matrix()
{
   cout <<"in generate rot pi 1m11 bfz matrix"<< endl;
   gf_k_rot_matrix_t rot_2pi3_1m11_bfz_matrix;
   double axis[3]={1./sqrt(3),-1./sqrt(3),1./sqrt(3)};
   for(unsigned int i= 0; i< PATCH_COUNT; ++i){
      rot_2pi3_1m11_bfz_matrix[i] = get_rot_k_bfz(axis,3,i);
   }  
   return rot_2pi3_1m11_bfz_matrix;
}

gf_k_rot_matrix_t K_Grid::generate_rot_2pi3_m111_matrix()
{
   cout <<"in generate rot pi m111 matrix"<< endl;
   gf_k_rot_matrix_t rot_2pi3_m111_matrix;
   double axis[3]={-1./sqrt(3),1./sqrt(3),1./sqrt(3)};
   for(unsigned int i= 0; i< PATCH_COUNT; ++i){
      rot_2pi3_m111_matrix[i] = get_rot_k(axis,3,i);
   }  
   return rot_2pi3_m111_matrix;
}
gf_k_rot_matrix_t K_Grid::generate_rot_2pi3_m111_bfx_matrix()
{
   cout <<"in generate rot pi m111 bfx matrix"<< endl;
   gf_k_rot_matrix_t rot_2pi3_m111_bfx_matrix;
   double axis[3]={-1./sqrt(3),1./sqrt(3),1./sqrt(3)};
   for(unsigned int i= 0; i< PATCH_COUNT; ++i){
      rot_2pi3_m111_bfx_matrix[i] = get_rot_k_bfx(axis,3,i);
   }  
   return rot_2pi3_m111_bfx_matrix;
}
gf_k_rot_matrix_t K_Grid::generate_rot_2pi3_m111_bfy_matrix()
{
   cout <<"in generate rot pi m111 bfy matrix"<< endl;
   gf_k_rot_matrix_t rot_2pi3_m111_bfy_matrix;
   double axis[3]={-1./sqrt(3),1./sqrt(3),1./sqrt(3)};
   for(unsigned int i= 0; i< PATCH_COUNT; ++i){
      rot_2pi3_m111_bfy_matrix[i] = get_rot_k_bfy(axis,3,i);
   }  
   return rot_2pi3_m111_bfy_matrix;
}
gf_k_rot_matrix_t K_Grid::generate_rot_2pi3_m111_bfz_matrix()
{
   cout <<"in generate rot pi m111 bfz matrix"<< endl;
   gf_k_rot_matrix_t rot_2pi3_m111_bfz_matrix;
   double axis[3]={-1./sqrt(3),1./sqrt(3),1./sqrt(3)};
   for(unsigned int i= 0; i< PATCH_COUNT; ++i){
      rot_2pi3_m111_bfz_matrix[i] = get_rot_k_bfz(axis,3,i);
   }  
   return rot_2pi3_m111_bfz_matrix;
}

gf_k_rot_matrix_t K_Grid::generate_rot_2pi3_m1m11_matrix()
{
   cout <<"in generate rot pi m1m11 matrix"<< endl;
   gf_k_rot_matrix_t rot_2pi3_m1m11_matrix;
   double axis[3]={-1./sqrt(3),-1./sqrt(3),1./sqrt(3)};
   for(unsigned int i= 0; i< PATCH_COUNT; ++i){
      rot_2pi3_m1m11_matrix[i] = get_rot_k(axis,3,i);
   }  
   return rot_2pi3_m1m11_matrix;
}
gf_k_rot_matrix_t K_Grid::generate_rot_2pi3_m1m11_bfx_matrix()
{
   cout <<"in generate rot pi m1m11 bfx matrix"<< endl;
   gf_k_rot_matrix_t rot_2pi3_m1m11_bfx_matrix;
   double axis[3]={-1./sqrt(3),-1./sqrt(3),1./sqrt(3)};
   for(unsigned int i= 0; i< PATCH_COUNT; ++i){
      rot_2pi3_m1m11_bfx_matrix[i] = get_rot_k_bfx(axis,3,i);
   }  
   return rot_2pi3_m1m11_bfx_matrix;
}
gf_k_rot_matrix_t K_Grid::generate_rot_2pi3_m1m11_bfy_matrix()
{
   cout <<"in generate rot pi m1m11 bfy matrix"<< endl;
   gf_k_rot_matrix_t rot_2pi3_m1m11_bfy_matrix;
   double axis[3]={-1./sqrt(3),-1./sqrt(3),1./sqrt(3)};
   for(unsigned int i= 0; i< PATCH_COUNT; ++i){
      rot_2pi3_m1m11_bfy_matrix[i] = get_rot_k_bfy(axis,3,i);
   }  
   return rot_2pi3_m1m11_bfy_matrix;
}
gf_k_rot_matrix_t K_Grid::generate_rot_2pi3_m1m11_bfz_matrix()
{
   cout <<"in generate rot pi m1m11 bfz matrix"<< endl;
   gf_k_rot_matrix_t rot_2pi3_m1m11_bfz_matrix;
   double axis[3]={-1./sqrt(3),-1./sqrt(3),1./sqrt(3)};
   for(unsigned int i= 0; i< PATCH_COUNT; ++i){
      rot_2pi3_m1m11_bfz_matrix[i] = get_rot_k_bfz(axis,3,i);
   }  
   return rot_2pi3_m1m11_bfz_matrix;
}

// Rotations around the plane diagonals
gf_k_rot_matrix_t K_Grid::generate_rot_pi_110_matrix()
{
   cout <<"in generate rot pi 110 matrix"<< endl;
   gf_k_rot_matrix_t rot_pi_110_matrix;
   double axis[3]={1./sqrt(2),1./sqrt(2),0.};
   for(unsigned int i= 0; i< PATCH_COUNT; ++i){
      rot_pi_110_matrix[i] = get_rot_k(axis,2,i);
   }  
   return rot_pi_110_matrix;
}
gf_k_rot_matrix_t K_Grid::generate_rot_pi_110_bfx_matrix()
{
   cout <<"in generate rot pi 110 bfx matrix"<< endl;
   gf_k_rot_matrix_t rot_pi_110_bfx_matrix;
   double axis[3]={1./sqrt(2),1./sqrt(2),0.};
   for(unsigned int i= 0; i< PATCH_COUNT; ++i){
      rot_pi_110_bfx_matrix[i] = get_rot_k_bfx(axis,2,i);
   }  
   return rot_pi_110_bfx_matrix;
}
gf_k_rot_matrix_t K_Grid::generate_rot_pi_110_bfy_matrix()
{
   cout <<"in generate rot pi 110 bfy matrix"<< endl;
   gf_k_rot_matrix_t rot_pi_110_bfy_matrix;
   double axis[3]={1./sqrt(2),1./sqrt(2),0.};
   for(unsigned int i= 0; i< PATCH_COUNT; ++i){
      rot_pi_110_bfy_matrix[i] = get_rot_k_bfy(axis,2,i);
   }  
   return rot_pi_110_bfy_matrix;
}
gf_k_rot_matrix_t K_Grid::generate_rot_pi_110_bfz_matrix()
{
   cout <<"in generate rot pi 110 bfz matrix"<< endl;
   gf_k_rot_matrix_t rot_pi_110_bfz_matrix;
   double axis[3]={1./sqrt(2),1./sqrt(2),0.};
   for(unsigned int i= 0; i< PATCH_COUNT; ++i){
      rot_pi_110_bfz_matrix[i] = get_rot_k_bfz(axis,2,i);
   }  
   return rot_pi_110_bfz_matrix;
}

gf_k_rot_matrix_t K_Grid::generate_rot_pi_1m10_matrix()
{
   cout <<"in generate rot pi 1m10 matrix"<< endl;
   gf_k_rot_matrix_t rot_pi_1m10_matrix;
   double axis[3]={1./sqrt(2),-1./sqrt(2),0.};
   for(unsigned int i= 0; i< PATCH_COUNT; ++i){
      rot_pi_1m10_matrix[i] = get_rot_k(axis,2,i);
   }  
   return rot_pi_1m10_matrix;
}
gf_k_rot_matrix_t K_Grid::generate_rot_pi_1m10_bfx_matrix()
{
   cout <<"in generate rot pi 1m10 bfx matrix"<< endl;
   gf_k_rot_matrix_t rot_pi_1m10_bfx_matrix;
   double axis[3]={1./sqrt(2),-1./sqrt(2),0.};
   for(unsigned int i= 0; i< PATCH_COUNT; ++i){
      rot_pi_1m10_bfx_matrix[i] = get_rot_k_bfx(axis,2,i);
   }  
   return rot_pi_1m10_bfx_matrix;
}
gf_k_rot_matrix_t K_Grid::generate_rot_pi_1m10_bfy_matrix()
{
   cout <<"in generate rot pi 1m10 bfy matrix"<< endl;
   gf_k_rot_matrix_t rot_pi_1m10_bfy_matrix;
   double axis[3]={1./sqrt(2),-1./sqrt(2),0.};
   for(unsigned int i= 0; i< PATCH_COUNT; ++i){
      rot_pi_1m10_bfy_matrix[i] = get_rot_k_bfy(axis,2,i);
   }  
   return rot_pi_1m10_bfy_matrix;
}
gf_k_rot_matrix_t K_Grid::generate_rot_pi_1m10_bfz_matrix()
{
   cout <<"in generate rot pi 1m10 bfz matrix"<< endl;
   gf_k_rot_matrix_t rot_pi_1m10_bfz_matrix;
   double axis[3]={1./sqrt(2),-1./sqrt(2),0.};
   for(unsigned int i= 0; i< PATCH_COUNT; ++i){
      rot_pi_1m10_bfz_matrix[i] = get_rot_k_bfz(axis,2,i);
   }  
   return rot_pi_1m10_bfz_matrix;
}

gf_k_rot_matrix_t K_Grid::generate_rot_pi_011_matrix()
{
   cout <<"in generate rot pi 011 matrix"<< endl;
   gf_k_rot_matrix_t rot_pi_011_matrix;
   double axis[3]={0.,1./sqrt(2),1./sqrt(2)};
   for(unsigned int i= 0; i< PATCH_COUNT; ++i){
      rot_pi_011_matrix[i] = get_rot_k(axis,2,i);
   }  
   return rot_pi_011_matrix;
}
gf_k_rot_matrix_t K_Grid::generate_rot_pi_011_bfx_matrix()
{
   cout <<"in generate rot pi 011 bfx matrix"<< endl;
   gf_k_rot_matrix_t rot_pi_011_bfx_matrix;
   double axis[3]={0.,1./sqrt(2),1./sqrt(2)};
   for(unsigned int i= 0; i< PATCH_COUNT; ++i){
      rot_pi_011_bfx_matrix[i] = get_rot_k_bfx(axis,2,i);
   }  
   return rot_pi_011_bfx_matrix;
}
gf_k_rot_matrix_t K_Grid::generate_rot_pi_011_bfy_matrix()
{
   cout <<"in generate rot pi 011 bfy matrix"<< endl;
   gf_k_rot_matrix_t rot_pi_011_bfy_matrix;
   double axis[3]={0.,1./sqrt(2),1./sqrt(2)};
   for(unsigned int i= 0; i< PATCH_COUNT; ++i){
      rot_pi_011_bfy_matrix[i] = get_rot_k_bfy(axis,2,i);
   }  
   return rot_pi_011_bfy_matrix;
}
gf_k_rot_matrix_t K_Grid::generate_rot_pi_011_bfz_matrix()
{
   cout <<"in generate rot pi 011 bfz matrix"<< endl;
   gf_k_rot_matrix_t rot_pi_011_bfz_matrix;
   double axis[3]={0.,1./sqrt(2),1./sqrt(2)};
   for(unsigned int i= 0; i< PATCH_COUNT; ++i){
      rot_pi_011_bfz_matrix[i] = get_rot_k_bfz(axis,2,i);
   }  
   return rot_pi_011_bfz_matrix;
}

gf_k_rot_matrix_t K_Grid::generate_rot_pi_01m1_matrix()
{
   cout <<"in generate rot pi 01m1 matrix"<< endl;
   gf_k_rot_matrix_t rot_pi_01m1_matrix;
   double axis[3]={0.,1./sqrt(2),-1./sqrt(2)};
   for(unsigned int i= 0; i< PATCH_COUNT; ++i){
      rot_pi_01m1_matrix[i] = get_rot_k(axis,2,i);
   }  
   return rot_pi_01m1_matrix;
}
gf_k_rot_matrix_t K_Grid::generate_rot_pi_01m1_bfx_matrix()
{
   cout <<"in generate rot pi 01m1 bfx matrix"<< endl;
   gf_k_rot_matrix_t rot_pi_01m1_bfx_matrix;
   double axis[3]={0.,1./sqrt(2),-1./sqrt(2)};
   for(unsigned int i= 0; i< PATCH_COUNT; ++i){
      rot_pi_01m1_bfx_matrix[i] = get_rot_k_bfx(axis,2,i);
   }  
   return rot_pi_01m1_bfx_matrix;
}
gf_k_rot_matrix_t K_Grid::generate_rot_pi_01m1_bfy_matrix()
{
   cout <<"in generate rot pi 01m1 bfy matrix"<< endl;
   gf_k_rot_matrix_t rot_pi_01m1_bfy_matrix;
   double axis[3]={0.,1./sqrt(2),-1./sqrt(2)};
   for(unsigned int i= 0; i< PATCH_COUNT; ++i){
      rot_pi_01m1_bfy_matrix[i] = get_rot_k_bfy(axis,2,i);
   }  
   return rot_pi_01m1_bfy_matrix;
}
gf_k_rot_matrix_t K_Grid::generate_rot_pi_01m1_bfz_matrix()
{
   cout <<"in generate rot pi 01m1 bfz matrix"<< endl;
   gf_k_rot_matrix_t rot_pi_01m1_bfz_matrix;
   double axis[3]={0.,1./sqrt(2),-1./sqrt(2)};
   for(unsigned int i= 0; i< PATCH_COUNT; ++i){
      rot_pi_01m1_bfz_matrix[i] = get_rot_k_bfz(axis,2,i);
   }  
   return rot_pi_01m1_bfz_matrix;
}

gf_k_rot_matrix_t K_Grid::generate_rot_pi_101_matrix()
{
   cout <<"in generate rot pi 101 matrix"<< endl;
   gf_k_rot_matrix_t rot_pi_101_matrix;
   double axis[3]={1./sqrt(2),0.,1./sqrt(2)};
   for(unsigned int i= 0; i< PATCH_COUNT; ++i){
      rot_pi_101_matrix[i] = get_rot_k(axis,2,i);
   }  
   return rot_pi_101_matrix;
}
gf_k_rot_matrix_t K_Grid::generate_rot_pi_101_bfx_matrix()
{
   cout <<"in generate rot pi 101 bfx matrix"<< endl;
   gf_k_rot_matrix_t rot_pi_101_bfx_matrix;
   double axis[3]={1./sqrt(2),0.,1./sqrt(2)};
   for(unsigned int i= 0; i< PATCH_COUNT; ++i){
      rot_pi_101_bfx_matrix[i] = get_rot_k_bfx(axis,2,i);
   }  
   return rot_pi_101_bfx_matrix;
}
gf_k_rot_matrix_t K_Grid::generate_rot_pi_101_bfy_matrix()
{
   cout <<"in generate rot pi 101 bfy matrix"<< endl;
   gf_k_rot_matrix_t rot_pi_101_bfy_matrix;
   double axis[3]={1./sqrt(2),0.,1./sqrt(2)};
   for(unsigned int i= 0; i< PATCH_COUNT; ++i){
      rot_pi_101_bfy_matrix[i] = get_rot_k_bfy(axis,2,i);
   }  
   return rot_pi_101_bfy_matrix;
}
gf_k_rot_matrix_t K_Grid::generate_rot_pi_101_bfz_matrix()
{
   cout <<"in generate rot pi 101 bfz matrix"<< endl;
   gf_k_rot_matrix_t rot_pi_101_bfz_matrix;
   double axis[3]={1./sqrt(2),0.,1./sqrt(2)};
   for(unsigned int i= 0; i< PATCH_COUNT; ++i){
      rot_pi_101_bfz_matrix[i] = get_rot_k_bfz(axis,2,i);
   }  
   return rot_pi_101_bfz_matrix;
}

gf_k_rot_matrix_t K_Grid::generate_rot_pi_10m1_matrix()
{
   cout <<"in generate rot pi 10m1 matrix"<< endl;
   gf_k_rot_matrix_t rot_pi_10m1_matrix;
   double axis[3]={1./sqrt(2),0.,-1./sqrt(2)};
   for(unsigned int i= 0; i< PATCH_COUNT; ++i){
      rot_pi_10m1_matrix[i] = get_rot_k(axis,2,i);
   }  
   return rot_pi_10m1_matrix;
}
gf_k_rot_matrix_t K_Grid::generate_rot_pi_10m1_bfx_matrix()
{
   cout <<"in generate rot pi 10m1 bfx matrix"<< endl;
   gf_k_rot_matrix_t rot_pi_10m1_bfx_matrix;
   double axis[3]={1./sqrt(2),0.,-1./sqrt(2)};
   for(unsigned int i= 0; i< PATCH_COUNT; ++i){
      rot_pi_10m1_bfx_matrix[i] = get_rot_k_bfx(axis,2,i);
   }  
   return rot_pi_10m1_bfx_matrix;
}
gf_k_rot_matrix_t K_Grid::generate_rot_pi_10m1_bfy_matrix()
{
   cout <<"in generate rot pi 10m1 bfy matrix"<< endl;
   gf_k_rot_matrix_t rot_pi_10m1_bfy_matrix;
   double axis[3]={1./sqrt(2),0.,-1./sqrt(2)};
   for(unsigned int i= 0; i< PATCH_COUNT; ++i){
      rot_pi_10m1_bfy_matrix[i] = get_rot_k_bfy(axis,2,i);
   }  
   return rot_pi_10m1_bfy_matrix;
}
gf_k_rot_matrix_t K_Grid::generate_rot_pi_10m1_bfz_matrix()
{
   cout <<"in generate rot pi 10m1 bfz matrix"<< endl;
   gf_k_rot_matrix_t rot_pi_10m1_bfz_matrix;
   double axis[3]={1./sqrt(2),0.,-1./sqrt(2)};
   for(unsigned int i= 0; i< PATCH_COUNT; ++i){
      rot_pi_10m1_bfz_matrix[i] = get_rot_k_bfz(axis,2,i);
   }  
   return rot_pi_10m1_bfz_matrix;
}

gf_k_rot_matrix_t K_Grid::generate_mirror_y_matrix()
{
   cout <<"in generate mirror y matrix"<< endl;
   gf_k_rot_matrix_t mirror_y_matrix;
   for(unsigned int i= 0; i< PATCH_COUNT; ++i){
      mirror_y_matrix[i] = get_mirror_y(i);
   }  
   return mirror_y_matrix;
}
gf_k_rot_matrix_t K_Grid::generate_mirror_y_bfx_matrix()
{
   cout <<"in generate mirror y bfx matrix"<< endl;
   gf_k_rot_matrix_t mirror_y_matrix;
   for(unsigned int i= 0; i< PATCH_COUNT; ++i){
      mirror_y_matrix[i] = get_mirror_y_bfx(i);
   }  
   return mirror_y_bfx_matrix;
}
gf_k_rot_matrix_t K_Grid::generate_mirror_y_bfy_matrix()
{
   cout <<"in generate rot pi 10m1 bfy matrix"<< endl;
   gf_k_rot_matrix_t mirror_y_bfy_matrix;
   for(unsigned int i= 0; i< PATCH_COUNT; ++i){
      mirror_y_bfy_matrix[i] = get_mirror_y_bfy(i);
   }  
   return mirror_y_bfy_matrix;
}
gf_k_rot_matrix_t K_Grid::generate_mirror_y_bfz_matrix()
{
   cout <<"in generate rot pi 10m1 bfz matrix"<< endl;
   gf_k_rot_matrix_t mirror_y_bfz_matrix;
   for(unsigned int i= 0; i< PATCH_COUNT; ++i){
      mirror_y_bfz_matrix[i] = get_mirror_y_bfz(i);
   }  
   return mirror_y_bfz_matrix;
}

#endif
#ifdef SIMPLE_PATCH
K_Grid::K_Grid()
{
   push_back( tuple<double, double, double>( 0.0, 	0.0 ,		0.0));		// PATCH 0
   push_back( tuple<double, double, double>( PI/2.0, 	PI/2.0, 	PI/2.0 ));	// PATCH 1
   push_back( tuple<double, double, double>( -PI/2.0, 	PI/2.0,		PI/2.0 ));	// PATCH 2
   push_back( tuple<double, double, double>( -PI/2.0,    -PI/2.0,	PI/2.0 ));	// PATCH 3
   push_back( tuple<double, double, double>( PI/2.0,     -PI/2.0,        PI/2.0 ));      // PATCH 4
   push_back( tuple<double, double, double>( PI/2.0,     PI/2.0,         -PI/2.0 ));      // PATCH 5  
   push_back( tuple<double, double, double>( -PI/2.0,    PI/2.0,         -PI/2.0 ));      // PATCH 6
   push_back( tuple<double, double, double>( -PI/2.0,    -PI/2.0,        -PI/2.0 ));      // PATCH 7  
   push_back( tuple<double, double, double>( PI/2.0,     -PI/2.0,        -PI/2.0 ));      // PATCH 8 
   push_back( tuple<double, double, double>( PI,         0.0,       	0.0 ));      // PATCH 9
   push_back( tuple<double, double, double>( PI,    	PI,       	0.0 ));      // PATCH 10
   push_back( tuple<double, double, double>( 0.0,   	PI,       	0.0 ));      // PATCH 11 
   push_back( tuple<double, double, double>( 0.0,   	PI,       	PI ));      // PATCH 12
   push_back( tuple<double, double, double>( 0.0,    	0.0,       	PI ));      // PATCH 13
   push_back( tuple<double, double, double>( PI,    	0.0,       	PI ));      // PATCH 14
   push_back( tuple<double, double, double>( PI,    	PI,       	PI ));      // PATCH 15
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
	    data << "# kx\t\tky\t\tkz" << endl;         
	       for( tuple<double, double, double> mom : *this )
		        data << std::get<0>(mom) << "\t\t" << std::get<1>(mom) << "\t\t" << std::get<2>(mom) << endl;
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
	   push_back([](double kx, double ky, double kz){return 1./sqrt(8*PI*PI*PI);}); //s-wave
	   //parity_arr[0] = 1.0;
	}
       //FIRST NEAREST NEIGHBOR
	if(NN_COUNT == 1){
	   //cout << " s-wave and p-wave" << endl;
	   push_back([](double kx, double ky, double kz){return 1./sqrt(8*PI*PI*PI);}); //s-wave}      
           push_back([](double kx, double ky, double kz){return 1./sqrt(4*PI*PI*PI)*cos(kx);}); //A1g
           push_back([](double kx, double ky, double kz){return 1./sqrt(4*PI*PI*PI)*cos(ky);}); //A1g
           push_back([](double kx, double ky, double kz){return 1./sqrt(4*PI*PI*PI)*cos(kz);}); //A1g
           push_back([](double kx, double ky, double kz){return 1./sqrt(4*PI*PI*PI)*sin(kx);}); //T1u x
           push_back([](double kx, double ky, double kz){return 1./sqrt(4*PI*PI*PI)*sin(ky);}); //T1u y
           push_back([](double kx, double ky, double kz){return 1./sqrt(4*PI*PI*PI)*sin(kz);}); //T1u z
    }
//	if(NN_COUNT == 1){
//	   //cout << " s-wave and p-wave" << endl;
//	   push_back([](double kx, double ky, double kz){return 1./sqrt(8*PI*PI*PI);}); //s-wave}      
//           push_back([](double kx, double ky, double kz){return 1./sqrt(12*PI*PI*PI)*(cos(kx)+cos(ky)+cos(kz));}); //A1g
//           push_back([](double kx, double ky, double kz){return 1./sqrt(4*PI*PI*PI)*sin(kx);}); //T1u x
//           push_back([](double kx, double ky, double kz){return 1./sqrt(4*PI*PI*PI)*sin(ky);}); //T1u y
//           push_back([](double kx, double ky, double kz){return 1./sqrt(4*PI*PI*PI)*sin(kz);}); //T1u z
//           push_back([](double kx, double ky, double kz){return 1./sqrt(24*PI*PI*PI)*(2*cos(kx)-cos(ky)-cos(kz));}); //Eg x
//           push_back([](double kx, double ky, double kz){return 1./sqrt(24*PI*PI*PI)*(2*cos(ky)-cos(kx)-cos(kz));}); //Eg y
//       	   }
        ffactor_count = size();
    }


unsigned int F_factors::get_ffactor_count()
{
   return ffactor_count;
}

#ifdef LOCAL
int F_factors::parity_arr[FFACTOR_COUNT] = {1};
int F_factors::translate_2pi_arr[FFACTOR_COUNT][3]= { 1, 1, 1} ;
#elif defined FIRST_NN
int F_factors::parity_arr[FFACTOR_COUNT] = {1,1,1,1,-1,-1,-1};
int F_factors::translate_2pi_arr[FFACTOR_COUNT][3] = { { 1, 1, 1}, {-1, 1, 1}, { 1,-1, 1}, { 1, 1,-1}, {-1, 1, 1}, { 1,-1, 1}, { 1, 1,-1}}; 
#endif

constexpr int F_factors::FF_rot_pi2_x_arr[FFACTOR_COUNT]; 
constexpr int F_factors::FF_rot_pi2_y_arr[FFACTOR_COUNT]; 
constexpr int F_factors::FF_rot_pi2_z_arr[FFACTOR_COUNT]; 
constexpr int F_factors::FF_rot_2pi3_111_arr[FFACTOR_COUNT]; 
constexpr int F_factors::FF_rot_2pi3_1m11_arr[FFACTOR_COUNT]; 
constexpr int F_factors::FF_rot_2pi3_m111_arr[FFACTOR_COUNT]; 
constexpr int F_factors::FF_rot_2pi3_m1m11_arr[FFACTOR_COUNT]; 
constexpr int F_factors::FF_rot_pi_110_arr[FFACTOR_COUNT]; 
constexpr int F_factors::FF_rot_pi_1m10_arr[FFACTOR_COUNT]; 
constexpr int F_factors::FF_rot_pi_011_arr[FFACTOR_COUNT]; 
constexpr int F_factors::FF_rot_pi_01m1_arr[FFACTOR_COUNT]; 
constexpr int F_factors::FF_rot_pi_101_arr[FFACTOR_COUNT]; 
constexpr int F_factors::FF_rot_pi_10m1_arr[FFACTOR_COUNT]; 
constexpr int F_factors::FF_mirror_y_arr[FFACTOR_COUNT]; 
constexpr int F_factors::FF_rot_sign_pi2_x_arr[FFACTOR_COUNT]; 
constexpr int F_factors::FF_rot_sign_pi2_y_arr[FFACTOR_COUNT]; 
constexpr int F_factors::FF_rot_sign_pi2_z_arr[FFACTOR_COUNT]; 
constexpr int F_factors::FF_rot_sign_2pi3_111_arr[FFACTOR_COUNT]; 
constexpr int F_factors::FF_rot_sign_2pi3_1m11_arr[FFACTOR_COUNT]; 
constexpr int F_factors::FF_rot_sign_2pi3_m111_arr[FFACTOR_COUNT]; 
constexpr int F_factors::FF_rot_sign_2pi3_m1m11_arr[FFACTOR_COUNT]; 
constexpr int F_factors::FF_rot_sign_pi_110_arr[FFACTOR_COUNT]; 
constexpr int F_factors::FF_rot_sign_pi_1m10_arr[FFACTOR_COUNT]; 
constexpr int F_factors::FF_rot_sign_pi_011_arr[FFACTOR_COUNT]; 
constexpr int F_factors::FF_rot_sign_pi_01m1_arr[FFACTOR_COUNT]; 
constexpr int F_factors::FF_rot_sign_pi_101_arr[FFACTOR_COUNT]; 
constexpr int F_factors::FF_rot_sign_pi_10m1_arr[FFACTOR_COUNT]; 
constexpr int F_factors::FF_mirror_y_sign_arr[FFACTOR_COUNT]; 

F_factors_real::F_factors_real()
    {
#ifdef LOCAL
       //LOCAL FORM FACTORS
        if(NN_COUNT == 0){
//	   cout << "Just s-wave" << endl;
	   push_back(std::array<std::array<std::array<dcomplex,FFREAL_DIM>,FFREAL_DIM>,FFREAL_DIM>{{{{1}}}}); //s-wave
	}
       //FIRST NEAREST NEIGHBOR
#elif defined FIRST_NN
	if(NN_COUNT == 1){
//	   cout << " s-wave and p-wave" << endl;
	   push_back(array<array<array<dcomplex,FFREAL_DIM>,FFREAL_DIM>,FFREAL_DIM>{{{{{0,0,0},{0,0,0},{0,0,0}}},{{{0,0,0},{0,1,0},{0,0,0}}},{{{0,0,0},{0,0,0},{0,0,0}}}}}); // s-wave
	   push_back(array<array<array<dcomplex,FFREAL_DIM>,FFREAL_DIM>,FFREAL_DIM>{{{{{0,0,0},{0,1./sqrt(2.),0},{0,0,0}}},{{{0,0,0},{0,0,0},{0,0,0}}},{{{0,0,0},{0,1./sqrt(2.),0},{0,0,0}}}}}); // cos x
	   push_back(array<array<array<dcomplex,FFREAL_DIM>,FFREAL_DIM>,FFREAL_DIM>{{{{{0,0,0},{0,0,0},{0,0,0}}},{{{0,1./sqrt(2.),0},{0,0,0},{0,1./sqrt(2.),0}}},{{{0,0,0},{0,0,0},{0,0,0}}}}}); // cos y
	   push_back(array<array<array<dcomplex,FFREAL_DIM>,FFREAL_DIM>,FFREAL_DIM>{{{{{0,0,0},{0,0,0},{0,0,0}}},{{{0,0,0},{1./sqrt(2.),0,1./sqrt(2.)},{0,0,0}}},{{{0,0,0},{0,0,0},{0,0,0}}}}}); // cos z
	   push_back(array<array<array<dcomplex,FFREAL_DIM>,FFREAL_DIM>,FFREAL_DIM>{{{{{0,0,0},{0,-I/sqrt(2.),0},{0,0,0}}},{{{0,0,0},{0,0,0},{0,0,0}}},{{{0,0,0},{0,I/sqrt(2.),0},{0,0,0}}}}}); // sin x
	   push_back(array<array<array<dcomplex,FFREAL_DIM>,FFREAL_DIM>,FFREAL_DIM>{{{{{0,0,0},{0,0,0},{0,0,0}}},{{{0,-I/sqrt(2.),0},{0,0,0},{0,I/sqrt(2.),0}}},{{{0,0,0},{0,0,0},{0,0,0}}}}}); // sin y
	   push_back(array<array<array<dcomplex,FFREAL_DIM>,FFREAL_DIM>,FFREAL_DIM>{{{{{0,0,0},{0,0,0},{0,0,0}}},{{{0,0,0},{-I/sqrt(2.),0,I/sqrt(2.)},{0,0,0}}},{{{0,0,0},{0,0,0},{0,0,0}}}}}); // sin z
    }
//	if(NN_COUNT == 1){
////	   cout << " s-wave and p-wave" << endl;
//	   push_back(array<array<array<dcomplex,FFREAL_DIM>,FFREAL_DIM>,FFREAL_DIM>{{{{{0,0,0},{0,0,0},{0,0,0}}},{{{0,0,0},{0,1,0},{0,0,0}}},{{{0,0,0},{0,0,0},{0,0,0}}}}}); // s-wave
//	   push_back(array<array<array<dcomplex,FFREAL_DIM>,FFREAL_DIM>,FFREAL_DIM>{{{{{0,0,0},{0,1./sqrt(6.),0},{0,0,0}}},{{{0,1./sqrt(6.),0},{1./sqrt(6.),0,1./sqrt(6.)},{0,1./sqrt(6.),0}}},{{{0,0,0},{0,1./sqrt(6.),0},{0,0,0}}}}}); // A1g
//	   push_back(array<array<array<dcomplex,FFREAL_DIM>,FFREAL_DIM>,FFREAL_DIM>{{{{{0,0,0},{0,I/sqrt(2.),0},{0,0,0}}},{{{0,0,0},{0,0,0},{0,0,0}}},{{{0,0,0},{0,-I/sqrt(2.),0},{0,0,0}}}}}); // T1u x
//	   push_back(array<array<array<dcomplex,FFREAL_DIM>,FFREAL_DIM>,FFREAL_DIM>{{{{{0,0,0},{0,0,0},{0,0,0}}},{{{0,I/sqrt(2.),0},{0,0,0},{0,-I/sqrt(2.),0}}},{{{0,0,0},{0,0,0},{0,0,0}}}}}); // T1u y
//	   push_back(array<array<array<dcomplex,FFREAL_DIM>,FFREAL_DIM>,FFREAL_DIM>{{{{{0,0,0},{0,0,0},{0,0,0}}},{{{0,0,0},{I/sqrt(2.),0,-I/sqrt(2.)},{0,0,0}}},{{{0,0,0},{0,0,0},{0,0,0}}}}}); // T1u z
//	   push_back(array<array<array<dcomplex,FFREAL_DIM>,FFREAL_DIM>,FFREAL_DIM>{{{{{0,0,0},{0,2./sqrt(12.),0},{0,0,0}}},{{{0,-1./sqrt(12.),0},{-1./sqrt(12.),0,-1./sqrt(12.)},{0,-1./sqrt(12.),0}}},{{{0,0,0},{0,2./sqrt(12.),0},{0,0,0}}}}}); // Eg x
//	   push_back(array<array<array<dcomplex,FFREAL_DIM>,FFREAL_DIM>,FFREAL_DIM>{{{{{0,0,0},{0,-1./sqrt(12.),0},{0,0,0}}},{{{0,2./sqrt(12.),0},{-1./sqrt(12.),0,-1./sqrt(12.)},{0,2./sqrt(12.),0}}},{{{0,0,0},{0,-1./sqrt(12.),0},{0,0,0}}}}}); // Eg y
//    }
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
      int rz = std::get<2>(R_Grid()[R]);
    for(int m=0; m< FFACTOR_COUNT; ++m)
       for(int n=0; n< FFACTOR_COUNT; ++n){
	  dcomplex weight(0.0,0.0);
	  for(int rp=0; rp< REAL_GRID_FF_SHELL; ++rp){
	     int rpx = std::get<0>(R_Grid()[rp]);
	     int rpy = std::get<1>(R_Grid()[rp]);
	     int rpz = std::get<2>(R_Grid()[rp]);
	     int diffx = rpx-rx;
   	     int diffy = rpy-ry;
   	     int diffz = rpz-rz;	   
	     if(abs(diffx)<= (FFREAL_DIM-1)/2 && abs(diffy)<= (FFREAL_DIM-1)/2 && abs(diffz)<= (FFREAL_DIM-1)/2 ){
		weight += conj(ffactor_real(m, diffx, diffy, diffz))*ffactor_real(n, rpx, rpy,rpz);
	     }
	  }

	  fft_fxf_weight[R][m][n] = weight;
       }
 }
 return fft_fxf_weight;
}
