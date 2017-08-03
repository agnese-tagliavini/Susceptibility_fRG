
/******************************************************************************************//** @file
 *  		
 * 	file: 		grid.h
 * 	contents:	Classes for frequency and momentum grid  
 * 
 ****************************************************************************************************/


#pragma once

#include <vector>
#include <string>
#include <iostream>
#include <const.h>
#include <functional>
#include <def.h>
#include <mymath.h>


class F_Grid
{	
   public: 	
      /********************* Constructor ********************/
      F_Grid( unsigned int pos_freq_count, double step_size );	///< Initialing constructor for equidistant frequency grid
      //F_Grid( F_Grid&& ) = default;			///< Default move constructor

      /***************************************************/
      void print( std::string fname = "f_grid.dat");	///< Prints frequency grid to file "f_grid.dat"

      inline const double& operator[]( int w )
      {
	 return grid_points[w + pos_freq_count]; 
      }
      
      inline const double * data() const

      {
	 return grid_points.data(); 
      }

      unsigned int get_pos_freq_count() const;
      double get_step_size() const;

   private:

      std::vector<double> grid_points;  		///< Values of grid
      const int 	 pos_freq_count;		///< Number of elements in grid
      const double 	 step_size;			///< Step size of equidistant grid 

};

class Bos_Grid
{	
   public:
      /********************* Constructor ********************/
      Bos_Grid( unsigned int pos_freq_count, double step_size );	///< Initialing constructor for equidistant frequency grid
      //Bos_Grid( Bos_Grid&& ) = default;			///< Default move constructor

      /***************************************************/
      void print( std::string fname = "bos_grid.dat");	///< Prints frequency grid to file "f_grid.dat"

      inline const double& operator[]( int w )
      {
	 return grid_points[w + pos_freq_count]; 

      }

      inline const double * data() const
      {
	 return grid_points.data(); 

      }

      unsigned int get_pos_freq_count() const;
      double get_step_size() const;

   private:

      std::vector<double> grid_points; 			///< Values of grid
      const int 	 pos_freq_count;		///< Number of elements in grid
      const double 	 step_size;			///< Step size of equidistant grid 

};

/********************* R_Grid Class ********************/

class R_Grid : public std::vector<std::tuple<int, int>>
{	
   public:

      /********************* Constructor ********************/
      R_Grid();						///< Constructor for real lattice grid accordingly with NN_COUNT
      
      /***************************************************/
};

/********************* K_Grid Class ********************/

class K_Grid : public std::vector<std::tuple<double, double>>
{	
   public:

      /********************* Constructor ********************/
      K_Grid();						///< Constructor for 8-patch momentum grid
      
      /***************************************************/

      void print( std::string fname = "k_grid.dat"); 			///< Prints frequency grid to file "f_grid.dat"
      unsigned int get_patch_count();

#ifdef UNIFORM_GRID 
      static unsigned int get_patch(int ind[2] );
      static unsigned int get_bf(int indx );
      static std::tuple<int, int> get_indices(double kx, double ky);
      static std::tuple<int, int> get_indices(unsigned int patch);     //kx and ky starting from -MAX_KPOS+1 (corresponding to physical -pi+delta)
      static std::tuple<int, int> get_indices_enum(unsigned int patch);//kx and ky starting from 0 (but corresponding to -pi+delta)
      gf_k_4MAXKPOS_matrix_t generate_cos_matrix();
      gf_k_4MAXKPOS_matrix_t generate_sin_matrix();

      static gf_k_4MAXKPOS_matrix_t COS;  // cosine( 1/2 * Q_x * PI/MAX_KPOS ) starting from 0 to 4*MAX_KPOS-1
      static gf_k_4MAXKPOS_matrix_t SIN;  // sine( 1/2 * Q_x * PI/MAX_KPOS )     !pay attention if input becomes minus!

      
      gf_k_sum_matrix_t generate_add_k_matrix();
      gf_k_sum_matrix_t generate_add_k_bfx_matrix();
      gf_k_sum_matrix_t generate_add_k_bfy_matrix();
      gf_k_sum_matrix_t generate_dif_k_matrix();
      gf_k_sum_matrix_t generate_dif_k_bfx_matrix();
      gf_k_sum_matrix_t generate_dif_k_bfy_matrix();
      
      gf_k_rot_matrix_t generate_rot_pi2_z_matrix();
      gf_k_rot_matrix_t generate_rot_pi2_z_bfx_matrix();
      gf_k_rot_matrix_t generate_rot_pi2_z_bfy_matrix();
     
      gf_k_rot_matrix_t generate_mirror_y_matrix();
      gf_k_rot_matrix_t generate_mirror_y_bfx_matrix();
      gf_k_rot_matrix_t generate_mirror_y_bfy_matrix();
      
      gf_k_rot_matrix_t generate_mirror_diagonal_matrix();
      
      static int k_origin[2]; //= {0,0,0};
      static unsigned int patch_origin;//  = K_Grid().get_patch(k_origin);
      static gf_k_sum_matrix_t add_k_matrix;
      static gf_k_sum_matrix_t add_k_bfx_matrix;
      static gf_k_sum_matrix_t add_k_bfy_matrix;
      static gf_k_sum_matrix_t dif_k_matrix;
      static gf_k_sum_matrix_t dif_k_bfx_matrix;
      static gf_k_sum_matrix_t dif_k_bfy_matrix;

      static gf_k_rot_matrix_t rot_pi2_z_matrix;
      static gf_k_rot_matrix_t rot_pi2_z_bfx_matrix;
      static gf_k_rot_matrix_t rot_pi2_z_bfy_matrix;
 
      static gf_k_rot_matrix_t mirror_y_matrix;
      static gf_k_rot_matrix_t mirror_y_bfx_matrix;
      static gf_k_rot_matrix_t mirror_y_bfy_matrix;
      
      static gf_k_rot_matrix_t mirror_diagonal_matrix;

      int get_add_k(unsigned int i, unsigned int j);
      int get_add_k_bfx(unsigned int i, unsigned int j);
      int get_add_k_bfy(unsigned int i, unsigned int j);
      int get_dif_k(unsigned int i, unsigned int j); 
      int get_dif_k_bfx(unsigned int i, unsigned int j); 
      int get_dif_k_bfy(unsigned int i, unsigned int j); 
      unsigned int get_rot_k(int n, unsigned int p); // axis of the rotation, n indicates "angle" (n=2 is pi, n=4 is pi/2, p patch to be rotated )
      unsigned int get_rot_k_bfx(int n, unsigned int p); // axis of the rotation, n indicates "angle" (n=2 is pi, n=4 is pi/2, p patch to be rotated )
      unsigned int get_rot_k_bfy(int n, unsigned int p); // axis of the rotation, n indicates "angle" (n=2 is pi, n=4 is pi/2, p patch to be rotated )
      unsigned int get_mirror_k(unsigned int p); 	// axis of the rotation, n indicates "angle" (n=2 is pi, n=4 is pi/2, p patch to be rotated )
      unsigned int get_mirror_k_bfx(unsigned int p); 	// axis of the rotation, n indicates "angle" (n=2 is pi, n=4 is pi/2, p patch to be rotated )
      unsigned int get_mirror_k_bfy(unsigned int p); 	// axis of the rotation, n indicates "angle" (n=2 is pi, n=4 is pi/2, p patch to be rotated )
      unsigned int get_mirror_diagonal_k(unsigned int p); 	// (x<->y )
//      
#endif
      unsigned int patch_count;
};    

#ifdef UNIFORM_GRID
inline int add_k( int k1, int k2, int bf[2] )
{
   bf[0] = abs(K_Grid::add_k_bfx_matrix[k1][k2])%2;
   bf[1] = abs(K_Grid::add_k_bfy_matrix[k1][k2])%2;
   return K_Grid::add_k_matrix[k1][k2]; 
}

inline int dif_k( int k1, int k2, int bf[2] )
{
   bf[0] = abs(K_Grid::dif_k_bfx_matrix[k1][k2])%2;
   bf[1] = abs(K_Grid::dif_k_bfy_matrix[k1][k2])%2;
   return K_Grid::dif_k_matrix[k1][k2];
}

inline int neg_k( int k, int bf[2])
{
   bf[0]= abs(K_Grid::dif_k_bfx_matrix[K_Grid::patch_origin][k])%2 ;
   bf[1]= abs(K_Grid::dif_k_bfy_matrix[K_Grid::patch_origin][k])%2 ;
   return K_Grid::dif_k_matrix[K_Grid::patch_origin][k];
}

inline int get_k_patch(int x, int y)
{
   int vec[2] = {x,y};
   return K_Grid::get_patch(vec);
}

inline double get_kx( unsigned int k)
{
   return std::get<0>(K_Grid::get_indices(k))*PI/MAX_KPOS;
}

inline double get_ky( unsigned int k)
{
   return std::get<1>(K_Grid::get_indices(k))*PI/MAX_KPOS;
}


// Inline functions for the finer k grid -> used in FFT
// range is going from px={0,...,pi-delta,-pi,...0-delta}


inline double get_px(int p)
{   
   int px = (((p/FFT_DIM) % FFT_DIM)+FFT_DIM/2) % FFT_DIM - FFT_DIM/2;
   return px * 2.0 * PI/FFT_DIM;
}

inline double get_py(int p)
{   
   int py = (p % FFT_DIM+FFT_DIM/2) % FFT_DIM - FFT_DIM/2;
   return py * 2.0 * PI/FFT_DIM;
}
// inline functionto get back x,y,z index from patch in finer grid
// range is going from px_idx ={0,...FFT_DIM}

inline int get_px_idx(int p)
{   
   return (p/FFT_DIM) % FFT_DIM;
}

inline int get_py_idx(int p)
{   
   return (p % FFT_DIM);
}

inline int get_p_patch(int x, int y)
{
   int x_patch, y_patch; 
   x_patch = (x+1000*FFT_DIM)%(FFT_DIM);
   y_patch = (y+1000*FFT_DIM)%(FFT_DIM);

   return y_patch + FFT_DIM*x_patch;      
}

inline int p_to_k_patch(int p)
{
   int px, py;
   double kx, ky;
   px = (p/FFT_DIM) % FFT_DIM;
   py = p % FFT_DIM;

   kx = (double)( px * 2. * MAX_KPOS/FFT_DIM );
   ky = (double)( py * 2. * MAX_KPOS/FFT_DIM );
   
   return get_k_patch(round(kx), round(ky));
 

}
inline int k_to_p_patch(int k)
{
   double px, py;
   int kx, ky;
   kx = (k/2/MAX_KPOS) % (2 * MAX_KPOS);
   ky = k % (2 * MAX_KPOS);

   px = (double)(kx*FFT_DIM)/2./(double)(MAX_KPOS)-(double)(FFT_DIM*(MAX_KPOS-1))/2./(double)(MAX_KPOS);
   py = (double)(ky*FFT_DIM)/2./(double)(MAX_KPOS)-(double)(FFT_DIM*(MAX_KPOS-1))/2./(double)(MAX_KPOS);
   
   return get_p_patch(round(px), round(py));
 

}
// returns "ceiling" of a patching point divided by 2. Needed for the bubble integration
inline unsigned int div2_ceil_patch( int k )
{
   std::tuple<int,int> k_vec = K_Grid::get_indices(k);
   int div2_ceil_x = div2_ceil(std::get<0>(k_vec));
   int div2_ceil_y = div2_ceil(std::get<1>(k_vec));
   int div2_ceil[2]={div2_ceil_x,div2_ceil_y};
   return K_Grid::get_patch(div2_ceil);
}

// returns "floor" of a patching point divided by 2. Needed for the bubble integration
inline unsigned int div2_floor_patch( int k )
{
   std::tuple<int,int> k_vec = K_Grid::get_indices(k);
   int div2_floor_x = div2_floor(std::get<0>(k_vec));
   int div2_floor_y = div2_floor(std::get<1>(k_vec));
   int div2_floor[2]={div2_floor_x,div2_floor_y};
   return K_Grid::get_patch(div2_floor);
}

// returns the half a step in the k-grid if k is "odd" and 0 if k is "eveevenn"

inline double mod_patch_x( int k )
{
   return  ( (std::get<0>(K_Grid::get_indices(k))) + 100000) %2 ;
}

inline double mod_patch_y( int k )
{
   return  ( (std::get<1>(K_Grid::get_indices(k))) + 100000) %2 ;
}

inline double cos_k(int kx)
{
   return K_Grid::COS[(kx+40*MAX_KPOS)%(4*MAX_KPOS)];
}

inline double sin_k(int kx)
{
   return K_Grid::SIN[(kx+40*MAX_KPOS)%(4*MAX_KPOS)];
}

#endif

//TODO: Now only local ffactors are considered(s-wave) so we have just one ffactor function
//      When we wat to extend to next nearest neighbors the ffactor_mom list will be extended accordignly with the next-nearest neighbor cosidered and the new ffcator function in ka space must be defined in private and pushed back to fill fffact_mom (maybe with if conditions on the number of next-nearest neighbors)...in that case a fuction to get the size of ffactor_mom has to be included in the class
//TODO: Include the all the rotation matrices for the form factors considered

class F_factors : public std::vector<std::function<double(double,double)>>
{
       public:
	      F_factors();
	      static int parity_arr[FFACTOR_COUNT];
      static int translate_2pi_arr[FFACTOR_COUNT][2];
#ifdef LOCAL
      static constexpr int FF_rot_pi2_z_arr[FFACTOR_COUNT]={ 0 };
      static constexpr int FF_mirror_y_arr[FFACTOR_COUNT]={ 0 };
      static constexpr int FF_mirror_diagonal_arr[FFACTOR_COUNT]={ 0 }; 
      static constexpr int FF_rot_sign_pi2_z_arr[FFACTOR_COUNT]={ 1 };
      static constexpr int FF_mirror_sign_y_arr[FFACTOR_COUNT]={ 1 };
      static constexpr int FF_mirror_sign_diagonal_arr[FFACTOR_COUNT]={ 1 };
#elif defined FIRST_NN
      static constexpr int FF_rot_pi2_z_arr[FFACTOR_COUNT] 		= { 0, 2, 1, 4, 3};
      static constexpr int FF_mirror_y_arr[FFACTOR_COUNT] 		= { 0, 1, 2, 3, 4};
      static constexpr int FF_mirror_diagonal_arr[FFACTOR_COUNT] 	= { 0, 2, 1, 4, 3};
      static constexpr int FF_rot_sign_pi2_z_arr[FFACTOR_COUNT] 	= { 1, 1, 1,-1, 1};
      static constexpr int FF_mirror_sign_y_arr[FFACTOR_COUNT] 		= { 1, 1, 1,-1, 1};
      static constexpr int FF_mirror_sign_diagonal_arr[FFACTOR_COUNT] 	= { 1, 1, 1, 1, 1};
#endif
	      unsigned int get_ffactor_count();
       private:
	      int ffactor_count;
              	      
};

class F_factors_real : public std::vector<std::array<std::array<dcomplex,FFREAL_DIM>,FFREAL_DIM>>
{
       public:
	      F_factors_real();
	      static gf_fft_fxf_weight_t fft_fxf_weight;
	      gf_fft_fxf_weight_t generate_fft_fxf_weight();
	      unsigned int get_ffactor_real_count();
       private:
	      int ffactor_real_count;
              	      
};

inline dcomplex ffactor_real( unsigned int m, const int x, const int y )
{
   std::array<std::array<dcomplex,FFREAL_DIM>,FFREAL_DIM> Arr = F_factors_real()[m];
   return  Arr[x+(FFREAL_DIM-1)/2][y+(FFREAL_DIM-1)/2];   // The element [0][0][0] correspomds to the cube point (-1,-1,-1)
}

inline int Parity(int n)
{
   return F_factors::parity_arr[n];
}
 
// rotations of the Formfactors

inline int FF_rot_pi2_z(int n)
{
   return F_factors::FF_rot_pi2_z_arr[n];
}
inline int FF_mirror_y(int n)
{
   return F_factors::FF_mirror_y_arr[n];
}
inline int FF_mirror_diagonal(int n)
{
   return F_factors::FF_mirror_diagonal_arr[n];
}

// sign for rotation

inline int FF_rot_sign_pi2_z(int n)
{
   return F_factors::FF_rot_sign_pi2_z_arr[n];
}
inline int FF_mirror_sign_y(int n)
{
   return F_factors::FF_mirror_sign_y_arr[n];
}
inline int FF_mirror_sign_diagonal(int n)
{
   return F_factors::FF_mirror_sign_diagonal_arr[n];
}

// gives -sign (or +sign) for folding back by 2pi 

inline int translate_2pi_x(int m)
{
   return F_factors::translate_2pi_arr[m][0];
}                                           
                                            
inline int translate_2pi_y(int m)           
{                                           
   return F_factors::translate_2pi_arr[m][1];
}
