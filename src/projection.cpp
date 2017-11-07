
/******************************************************************************************//** @file
 *  		
 * 	file: 		projection.cpp
 * 	contents:  	See projection.h
 * 
 ****************************************************************************************************/


#include <projection.h>
//#include <mymath.h>
#include <frg.h>
#include <cmath>
//#include <cubature.h>

using namespace std; 

F_factors ffactor;
R_Grid rgrid;


dcomplex generate_proj_matrix_ph_to_pp( const idx_proj_matrix_t& idx  )
{
   int K=idx(0);
   int Q = idx(1);
   int m = idx(2);
   int n = idx(3);
   int mp = idx(4);
   int np = idx(5);
   dcomplex val( 0.0, 0.0 );
   
   int Qx_idx = std::get<0>(K_Grid::get_indices(Q));  
   int Qy_idx = std::get<1>(K_Grid::get_indices(Q));  
   int Kx_idx = std::get<0>(K_Grid::get_indices(K));  
   int Ky_idx = std::get<1>(K_Grid::get_indices(K));  

  for(int r1 = 0; r1< REAL_GRID_FF_SHELL;++r1){
     int r1x = std::get<0>(rgrid[r1]);
     int r1y = std::get<1>(rgrid[r1]);
     for(int r2 = 0; r2< REAL_GRID_FF_SHELL;++r2){
	int r2x = std::get<0>(rgrid[r2]);
	int r2y = std::get<1>(rgrid[r2]);
	int weight_x = (1 + (translate_2pi_x(mp) * translate_2pi_x(np) * (1 - 2 *((r1x+r2x)%2)))); 
	int weight_y = (1 + (translate_2pi_y(mp) * translate_2pi_y(np) * (1 - 2 *((r1y+r2y)%2)))); 
	int weight = weight_x * weight_y / 4;  //in practice weight is either 0 or 1
	if(weight==1){
	for(int r3x = - MAX_PROJ_R_GRID+abs((r1x+r2x)%2) ; r3x< MAX_PROJ_R_GRID+1-abs((r1x+r2x)%2) ;r3x+=2)
   	   for(int r3y = - MAX_PROJ_R_GRID+abs((r1y+r2y)%2) ; r3y< MAX_PROJ_R_GRID+1-abs((r1y+r2y)%2) ;r3y+=2){
   		 if(abs(r1x+r2x-r3x) <= (FFREAL_DIM-1) && abs(r1y+r2y-r3y) <= (FFREAL_DIM-1))		
		    if(abs(-r1x-r2x-r3x) <= (FFREAL_DIM-1) && abs(-r1y-r2y-r3y)<= (FFREAL_DIM-1) )
   		       val+= conj(ffactor_real(m,(r1x+r2x-r3x)/2 ,(r1y+r2y-r3y)/2))*
   			     ffactor_real(n,-(r1x+r2x+r3x)/2 ,-(r1y+r2y+r3y)/2)*
   			     ffactor_real(mp, r1x ,r1y)*conj(ffactor_real(np, r2x ,r2y))*
			     (cos_k(Qx_idx * r3x + Qy_idx * r3y + Kx_idx * (r2x-r1x) + Ky_idx * (r2y-r1y)) + I * sin_k(Qx_idx * r3x + Qy_idx * r3y + Kx_idx * (r2x-r1x) + Ky_idx * (r2y-r1y) ));
	    		     }
	   }
     }}
   return val;
}


dcomplex generate_proj_matrix_pp_to_ph( const idx_proj_matrix_t& idx )
{
   int K=idx(0);
   int Q = idx(1);
   int m = idx(2);
   int n = idx(3);
   int mp = idx(4);
   int np = idx(5);
   dcomplex val( 0.0, 0.0 );
   
   int Qx_idx = std::get<0>(K_Grid::get_indices(Q));  
   int Qy_idx = std::get<1>(K_Grid::get_indices(Q));  
   int Kx_idx = std::get<0>(K_Grid::get_indices(K));  
   int Ky_idx = std::get<1>(K_Grid::get_indices(K));  
  
  for(int r1 = 0; r1< REAL_GRID_FF_SHELL;++r1){
     int r1x = std::get<0>(rgrid[r1]);
     int r1y = std::get<1>(rgrid[r1]);
     for(int r2 = 0; r2< REAL_GRID_FF_SHELL;++r2){
	int r2x = std::get<0>(rgrid[r2]);
	int r2y = std::get<1>(rgrid[r2]);
	int weight_x = (1 + (translate_2pi_x(mp) * translate_2pi_x(np) * (1 - 2 *((r1x+r2x)%2)))); 
	int weight_y = (1 + (translate_2pi_y(mp) * translate_2pi_y(np) * (1 - 2 *((r1y+r2y)%2)))); 
	int weight = weight_x * weight_y / 4;  //in practice weight is either 0 or 1
	if(weight==1){
	for(int r3x = - MAX_PROJ_R_GRID+abs((r1x+r2x)%2) ; r3x< MAX_PROJ_R_GRID+1-abs((r1x+r2x)%2) ;r3x+=2)
   	   for(int r3y = - MAX_PROJ_R_GRID+abs((r1y+r2y)%2) ; r3y< MAX_PROJ_R_GRID+1-abs((r1y+r2y)%2) ;r3y+=2){
   		 if(abs(r1x-r2x+r3x) <= (FFREAL_DIM-1) && abs(r1y-r2y+r3y) <= (FFREAL_DIM-1) )
   		    if(abs(r1x-r2x-r3x) <= (FFREAL_DIM-1) && abs(r1y-r2y-r3y)<= (FFREAL_DIM-1) )
   		       val+= conj(ffactor_real(m,(r1x-r2x+r3x)/2 ,(r1y-r2y+r3y)/2))*
			  ffactor_real(n,(r1x-r2x-r3x)/2 ,(r1y-r2y-r3y)/2)*
			  ffactor_real(mp, r1x ,r1y)*conj(ffactor_real(np, r2x ,r2y))*
			  (cos_k(Qx_idx * r3x + Qy_idx * r3y + Kx_idx * (r2x+r1x) + Ky_idx * (r2y+r1y) ) + I * sin_k(Qx_idx * r3x + Qy_idx * r3y + Kx_idx * (r2x+r1x) + Ky_idx * (r2y+r1y) ));
   		    }
	}
     }}

   return val;
}



dcomplex generate_proj_matrix_ph_to_xph( const idx_proj_matrix_t& idx )
{
   int K=idx(0);
   int Q = idx(1);
   int m = idx(2);
   int n = idx(3);
   int mp = idx(4);
   int np = idx(5);
   dcomplex val( 0.0, 0.0 );
   
   int Qx_idx = std::get<0>(K_Grid::get_indices(Q));  
   int Qy_idx = std::get<1>(K_Grid::get_indices(Q));  
   int Kx_idx = std::get<0>(K_Grid::get_indices(K));  
   int Ky_idx = std::get<1>(K_Grid::get_indices(K));  
  
  for(int r1 = 0; r1< REAL_GRID_FF_SHELL;++r1){
     int r1x = std::get<0>(rgrid[r1]);
     int r1y = std::get<1>(rgrid[r1]);
     for(int r2 = 0; r2< REAL_GRID_FF_SHELL;++r2){
	int r2x = std::get<0>(rgrid[r2]);
	int r2y = std::get<1>(rgrid[r2]);
	int weight_x = (1 + (translate_2pi_x(mp) * translate_2pi_x(np) * (1 - 2 *((r1x+r2x)%2)))); 
	int weight_y = (1 + (translate_2pi_y(mp) * translate_2pi_y(np) * (1 - 2 *((r1y+r2y)%2)))); 
	int weight = weight_x * weight_y / 4;  //in practice weight is either 0 or 1
	if(weight==1){
	for(int r3x = - MAX_PROJ_R_GRID+abs((r1x+r2x)%2) ; r3x< MAX_PROJ_R_GRID+1-abs((r1x+r2x)%2) ;r3x+=2)
   	   for(int r3y = - MAX_PROJ_R_GRID+abs((r1y+r2y)%2) ; r3y< MAX_PROJ_R_GRID+1-abs((r1y+r2y)%2) ;r3y+=2){
   		 if(abs(r1x-r2x-r3x) <= (FFREAL_DIM-1) && abs(r1y-r2y-r3y) <= (FFREAL_DIM-1) )
   		    if(abs(-r1x+r2x-r3x) <= (FFREAL_DIM-1) && abs(-r1y+r2y-r3y)<= (FFREAL_DIM-1) )
   		       val+= conj(ffactor_real(m,(r1x-r2x-r3x)/2 ,(r1y-r2y-r3y)/2))*
			  ffactor_real(n,(-r1x+r2x-r3x)/2 ,(-r1y+r2y-r3y)/2)*
			  ffactor_real(mp, r1x ,r1y)*conj(ffactor_real(np, r2x ,r2y))*
			  (cos_k(Qx_idx * r3x + Qy_idx * r3y + Kx_idx * (r2x+r1x) + Ky_idx * (r2y+r1y) ) + I * sin_k(Qx_idx * r3x + Qy_idx * r3y + Kx_idx * (r2x+r1x) + Ky_idx * (r2y+r1y) ));
   		    }
	}
     }}
   
  return val;
} 


//BARE VERTEX PROJECTION

dcomplex generate_proj_vert_bare(const idx_vert_bare_ff_t& idx )
{
   int n_in = idx(0);
   int n_out = idx(1);
   int s1_in = idx(2);
   int s2_in = idx(3);
   int s1_out = idx(4);
   int s2_out = idx(5);

   if( n_in==0 && n_out==0)
      return 4.0 * PI * PI * vert_bare(s1_in,s2_in,s1_out,s2_out);
   else
      return dcomplex(0.0,0.0);
}

