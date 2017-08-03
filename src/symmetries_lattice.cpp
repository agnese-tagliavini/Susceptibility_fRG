
/************************************************************************************************//**
 *  		
 * 	file: 		symmetries_lattice.cpp
 * 	contents:   	see symmetries.h
 * 
 ****************************************************************************************************/

#include <symmetries.h>
#include <grid.h>
// ROTATION SYMMETRIES IMPLEMENTED JUST FOR THE MOMENTUM TRANSFER Q (S-WAVE FFACTOR ASSUMED)-> TODO: INCLUDE THE FFACTOR ROTATION MATRICES IN GRID.H
// -- Phi

operation rot_pi2_z_phi( idx_phi_t& idx )
{
   bool sign = (1- FF_rot_sign_pi2_z( idx(4) ) * FF_rot_sign_pi2_z( idx(5) ))/2;
   bool bf_sign_x = (1 - (translate_2pi_x(idx(4))*translate_2pi_x(idx(5))))/2;
   bool bf_sign_y = (1 - (translate_2pi_y(idx(4))*translate_2pi_y(idx(5))))/2;
   int bf[2] = {0,0};
   rot_pi2_z( idx(3), bf );
   idx(4)= FF_rot_pi2_z( idx(4) );
   idx(5)= FF_rot_pi2_z( idx(5) );
   bf_sign_x &= bool(bf[0]);
   bf_sign_y &= bool(bf[1]);

   return operation(sign^(bf_sign_x^bf_sign_y), false );
}
   

operation mirror_y_phi( idx_phi_t& idx )
{
   bool sign = (1- FF_mirror_sign_y( idx(4) ) * FF_mirror_sign_y( idx(5) ))/2;
   bool bf_sign_x = (1 - (translate_2pi_x(idx(4))*translate_2pi_x(idx(5))))/2;
   bool bf_sign_y = (1 - (translate_2pi_y(idx(4))*translate_2pi_y(idx(5))))/2;
   int bf[2] = {0,0};
   mirror_y( idx(3), bf );
   idx(4)= FF_mirror_y( idx(4) );
   idx(5)= FF_mirror_y( idx(5) );
   bf_sign_x &= bool(bf[0]);
   bf_sign_y &= bool(bf[1]);

   return operation(sign^(bf_sign_x^bf_sign_y), false );
}
// ---- P ----

operation rot_pi2_z_P( idx_P_t& idx )
{
   bool sign =(1- FF_rot_sign_pi2_z( idx(3) ))/2;
   bool bf_sign_x = (1 - translate_2pi_x(idx(3)))/2;
   bool bf_sign_y = (1 - translate_2pi_y(idx(3)))/2;
   int bf[2] = {0,0};
   rot_pi2_z( idx(2), bf );
   idx(3) = FF_rot_pi2_z( idx(3) );
   bf_sign_x &= bool(bf[0]);
   bf_sign_y &= bool(bf[1]);

   return operation(sign^(bf_sign_x^bf_sign_y), false );
}
   


operation mirror_y_P( idx_P_t& idx )
{
   bool sign = (1- FF_mirror_sign_y( idx(3) ) )/2;
   bool bf_sign_x = (1 - translate_2pi_x(idx(3)))/2;
   bool bf_sign_y = (1 - translate_2pi_y(idx(3)))/2;
   int bf[2] = {0,0};
   mirror_y( idx(2), bf );
   idx(3) = FF_mirror_y( idx(3) );
   bf_sign_x &= bool(bf[0]);
   bf_sign_y &= bool(bf[1]);

   return operation(sign^(bf_sign_x^bf_sign_y), false );
}

// ---- chi ----

operation rot_pi2_z_chi( idx_chi_t& idx )
{
   int bf[2] = {0,0};
   rot_pi2_z( idx(1), bf );

   return operation( false, false );
}
   
operation mirror_y_chi( idx_chi_t& idx )
{
   int bf[2] = {0,0};
   mirror_y( idx(1), bf );

   return operation( false, false );
}

// ---- suscept ---
operation rot_pi2_z_suscept( idx_suscept_t& idx )
{
   bool sign = (1- FF_rot_sign_pi2_z( idx(2) ) * FF_rot_sign_pi2_z( idx(3) ))/2;
   bool bf_sign_x = (1 - (translate_2pi_x(idx(2))*translate_2pi_x(idx(3))))/2;
   bool bf_sign_y = (1 - (translate_2pi_y(idx(2))*translate_2pi_y(idx(3))))/2;
   int bf[2] = {0,0};
   rot_pi2_z( idx(1), bf );
   idx(2)= FF_rot_pi2_z( idx(2) );
   idx(3)= FF_rot_pi2_z( idx(3) );
   bf_sign_x &= bool(bf[0]);
   bf_sign_y &= bool(bf[1]);

   return operation(sign^(bf_sign_x^bf_sign_y), false );
}
   

operation mirror_y_suscept( idx_suscept_t& idx )
{
   bool sign = (1- FF_mirror_sign_y( idx(2) ) * FF_mirror_sign_y( idx(3) ))/2;
   bool bf_sign_x = (1 - (translate_2pi_x(idx(2))*translate_2pi_x(idx(3))))/2;
   bool bf_sign_y = (1 - (translate_2pi_y(idx(2))*translate_2pi_y(idx(3))))/2;
   int bf[2] = {0,0};
   mirror_y( idx(1), bf );
   idx(2)= FF_mirror_y( idx(2) );
   idx(3)= FF_mirror_y( idx(3) );
   bf_sign_x &= bool(bf[0]);
   bf_sign_y &= bool(bf[1]);

   return operation(sign^(bf_sign_x^bf_sign_y), false );
}

// ---- Sig ----

operation rot_pi2_z_sig( idx_1p_t& idx )
{
   int bf[2] = {0,0};
   rot_pi2_z( idx(1), bf );

   return operation( false, false );
}
   


operation mirror_y_sig( idx_1p_t& idx )
{
   int bf[2] = {0,0};
   mirror_y( idx(1), bf );

   return operation( false, false );
}

// -- BARE VERTEX PROJECTED
// not used because local bare interaction is really easy in FF
/*******************************************************************
 *
 * 			SYMMETRIES PROJECTION MATRICES
 *
 ******************************************************************/

operation rot_pi2_z_projmat( idx_proj_matrix_t& idx )
{
   bool sign= ((1- (FF_rot_sign_pi2_z( idx(4) ) * FF_rot_sign_pi2_z( idx(5) ) * FF_rot_sign_pi2_z( idx(2) ) * FF_rot_sign_pi2_z( idx(3) )))/2);
   bool bf_sign_x = (1 - (translate_2pi_x(idx(4))*translate_2pi_x(idx(5))))/2;
   bool bf_sign_y = (1 - (translate_2pi_y(idx(4))*translate_2pi_y(idx(5))))/2;
   bool bf_sign_p_x = (1 - (translate_2pi_x(idx(2))*translate_2pi_x(idx(3))))/2;
   bool bf_sign_p_y = (1 - (translate_2pi_y(idx(2))*translate_2pi_y(idx(3))))/2;
   int bf[2] = {0,0};
   int bf_p[2] = {0,0};
   rot_pi2_z( idx(0), bf_p );
   rot_pi2_z( idx(1), bf );
   idx(2) = FF_rot_pi2_z( idx(2) );
   idx(3) = FF_rot_pi2_z( idx(3) );
   idx(4) = FF_rot_pi2_z( idx(4) );
   idx(5) = FF_rot_pi2_z( idx(5) );
   bf_sign_x &= bool(bf[0]);
   bf_sign_y &= bool(bf[1]);
   bf_sign_p_x &= bool(bf_p[0]);
   bf_sign_p_y &= bool(bf_p[1]);

   return operation( sign ^((bf_sign_x^bf_sign_y)^(bf_sign_p_x^bf_sign_p_y)), false );
}
   


operation mirror_y_projmat( idx_proj_matrix_t& idx )
{
   bool sign= ((1- (FF_mirror_sign_y( idx(4) ) * FF_mirror_sign_y( idx(5) ) * FF_mirror_sign_y( idx(2) ) * FF_mirror_sign_y( idx(3) )))/2);
   bool bf_sign_x = (1 - (translate_2pi_x(idx(4))*translate_2pi_x(idx(5))))/2;
   bool bf_sign_y = (1 - (translate_2pi_y(idx(4))*translate_2pi_y(idx(5))))/2;
   bool bf_sign_p_x = (1 - (translate_2pi_x(idx(2))*translate_2pi_x(idx(3))))/2;
   bool bf_sign_p_y = (1 - (translate_2pi_y(idx(2))*translate_2pi_y(idx(3))))/2;
   int bf[2] = {0,0};
   int bf_p[2] = {0,0};
   mirror_y( idx(0), bf_p );
   mirror_y( idx(1), bf );
   idx(2) = FF_mirror_y( idx(2) );
   idx(3) = FF_mirror_y( idx(3) );
   idx(4) = FF_mirror_y( idx(4) );
   idx(5) = FF_mirror_y( idx(5) );
   bf_sign_x &= bool(bf[0]);
   bf_sign_y &= bool(bf[1]);
   bf_sign_p_x &= bool(bf_p[0]);
   bf_sign_p_y &= bool(bf_p[1]);

   return operation( sign ^((bf_sign_x^bf_sign_y)^(bf_sign_p_x^bf_sign_p_y)), false );
}

operation mirror_diagonal_projmat( idx_proj_matrix_t& idx )
{
   bool sign= false;
   bool bf_sign_x = false;
   bool bf_sign_y = false;
   bool bf_sign_p_x = false;
   bool bf_sign_p_y = false;
   mirror_diagonal( idx(0) );
   mirror_diagonal( idx(1) );
   idx(2) = FF_mirror_diagonal( idx(2) ); // x <-> y is the same as the rotation of form factors around z 
   idx(3) = FF_mirror_diagonal( idx(3) ); // except that the signs multiplied are ALWAYS +
   idx(4) = FF_mirror_diagonal( idx(4) );
   idx(5) = FF_mirror_diagonal( idx(5) );

   return operation( sign ^((bf_sign_x^bf_sign_y)^(bf_sign_p_x^bf_sign_p_y)), false );
}

