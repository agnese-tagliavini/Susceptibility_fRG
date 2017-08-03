
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

operation rot_pi2_x_phi( idx_phi_t& idx )
{
   bool sign = (1- FF_rot_sign_pi2_x( idx(4) ) * FF_rot_sign_pi2_x( idx(5) ))/2;
   bool bf_sign_x = (1 - (translate_2pi_x(idx(4))*translate_2pi_x(idx(5))))/2;
   bool bf_sign_y = (1 - (translate_2pi_y(idx(4))*translate_2pi_y(idx(5))))/2;
   bool bf_sign_z = (1 - (translate_2pi_z(idx(4))*translate_2pi_z(idx(5))))/2;
   int bf[3] = {0,0,0};
   rot_pi2_x( idx(3), bf );
   idx(4)= FF_rot_pi2_x( idx(4) );
   idx(5)= FF_rot_pi2_x( idx(5) );
   bf_sign_x &= bool(bf[0]);
   bf_sign_y &= bool(bf[1]);
   bf_sign_z &= bool(bf[2]);

   return operation(sign^(bf_sign_x^bf_sign_y^bf_sign_z), false );
}

operation rot_pi2_y_phi( idx_phi_t& idx )
{
   bool sign = (1- FF_rot_sign_pi2_y( idx(4) ) * FF_rot_sign_pi2_y( idx(5) ))/2;
   bool bf_sign_x = (1 - (translate_2pi_x(idx(4))*translate_2pi_x(idx(5))))/2;
   bool bf_sign_y = (1 - (translate_2pi_y(idx(4))*translate_2pi_y(idx(5))))/2;
   bool bf_sign_z = (1 - (translate_2pi_z(idx(4))*translate_2pi_z(idx(5))))/2;
   int bf[3] = {0,0,0};
   rot_pi2_y( idx(3), bf );
   idx(4)= FF_rot_pi2_y( idx(4) );
   idx(5)= FF_rot_pi2_y( idx(5) );
   bf_sign_x &= bool(bf[0]);
   bf_sign_y &= bool(bf[1]);
   bf_sign_z &= bool(bf[2]);

   return operation(sign^(bf_sign_x^bf_sign_y^bf_sign_z), false );
}
operation rot_pi2_z_phi( idx_phi_t& idx )
{
   bool sign = (1- FF_rot_sign_pi2_z( idx(4) ) * FF_rot_sign_pi2_z( idx(5) ))/2;
   bool bf_sign_x = (1 - (translate_2pi_x(idx(4))*translate_2pi_x(idx(5))))/2;
   bool bf_sign_y = (1 - (translate_2pi_y(idx(4))*translate_2pi_y(idx(5))))/2;
   bool bf_sign_z = (1 - (translate_2pi_z(idx(4))*translate_2pi_z(idx(5))))/2;
   int bf[3] = {0,0,0};
   rot_pi2_z( idx(3), bf );
   idx(4)= FF_rot_pi2_z( idx(4) );
   idx(5)= FF_rot_pi2_z( idx(5) );
   bf_sign_x &= bool(bf[0]);
   bf_sign_y &= bool(bf[1]);
   bf_sign_z &= bool(bf[2]);

   return operation(sign^(bf_sign_x^bf_sign_y^bf_sign_z), false );
}
   

operation rot_2pi3_111_phi( idx_phi_t& idx )
{
   bool sign = (1- FF_rot_sign_2pi3_111( idx(4) ) * FF_rot_sign_2pi3_111( idx(5) ))/2;
   bool bf_sign_x = (1 - (translate_2pi_x(idx(4))*translate_2pi_x(idx(5))))/2;
   bool bf_sign_y = (1 - (translate_2pi_y(idx(4))*translate_2pi_y(idx(5))))/2;
   bool bf_sign_z = (1 - (translate_2pi_z(idx(4))*translate_2pi_z(idx(5))))/2;
   int bf[3] = {0,0,0};
   rot_2pi3_111( idx(3), bf );
   idx(4)= FF_rot_2pi3_111( idx(4) );
   idx(5)= FF_rot_2pi3_111( idx(5) );
   bf_sign_x &= bool(bf[0]);
   bf_sign_y &= bool(bf[1]);
   bf_sign_z &= bool(bf[2]);

   return operation(sign^(bf_sign_x^bf_sign_y^bf_sign_z), false );
}

   
operation rot_2pi3_1m11_phi( idx_phi_t& idx )
{
   bool sign = (1- FF_rot_sign_2pi3_1m11( idx(4) ) * FF_rot_sign_2pi3_1m11( idx(5) ))/2;
   bool bf_sign_x = (1 - (translate_2pi_x(idx(4))*translate_2pi_x(idx(5))))/2;
   bool bf_sign_y = (1 - (translate_2pi_y(idx(4))*translate_2pi_y(idx(5))))/2;
   bool bf_sign_z = (1 - (translate_2pi_z(idx(4))*translate_2pi_z(idx(5))))/2;
   int bf[3] = {0,0,0};
   rot_2pi3_1m11( idx(3), bf );
   idx(4)= FF_rot_2pi3_1m11( idx(4) );
   idx(5)= FF_rot_2pi3_1m11( idx(5) );
   bf_sign_x &= bool(bf[0]);
   bf_sign_y &= bool(bf[1]);
   bf_sign_z &= bool(bf[2]);

   return operation(sign^(bf_sign_x^bf_sign_y^bf_sign_z), false );
}


operation rot_2pi3_m111_phi( idx_phi_t& idx )
{
   bool sign =(1- FF_rot_sign_2pi3_m111( idx(4) ) * FF_rot_sign_2pi3_m111( idx(5) ))/2;
   bool bf_sign_x = (1 - (translate_2pi_x(idx(4))*translate_2pi_x(idx(5))))/2;
   bool bf_sign_y = (1 - (translate_2pi_y(idx(4))*translate_2pi_y(idx(5))))/2;
   bool bf_sign_z = (1 - (translate_2pi_z(idx(4))*translate_2pi_z(idx(5))))/2;
   int bf[3] = {0,0,0};
   rot_2pi3_m111( idx(3), bf );
   idx(4)= FF_rot_2pi3_m111( idx(4) );
   idx(5)= FF_rot_2pi3_m111( idx(5) );
   bf_sign_x &= bool(bf[0]);
   bf_sign_y &= bool(bf[1]);
   bf_sign_z &= bool(bf[2]);

   return operation(sign^(bf_sign_x^bf_sign_y^bf_sign_z), false );
}


operation rot_2pi3_m1m11_phi( idx_phi_t& idx )
{
   bool sign =(1- FF_rot_sign_2pi3_m1m11( idx(4) ) * FF_rot_sign_2pi3_m1m11( idx(5) ))/2;
   bool bf_sign_x = (1 - (translate_2pi_x(idx(4))*translate_2pi_x(idx(5))))/2;
   bool bf_sign_y = (1 - (translate_2pi_y(idx(4))*translate_2pi_y(idx(5))))/2;
   bool bf_sign_z = (1 - (translate_2pi_z(idx(4))*translate_2pi_z(idx(5))))/2;
   int bf[3] = {0,0,0};
   rot_2pi3_m1m11( idx(3), bf );
   idx(4)= FF_rot_2pi3_m1m11( idx(4) );
   idx(5)= FF_rot_2pi3_m1m11( idx(5) );
   bf_sign_x &= bool(bf[0]);
   bf_sign_y &= bool(bf[1]);
   bf_sign_z &= bool(bf[2]);

   return operation(sign^(bf_sign_x^bf_sign_y^bf_sign_z), false );
}


operation rot_pi_110_phi( idx_phi_t& idx )
{
   bool sign =(1- FF_rot_sign_pi_110( idx(4) ) * FF_rot_sign_pi_110( idx(5) ))/2;
   bool bf_sign_x = (1 - (translate_2pi_x(idx(4))*translate_2pi_x(idx(5))))/2;
   bool bf_sign_y = (1 - (translate_2pi_y(idx(4))*translate_2pi_y(idx(5))))/2;
   bool bf_sign_z = (1 - (translate_2pi_z(idx(4))*translate_2pi_z(idx(5))))/2;
   int bf[3] = {0,0,0};
   rot_pi_110( idx(3), bf );
   idx(4)= FF_rot_pi_110( idx(4) );
   idx(5)= FF_rot_pi_110( idx(5) );
   bf_sign_x &= bool(bf[0]);
   bf_sign_y &= bool(bf[1]);
   bf_sign_z &= bool(bf[2]);

   return operation(sign^(bf_sign_x^bf_sign_y^bf_sign_z), false );
}

operation rot_pi_1m10_phi( idx_phi_t& idx )
{
   bool sign = (1- FF_rot_sign_pi_1m10( idx(4) ) * FF_rot_sign_pi_1m10( idx(5) ))/2;
   bool bf_sign_x = (1 - (translate_2pi_x(idx(4))*translate_2pi_x(idx(5))))/2;
   bool bf_sign_y = (1 - (translate_2pi_y(idx(4))*translate_2pi_y(idx(5))))/2;
   bool bf_sign_z = (1 - (translate_2pi_z(idx(4))*translate_2pi_z(idx(5))))/2;
   int bf[3] = {0,0,0};
   rot_pi_1m10( idx(3), bf );
   idx(4)= FF_rot_pi_1m10( idx(4) );
   idx(5)= FF_rot_pi_1m10( idx(5) );
   bf_sign_x &= bool(bf[0]);
   bf_sign_y &= bool(bf[1]);
   bf_sign_z &= bool(bf[2]);

   return operation(sign^(bf_sign_x^bf_sign_y^bf_sign_z), false );
}

operation rot_pi_011_phi( idx_phi_t& idx )
{
   bool sign =(1- FF_rot_sign_pi_011( idx(4) ) * FF_rot_sign_pi_011( idx(5) ))/2;
   bool bf_sign_x = (1 - (translate_2pi_x(idx(4))*translate_2pi_x(idx(5))))/2;
   bool bf_sign_y = (1 - (translate_2pi_y(idx(4))*translate_2pi_y(idx(5))))/2;
   bool bf_sign_z = (1 - (translate_2pi_z(idx(4))*translate_2pi_z(idx(5))))/2;
   int bf[3] = {0,0,0};
   rot_pi_011( idx(3), bf );
   idx(4)= FF_rot_pi_011( idx(4) );
   idx(5)= FF_rot_pi_011( idx(5) );
   bf_sign_x &= bool(bf[0]);
   bf_sign_y &= bool(bf[1]);
   bf_sign_z &= bool(bf[2]);

   return operation(sign^(bf_sign_x^bf_sign_y^bf_sign_z), false );
}


operation rot_pi_01m1_phi( idx_phi_t& idx )
{
   bool sign = (1- FF_rot_sign_pi_01m1( idx(4) ) * FF_rot_sign_pi_01m1( idx(5) ))/2;
   bool bf_sign_x = (1 - (translate_2pi_x(idx(4))*translate_2pi_x(idx(5))))/2;
   bool bf_sign_y = (1 - (translate_2pi_y(idx(4))*translate_2pi_y(idx(5))))/2;
   bool bf_sign_z = (1 - (translate_2pi_z(idx(4))*translate_2pi_z(idx(5))))/2;
   int bf[3] = {0,0,0};
   rot_pi_01m1( idx(3), bf );
   idx(4)= FF_rot_pi_01m1( idx(4) );
   idx(5)= FF_rot_pi_01m1( idx(5) );
   bf_sign_x &= bool(bf[0]);
   bf_sign_y &= bool(bf[1]);
   bf_sign_z &= bool(bf[2]);

   return operation(sign^(bf_sign_x^bf_sign_y^bf_sign_z), false );
}

operation rot_pi_101_phi( idx_phi_t& idx )
{
   bool sign =(1- FF_rot_sign_pi_101( idx(4) ) * FF_rot_sign_pi_101( idx(5) ))/2;
   bool bf_sign_x = (1 - (translate_2pi_x(idx(4))*translate_2pi_x(idx(5))))/2;
   bool bf_sign_y = (1 - (translate_2pi_y(idx(4))*translate_2pi_y(idx(5))))/2;
   bool bf_sign_z = (1 - (translate_2pi_z(idx(4))*translate_2pi_z(idx(5))))/2;
   int bf[3] = {0,0,0};
   rot_pi_101( idx(3), bf );
   idx(4)= FF_rot_pi_101( idx(4) );
   idx(5)= FF_rot_pi_101( idx(5) );
   bf_sign_x &= bool(bf[0]);
   bf_sign_y &= bool(bf[1]);
   bf_sign_z &= bool(bf[2]);

   return operation(sign^(bf_sign_x^bf_sign_y^bf_sign_z), false );
}

operation rot_pi_10m1_phi( idx_phi_t& idx )
{
   bool sign = (1- FF_rot_sign_pi_10m1( idx(4) ) * FF_rot_sign_pi_10m1( idx(5) ))/2;
   bool bf_sign_x = (1 - (translate_2pi_x(idx(4))*translate_2pi_x(idx(5))))/2;
   bool bf_sign_y = (1 - (translate_2pi_y(idx(4))*translate_2pi_y(idx(5))))/2;
   bool bf_sign_z = (1 - (translate_2pi_z(idx(4))*translate_2pi_z(idx(5))))/2;
   int bf[3] = {0,0,0};
   rot_pi_10m1( idx(3), bf );
   idx(4)= FF_rot_pi_10m1( idx(4) );
   idx(5)= FF_rot_pi_10m1( idx(5) );
   bf_sign_x &= bool(bf[0]);
   bf_sign_y &= bool(bf[1]);
   bf_sign_z &= bool(bf[2]);

   return operation(sign^(bf_sign_x^bf_sign_y^bf_sign_z), false );
}

operation mirror_y_phi( idx_phi_t& idx )
{
   bool sign = (1- FF_mirror_y_sign( idx(4) ) * FF_mirror_y_sign( idx(5) ))/2;
   bool bf_sign_x = (1 - (translate_2pi_x(idx(4))*translate_2pi_x(idx(5))))/2;
   bool bf_sign_y = (1 - (translate_2pi_y(idx(4))*translate_2pi_y(idx(5))))/2;
   bool bf_sign_z = (1 - (translate_2pi_z(idx(4))*translate_2pi_z(idx(5))))/2;
   int bf[3] = {0,0,0};
   mirror_y( idx(3), bf );
   idx(4)= FF_mirror_y( idx(4) );
   idx(5)= FF_mirror_y( idx(5) );
   bf_sign_x &= bool(bf[0]);
   bf_sign_y &= bool(bf[1]);
   bf_sign_z &= bool(bf[2]);

   return operation(sign^(bf_sign_x^bf_sign_y^bf_sign_z), false );
}
// ---- P ----

operation rot_pi2_x_P( idx_P_t& idx )
{
   bool sign =(1- FF_rot_sign_pi2_x( idx(3) ) )/2;
   bool bf_sign_x = (1 - (translate_2pi_x(idx(3))))/2;
   bool bf_sign_y = (1 - (translate_2pi_y(idx(3))))/2;
   bool bf_sign_z = (1 - (translate_2pi_z(idx(3))))/2;
   int bf[3] = {0,0,0};
   rot_pi2_x( idx(2), bf );
   idx(3) = FF_rot_pi2_x( idx(3) );
   bf_sign_x &= bool(bf[0]);
   bf_sign_y &= bool(bf[1]);
   bf_sign_z &= bool(bf[2]);

   return operation(sign^(bf_sign_x^bf_sign_y^bf_sign_z), false );
}

operation rot_pi2_y_P( idx_P_t& idx )
{
   bool sign =(1- FF_rot_sign_pi2_y( idx(3) ) )/2;
   bool bf_sign_x = (1 - (translate_2pi_x(idx(3))))/2;
   bool bf_sign_y = (1 - (translate_2pi_y(idx(3))))/2;
   bool bf_sign_z = (1 - (translate_2pi_z(idx(3))))/2;
   int bf[3] = {0,0,0};
   rot_pi2_y( idx(2), bf );
   idx(3) = FF_rot_pi2_y( idx(3) );
   bf_sign_x &= bool(bf[0]);
   bf_sign_y &= bool(bf[1]);
   bf_sign_z &= bool(bf[2]);

   return operation(sign^(bf_sign_x^bf_sign_y^bf_sign_z), false );
}
operation rot_pi2_z_P( idx_P_t& idx )
{
   bool sign =(1- FF_rot_sign_pi2_z( idx(3) ))/2;
   bool bf_sign_x = (1 - (translate_2pi_x(idx(3))))/2;
   bool bf_sign_y = (1 - (translate_2pi_y(idx(3))))/2;
   bool bf_sign_z = (1 - (translate_2pi_z(idx(3))))/2;
   int bf[3] = {0,0,0};
   rot_pi2_z( idx(2), bf );
   bf_sign_x &= bool(bf[0]);
   bf_sign_y &= bool(bf[1]);
   bf_sign_z &= bool(bf[2]);

   return operation(sign^(bf_sign_x^bf_sign_y^bf_sign_z), false );
}   

operation rot_2pi3_111_P( idx_P_t& idx )
{
   bool sign =(1- FF_rot_sign_2pi3_111( idx(3) ) )/2;
   bool bf_sign_x = (1 - (translate_2pi_x(idx(3))))/2;
   bool bf_sign_y = (1 - (translate_2pi_y(idx(3))))/2;
   bool bf_sign_z = (1 - (translate_2pi_z(idx(3))))/2;
   int bf[3] = {0,0,0};
   rot_2pi3_111( idx(2), bf );
   idx(3) = FF_rot_2pi3_111( idx(3) );
   bf_sign_x &= bool(bf[0]);
   bf_sign_y &= bool(bf[1]);
   bf_sign_z &= bool(bf[2]);

   return operation(sign^(bf_sign_x^bf_sign_y^bf_sign_z), false );
}

   
operation rot_2pi3_1m11_P( idx_P_t& idx )
{
   bool sign =(1- FF_rot_sign_2pi3_1m11( idx(3) ) )/2;
   bool bf_sign_x = (1 - (translate_2pi_x(idx(3))))/2;
   bool bf_sign_y = (1 - (translate_2pi_y(idx(3))))/2;
   bool bf_sign_z = (1 - (translate_2pi_z(idx(3))))/2;
   int bf[3] = {0,0,0};
   rot_2pi3_1m11( idx(2), bf );
   idx(3) = FF_rot_2pi3_1m11( idx(3) );
   bf_sign_x &= bool(bf[0]);
   bf_sign_y &= bool(bf[1]);
   bf_sign_z &= bool(bf[2]);

   return operation(sign^(bf_sign_x^bf_sign_y^bf_sign_z), false );
}


operation rot_2pi3_m111_P( idx_P_t& idx )
{
   bool sign = (1- FF_rot_sign_2pi3_m111( idx(3) ) )/2;
   bool bf_sign_x = (1 - (translate_2pi_x(idx(3))))/2;
   bool bf_sign_y = (1 - (translate_2pi_y(idx(3))))/2;
   bool bf_sign_z = (1 - (translate_2pi_z(idx(3))))/2;
   int bf[3] = {0,0,0};
   rot_2pi3_m111( idx(2), bf );
   idx(3) = FF_rot_2pi3_m111( idx(3) );
   bf_sign_x &= bool(bf[0]);
   bf_sign_y &= bool(bf[1]);
   bf_sign_z &= bool(bf[2]);

   return operation(sign^(bf_sign_x^bf_sign_y^bf_sign_z), false );
}


operation rot_2pi3_m1m11_P( idx_P_t& idx )
{
   bool sign =(1- FF_rot_sign_2pi3_m1m11( idx(3) ) )/2;
   bool bf_sign_x = (1 - (translate_2pi_x(idx(3))))/2;
   bool bf_sign_y = (1 - (translate_2pi_y(idx(3))))/2;
   bool bf_sign_z = (1 - (translate_2pi_z(idx(3))))/2;
   int bf[3] = {0,0,0};
   rot_2pi3_m1m11( idx(2), bf );
   idx(3) = FF_rot_2pi3_m1m11( idx(3) );
   bf_sign_x &= bool(bf[0]);
   bf_sign_y &= bool(bf[1]);
   bf_sign_z &= bool(bf[2]);

   return operation(sign^(bf_sign_x^bf_sign_y^bf_sign_z), false );
}


operation rot_pi_110_P( idx_P_t& idx )
{
   bool sign =(1- FF_rot_sign_pi_110( idx(3) ) )/2;
   bool bf_sign_x = (1 - (translate_2pi_x(idx(3))))/2;
   bool bf_sign_y = (1 - (translate_2pi_y(idx(3))))/2;
   bool bf_sign_z = (1 - (translate_2pi_z(idx(3))))/2;
   int bf[3] = {0,0,0};
   rot_pi_110( idx(2), bf );
   idx(3) = FF_rot_pi_110( idx(3) );
   bf_sign_x &= bool(bf[0]);
   bf_sign_y &= bool(bf[1]);
   bf_sign_z &= bool(bf[2]);

   return operation(sign^(bf_sign_x^bf_sign_y^bf_sign_z), false );
}

operation rot_pi_1m10_P( idx_P_t& idx )
{
   bool sign =(1- FF_rot_sign_pi_1m10( idx(3) ) )/2;
   bool bf_sign_x = (1 - (translate_2pi_x(idx(3))))/2;
   bool bf_sign_y = (1 - (translate_2pi_y(idx(3))))/2;
   bool bf_sign_z = (1 - (translate_2pi_z(idx(3))))/2;
   int bf[3] = {0,0,0};
   rot_pi_1m10( idx(2), bf );
   idx(3) = FF_rot_pi_1m10( idx(3) );
   bf_sign_x &= bool(bf[0]);
   bf_sign_y &= bool(bf[1]);
   bf_sign_z &= bool(bf[2]);

   return operation(sign^(bf_sign_x^bf_sign_y^bf_sign_z), false );
}

operation rot_pi_011_P( idx_P_t& idx )
{
   bool sign =(1- FF_rot_sign_pi_011( idx(3) ) )/2;
   bool bf_sign_x = (1 - (translate_2pi_x(idx(3))))/2;
   bool bf_sign_y = (1 - (translate_2pi_y(idx(3))))/2;
   bool bf_sign_z = (1 - (translate_2pi_z(idx(3))))/2;
   int bf[3] = {0,0,0};
   rot_pi_011( idx(2), bf );
   idx(3) = FF_rot_pi_011( idx(3) );
   bf_sign_x &= bool(bf[0]);
   bf_sign_y &= bool(bf[1]);
   bf_sign_z &= bool(bf[2]);

   return operation(sign^(bf_sign_x^bf_sign_y^bf_sign_z), false );
}


operation rot_pi_01m1_P( idx_P_t& idx )
{
   bool sign =(1- FF_rot_sign_pi_01m1( idx(3) ) )/2;
   bool bf_sign_x = (1 - (translate_2pi_x(idx(3))))/2;
   bool bf_sign_y = (1 - (translate_2pi_y(idx(3))))/2;
   bool bf_sign_z = (1 - (translate_2pi_z(idx(3))))/2;
   int bf[3] = {0,0,0};
   rot_pi_01m1( idx(2), bf );
   idx(3) = FF_rot_pi_01m1( idx(3) );
   bf_sign_x &= bool(bf[0]);
   bf_sign_y &= bool(bf[1]);
   bf_sign_z &= bool(bf[2]);

   return operation(sign^(bf_sign_x^bf_sign_y^bf_sign_z), false );
}

operation rot_pi_101_P( idx_P_t& idx )
{
   bool sign =(1- FF_rot_sign_pi_101( idx(3) ) )/2;
   bool bf_sign_x = (1 - (translate_2pi_x(idx(3))))/2;
   bool bf_sign_y = (1 - (translate_2pi_y(idx(3))))/2;
   bool bf_sign_z = (1 - (translate_2pi_z(idx(3))))/2;
   int bf[3] = {0,0,0};
   rot_pi_101( idx(2), bf );
   idx(3) = FF_rot_pi_101( idx(3) );
   bf_sign_x &= bool(bf[0]);
   bf_sign_y &= bool(bf[1]);
   bf_sign_z &= bool(bf[2]);

   return operation(sign^(bf_sign_x^bf_sign_y^bf_sign_z), false );
}

operation rot_pi_10m1_P( idx_P_t& idx )
{
   bool sign = (1- FF_rot_sign_pi_10m1( idx(3) ) )/2;
   bool bf_sign_x = (1 - (translate_2pi_x(idx(3))))/2;
   bool bf_sign_y = (1 - (translate_2pi_y(idx(3))))/2;
   bool bf_sign_z = (1 - (translate_2pi_z(idx(3))))/2;
   int bf[3] = {0,0,0};
   rot_pi_10m1( idx(2), bf );
   idx(3) = FF_rot_pi_10m1( idx(3) );
   bf_sign_x &= bool(bf[0]);
   bf_sign_y &= bool(bf[1]);
   bf_sign_z &= bool(bf[2]);

   return operation(sign^(bf_sign_x^bf_sign_y^bf_sign_z), false );
}

operation mirror_y_P( idx_P_t& idx )
{
   bool sign = (1- FF_mirror_y_sign( idx(3) ) )/2;
   bool bf_sign_x = (1 - (translate_2pi_x(idx(3))))/2;
   bool bf_sign_y = (1 - (translate_2pi_y(idx(3))))/2;
   bool bf_sign_z = (1 - (translate_2pi_z(idx(3))))/2;
   int bf[3] = {0,0,0};
   mirror_y( idx(2), bf );
   idx(3)= FF_mirror_y( idx(3) );
   bf_sign_x &= bool(bf[0]);
   bf_sign_y &= bool(bf[1]);
   bf_sign_z &= bool(bf[2]);

   return operation(sign^(bf_sign_x^bf_sign_y^bf_sign_z), false );
}

// ---- chi ----

operation rot_pi2_x_chi( idx_chi_t& idx )
{
   int bf[3] = {0,0,0};
   rot_pi2_x( idx(1), bf );

   return operation( false, false );
}

operation rot_pi2_y_chi( idx_chi_t& idx )
{
   int bf[3] = {0,0,0};
   rot_pi2_y( idx(1), bf );

   return operation( false, false );
}
operation rot_pi2_z_chi( idx_chi_t& idx )
{
   int bf[3] = {0,0,0};
   rot_pi2_z( idx(1), bf );

   return operation( false, false );
}
   

operation rot_2pi3_111_chi( idx_chi_t& idx )
{
   int bf[3] = {0,0,0};
   rot_2pi3_111( idx(1), bf );

   return operation( false, false );
}

   
operation rot_2pi3_1m11_chi( idx_chi_t& idx )
{
   int bf[3] = {0,0,0};
   rot_2pi3_1m11( idx(1), bf );

   return operation( false, false );
}


operation rot_2pi3_m111_chi( idx_chi_t& idx )
{
   int bf[3] = {0,0,0};
   rot_2pi3_m111( idx(1), bf );

   return operation( false, false );
}


operation rot_2pi3_m1m11_chi( idx_chi_t& idx )
{
   int bf[3] = {0,0,0};
   rot_2pi3_m1m11( idx(1), bf );

   return operation( false, false );
}


operation rot_pi_110_chi( idx_chi_t& idx )
{
   int bf[3] = {0,0,0};
   rot_pi_110( idx(1), bf );

   return operation( false, false );
}

operation rot_pi_1m10_chi( idx_chi_t& idx )
{
   int bf[3] = {0,0,0};
   rot_pi_1m10( idx(1), bf );

   return operation( false, false );
}

operation rot_pi_011_chi( idx_chi_t& idx )
{
   int bf[3] = {0,0,0};
   rot_pi_011( idx(1), bf );

   return operation( false, false );
}


operation rot_pi_01m1_chi( idx_chi_t& idx )
{
   int bf[3] = {0,0,0};
   rot_pi_01m1( idx(1), bf );

   return operation( false, false );
}

operation rot_pi_101_chi( idx_chi_t& idx )
{
   int bf[3] = {0,0,0};
   rot_pi_101( idx(1), bf );

   return operation( false, false );
}

operation rot_pi_10m1_chi( idx_chi_t& idx )
{
   int bf[3] = {0,0,0};
   rot_pi_10m1( idx(1), bf );

   return operation( false, false );
}

operation mirror_y_chi( idx_chi_t& idx )
{
   int bf[3] = {0,0,0};
   mirror_y( idx(1), bf );

   return operation(false, false );
}

// ---- suscept ---
operation rot_pi2_x_suscept( idx_suscept_t& idx )
{
   bool sign = (1- FF_rot_sign_pi2_x( idx(2) ) * FF_rot_sign_pi2_x( idx(3) ))/2;
   bool bf_sign_x = (1 - (translate_2pi_x(idx(2))*translate_2pi_x(idx(3))))/2;
   bool bf_sign_y = (1 - (translate_2pi_y(idx(2))*translate_2pi_y(idx(3))))/2;
   bool bf_sign_z = (1 - (translate_2pi_z(idx(2))*translate_2pi_z(idx(3))))/2;
   int bf[3] = {0,0,0};
   rot_pi2_x( idx(3), bf );
   idx(2)= FF_rot_pi2_x( idx(2) );
   idx(3)= FF_rot_pi2_x( idx(3) );
   bf_sign_x &= bool(bf[0]);
   bf_sign_y &= bool(bf[1]);
   bf_sign_z &= bool(bf[2]);

   return operation(sign^(bf_sign_x^bf_sign_y^bf_sign_z), false );
}

operation rot_pi2_y_suscept( idx_suscept_t& idx )
{
   bool sign = (1- FF_rot_sign_pi2_y( idx(2) ) * FF_rot_sign_pi2_y( idx(3) ))/2;
   bool bf_sign_x = (1 - (translate_2pi_x(idx(2))*translate_2pi_x(idx(3))))/2;
   bool bf_sign_y = (1 - (translate_2pi_y(idx(2))*translate_2pi_y(idx(3))))/2;
   bool bf_sign_z = (1 - (translate_2pi_z(idx(2))*translate_2pi_z(idx(3))))/2;
   int bf[3] = {0,0,0};
   rot_pi2_y( idx(3), bf );
   idx(2)= FF_rot_pi2_y( idx(2) );
   idx(3)= FF_rot_pi2_y( idx(3) );
   bf_sign_x &= bool(bf[0]);
   bf_sign_y &= bool(bf[1]);
   bf_sign_z &= bool(bf[2]);

   return operation(sign^(bf_sign_x^bf_sign_y^bf_sign_z), false );
}
operation rot_pi2_z_suscept( idx_suscept_t& idx )
{
   bool sign = (1- FF_rot_sign_pi2_z( idx(2) ) * FF_rot_sign_pi2_z( idx(3) ))/2;
   bool bf_sign_x = (1 - (translate_2pi_x(idx(2))*translate_2pi_x(idx(3))))/2;
   bool bf_sign_y = (1 - (translate_2pi_y(idx(2))*translate_2pi_y(idx(3))))/2;
   bool bf_sign_z = (1 - (translate_2pi_z(idx(2))*translate_2pi_z(idx(3))))/2;
   int bf[3] = {0,0,0};
   rot_pi2_z( idx(3), bf );
   idx(2)= FF_rot_pi2_z( idx(2) );
   idx(3)= FF_rot_pi2_z( idx(3) );
   bf_sign_x &= bool(bf[0]);
   bf_sign_y &= bool(bf[1]);
   bf_sign_z &= bool(bf[2]);

   return operation(sign^(bf_sign_x^bf_sign_y^bf_sign_z), false );
}
   

operation rot_2pi3_111_suscept( idx_suscept_t& idx )
{
   bool sign = (1- FF_rot_sign_2pi3_111( idx(2) ) * FF_rot_sign_2pi3_111( idx(3) ))/2;
   bool bf_sign_x = (1 - (translate_2pi_x(idx(2))*translate_2pi_x(idx(3))))/2;
   bool bf_sign_y = (1 - (translate_2pi_y(idx(2))*translate_2pi_y(idx(3))))/2;
   bool bf_sign_z = (1 - (translate_2pi_z(idx(2))*translate_2pi_z(idx(3))))/2;
   int bf[3] = {0,0,0};
   rot_2pi3_111( idx(3), bf );
   idx(2)= FF_rot_2pi3_111( idx(2) );
   idx(3)= FF_rot_2pi3_111( idx(3) );
   bf_sign_x &= bool(bf[0]);
   bf_sign_y &= bool(bf[1]);
   bf_sign_z &= bool(bf[2]);

   return operation(sign^(bf_sign_x^bf_sign_y^bf_sign_z), false );
}

   
operation rot_2pi3_1m11_suscept( idx_suscept_t& idx )
{
   bool sign = (1- FF_rot_sign_2pi3_1m11( idx(2) ) * FF_rot_sign_2pi3_1m11( idx(3) ))/2;
   bool bf_sign_x = (1 - (translate_2pi_x(idx(2))*translate_2pi_x(idx(3))))/2;
   bool bf_sign_y = (1 - (translate_2pi_y(idx(2))*translate_2pi_y(idx(3))))/2;
   bool bf_sign_z = (1 - (translate_2pi_z(idx(2))*translate_2pi_z(idx(3))))/2;
   int bf[3] = {0,0,0};
   rot_2pi3_1m11( idx(3), bf );
   idx(2)= FF_rot_2pi3_1m11( idx(2) );
   idx(3)= FF_rot_2pi3_1m11( idx(3) );
   bf_sign_x &= bool(bf[0]);
   bf_sign_y &= bool(bf[1]);
   bf_sign_z &= bool(bf[2]);

   return operation(sign^(bf_sign_x^bf_sign_y^bf_sign_z), false );
}


operation rot_2pi3_m111_suscept( idx_suscept_t& idx )
{
   bool sign =(1- FF_rot_sign_2pi3_m111( idx(2) ) * FF_rot_sign_2pi3_m111( idx(3) ))/2;
   bool bf_sign_x = (1 - (translate_2pi_x(idx(2))*translate_2pi_x(idx(3))))/2;
   bool bf_sign_y = (1 - (translate_2pi_y(idx(2))*translate_2pi_y(idx(3))))/2;
   bool bf_sign_z = (1 - (translate_2pi_z(idx(2))*translate_2pi_z(idx(3))))/2;
   int bf[3] = {0,0,0};
   rot_2pi3_m111( idx(3), bf );
   idx(2)= FF_rot_2pi3_m111( idx(2) );
   idx(3)= FF_rot_2pi3_m111( idx(3) );
   bf_sign_x &= bool(bf[0]);
   bf_sign_y &= bool(bf[1]);
   bf_sign_z &= bool(bf[2]);

   return operation(sign^(bf_sign_x^bf_sign_y^bf_sign_z), false );
}


operation rot_2pi3_m1m11_suscept( idx_suscept_t& idx )
{
   bool sign =(1- FF_rot_sign_2pi3_m1m11( idx(2) ) * FF_rot_sign_2pi3_m1m11( idx(3) ))/2;
   bool bf_sign_x = (1 - (translate_2pi_x(idx(2))*translate_2pi_x(idx(3))))/2;
   bool bf_sign_y = (1 - (translate_2pi_y(idx(2))*translate_2pi_y(idx(3))))/2;
   bool bf_sign_z = (1 - (translate_2pi_z(idx(2))*translate_2pi_z(idx(3))))/2;
   int bf[3] = {0,0,0};
   rot_2pi3_m1m11( idx(3), bf );
   idx(2)= FF_rot_2pi3_m1m11( idx(2) );
   idx(3)= FF_rot_2pi3_m1m11( idx(3) );
   bf_sign_x &= bool(bf[0]);
   bf_sign_y &= bool(bf[1]);
   bf_sign_z &= bool(bf[2]);

   return operation(sign^(bf_sign_x^bf_sign_y^bf_sign_z), false );
}


operation rot_pi_110_suscept( idx_suscept_t& idx )
{
   bool sign =(1- FF_rot_sign_pi_110( idx(2) ) * FF_rot_sign_pi_110( idx(3) ))/2;
   bool bf_sign_x = (1 - (translate_2pi_x(idx(2))*translate_2pi_x(idx(3))))/2;
   bool bf_sign_y = (1 - (translate_2pi_y(idx(2))*translate_2pi_y(idx(3))))/2;
   bool bf_sign_z = (1 - (translate_2pi_z(idx(2))*translate_2pi_z(idx(3))))/2;
   int bf[3] = {0,0,0};
   rot_pi_110( idx(3), bf );
   idx(2)= FF_rot_pi_110( idx(2) );
   idx(3)= FF_rot_pi_110( idx(3) );
   bf_sign_x &= bool(bf[0]);
   bf_sign_y &= bool(bf[1]);
   bf_sign_z &= bool(bf[2]);

   return operation(sign^(bf_sign_x^bf_sign_y^bf_sign_z), false );
}

operation rot_pi_1m10_suscept( idx_suscept_t& idx )
{
   bool sign = (1- FF_rot_sign_pi_1m10( idx(2) ) * FF_rot_sign_pi_1m10( idx(3) ))/2;
   bool bf_sign_x = (1 - (translate_2pi_x(idx(2))*translate_2pi_x(idx(3))))/2;
   bool bf_sign_y = (1 - (translate_2pi_y(idx(2))*translate_2pi_y(idx(3))))/2;
   bool bf_sign_z = (1 - (translate_2pi_z(idx(2))*translate_2pi_z(idx(3))))/2;
   int bf[3] = {0,0,0};
   rot_pi_1m10( idx(3), bf );
   idx(2)= FF_rot_pi_1m10( idx(2) );
   idx(3)= FF_rot_pi_1m10( idx(3) );
   bf_sign_x &= bool(bf[0]);
   bf_sign_y &= bool(bf[1]);
   bf_sign_z &= bool(bf[2]);

   return operation(sign^(bf_sign_x^bf_sign_y^bf_sign_z), false );
}

operation rot_pi_011_suscept( idx_suscept_t& idx )
{
   bool sign =(1- FF_rot_sign_pi_011( idx(2) ) * FF_rot_sign_pi_011( idx(3) ))/2;
   bool bf_sign_x = (1 - (translate_2pi_x(idx(2))*translate_2pi_x(idx(3))))/2;
   bool bf_sign_y = (1 - (translate_2pi_y(idx(2))*translate_2pi_y(idx(3))))/2;
   bool bf_sign_z = (1 - (translate_2pi_z(idx(2))*translate_2pi_z(idx(3))))/2;
   int bf[3] = {0,0,0};
   rot_pi_011( idx(3), bf );
   idx(2)= FF_rot_pi_011( idx(2) );
   idx(3)= FF_rot_pi_011( idx(3) );
   bf_sign_x &= bool(bf[0]);
   bf_sign_y &= bool(bf[1]);
   bf_sign_z &= bool(bf[2]);

   return operation(sign^(bf_sign_x^bf_sign_y^bf_sign_z), false );
}


operation rot_pi_01m1_suscept( idx_suscept_t& idx )
{
   bool sign = (1- FF_rot_sign_pi_01m1( idx(2) ) * FF_rot_sign_pi_01m1( idx(3) ))/2;
   bool bf_sign_x = (1 - (translate_2pi_x(idx(2))*translate_2pi_x(idx(3))))/2;
   bool bf_sign_y = (1 - (translate_2pi_y(idx(2))*translate_2pi_y(idx(3))))/2;
   bool bf_sign_z = (1 - (translate_2pi_z(idx(2))*translate_2pi_z(idx(3))))/2;
   int bf[3] = {0,0,0};
   rot_pi_01m1( idx(3), bf );
   idx(2)= FF_rot_pi_01m1( idx(2) );
   idx(3)= FF_rot_pi_01m1( idx(3) );
   bf_sign_x &= bool(bf[0]);
   bf_sign_y &= bool(bf[1]);
   bf_sign_z &= bool(bf[2]);

   return operation(sign^(bf_sign_x^bf_sign_y^bf_sign_z), false );
}

operation rot_pi_101_suscept( idx_suscept_t& idx )
{
   bool sign =(1- FF_rot_sign_pi_101( idx(2) ) * FF_rot_sign_pi_101( idx(3) ))/2;
   bool bf_sign_x = (1 - (translate_2pi_x(idx(2))*translate_2pi_x(idx(3))))/2;
   bool bf_sign_y = (1 - (translate_2pi_y(idx(2))*translate_2pi_y(idx(3))))/2;
   bool bf_sign_z = (1 - (translate_2pi_z(idx(2))*translate_2pi_z(idx(3))))/2;
   int bf[3] = {0,0,0};
   rot_pi_101( idx(3), bf );
   idx(2)= FF_rot_pi_101( idx(2) );
   idx(3)= FF_rot_pi_101( idx(3) );
   bf_sign_x &= bool(bf[0]);
   bf_sign_y &= bool(bf[1]);
   bf_sign_z &= bool(bf[2]);

   return operation(sign^(bf_sign_x^bf_sign_y^bf_sign_z), false );
}

operation rot_pi_10m1_suscept( idx_suscept_t& idx )
{
   bool sign = (1- FF_rot_sign_pi_10m1( idx(2) ) * FF_rot_sign_pi_10m1( idx(3) ))/2;
   bool bf_sign_x = (1 - (translate_2pi_x(idx(2))*translate_2pi_x(idx(3))))/2;
   bool bf_sign_y = (1 - (translate_2pi_y(idx(2))*translate_2pi_y(idx(3))))/2;
   bool bf_sign_z = (1 - (translate_2pi_z(idx(2))*translate_2pi_z(idx(3))))/2;
   int bf[3] = {0,0,0};
   rot_pi_10m1( idx(3), bf );
   idx(2)= FF_rot_pi_10m1( idx(2) );
   idx(3)= FF_rot_pi_10m1( idx(3) );
   bf_sign_x &= bool(bf[0]);
   bf_sign_y &= bool(bf[1]);
   bf_sign_z &= bool(bf[2]);

   return operation(sign^(bf_sign_x^bf_sign_y^bf_sign_z), false );
}

operation mirror_y_suscept( idx_suscept_t& idx )
{
   bool sign = (1- FF_mirror_y_sign( idx(2) ) * FF_mirror_y_sign( idx(3) ))/2;
   bool bf_sign_x = (1 - (translate_2pi_x(idx(2))*translate_2pi_x(idx(3))))/2;
   bool bf_sign_y = (1 - (translate_2pi_y(idx(2))*translate_2pi_y(idx(3))))/2;
   bool bf_sign_z = (1 - (translate_2pi_z(idx(2))*translate_2pi_z(idx(3))))/2;
   int bf[3] = {0,0,0};
   mirror_y( idx(3), bf );
   idx(2)= FF_mirror_y( idx(2) );
   idx(3)= FF_mirror_y( idx(3) );
   bf_sign_x &= bool(bf[0]);
   bf_sign_y &= bool(bf[1]);
   bf_sign_z &= bool(bf[2]);

   return operation(sign^(bf_sign_x^bf_sign_y^bf_sign_z), false );
}
// ---- Sig ----

operation rot_pi2_x_sig( idx_1p_t& idx )
{
   int bf[3] = {0,0,0};
   rot_pi2_x( idx(1), bf );

   return operation( false, false );
}

operation rot_pi2_y_sig( idx_1p_t& idx )
{
   int bf[3] = {0,0,0};
   rot_pi2_y( idx(1), bf );

   return operation( false, false );
}
operation rot_pi2_z_sig( idx_1p_t& idx )
{
   int bf[3] = {0,0,0};
   rot_pi2_z( idx(1), bf );

   return operation( false, false );
}
   

operation rot_2pi3_111_sig( idx_1p_t& idx )
{
   int bf[3] = {0,0,0};
   rot_2pi3_111( idx(1), bf );

   return operation( false, false );
}

   
operation rot_2pi3_1m11_sig( idx_1p_t& idx )
{
   int bf[3] = {0,0,0};
   rot_2pi3_1m11( idx(1), bf );

   return operation( false, false );
}


operation rot_2pi3_m111_sig( idx_1p_t& idx )
{
   int bf[3] = {0,0,0};
   rot_2pi3_m111( idx(1), bf );

   return operation( false, false );
}


operation rot_2pi3_m1m11_sig( idx_1p_t& idx )
{
   int bf[3] = {0,0,0};
   rot_2pi3_m1m11( idx(1), bf );

   return operation( false, false );
}


operation rot_pi_110_sig( idx_1p_t& idx )
{
   int bf[3] = {0,0,0};
   rot_pi_110( idx(1), bf );

   return operation( false, false );
}

operation rot_pi_1m10_sig( idx_1p_t& idx )
{
   int bf[3] = {0,0,0};
   rot_pi_1m10( idx(1), bf );

   return operation( false, false );
}

operation rot_pi_011_sig( idx_1p_t& idx )
{
   int bf[3] = {0,0,0};
   rot_pi_011( idx(1), bf );

   return operation( false, false );
}


operation rot_pi_01m1_sig( idx_1p_t& idx )
{
   int bf[3] = {0,0,0};
   rot_pi_01m1( idx(1), bf );

   return operation( false, false );
}

operation rot_pi_101_sig( idx_1p_t& idx )
{
   int bf[3] = {0,0,0};
   rot_pi_101( idx(1), bf );

   return operation( false, false );
}

operation rot_pi_10m1_sig( idx_1p_t& idx )
{
   int bf[3] = {0,0,0};
   rot_pi_10m1( idx(1), bf );

   return operation( false, false );
}

operation mirror_y_sig( idx_1p_t& idx )
{
   int bf[3] = {0,0,0};
   mirror_y( idx(1), bf );

   return operation(false, false );
}

// -- BARE VERTEX PROJECTED
// not used because local bare interaction is really easy in FF
/*******************************************************************
 *
 * 			SYMMETRIES PROJECTION MATRICES
 *
 ******************************************************************/

operation rot_pi2_x_projmat( idx_proj_matrix_t& idx )
{
   bool sign= ((1- (FF_rot_sign_pi2_x( idx(4) ) * FF_rot_sign_pi2_x( idx(5) ) * FF_rot_sign_pi2_x( idx(2) ) * FF_rot_sign_pi2_x( idx(3) )))/2);
   bool bf_sign_x = (1 - (translate_2pi_x(idx(4))*translate_2pi_x(idx(5))))/2;
   bool bf_sign_y = (1 - (translate_2pi_y(idx(4))*translate_2pi_y(idx(5))))/2;
   bool bf_sign_z = (1 - (translate_2pi_z(idx(4))*translate_2pi_z(idx(5))))/2;
   bool bf_sign_p_x = (1 - (translate_2pi_x(idx(2))*translate_2pi_x(idx(3))))/2;
   bool bf_sign_p_y = (1 - (translate_2pi_y(idx(2))*translate_2pi_y(idx(3))))/2;
   bool bf_sign_p_z = (1 - (translate_2pi_z(idx(2))*translate_2pi_z(idx(3))))/2;
   int bf[3] = {0,0,0};
   int bf_p[3] = {0,0,0};
   rot_pi2_x( idx(0), bf_p );
   rot_pi2_x( idx(1), bf );
   idx(2) = FF_rot_pi2_x( idx(2) );
   idx(3) = FF_rot_pi2_x( idx(3) );
   idx(4) = FF_rot_pi2_x( idx(4) );
   idx(5) = FF_rot_pi2_x( idx(5) );
   bf_sign_x &= bool(bf[0]);
   bf_sign_y &= bool(bf[1]);
   bf_sign_p_x &= bool(bf_p[0]);
   bf_sign_p_y &= bool(bf_p[1]);

   return operation( sign ^((bf_sign_x^bf_sign_y^bf_sign_z)^(bf_sign_p_x^bf_sign_p_y^bf_sign_p_z)), false );
}

operation rot_pi2_y_projmat( idx_proj_matrix_t& idx )
{
   bool sign= ((1- (FF_rot_sign_pi2_y( idx(4) ) * FF_rot_sign_pi2_y( idx(5) ) * FF_rot_sign_pi2_y( idx(2) ) * FF_rot_sign_pi2_y( idx(3) )))/2);
   bool bf_sign_x = (1 - (translate_2pi_x(idx(4))*translate_2pi_x(idx(5))))/2;
   bool bf_sign_y = (1 - (translate_2pi_y(idx(4))*translate_2pi_y(idx(5))))/2;
   bool bf_sign_z = (1 - (translate_2pi_z(idx(4))*translate_2pi_z(idx(5))))/2;
   bool bf_sign_p_x = (1 - (translate_2pi_x(idx(2))*translate_2pi_x(idx(3))))/2;
   bool bf_sign_p_y = (1 - (translate_2pi_y(idx(2))*translate_2pi_y(idx(3))))/2;
   bool bf_sign_p_z = (1 - (translate_2pi_z(idx(2))*translate_2pi_z(idx(3))))/2;
   int bf[3] = {0,0,0};
   int bf_p[3] = {0,0,0};
   rot_pi2_y( idx(0), bf_p );
   rot_pi2_y( idx(1), bf );
   idx(2) = FF_rot_pi2_y( idx(2) );
   idx(3) = FF_rot_pi2_y( idx(3) );
   idx(4) = FF_rot_pi2_y( idx(4) );
   idx(5) = FF_rot_pi2_y( idx(5) );
   bf_sign_x &= bool(bf[0]);
   bf_sign_y &= bool(bf[1]);
   bf_sign_p_x &= bool(bf_p[0]);
   bf_sign_p_y &= bool(bf_p[1]);

   return operation( sign ^((bf_sign_x^bf_sign_y^bf_sign_z)^(bf_sign_p_x^bf_sign_p_y^bf_sign_p_z)), false );
}
operation rot_pi2_z_projmat( idx_proj_matrix_t& idx )
{
   bool sign= ((1- (FF_rot_sign_pi2_z( idx(4) ) * FF_rot_sign_pi2_z( idx(5) ) * FF_rot_sign_pi2_z( idx(2) ) * FF_rot_sign_pi2_z( idx(3) )))/2);
   bool bf_sign_x = (1 - (translate_2pi_x(idx(4))*translate_2pi_x(idx(5))))/2;
   bool bf_sign_y = (1 - (translate_2pi_y(idx(4))*translate_2pi_y(idx(5))))/2;
   bool bf_sign_z = (1 - (translate_2pi_z(idx(4))*translate_2pi_z(idx(5))))/2;
   bool bf_sign_p_x = (1 - (translate_2pi_x(idx(2))*translate_2pi_x(idx(3))))/2;
   bool bf_sign_p_y = (1 - (translate_2pi_y(idx(2))*translate_2pi_y(idx(3))))/2;
   bool bf_sign_p_z = (1 - (translate_2pi_z(idx(2))*translate_2pi_z(idx(3))))/2;
   int bf[3] = {0,0,0};
   int bf_p[3] = {0,0,0};
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

   return operation( sign ^((bf_sign_x^bf_sign_y^bf_sign_z)^(bf_sign_p_x^bf_sign_p_y^bf_sign_p_z)), false );
}
   

operation rot_2pi3_111_projmat( idx_proj_matrix_t& idx )
{
   bool sign= ((1- (FF_rot_sign_2pi3_111( idx(4) ) * FF_rot_sign_2pi3_111( idx(5) ) * FF_rot_sign_2pi3_111( idx(2) ) * FF_rot_sign_2pi3_111( idx(3) )))/2);
   bool bf_sign_x = (1 - (translate_2pi_x(idx(4))*translate_2pi_x(idx(5))))/2;
   bool bf_sign_y = (1 - (translate_2pi_y(idx(4))*translate_2pi_y(idx(5))))/2;
   bool bf_sign_z = (1 - (translate_2pi_z(idx(4))*translate_2pi_z(idx(5))))/2;
   bool bf_sign_p_x = (1 - (translate_2pi_x(idx(2))*translate_2pi_x(idx(3))))/2;
   bool bf_sign_p_y = (1 - (translate_2pi_y(idx(2))*translate_2pi_y(idx(3))))/2;
   bool bf_sign_p_z = (1 - (translate_2pi_z(idx(2))*translate_2pi_z(idx(3))))/2;
   int bf[3] = {0,0,0};
   int bf_p[3] = {0,0,0};
   rot_2pi3_111( idx(0), bf_p );
   rot_2pi3_111( idx(1), bf );
   idx(2) = FF_rot_2pi3_111( idx(2) );
   idx(3) = FF_rot_2pi3_111( idx(3) );
   idx(4) = FF_rot_2pi3_111( idx(4) );
   idx(5) = FF_rot_2pi3_111( idx(5) );
   bf_sign_x &= bool(bf[0]);
   bf_sign_y &= bool(bf[1]);
   bf_sign_p_x &= bool(bf_p[0]);
   bf_sign_p_y &= bool(bf_p[1]);

   return operation( sign ^((bf_sign_x^bf_sign_y^bf_sign_z)^(bf_sign_p_x^bf_sign_p_y^bf_sign_p_z)), false );
}

   
operation rot_2pi3_1m11_projmat( idx_proj_matrix_t& idx )
{
   bool sign= ((1- (FF_rot_sign_2pi3_1m11( idx(4) ) * FF_rot_sign_2pi3_1m11( idx(5) ) * FF_rot_sign_2pi3_1m11( idx(2) ) * FF_rot_sign_2pi3_1m11( idx(3) )))/2);
   bool bf_sign_x = (1 - (translate_2pi_x(idx(4))*translate_2pi_x(idx(5))))/2;
   bool bf_sign_y = (1 - (translate_2pi_y(idx(4))*translate_2pi_y(idx(5))))/2;
   bool bf_sign_z = (1 - (translate_2pi_z(idx(4))*translate_2pi_z(idx(5))))/2;
   bool bf_sign_p_x = (1 - (translate_2pi_x(idx(2))*translate_2pi_x(idx(3))))/2;
   bool bf_sign_p_y = (1 - (translate_2pi_y(idx(2))*translate_2pi_y(idx(3))))/2;
   bool bf_sign_p_z = (1 - (translate_2pi_z(idx(2))*translate_2pi_z(idx(3))))/2;
   int bf[3] = {0,0,0};
   int bf_p[3] = {0,0,0};
   rot_2pi3_1m11( idx(0), bf_p );
   rot_2pi3_1m11( idx(1), bf );
   idx(2) = FF_rot_2pi3_1m11( idx(2) );
   idx(3) = FF_rot_2pi3_1m11( idx(3) );
   idx(4) = FF_rot_2pi3_1m11( idx(4) );
   idx(5) = FF_rot_2pi3_1m11( idx(5) );
   bf_sign_x &= bool(bf[0]);
   bf_sign_y &= bool(bf[1]);
   bf_sign_p_x &= bool(bf_p[0]);
   bf_sign_p_y &= bool(bf_p[1]);

   return operation( sign ^((bf_sign_x^bf_sign_y^bf_sign_z)^(bf_sign_p_x^bf_sign_p_y^bf_sign_p_z)), false );
}


operation rot_2pi3_m111_projmat( idx_proj_matrix_t& idx )
{
   bool sign= ((1- (FF_rot_sign_2pi3_m111( idx(4) ) * FF_rot_sign_2pi3_m111( idx(5) ) * FF_rot_sign_2pi3_m111( idx(2) ) * FF_rot_sign_2pi3_m111( idx(3) )))/2);
   bool bf_sign_x = (1 - (translate_2pi_x(idx(4))*translate_2pi_x(idx(5))))/2;
   bool bf_sign_y = (1 - (translate_2pi_y(idx(4))*translate_2pi_y(idx(5))))/2;
   bool bf_sign_z = (1 - (translate_2pi_z(idx(4))*translate_2pi_z(idx(5))))/2;
   bool bf_sign_p_x = (1 - (translate_2pi_x(idx(2))*translate_2pi_x(idx(3))))/2;
   bool bf_sign_p_y = (1 - (translate_2pi_y(idx(2))*translate_2pi_y(idx(3))))/2;
   bool bf_sign_p_z = (1 - (translate_2pi_z(idx(2))*translate_2pi_z(idx(3))))/2;
   int bf[3] = {0,0,0};
   int bf_p[3] = {0,0,0};
   rot_2pi3_m111( idx(0), bf_p );
   rot_2pi3_m111( idx(1), bf );
   idx(2) = FF_rot_2pi3_m111( idx(2) );
   idx(3) = FF_rot_2pi3_m111( idx(3) );
   idx(4) = FF_rot_2pi3_m111( idx(4) );
   idx(5) = FF_rot_2pi3_m111( idx(5) );
   bf_sign_x &= bool(bf[0]);
   bf_sign_y &= bool(bf[1]);
   bf_sign_p_x &= bool(bf_p[0]);
   bf_sign_p_y &= bool(bf_p[1]);

   return operation( sign ^((bf_sign_x^bf_sign_y^bf_sign_z)^(bf_sign_p_x^bf_sign_p_y^bf_sign_p_z)), false );
}


operation rot_2pi3_m1m11_projmat( idx_proj_matrix_t& idx )
{
   bool sign= ((1- (FF_rot_sign_2pi3_m1m11( idx(4) ) * FF_rot_sign_2pi3_m1m11( idx(5) ) * FF_rot_sign_2pi3_m1m11( idx(2) ) * FF_rot_sign_2pi3_m1m11( idx(3) )))/2);
   bool bf_sign_x = (1 - (translate_2pi_x(idx(4))*translate_2pi_x(idx(5))))/2;
   bool bf_sign_y = (1 - (translate_2pi_y(idx(4))*translate_2pi_y(idx(5))))/2;
   bool bf_sign_z = (1 - (translate_2pi_z(idx(4))*translate_2pi_z(idx(5))))/2;
   bool bf_sign_p_x = (1 - (translate_2pi_x(idx(2))*translate_2pi_x(idx(3))))/2;
   bool bf_sign_p_y = (1 - (translate_2pi_y(idx(2))*translate_2pi_y(idx(3))))/2;
   bool bf_sign_p_z = (1 - (translate_2pi_z(idx(2))*translate_2pi_z(idx(3))))/2;
   int bf[3] = {0,0,0};
   int bf_p[3] = {0,0,0};
   rot_2pi3_m1m11( idx(0), bf_p );
   rot_2pi3_m1m11( idx(1), bf );
   idx(2) = FF_rot_2pi3_m1m11( idx(2) );
   idx(3) = FF_rot_2pi3_m1m11( idx(3) );
   idx(4) = FF_rot_2pi3_m1m11( idx(4) );
   idx(5) = FF_rot_2pi3_m1m11( idx(5) );
   bf_sign_x &= bool(bf[0]);
   bf_sign_y &= bool(bf[1]);
   bf_sign_p_x &= bool(bf_p[0]);
   bf_sign_p_y &= bool(bf_p[1]);

   return operation( sign ^((bf_sign_x^bf_sign_y^bf_sign_z)^(bf_sign_p_x^bf_sign_p_y^bf_sign_p_z)), false );
}


operation rot_pi_110_projmat( idx_proj_matrix_t& idx )
{
   bool sign= ((1- (FF_rot_sign_pi_110( idx(4) ) * FF_rot_sign_pi_110( idx(5) ) * FF_rot_sign_pi_110( idx(2) ) * FF_rot_sign_pi_110( idx(3) )))/2);
   bool bf_sign_x = (1 - (translate_2pi_x(idx(4))*translate_2pi_x(idx(5))))/2;
   bool bf_sign_y = (1 - (translate_2pi_y(idx(4))*translate_2pi_y(idx(5))))/2;
   bool bf_sign_z = (1 - (translate_2pi_z(idx(4))*translate_2pi_z(idx(5))))/2;
   bool bf_sign_p_x = (1 - (translate_2pi_x(idx(2))*translate_2pi_x(idx(3))))/2;
   bool bf_sign_p_y = (1 - (translate_2pi_y(idx(2))*translate_2pi_y(idx(3))))/2;
   bool bf_sign_p_z = (1 - (translate_2pi_z(idx(2))*translate_2pi_z(idx(3))))/2;
   int bf[3] = {0,0,0};
   int bf_p[3] = {0,0,0};
   rot_pi_110( idx(0), bf_p );
   rot_pi_110( idx(1), bf );
   idx(2) = FF_rot_pi_110( idx(2) );
   idx(3) = FF_rot_pi_110( idx(3) );
   idx(4) = FF_rot_pi_110( idx(4) );
   idx(5) = FF_rot_pi_110( idx(5) );
   bf_sign_x &= bool(bf[0]);
   bf_sign_y &= bool(bf[1]);
   bf_sign_p_x &= bool(bf_p[0]);
   bf_sign_p_y &= bool(bf_p[1]);

   return operation( sign ^((bf_sign_x^bf_sign_y^bf_sign_z)^(bf_sign_p_x^bf_sign_p_y^bf_sign_p_z)), false );
}

operation rot_pi_1m10_projmat( idx_proj_matrix_t& idx )
{
   bool sign= ((1- (FF_rot_sign_pi_1m10( idx(4) ) * FF_rot_sign_pi_1m10( idx(5) ) * FF_rot_sign_pi_1m10( idx(2) ) * FF_rot_sign_pi_1m10( idx(3) )))/2);
   bool bf_sign_x = (1 - (translate_2pi_x(idx(4))*translate_2pi_x(idx(5))))/2;
   bool bf_sign_y = (1 - (translate_2pi_y(idx(4))*translate_2pi_y(idx(5))))/2;
   bool bf_sign_z = (1 - (translate_2pi_z(idx(4))*translate_2pi_z(idx(5))))/2;
   bool bf_sign_p_x = (1 - (translate_2pi_x(idx(2))*translate_2pi_x(idx(3))))/2;
   bool bf_sign_p_y = (1 - (translate_2pi_y(idx(2))*translate_2pi_y(idx(3))))/2;
   bool bf_sign_p_z = (1 - (translate_2pi_z(idx(2))*translate_2pi_z(idx(3))))/2;
   int bf[3] = {0,0,0};
   int bf_p[3] = {0,0,0};
   rot_pi_1m10( idx(0), bf_p );
   rot_pi_1m10( idx(1), bf );
   idx(2) = FF_rot_pi_1m10( idx(2) );
   idx(3) = FF_rot_pi_1m10( idx(3) );
   idx(4) = FF_rot_pi_1m10( idx(4) );
   idx(5) = FF_rot_pi_1m10( idx(5) );
   bf_sign_x &= bool(bf[0]);
   bf_sign_y &= bool(bf[1]);
   bf_sign_p_x &= bool(bf_p[0]);
   bf_sign_p_y &= bool(bf_p[1]);

   return operation( sign ^((bf_sign_x^bf_sign_y^bf_sign_z)^(bf_sign_p_x^bf_sign_p_y^bf_sign_p_z)), false );
}

operation rot_pi_011_projmat( idx_proj_matrix_t& idx )
{
   bool sign= ((1- (FF_rot_sign_pi_011( idx(4) ) * FF_rot_sign_pi_011( idx(5) ) * FF_rot_sign_pi_011( idx(2) ) * FF_rot_sign_pi_011( idx(3) )))/2);
   bool bf_sign_x = (1 - (translate_2pi_x(idx(4))*translate_2pi_x(idx(5))))/2;
   bool bf_sign_y = (1 - (translate_2pi_y(idx(4))*translate_2pi_y(idx(5))))/2;
   bool bf_sign_z = (1 - (translate_2pi_z(idx(4))*translate_2pi_z(idx(5))))/2;
   bool bf_sign_p_x = (1 - (translate_2pi_x(idx(2))*translate_2pi_x(idx(3))))/2;
   bool bf_sign_p_y = (1 - (translate_2pi_y(idx(2))*translate_2pi_y(idx(3))))/2;
   bool bf_sign_p_z = (1 - (translate_2pi_z(idx(2))*translate_2pi_z(idx(3))))/2;
   int bf[3] = {0,0,0};
   int bf_p[3] = {0,0,0};
   rot_pi_011( idx(0), bf_p );
   rot_pi_011( idx(1), bf );
   idx(2) = FF_rot_pi_011( idx(2) );
   idx(3) = FF_rot_pi_011( idx(3) );
   idx(4) = FF_rot_pi_011( idx(4) );
   idx(5) = FF_rot_pi_011( idx(5) );
   bf_sign_x &= bool(bf[0]);
   bf_sign_y &= bool(bf[1]);
   bf_sign_p_x &= bool(bf_p[0]);
   bf_sign_p_y &= bool(bf_p[1]);

   return operation( sign ^((bf_sign_x^bf_sign_y^bf_sign_z)^(bf_sign_p_x^bf_sign_p_y^bf_sign_p_z)), false );
}


operation rot_pi_01m1_projmat( idx_proj_matrix_t& idx )
{
   bool sign= ((1- (FF_rot_sign_pi_01m1( idx(4) ) * FF_rot_sign_pi_01m1( idx(5) ) * FF_rot_sign_pi_01m1( idx(2) ) * FF_rot_sign_pi_01m1( idx(3) )))/2);
   bool bf_sign_x = (1 - (translate_2pi_x(idx(4))*translate_2pi_x(idx(5))))/2;
   bool bf_sign_y = (1 - (translate_2pi_y(idx(4))*translate_2pi_y(idx(5))))/2;
   bool bf_sign_z = (1 - (translate_2pi_z(idx(4))*translate_2pi_z(idx(5))))/2;
   bool bf_sign_p_x = (1 - (translate_2pi_x(idx(2))*translate_2pi_x(idx(3))))/2;
   bool bf_sign_p_y = (1 - (translate_2pi_y(idx(2))*translate_2pi_y(idx(3))))/2;
   bool bf_sign_p_z = (1 - (translate_2pi_z(idx(2))*translate_2pi_z(idx(3))))/2;
   int bf[3] = {0,0,0};
   int bf_p[3] = {0,0,0};
   rot_pi_01m1( idx(0), bf_p );
   rot_pi_01m1( idx(1), bf );
   idx(2) = FF_rot_pi_01m1( idx(2) );
   idx(3) = FF_rot_pi_01m1( idx(3) );
   idx(4) = FF_rot_pi_01m1( idx(4) );
   idx(5) = FF_rot_pi_01m1( idx(5) );
   bf_sign_x &= bool(bf[0]);
   bf_sign_y &= bool(bf[1]);
   bf_sign_p_x &= bool(bf_p[0]);
   bf_sign_p_y &= bool(bf_p[1]);

   return operation( sign ^((bf_sign_x^bf_sign_y^bf_sign_z)^(bf_sign_p_x^bf_sign_p_y^bf_sign_p_z)), false );
}

operation rot_pi_101_projmat( idx_proj_matrix_t& idx )
{
   bool sign= ((1- (FF_rot_sign_pi_101( idx(4) ) * FF_rot_sign_pi_101( idx(5) ) * FF_rot_sign_pi_101( idx(2) ) * FF_rot_sign_pi_101( idx(3) )))/2);
   bool bf_sign_x = (1 - (translate_2pi_x(idx(4))*translate_2pi_x(idx(5))))/2;
   bool bf_sign_y = (1 - (translate_2pi_y(idx(4))*translate_2pi_y(idx(5))))/2;
   bool bf_sign_z = (1 - (translate_2pi_z(idx(4))*translate_2pi_z(idx(5))))/2;
   bool bf_sign_p_x = (1 - (translate_2pi_x(idx(2))*translate_2pi_x(idx(3))))/2;
   bool bf_sign_p_y = (1 - (translate_2pi_y(idx(2))*translate_2pi_y(idx(3))))/2;
   bool bf_sign_p_z = (1 - (translate_2pi_z(idx(2))*translate_2pi_z(idx(3))))/2;
   int bf[3] = {0,0,0};
   int bf_p[3] = {0,0,0};
   rot_pi_101( idx(0), bf_p );
   rot_pi_101( idx(1), bf );
   idx(2) = FF_rot_pi_101( idx(2) );
   idx(3) = FF_rot_pi_101( idx(3) );
   idx(4) = FF_rot_pi_101( idx(4) );
   idx(5) = FF_rot_pi_101( idx(5) );
   bf_sign_x &= bool(bf[0]);
   bf_sign_y &= bool(bf[1]);
   bf_sign_p_x &= bool(bf_p[0]);
   bf_sign_p_y &= bool(bf_p[1]);

   return operation( sign ^((bf_sign_x^bf_sign_y^bf_sign_z)^(bf_sign_p_x^bf_sign_p_y^bf_sign_p_z)), false );
}

operation rot_pi_10m1_projmat( idx_proj_matrix_t& idx )
{
   bool sign= ((1- (FF_rot_sign_pi_10m1( idx(4) ) * FF_rot_sign_pi_10m1( idx(5) ) * FF_rot_sign_pi_10m1( idx(2) ) * FF_rot_sign_pi_10m1( idx(3) )))/2);
   bool bf_sign_x = (1 - (translate_2pi_x(idx(4))*translate_2pi_x(idx(5))))/2;
   bool bf_sign_y = (1 - (translate_2pi_y(idx(4))*translate_2pi_y(idx(5))))/2;
   bool bf_sign_z = (1 - (translate_2pi_z(idx(4))*translate_2pi_z(idx(5))))/2;
   bool bf_sign_p_x = (1 - (translate_2pi_x(idx(2))*translate_2pi_x(idx(3))))/2;
   bool bf_sign_p_y = (1 - (translate_2pi_y(idx(2))*translate_2pi_y(idx(3))))/2;
   bool bf_sign_p_z = (1 - (translate_2pi_z(idx(2))*translate_2pi_z(idx(3))))/2;
   int bf[3] = {0,0,0};
   int bf_p[3] = {0,0,0};
   rot_pi_10m1( idx(0), bf_p );
   rot_pi_10m1( idx(1), bf );
   idx(2) = FF_rot_pi_10m1( idx(2) );
   idx(3) = FF_rot_pi_10m1( idx(3) );
   idx(4) = FF_rot_pi_10m1( idx(4) );
   idx(5) = FF_rot_pi_10m1( idx(5) );
   bf_sign_x &= bool(bf[0]);
   bf_sign_y &= bool(bf[1]);
   bf_sign_p_x &= bool(bf_p[0]);
   bf_sign_p_y &= bool(bf_p[1]);

   return operation( sign ^((bf_sign_x^bf_sign_y^bf_sign_z)^(bf_sign_p_x^bf_sign_p_y^bf_sign_p_z)), false );
}

operation mirror_y_projmat( idx_proj_matrix_t& idx )
{
   bool sign= ((1- (FF_mirror_y_sign( idx(4) ) * FF_mirror_y_sign( idx(5) ) * FF_mirror_y_sign( idx(2) ) * FF_mirror_y_sign( idx(3) )))/2);
   bool bf_sign_x = (1 - (translate_2pi_x(idx(4))*translate_2pi_x(idx(5))))/2;
   bool bf_sign_y = (1 - (translate_2pi_y(idx(4))*translate_2pi_y(idx(5))))/2;
   bool bf_sign_z = (1 - (translate_2pi_z(idx(4))*translate_2pi_z(idx(5))))/2;
   bool bf_sign_p_x = (1 - (translate_2pi_x(idx(2))*translate_2pi_x(idx(3))))/2;
   bool bf_sign_p_y = (1 - (translate_2pi_y(idx(2))*translate_2pi_y(idx(3))))/2;
   bool bf_sign_p_z = (1 - (translate_2pi_z(idx(2))*translate_2pi_z(idx(3))))/2;
   int bf[3] = {0,0,0};
   int bf_p[3] = {0,0,0};
   mirror_y( idx(0), bf_p );
   mirror_y( idx(1), bf );
   idx(2) = FF_mirror_y( idx(2) );
   idx(3) = FF_mirror_y( idx(3) );
   idx(4) = FF_mirror_y( idx(4) );
   idx(5) = FF_mirror_y( idx(5) );
   bf_sign_x &= bool(bf[0]);
   bf_sign_y &= bool(bf[1]);
   bf_sign_z &= bool(bf[2]);
   bf_sign_p_x &= bool(bf_p[0]);
   bf_sign_p_y &= bool(bf_p[1]);
   bf_sign_p_z &= bool(bf_p[2]);

   return operation( sign ^((bf_sign_x^bf_sign_y^bf_sign_z)^(bf_sign_p_x^bf_sign_p_y^bf_sign_p_z)), false );
}
