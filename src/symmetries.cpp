
/************************************************************************************************//**
 *  		
 * 	file: 		symmetries.cpp
 * 	contents:   	see symmetries.h
 * 
 ****************************************************************************************************/

#include <symmetries.h>
#include <grid.h>

using namespace std;

// SYMMETRIES FOR THE FORM FACTORS IMPLEMENTATION -> FOR THE MOMENT THE FFFACTORS ARE ASSUMED TO BE SYMMETRIC UNDER INVERSION AND REAL
// Shifted notation adopted for frequency and momenta

// -- Phi

operation hmirror_phi_pp( idx_phi_t& idx )
{
   bool sign= (1 - (Parity(idx(4))*Parity(idx(5))))/2;
   bool bf_sign_x = (1 - (translate_2pi_x(idx(4))*translate_2pi_x(idx(5))))/2;
   bool bf_sign_y = (1 - (translate_2pi_y(idx(4))*translate_2pi_y(idx(5))))/2;
   int bf[2] = {0,0};
   idx(0) *= 1; 
   idx(1) =-idx(1)-1-abs(idx(0)%2);
   idx(2) =-idx(2)-1-abs(idx(0)%2);
   idx(3) *= 1;
   idx(4) *= 1;
   idx(5) *= 1;
   swap(idx(6),idx(7));
   swap(idx(8),idx(9));
   bf_sign_x &= bool(bf[0]);
   bf_sign_y &= bool(bf[1]);
   return operation( sign^(bf_sign_x^bf_sign_y) , false ); 	// true=1 is multiplying a factor -1 to the matrix element
   										// if Parity*Parity==-1 -> first = true -> factor -1
   										// elseif Parity*Parity==1 -> first = false -> no factor
}

operation hmirror_phi_ph( idx_phi_t& idx )
{
   bool sign = false;
   bool bf_sign_x = (1 - (translate_2pi_x(idx(4))*translate_2pi_x(idx(5))))/2;
   bool bf_sign_y = (1 - (translate_2pi_y(idx(4))*translate_2pi_y(idx(5))))/2;
   int bf[2] = {0,0};
   idx(0) *= -1;
   swap(idx(1),idx(2)); 
   idx(3) =neg_k(idx(3),bf);
   swap(idx(4),idx(5));
   swap(idx(6),idx(7));
   swap(idx(8),idx(9));
   bf_sign_x &= bool(bf[0]);
   bf_sign_y &= bool(bf[1]);
   return operation( sign^(bf_sign_x^bf_sign_y) , false ); 	// true=1 is multiplying a factor -1 to the matrix element
   										// if Parity*Parity==-1 -> first = true -> factor -1
   										// elseif Parity*Parity==1 -> first = false -> no factor
}

operation hmirror_phi_xph( idx_phi_t& idx )
{
   bool sign = false;
   bool bf_sign_x = (1 - (translate_2pi_x(idx(4))*translate_2pi_x(idx(5))))/2;
   bool bf_sign_y = (1 - (translate_2pi_y(idx(4))*translate_2pi_y(idx(5))))/2;
   int bf[2] = {0,0};
   idx(0) *= -1;
   swap(idx(1),idx(2)); 
   idx(3) = neg_k(idx(3),bf);
   swap(idx(4),idx(5));
   swap(idx(6),idx(7));
   swap(idx(8),idx(9));
   bf_sign_x &= bool(bf[0]);
   bf_sign_y &= bool(bf[1]);
   return operation( sign^(bf_sign_x^bf_sign_y), false ); 
}

operation cconj_phi_pp( idx_phi_t& idx ) //notation changed but fine = cconj*hmirror in 3D_patch
{
   bool sign = false;
   bool bf_sign_x = (1 - (translate_2pi_x(idx(4))*translate_2pi_x(idx(5))))/2;
   bool bf_sign_y = (1 - (translate_2pi_y(idx(4))*translate_2pi_y(idx(5))))/2;
   int bf[2] = {0,0};
   idx(0) *= -1; 
   swap(idx(1),idx(2));
   idx(1) =-idx(1)-1-abs(idx(0)%2);
   idx(2) =-idx(2)-1-abs(idx(0)%2);
   idx(3) *=1;
   swap(idx(4),idx(5));
   swap(idx(6),idx(8));
   swap(idx(7),idx(9));
   bf_sign_x &= bool(bf[0]);
   bf_sign_y &= bool(bf[1]);
   return operation( sign^(bf_sign_x^bf_sign_y), true ); 
}

operation cconj_phi_ph( idx_phi_t& idx )
{
   bool sign = false;
   bool bf_sign_x = (1 - (translate_2pi_x(idx(4))*translate_2pi_x(idx(5))))/2;
   bool bf_sign_y = (1 - (translate_2pi_y(idx(4))*translate_2pi_y(idx(5))))/2;
   int bf[2] = {0,0};
   idx(0) *=  1;
   idx(1) =-idx(1)-1-abs(idx(0)%2);
   idx(2) =-idx(2)-1-abs(idx(0)%2);
   idx(3) = neg_k(idx(3),bf);
   swap(idx(6),idx(8));
   swap(idx(7),idx(9));
   bf_sign_x &= bool(bf[0]);
   bf_sign_y &= bool(bf[1]);
   return operation( sign^(bf_sign_x^bf_sign_y), true ); 
}

operation cconj_phi_xph( idx_phi_t& idx )
{
   bool sign = false;
   bool bf_sign_x = (1 - (translate_2pi_x(idx(4))*translate_2pi_x(idx(5))))/2;
   bool bf_sign_y = (1 - (translate_2pi_y(idx(4))*translate_2pi_y(idx(5))))/2;
   int bf[2] = {0,0};
   idx(0) *= -1;
   swap(idx(1),idx(2));
   idx(1) =-idx(1)-1-abs(idx(0)%2);
   idx(2) =-idx(2)-1-abs(idx(0)%2);
   swap(idx(4),idx(5));  
   swap(idx(6),idx(8));
   swap(idx(7),idx(9));
   bf_sign_x &= bool(bf[0]);
   bf_sign_y &= bool(bf[1]);
   return operation( sign^(bf_sign_x^bf_sign_y), true ); 
}

operation timerev_phi_pp( idx_phi_t& idx )
{
   bool sign = false;
   bool bf_sign_x = (1 - (translate_2pi_x(idx(4))*translate_2pi_x(idx(5))))/2;
   bool bf_sign_y = (1 - (translate_2pi_y(idx(4))*translate_2pi_y(idx(5))))/2;
   int bf[2] = {0,0};
   idx(0) *= 1; 
   swap(idx(1),idx(2));
   idx(3) *=1;
   swap(idx(4),idx(5));
   swap(idx(6),idx(8));
   swap(idx(7),idx(9));
   bf_sign_x &= bool(bf[0]);
   bf_sign_y &= bool(bf[1]);
   return operation( sign^(bf_sign_x^bf_sign_y), false ); 
}

operation timerev_phi_ph( idx_phi_t& idx )
{
   bool sign = false;
   bool bf_sign_x = (1 - (translate_2pi_x(idx(4))*translate_2pi_x(idx(5))))/2;
   bool bf_sign_y = (1 - (translate_2pi_y(idx(4))*translate_2pi_y(idx(5))))/2;
   int bf[2] = {0,0};
   idx(0) *= -1;
   idx(3) = neg_k(idx(3),bf);
   swap(idx(6),idx(8));
   swap(idx(7),idx(9));
   bf_sign_x &= bool(bf[0]);
   bf_sign_y &= bool(bf[1]);
   return operation( sign^(bf_sign_x^bf_sign_y), false ); 
}

operation timerev_phi_xph( idx_phi_t& idx )
{
   bool sign = false;
   bool bf_sign_x = (1 - (translate_2pi_x(idx(4))*translate_2pi_x(idx(5))))/2;
   bool bf_sign_y = (1 - (translate_2pi_y(idx(4))*translate_2pi_y(idx(5))))/2;
   int bf[2] = {0,0};
   idx(0) *= 1;
   swap(idx(1),idx(2));
   swap(idx(4),idx(5));
   swap(idx(6),idx(8));
   swap(idx(7),idx(9));
   bf_sign_x &= bool(bf[0]);
   bf_sign_y &= bool(bf[1]);
   return operation( sign^(bf_sign_x^bf_sign_y), false ); 
}

// --- P
operation hmirror_P_pp( idx_P_t& idx )
{
   bool sign= (1 - Parity(idx(3)))/2;
   bool bf_sign_x = (1 - translate_2pi_x(idx(3)))/2;
   bool bf_sign_y = (1 - translate_2pi_y(idx(3)))/2;
   int bf[2] = {0,0};
   idx(0) *= 1; 
   idx(1) =-idx(1)-1-abs(idx(0)%2);
   idx(2) *=1;
   idx(3) *= 1;
   swap(idx(4),idx(5));
   swap(idx(6),idx(7));
   bf_sign_x &= bool(bf[0]);
   bf_sign_y &= bool(bf[1]);
   return operation( sign^(bf_sign_x^bf_sign_y), false ); 
}

operation cconj_timerev_P_pp( idx_P_t& idx )
{
   bool sign = false;
   bool bf_sign_x = (1 - translate_2pi_x(idx(3)))/2;
   bool bf_sign_y = (1 - translate_2pi_y(idx(3)))/2;
   int bf[2] = {0,0};
   idx(0) *= -1; 
   idx(1) =-idx(1)-1-abs(idx(0)%2);
   idx(2) *= 1;
   idx(3) *= 1;
   bf_sign_x &= bool(bf[0]);
   bf_sign_y &= bool(bf[1]);
   return operation( sign^(bf_sign_x^bf_sign_y), true ); // Quantum numbers don't change after cconj and time reversal
}

operation timerev_P_ph( idx_P_t& idx )
{
   bool sign = false;
   bool bf_sign_x = (1 - translate_2pi_x(idx(3)))/2;
   bool bf_sign_y = (1 - translate_2pi_y(idx(3)))/2;
   int bf[2] = {0,0};
   idx(0) *= -1;
   idx(1) *= 1;
   idx(2) =neg_k(idx(2),bf);
   idx(3) *= 1;
   swap(idx(4),idx(6));
   swap(idx(5),idx(7));
   bf_sign_x &= bool(bf[0]);
   bf_sign_y &= bool(bf[1]);
   return operation( sign^(bf_sign_x^bf_sign_y), false ); 
}

operation cconj_P_ph( idx_P_t& idx )
{
   bool sign = false;
   bool bf_sign_x = (1 - translate_2pi_x(idx(3)))/2;
   bool bf_sign_y = (1 - translate_2pi_y(idx(3)))/2;
   int bf[2] = {0,0};
   idx(0) *=  1;
   idx(1) =-idx(1)-1-abs(idx(0)%2);
   idx(2) = neg_k(idx(2),bf);
   idx(3) *= 1;
   swap(idx(4),idx(6));
   swap(idx(5),idx(7));
   bf_sign_x &= bool(bf[0]);
   bf_sign_y &= bool(bf[1]);
   return operation( sign^(bf_sign_x^bf_sign_y), true ); 
}

operation hmirror_cconj_P_xph( idx_P_t& idx )
{
   bool sign = false;
   bool bf_sign_x = (1 - translate_2pi_x(idx(3)))/2;
   bool bf_sign_y = (1 - translate_2pi_y(idx(3)))/2;
   int bf[2] = {0,0};
   idx(0) *= 1;
   idx(1) =-idx(1)-1-abs(idx(0)%2);
   idx(2) = neg_k(idx(2),bf);
   idx(3) *= 1;
   swap(idx(4),idx(7));
   swap(idx(5),idx(6));
   bf_sign_x &= bool(bf[0]);
   bf_sign_y &= bool(bf[1]);
   return operation( sign^(bf_sign_x^bf_sign_y), true ); 
}

operation hmirror_timerev_P_xph( idx_P_t& idx )
{
   bool sign = false;
   bool bf_sign_x = (1 - translate_2pi_x(idx(3)))/2;
   bool bf_sign_y = (1 - translate_2pi_y(idx(3)))/2;
   int bf[2] = {0,0};
   idx(0) *= -1;
   idx(1) *= 1;
   idx(2) = neg_k(idx(2),bf);
   idx(3) *= 1;
   swap(idx(4),idx(7));
   swap(idx(5),idx(6));
   bf_sign_x &= bool(bf[0]);
   bf_sign_y &= bool(bf[1]);
   return operation( sign^(bf_sign_x^bf_sign_y), false ); 
}

//// --- Chi
//
operation hmirror_chi_pp( idx_chi_t& idx )
{
   bool sign= false;
   bool bf_sign_x = false;
   bool bf_sign_y = false;
   int bf[2] = {0,0};
   idx(0) *= 1; 
   idx(1) *= 1;
   swap(idx(2),idx(3));
   swap(idx(4),idx(5));
   return operation( sign^(bf_sign_x^bf_sign_y), false ); 
}
operation hmirror_chi_ph( idx_chi_t& idx )
{
   bool sign= false;
   bool bf_sign_x = false;
   bool bf_sign_y = false;
   int bf[2] = {0,0};
   idx(0) *= -1;
   idx(1) =neg_k(idx(1), bf);
   swap(idx(2),idx(3));
   swap(idx(4),idx(5));
   return operation( sign^(bf_sign_x^bf_sign_y), false ); 
}

operation hmirror_chi_xph( idx_chi_t& idx )
{
   bool sign= false;
   bool bf_sign_x = false;
   bool bf_sign_y = false;
   int bf[2] = {0,0};
   idx(0) *= -1;
   idx(1) = neg_k(idx(1), bf);
   swap(idx(2),idx(3));
   swap(idx(4),idx(5));
   return operation( sign^(bf_sign_x^bf_sign_y), false ); 
}

operation cconj_chi_pp( idx_chi_t& idx )
{
   bool sign= false;
   bool bf_sign_x = false;
   bool bf_sign_y = false;
   int bf[2] = {0,0};
   idx(0) *= -1; 
   idx(1) *=  1;
   swap(idx(2),idx(4));
   swap(idx(3),idx(5));
   return operation( sign^(bf_sign_x^bf_sign_y), true ); 
}

operation cconj_chi_ph( idx_chi_t& idx )
{
   bool sign= false;
   bool bf_sign_x = false;
   bool bf_sign_y = false;
   int bf[2] = {0,0};
   idx(0) *=  1;
   idx(1)  =  neg_k(idx(1), bf);
   swap(idx(2),idx(4));
   swap(idx(3),idx(5));
   return operation( sign^(bf_sign_x^bf_sign_y), true ); 
}

operation cconj_chi_xph( idx_chi_t& idx )
{
   bool sign= false;
   bool bf_sign_x = false;
   bool bf_sign_y = false;
   int bf[2] = {0,0};
   idx(0) *= -1;
   swap(idx(2),idx(4));
   swap(idx(3),idx(5));
   return operation( sign^(bf_sign_x^bf_sign_y), true ); 
}

operation timerev_chi_pp( idx_chi_t& idx )
{
   bool sign= false;
   bool bf_sign_x = false;
   bool bf_sign_y = false;
   int bf[2] = {0,0};
   idx(0) *= 1; 
   idx(1) *= 1;
   swap(idx(2),idx(4));
   swap(idx(3),idx(5));
   return operation( sign^(bf_sign_x^bf_sign_y), false ); 
}

operation timerev_chi_ph( idx_chi_t& idx )
{
   bool sign= false;
   bool bf_sign_x = false;
   bool bf_sign_y = false;
   int bf[2] = {0,0};
   idx(0) *= -1;
   idx(1) = neg_k(idx(1), bf);
   swap(idx(2),idx(4));
   swap(idx(3),idx(5));
   return operation( sign^(bf_sign_x^bf_sign_y), false ); 
}
operation timerev_chi_xph( idx_chi_t& idx )
{
   bool sign= false;
   bool bf_sign_x = false;
   bool bf_sign_y = false;
   int bf[2] = {0,0};
   swap(idx(2),idx(4));
   swap(idx(3),idx(5));
   return operation( sign^(bf_sign_x^bf_sign_y), false ); 
}

// --- Sig
operation cconj_sig( idx_1p_t& idx )
{
   idx(0) = -idx(0) -1;
   idx(1) *= 1;  
   swap(idx(2),idx(3));
   return operation( false, true ); 
}


/***************************************************************************
 *
 * 		SYMMETRIES FOR THE BUBBLES
 *
 **************************************************************************/


operation swap_bubble( idx_bubble_mat_t& idx )
{
   swap(idx(2),idx(3));
   return operation( false, false ); 
}

operation cconj_bubble( idx_bubble_mat_t& idx )
{
   idx(0) *= -1; 
   idx(1) = -idx(1)-1-abs(idx(0)%2);
   swap(idx(4),idx(5));
   swap(idx(6),idx(7));
   return operation( false, true ); 
}

operation hmirror_bubble_pp( idx_bubble_mat_t& idx )
{
   idx(0) *=  1;
   idx(1)  =  -idx(1)-1-abs(idx(0)%2);
   swap(idx(4),idx(6));
   swap(idx(5),idx(7));
   return operation((1 - (Parity(idx(2))*Parity(idx(3))))/2 , false ); 
}

operation hmirror_timerev_ph( idx_bubble_mat_t& idx)
{
   idx(0) *= -1;
   idx(1) *=  1;
   swap(idx(4),idx(7));
   swap(idx(5),idx(6));
   return operation((1 - (Parity(idx(2))*Parity(idx(3))))/2 , false ); 
}

/**************************************************************
 *
 * 		SYMMETRIES SUSCEPTIBILITIES
 *
 *************************************************************/

operation hmirror_suscept_pp( idx_suscept_t& idx )
{
   bool sign= (1 - (Parity(idx(2))*Parity(idx(3))))/2;
   bool bf_sign_x = (1 - (translate_2pi_x(idx(2))*translate_2pi_x(idx(3))))/2;
   bool bf_sign_y = (1 - (translate_2pi_y(idx(2))*translate_2pi_y(idx(3))))/2;
   int bf[2] = {0,0};
   idx(0) *= 1; 
   idx(1) *= 1;
   idx(2) *= 1;
   idx(3) *= 1;
   swap(idx(4),idx(5));
   swap(idx(6),idx(7));
   bf_sign_x &= bool(bf[0]);
   bf_sign_y &= bool(bf[1]);
   return operation( sign^(bf_sign_x^bf_sign_y) , false ); 	// true=1 is multiplying a factor -1 to the matrix element
   										// if Parity*Parity==-1 -> first = true -> factor -1
   										// elseif Parity*Parity==1 -> first = false -> no factor
}

operation hmirror_suscept_ph( idx_suscept_t& idx )
{
   bool sign = false;
   bool bf_sign_x = (1 - (translate_2pi_x(idx(2))*translate_2pi_x(idx(3))))/2;
   bool bf_sign_y = (1 - (translate_2pi_y(idx(2))*translate_2pi_y(idx(3))))/2;
   int bf[2] = {0,0};
   idx(0) *= -1;
   idx(1) =neg_k(idx(1),bf);
   swap(idx(2),idx(3));
   swap(idx(4),idx(5));
   swap(idx(6),idx(7));
   bf_sign_x &= bool(bf[0]);
   bf_sign_y &= bool(bf[1]);
   return operation( sign^(bf_sign_x^bf_sign_y) , false ); 	// true=1 is multiplying a factor -1 to the matrix element
   										// if Parity*Parity==-1 -> first = true -> factor -1
   										// elseif Parity*Parity==1 -> first = false -> no factor
}

operation hmirror_suscept_xph( idx_suscept_t& idx )
{
   bool sign = false;
   bool bf_sign_x = (1 - (translate_2pi_x(idx(2))*translate_2pi_x(idx(3))))/2;
   bool bf_sign_y = (1 - (translate_2pi_y(idx(2))*translate_2pi_y(idx(3))))/2;
   int bf[2] = {0,0};
   idx(0) *= -1;
   idx(1) = neg_k(idx(1),bf);
   swap(idx(2),idx(3));
   swap(idx(4),idx(5));
   swap(idx(6),idx(7));
   bf_sign_x &= bool(bf[0]);
   bf_sign_y &= bool(bf[1]);
   return operation( sign^(bf_sign_x^bf_sign_y), false ); 
}

operation cconj_suscept_pp( idx_suscept_t& idx ) //notation changed but fine = cconj*hmirror in 3D_patch
{
   bool sign = false;
   bool bf_sign_x = (1 - (translate_2pi_x(idx(2))*translate_2pi_x(idx(3))))/2;
   bool bf_sign_y = (1 - (translate_2pi_y(idx(2))*translate_2pi_y(idx(3))))/2;
   int bf[2] = {0,0};
   idx(0) *= -1; 
   idx(1) *=1;
   swap(idx(2),idx(3));
   swap(idx(4),idx(6));
   swap(idx(5),idx(7));
   bf_sign_x &= bool(bf[0]);
   bf_sign_y &= bool(bf[1]);
   return operation( sign^(bf_sign_x^bf_sign_y), true ); 
}

operation cconj_suscept_ph( idx_suscept_t& idx )
{
   bool sign = false;
   bool bf_sign_x = (1 - (translate_2pi_x(idx(2))*translate_2pi_x(idx(3))))/2;
   bool bf_sign_y = (1 - (translate_2pi_y(idx(2))*translate_2pi_y(idx(3))))/2;
   int bf[2] = {0,0};
   idx(0) *=  1;
   idx(1) = neg_k(idx(1),bf);
   swap(idx(4),idx(6));
   swap(idx(5),idx(7));
   bf_sign_x &= bool(bf[0]);
   bf_sign_y &= bool(bf[1]);
   return operation( sign^(bf_sign_x^bf_sign_y), true ); 
}

operation cconj_suscept_xph( idx_suscept_t& idx )
{
   bool sign = false;
   bool bf_sign_x = (1 - (translate_2pi_x(idx(2))*translate_2pi_x(idx(3))))/2;
   bool bf_sign_y = (1 - (translate_2pi_y(idx(2))*translate_2pi_y(idx(3))))/2;
   int bf[2] = {0,0};
   idx(0) *= -1;
   swap(idx(2),idx(3));  
   swap(idx(4),idx(6));
   swap(idx(5),idx(7));
   bf_sign_x &= bool(bf[0]);
   bf_sign_y &= bool(bf[1]);
   return operation( sign^(bf_sign_x^bf_sign_y), true ); 
}

operation timerev_suscept_pp( idx_suscept_t& idx )
{
   bool sign = false;
   bool bf_sign_x = (1 - (translate_2pi_x(idx(2))*translate_2pi_x(idx(3))))/2;
   bool bf_sign_y = (1 - (translate_2pi_y(idx(2))*translate_2pi_y(idx(3))))/2;
   int bf[2] = {0,0};
   idx(0) *= 1; 
   idx(1) *=1;
   swap(idx(2),idx(3));
   swap(idx(4),idx(6));
   swap(idx(5),idx(7));
   bf_sign_x &= bool(bf[0]);
   bf_sign_y &= bool(bf[1]);
   return operation( sign^(bf_sign_x^bf_sign_y), false ); 
}

operation timerev_suscept_ph( idx_suscept_t& idx )
{
   bool sign = false;
   bool bf_sign_x = (1 - (translate_2pi_x(idx(2))*translate_2pi_x(idx(3))))/2;
   bool bf_sign_y = (1 - (translate_2pi_y(idx(2))*translate_2pi_y(idx(3))))/2;
   int bf[2] = {0,0};
   idx(0) *= -1;
   idx(1) = neg_k(idx(1),bf);
   swap(idx(4),idx(6));
   swap(idx(5),idx(7));
   bf_sign_x &= bool(bf[0]);
   bf_sign_y &= bool(bf[1]);
   return operation( sign^(bf_sign_x^bf_sign_y), false ); 
}


operation timerev_suscept_xph( idx_suscept_t& idx )
{
   bool sign = false;
   bool bf_sign_x = (1 - (translate_2pi_x(idx(2))*translate_2pi_x(idx(3))))/2;
   bool bf_sign_y = (1 - (translate_2pi_y(idx(2))*translate_2pi_y(idx(3))))/2;
   int bf[2] = {0,0};
   idx(0) *= 1;
   swap(idx(2),idx(3));
   swap(idx(4),idx(6));
   swap(idx(5),idx(7));
   bf_sign_x &= bool(bf[0]);
   bf_sign_y &= bool(bf[1]);
   return operation( sign^(bf_sign_x^bf_sign_y), false ); 
}


