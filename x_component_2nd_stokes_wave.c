/*second order stokes wave at inlet Boundary, wave velocity components are from equations 3.27 and 3.58, Pengzhi lin. numerical modeling of water waves. CRC press, 2008*/
#include "udf.h"
#define pi 3.14159265359 /*define  constants*/
#define U  0.6 /*free stream velocity*/
#define H 0.076 /*wave height*/
#define g 9.81 /*gravity acceleration*/
#define L 4.8 /*wave length*/
#define d 1.6 /*water depth*/
#define T 1.456/*effective wave period (include doppler effect)*/ 

/*DEFINE_PROFILE: define an inlet velocity profile  that varies as a function of z coordinates or t.*/
DEFINE_PROFILE(streamwise_velocity,ft,var) /* DEFINE Macros, ft : thread; var:index */
{ 
	/*define variables*/
	real r[ND_ND]; /* array, Coordinates are stored in array, r[0] mean x coordinates, r[1] means y coordinates*/
	real k; /*wave number*/
	real z; /*z(vertical) axis*/
	real omega; /*effective wave angular velocity*/
	real t; /*t*/
	k = 2.0*pi/L; 	/*assign values to variables*/
	omega=2.0*pi/T;   /* effective wave angular velocity*/
	t = CURRENT_TIME; /* Special Fluent macro, current running t */

	face_t f; /* "f" is a face index for each face on the boundary */ 
	
	begin_f_loop(f,ft)/* face loop macro ,loop over all faces in a given face thread,i.e. "ft" */ 
	{		
		 F_CENTROID(r,f,ft); /*F_CENTROID finds the coordinate position of the centroid of the face "f" and stores the coordinates in the "r" array */
		z =r[2]; /* r[1] is y coordinate,r[2] is z coordinate */
		
F_PROFILE(f,ft,var) = U + H*g*k*cosh(k*(-z-0.782+d))*cos(-omega*t)/(2.0*(omega-k*U)*cosh(k*d)) + 3.0*H*H*(omega-k*U)*k*cosh(2.0*k*(-z-0.782+d))*cos(-2.0*omega*t)/(16.0*pow(sinh(k*d),4.0));
/*x-velocity component (flow direction): u_r+U= U + H*g*k*cosh(k*(z+d))*cos(-omega*t)/(2.0*(omega-k*U)*cosh(k*d)) + 3.0*H*H*(omega-k*U)*k*cosh(2.0*k*(z+d))*cos(-2.0*omega*t)/(16.0*pow(sinh(k*d),4.0)); z=0 is the mean free surface levelin the theory model, however, free surface level is z=-0.782m in the Fluent geometry model*/

	}
	end_f_loop(f,ft)
}



