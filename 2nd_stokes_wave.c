
/*********************/
/*linear wave*/
/*  u_r = H*g*k/(2*(omega-k*U))*exp(k*z)*(1+exp(-2*k*(d+z)))/(1+exp(-2*k*d))*cos(omega*tim)     */
/*  w_r = H*g*k/(2*(omega-k*U))*exp(k*z)*(1-exp(-2*k*(d+z)))/(1+exp(-2*k*d))*cos(omega*time)    */
/**************************/


include <udf.h>             /* include head file */
# define pi 3.14159265359  
/*define variable*/
real r[ND_ND]; /* Coordinates, r[0] mean x coordinates, r[1] means y coordinates*/
real U; /*free stream velocity*/
real H; /*wave height*/
real g; /*gravity acceleration*/
real k; /*wave number*/
real L; /*wave length*/
real z; /*z(vertical) axis*/
real d; /*water depth*/
real theta_t; /*wave phase*/
real T;/*wave period*/ 
real omega; /*wave angular velocity*/
real t; /*time*/
real theta_t; /*theta_t=-omega*t*/
/*assign values to variables*/
U=0.6;
H=0.076;
g=9.81;
L=4.8;
k=2.0*pi/L;
d=1.6;
T=1.43;
omega=2.0*pi/T;
t = CURRENT_TIME; /* Special Fluent macro, current running time */ 
theta_t=-omega*t;
z=r[2]; /*r[2] is z coordinates*/

/*DEFINE_PROFILE: define a custom boundary profile  that varies as a function of z coordinates or time.*/
DEFINE_PROFILE(x_velocity, ft, var) /* DEFINE Macros, ft is a thread; var:index */
{ 
   face_t f; /* Face index */ 

  begin_f_loop(f,ft) /* face loop macro ,loop over all faces in a given face thread “ft”*/ 
    {
      F_CENTROID(r,f,ft); /*only in the pressure-based solver, F_CENTROID finds the coordinate position of the centroid of the face "f" and stores the coordinates in the "r" array */
      F_PROFILE(f,ft,var) = U +  H*g*k/(2*(omega-k*U))*exp(k*z)*(1+exp(-2*k*(d+z)))/(1+exp(-2*k*d))*cos(omega*tim)     
}
end_f_loop(f,ft)
}




