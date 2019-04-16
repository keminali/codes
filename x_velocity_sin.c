/* U= U_infity ( 1+ 0.05 sin(omega*t)*/

#include "udf.h"
#define pi 3.14159265359 /*define  constants*/
#define U  0.6 /*free stream velocity*/
#define T 0.456/*effective wave period (include doppler effect)*/ 
#define k 0.1 
#define m 1
/*DEFINE_PROFILE: define an inlet velocity profile  that varies as a function of z coordinates or t.*/
DEFINE_PROFILE(x_vel,ft,var) /* DEFINE Macros, ft is a thread; var:index */
{ 
	/*define variables*/
	real r[ND_ND]; /*Coordinates, r[0] mean x coordinates, r[1] means y coordinates*/
	real z; /*z(vertical) axis*/
	real omega; /*effective wave angular velocity*/
	real t; /*t*/
	omega=2.0*pi/T;
	t = CURRENT_TIME; /* Special Fluent macro, current running t */

	face_t f; /* "f" is a face index for each face on the boundary */ 
	
	begin_f_loop(f,ft)/* face loop macro ,loop over all faces in a given face thread,i.e. "ft" */ 
	{		
		 F_CENTROID(r,f,ft); /*F_CENTROID finds the coordinate position of the centroid of the face "f" and stores the coordinates in the "r" array */
		z =r[2]; /* r[1] is y coordinate,r[2] is z coordinate */
		
F_PROFILE(f,ft,var) = U*( m + k*sin(omega*t));

	}
	end_f_loop(f,ft)

** save residual 
#+BEGIN_SRC emacs-lisp
;save residuals
(display "Save the residual in a file") (newline)
  (let
   ((
     writefile (lambda (p)
     (define np (length (residual-history "iteration")))
     (let loop ((i 0))
    (if (not (= i np))
    (begin (define j (+ i 1))
    (display (list-ref (residual-history "iteration") (- np j)) p) (display " " p)
   (display (list-ref (residual-history "continuity") (- np j)) p) (display " " p)
   (display (list-ref (residual-history "x-velocity") (- np j)) p) (display " " p)
   (display (list-ref (residual-history "y-velocity") (- np j)) p) (display " " p)
   (display (list-ref (residual-history "z-velocity") (- np j)) p) (display " " p)
   (display (list-ref (residual-history "k") (- np j)) p) (display " " p)
   (display (list-ref (residual-history "omega") (- np j)) p)
   (newline p)
   (loop (+ i 1))
 )
)
)
) )
    
(output-port (open-output-file "residual_2000_e387_udf.dat")))
(writefile output-port)
(close-output-port output-port))
#+END_SRC 

   
}



