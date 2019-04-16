/* heat source*/

#include "udf.h"
DEFINE_PROFILE(temp_wall_profile,t,i)
{
real x[ND_ND]; /* this holds the position vector */
face_t f; /* Index for each face on the boundary */
/* Loop over all faces of the boundary */
begin_f_loop(f,t)
{
F_CENTROID(x,f,t); /* Get the global coordinates of the face */
/* Calculate the temperature and assign to boundary */
F_PROFILE(f,t,i) = -3200.0*x[0]*x[0] + 1600.0*x[0] + 300.0;
}
end_f_loop(f,t)
}

DEFINE_SOURCE(my_heat_source, c, t, ds, eqn)
{
const real C1 = 1.0, C2 = 0.03;
real source;
source = -con*fabs(C_U(c, t))*C_U(c,t);
ds[eqn] = -2.*con*fabs(C_U(c,t));
return source;
}
