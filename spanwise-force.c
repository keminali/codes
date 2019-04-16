


/* Used to obtain a spanwise force report on blade/wings */
/*******************************************************/

#include "udf.h"

/* Spanwise direction */
#define AXIS_X 0
#define AXIS_Y 0
#define AXIS_Z 1

/* The direction along which the force is to be reported */
#define F_X 1
#define F_Y 0
#define F_Z 0

/* Number of strips into which the wing is to be divided */
/* Increasing this number increases the resolution of reporting */
/* It is however limited by the number of cells present along the spanwise direction */
/* Recommended value:  number of cells along spanwise dir. divided by 3 */

#define N_STRIPS 20

/* ID of the wing boundary zone */
/* You can get this from the Boundary Conditions Panel */
#define TID 5


/* Make sure you increase no. of User Defined Memory to 1 */
/* The UDM stores the strip number for each cell */
/* You can see contour plots of cell values of UDM-1 to check  if the strips are marked correctly */

DEFINE_ON_DEMAND(force)  /*a general-purpose macro which specify a UDF that is executed on demand" */
{
	Domain *d=Get_Domain(1);
 int i;
 real A[ND_ND];
 real axis[3], force_dir[3];
 face_t f;
 Thread *t=Lookup_Thread(d,TID);
 cell_t c;
 Thread *tc;
 real x[ND_ND];
 real NV_VEC(presforce),NV_VEC(viscforce);
 real NV_VEC(dforce_pres),NV_VEC(dforce_visc);
 real pforce[N_STRIPS],vforce[N_STRIPS];
 real press=0.0;
 real pressure=0.0;
 real visc=0.0;
 real min_dist=1000000000;
 real max_dist=-1000000000;
 real dist=0;
 real strip_length;
 FILE *fp;   /* Pointer to the file you are reading or writing*/

 NV_S(A,=,0.0);
 NV_S(dforce_pres,=,0.0);
 NV_S(dforce_visc,=,0.0);
 NV_S(presforce,=,0.0);
 NV_S(viscforce,=,0.0);

axis[0]=AXIS_X;
axis[1]=AXIS_Y;
axis[2]=AXIS_Z;

force_dir[0]=F_X;
force_dir[1]=F_Y;
force_dir[2]=F_Z;

 fp=fopen("spanwise-force-report.txt", "w+");  /* open a file, write/updte */

begin_f_loop(f, t)
 {
 F_CENTROID(x, f, t);
 dist=NV_DOT(x, axis);
 if(dist<min_dist) min_dist=dist;
 if(dist>max_dist) max_dist=dist;
 }
end_f_loop(f, t)

strip_length=(max_dist-min_dist)/N_STRIPS;

Message("Marking strips...");
begin_f_loop(f,t)
{
F_CENTROID(x,f,t);
c=F_C0(f, t);
tc=THREAD_T0(t);
dist=NV_DOT(x, axis);
C_UDMI(c,tc,0)=(int)((dist-min_dist)/strip_length);
}
end_f_loop(f,t)

Message("Donen");

for(i=0;i<N_STRIPS;i++) 
{
pforce[i]=0;
vforce[i]=0;
}

begin_f_loop(f,t)
{
F_AREA(A,f,t);
c=F_C0(f, t);
tc=THREAD_T0(t);
pforce[(int)C_UDMI(c, tc, 0)]+=F_P(f, t) * NV_DOT(A, force_dir);
vforce[(int)C_UDMI(c, tc, 0)]+=(-1*(NV_DOT(F_STORAGE_R_N3V(f,t,SV_WALL_SHEAR), force_dir)));
}
end_f_loop(f,t)

Message("nnPressure forcet Viscous force");
for(i=0;i<N_STRIPS;i++)
{
press+=pforce[i];
visc+=vforce[i];
 Message("n %g t %g",pforce[i],vforce[i]); /* The "Message" macro is a utility that displays data to the console. */
 fprintf(fp, "%d t %f t %f \n", i, pforce[i], vforce[i]); /* "fprintf" is a standard C function for writing to a file; */
}
Message("n---------- -----------");
	Message("n %gt%g",press,visc);
	fclose(fp);
} 

