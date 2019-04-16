/*******************************************************/
/* Used to obtain a time history of spanwise force report on 3D blade/cylinder/wing */
/*******************************************************/
#include "udf.h"
#include "para.h"
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
/* It is however limited by the number of cells present along the spanwise direction
*/
/* Recommended value: */
/* For unstructured cells: no of cells along spanwise dir. divided by 3 */
/* For structured cells: no. of cells along spanwise dir. divided by 3 */
#define N_STRIPS 20
/* ID of the wing boundary zone */
/* You can get this from the Boundary Conditions Panel */
#define TID 5

/* Make sure you increase no. of User Defined Memory to 1 */
/* The UDM stores the strip number for each cell */
/* You can see contour plots of cell values of UDM-1 to check */
/* if the strips are marked correctly */
DEFINE_EXECUTE_AT_END(force) /* execute at every iteration (steady), or time step( unsteady)*/
{
  #if !RP_HOST
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
	real pforce[N_STRIPS+1],vforce[N_STRIPS+1];
	real press=0.0;
	real pressure=0.0;
	real visc=0.0;
	real min_dist=1000000000;
	real max_dist=-1000000000;
	real dist=0;
	real strip_length;
	FILE *fp;
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
	
  #if PARALLEL
  if (I_AM_NODE_ZERO_P)
	fp=fopen("spanwise-force-report.txt", "a+");  /* open file and write data; "a+": append and update*/
  #else
    fp=fopen("spanwise-force-report.txt", "a+");
  #endif
	begin_f_loop(f, t)
	if (PRINCIPAL_FACE_P(f,t))
	{
		F_CENTROID(x, f, t);
		dist=NV_DOT(x, axis);
		if(dist<min_dist) min_dist=dist;
		if(dist>max_dist) max_dist=dist;
	}
	end_f_loop(f, t)
	
  #if PARALLEL
  /* get global min and max on node 0*/ 
	min_dist=PRF_GRLOW1(min_dist);
	max_dist=PRF_GRHIGH1(max_dist);
  #endif	
	strip_length=(max_dist-min_dist)/N_STRIPS;
	Message0("Marking strips...");
	
	begin_f_loop(f,t)
	{
		F_CENTROID(x,f,t);
		c=F_C0(f, t);
		tc=THREAD_T0(t);
		dist=NV_DOT(x, axis);
		/* assign strip number to UDM-0 */
		C_UDMI(c,tc,0)=(int)((dist-min_dist)/strip_length);
	}
	end_f_loop(f,t)
	Message0("Done\n");
	
	for(i=0;i<=N_STRIPS;i++)
	{
		pforce[i]=0;
		vforce[i]=0;
	}
 
	begin_f_loop(f,t)
	if (PRINCIPAL_FACE_P(f,t))
	{
		F_AREA(A,f,t);
		c=F_C0(f, t);
		tc=THREAD_T0(t);
		i=(int)(MIN(N_STRIPS,C_UDMI(c, tc, 0)));
		/* report force if thread id TID corresponds to a wall */
     if (THREAD_TYPE(t)==THREAD_F_WALL)
	 {
  /* pressure force of the strip i along vector (F_X,F_Y,F_Z) */		
		pforce[i]+=F_P(f, t) * NV_DOT(A, force_dir);
  /* viscous force of the strip i along vector (F_X,F_Y,F_Z) */				
		vforce[i]+=(-1*(NV_DOT(F_STORAGE_R_N3V(f,t,SV_WALL_SHEAR),force_dir)));	
		/*Message("CHECK: %d\n",sizeof(vforce) / sizeof(vforce[0]));*/	
	 }
	 else
	 {
	 if (f==0) Message0("No forces will be reported since, zone of ID %d is not a wall",TID);
	 }
	}
	end_f_loop(f,t)
	Message0("\n\nPressure force\t Viscous force");
	
  #if PARALLEL
/* Get forces of each partition on node 0 */
   	for(i=0;i<=N_STRIPS;i++)
	{
		pforce[i]=PRF_GRSUM1(pforce[i]);
		vforce[i]=PRF_GRSUM1(vforce[i]);
	}
  /* report force on node 0 in parallel- Node0 should have a disk system otherwise UDF has to be rewritten to exchange data through the host */	
  if (I_AM_NODE_ZERO_P)
  {
	for(i=0;i<=N_STRIPS;i++)
	{
		press+=pforce[i];
		visc+=vforce[i];
		Message0("\n %g \t %g",pforce[i],vforce[i]);
		fprintf(fp, "%d \t %f \t %f \n", i, pforce[i], vforce[i]);
	}
	Message0("\n---------- -----------");
	Message0("\n %g\t%g",press,visc);

	fclose(fp);
  }
  #else
    for(i=0;i<=N_STRIPS;i++)
	{
		press+=pforce[i];
		visc+=vforce[i];
		Message("\n %g \t %g",pforce[i],vforce[i]);
		fprintf(fp, "%d \t %f \t %f \n", i, pforce[i], vforce[i]);
	}
	Message("\n---------- -----------");
	Message("\n %g\t%g",press,visc);

	fclose(fp);
  #endif
	Message0("\nUDF DONE\n");
  #endif	/* !RP_HOST */
}
 
