#include "common_main.h"

int main(int argc, char **argv){

        initialization_MPI(&argc, &argv);
	assign_parameters();
	processor_splitting();
	allocate_fields();
	read_mesh();
	compute_volumes();
#ifdef LB
	design_lb();
#endif	
	initial_conditions(); 
	/*
	boundary_conditions(); 
	*/
	hydro_fields();

	/*	
  for (itime=0; itime<max_step; itime++) {
    if(itime%100==0) fprintf(stderr,"time step %d\n",itime);
	*/

#ifdef LB_FLUID
	/*
	boundary_conditions();
	*/
	compute_advection(p);
	compute_collision(p);
	/*
	add_forcing();
	time_stepping();
	hydro_fields();
	 */
#endif
	/*
   }
	*/
	/* Shut down MPI */
	MPI_Finalize();
	return 0;
}/* end main */
