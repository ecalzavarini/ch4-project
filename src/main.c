#include "common_main.h"

int main(int argc, char **argv){

  int i;
 
        initialization_MPI(&argc, &argv);
	assign_parameters();
	processor_splitting();
	allocate_fields();
#ifdef LB
	design_lb();
#endif
	read_mesh();
	compute_volumes();
	initial_conditions(); 
	/*
	boundary_conditions(); 
	*/
	hydro_fields();

	///*
	//for (time_now=0.0; time_now<max_max; time_now += time_dt) {
	  for (i=0; i<2000; i++) {
	  if(i%10==0 && ROOT) fprintf(stderr,"time step %d\n",i);
    //*/

#ifdef LB_FLUID
	  sendrecv_borders_pop(p);
	  /*
	  boundary_conditions();
	  */
	  compute_advection(p,rhs_p);
	  add_collision(p,rhs_p);
	  /*
	  add_forcing();
	  */
	  time_stepping(p,rhs_p,old_rhs_p);
	 hydro_fields();	
	 dump_averages(i);

#endif
	 ///*
   }
  //*/
	/* Shut down MPI */
	MPI_Finalize();
	return 0;
}/* end main */
