#include "common_main.h"

int main(int argc, char **argv){
 
        initialization_MPI(&argc, &argv);
	assign_parameters();
	processor_splitting();
	allocate_fields();
#ifdef LB
	design_lb();
#endif
	read_mesh();
	compute_volumes();
#ifdef METHOD_CENTERED
	compute_interpolation_coefficients();
#endif
#ifdef LB_BC
	prepare_boundary_conditions();
#endif
	initial_conditions(); 
	hydro_fields();


	itime=0;
	for (time_now=property.time_start; time_now <= property.time_end; time_now += property.time_dt){
	  itime++;
	  if(itime%10==0 && ROOT) fprintf(stderr,"time step %d\n",itime);


#ifdef LB_FLUID
	  sendrecv_borders_pop(p);
#ifdef LB_BC	  
	  boundary_conditions(p);
#endif
	  
	  compute_advection(p,rhs_p);
	  add_collision(p,rhs_p);

#ifdef LB_FLUID_FORCING
	   build_forcing();
	   add_forcing(p,rhs_p);
#endif
	  
	 time_stepping(p,rhs_p,old_rhs_p);
	 hydro_fields();	
	 dump_averages();

#endif
   }

	/* Shut down MPI */
	MPI_Finalize();
	return 0;
}/* end main */
