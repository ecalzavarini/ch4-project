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
	compute_interpolation_coefficients();

	//#ifdef LB_FLUID_BC
	//	prepare_boundary_conditions();
	//#endif
	initial_conditions(); 
	hydro_fields();
	//dump_averages();
	//exit(1);

	itime=0;
	for (time_now=property.time_start; time_now <= property.time_end; time_now += property.time_dt){
	  itime++;
	  if(itime%10==0 && ROOT) fprintf(stderr,"time step %d\n",itime);

#ifdef LB_FLUID
	  sendrecv_borders_pop(p);
#endif
#ifdef LB_TEMPERATURE
	  sendrecv_borders_pop(g);
#endif
#if (defined LB_FLUID_BC || defined LB_TEMPERATURE_BC)
	  boundary_conditions();
#endif
	  
#ifdef LB_FLUID
	  compute_advection(p,rhs_p);
	  add_collision(p,rhs_p);
#endif
#ifdef LB_TEMPERATURE
	  compute_advection(g,rhs_g);
	  add_collision(g,rhs_g);
#endif

#ifdef LB_FLUID_FORCING
	  build_forcing();
	  add_forcing(p,rhs_p);
#endif

#ifdef LB_FLUID	  
	 time_stepping(p,rhs_p,old_rhs_p);
#endif
#ifdef LB_TEMPERATURE
	 time_stepping(g,rhs_g,old_rhs_g);
#endif
	
	 hydro_fields();	
	 dump_averages();

	}

	/* Shut down MPI */
	MPI_Finalize();
	return 0;
}/* end main */
