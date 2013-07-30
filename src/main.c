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
	initial_conditions(resume); 
	hydro_fields();
	//	dump_averages();
	//      exit(1);

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
#ifdef LB_SCALAR
	  sendrecv_borders_pop(h);
#endif


#if (defined LB_FLUID_BC || defined LB_TEMPERATURE_BC || defined LB_SCALAR_BC)
#ifndef METHOD_STREAMING
	  boundary_conditions();
#endif
#endif
	  
#ifdef LB_FLUID
	  compute_advection(p,rhs_p);
	  add_collision(p,rhs_p,property.tau_u);
#endif
#ifdef LB_TEMPERATURE
	  compute_advection(g,rhs_g);
	  add_collision(g,rhs_g,property.tau_t);
#endif
#ifdef LB_SCALAR
	  compute_advection(h,rhs_h);
	  add_collision(h,rhs_h,property.tau_s);
#endif 

#if (defined LB_FLUID_FORCING || defined LB_TEMPERATURE_FORCING)
	  build_forcing();
#ifdef LB_TEMPERATURE_MELTING
	  melting();
#endif
	  add_forcing();
#endif

 
#ifdef METHOD_STREAMING
         boundary_conditions_for_streaming();
#endif

#ifdef LB_FLUID	  
	  time_stepping(p,rhs_p,old_rhs_p,property.tau_u);
#endif
#ifdef LB_TEMPERATURE
	  time_stepping(g,rhs_g,old_rhs_g,property.tau_t);
#endif
#ifdef LB_SCALAR
	  time_stepping(h,rhs_h,old_rhs_h,property.tau_s);
#endif	
	 hydro_fields();	
	 dump_averages();
	}

#ifdef OUTPUT_H5
       	write_pop_h5();
#endif

	/* Shut down MPI */
	MPI_Finalize();
	return 0;
}/* end main */
