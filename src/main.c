#include "common_main.h"

int main(int argc, char **argv){
 
        initialization_MPI(&argc, &argv);
	assign_parameters();
	processor_splitting();
	allocate_fields();
#ifdef LAGRANGE
        allocate_particles();
#endif
#ifdef LB
	design_lb();
#endif
	read_mesh();  /* the mesh is read */
	compute_volumes();
	compute_interpolation_coefficients();
#ifdef LB_FLUID_FORCING_LANDSCAPE
	read_landscape(); /* is there any topography ? */
#endif
	initial_conditions(resume); 
	hydro_fields();
	//	dump_averages();
	//      exit(1);
#ifdef LAGRANGE
        initial_conditions_particles();
        interpolate_vector_at_particles(u);
	//output_particles();
	//exit(-1);
#endif

#ifdef TIMING
	t1 = MPI_Wtime();
#endif

	itime=0;
	for (time_now=property.time_start; time_now <= property.time_end; time_now += property.time_dt){
	  itime++;
	  if(itime%10==0 && ROOT) fprintf(stderr,"time step %d\n",itime);

#ifndef METHOD_STREAMING
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

	  boundary_conditions();
#endif
#endif
	  
#ifdef LB_FLUID
	  compute_advection(p,rhs_p,property.tau_u,p_eq,'p');
#ifdef METHOD_HEUN
          copy_pop(rhs_p,old_rhs_p);
#endif
	  add_collision(p,rhs_p,property.tau_u,p_eq);
#endif

#ifdef LB_TEMPERATURE
	  compute_advection(g,rhs_g,property.tau_t,g_eq,'g');
#ifdef METHOD_HEUN
	  copy_pop(rhs_g,old_rhs_g);
#endif
	  add_collision(g,rhs_g,property.tau_t,g_eq);
#endif

#ifdef LB_SCALAR
	  compute_advection(h,rhs_h,property.tau_s,h_eq,'h');
#ifdef METHOD_HEUN
	  copy_pop(rhs_h,old_rhs_h);
#endif
	  add_collision(h,rhs_h,property.tau_s,h_eq);
#endif 

#if (defined LB_FLUID_FORCING || defined LB_TEMPERATURE_FORCING)
	  build_forcing();
#ifdef LB_TEMPERATURE_MELTING
	  melting();
#endif
	  add_forcing();
#endif

#ifdef METHOD_STREAMING
         boundary_and_pbc_conditions_for_streaming();
#endif

#ifdef LB_FLUID	  
	 time_stepping(p,rhs_p,old_rhs_p,old_old_rhs_p,property.tau_u,p_eq,'p');
#endif
#ifdef LB_TEMPERATURE
	 time_stepping(g,rhs_g,old_rhs_g,old_old_rhs_g,property.tau_t,g_eq,'g');
#endif
#ifdef LB_SCALAR
	 time_stepping(h,rhs_h,old_rhs_h,old_old_rhs_h,property.tau_s,h_eq,'h');
#endif	
	 hydro_fields();	
	 dump_averages();

#ifdef LAGRANGE       
	move_particles();
        output_particles();
#endif
	}/* loop on time: time_now */

#ifdef TIMING
	t2 = MPI_Wtime();
	measure_time();
#endif

#ifdef OUTPUT_H5
       	write_pop_h5();
#endif

	/* Shut down */
       free_fields();
       if(ROOT)fprintf(stdout,"The End\n");
	MPI_Finalize();
	return 0;
}/* end main */
