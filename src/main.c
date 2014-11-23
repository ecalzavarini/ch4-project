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
#ifdef LB_FLUID
	hydro_fields('p');
#endif
#ifdef LB_TEMPERATURE
        hydro_fields('g');
#endif
#ifdef LB_SCALAR
        hydro_fields('h');
#endif
	//	dump_averages();
	//      exit(1);
#ifdef LAGRANGE
        initial_conditions_particles(resume);
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

	  /* COMPUTE & ADD  COLLISION */
#ifdef LB_FLUID
          copy_pop(p,rhs_p);
	  add_collision(p,rhs_p,property.tau_u,p_eq,'p');
#endif

#ifdef LB_TEMPERATURE
          copy_pop(g,rhs_g);
	  add_collision(g,rhs_g,property.tau_t,g_eq,'g');
#endif

#ifdef LB_SCALAR
          copy_pop(h,rhs_h);
	  add_collision(h,rhs_h,property.tau_s,h_eq,'h');
#endif
	  /* COMPUTE & ADD  FORCING */
#if (defined LB_FLUID_FORCING || defined LB_TEMPERATURE_FORCING || defined LB_SCALAR_FORCING )
	  build_forcing();
#ifdef LB_TEMPERATURE_MELTING
	  melting();
#endif
	  add_forcing();
#endif

	  /* FIX PERIODICITY */
#ifdef METHOD_STREAMING
 #ifdef LB_FLUID
 sendrecv_borders_pop(rhs_p);
 #endif

 #ifdef LB_TEMPERATURE
 sendrecv_borders_pop(rhs_g);
 #endif

 #ifdef LB_SCALAR
 sendrecv_borders_pop(rhs_h);
 #endif
#else /* if not streaming, than is finite volume */
 #ifdef LB_FLUID
	  sendrecv_borders_pop(p);	 
 #endif
 #ifdef LB_TEMPERATURE
	  sendrecv_borders_pop(g);
 #endif
 #ifdef LB_SCALAR
	  sendrecv_borders_pop(h);
 #endif
#endif

	  /* IMPLEMENT BOUNDARY CONDITIONS */
#if (defined LB_FLUID_BC || defined LB_TEMPERATURE_BC || defined LB_SCALAR_BC)
 #ifdef METHOD_STREAMING
         boundary_conditions_for_streaming();
 #else
      	 boundary_conditions();
	  /* if REDEFINED_POP is defined the BC are computed on f_aux in advection */
 #endif
#endif 

	 /* COMPUTE & ADD ADVECTION or PERFORM STREAMING */
#ifdef LB_FLUID
	  compute_advection(p,rhs_p,property.tau_u,p_eq,'p');
#ifdef METHOD_HEUN
          copy_pop(rhs_p,old_rhs_p);
#endif
#endif

#ifdef LB_TEMPERATURE
	  compute_advection(g,rhs_g,property.tau_t,g_eq,'g');
#ifdef METHOD_HEUN
	  copy_pop(rhs_g,old_rhs_g);
#endif
#endif

#ifdef LB_SCALAR
	  compute_advection(h,rhs_h,property.tau_s,h_eq,'h');
#ifdef METHOD_HEUN
	  copy_pop(rhs_h,old_rhs_h);
#endif
#endif 

	  /* TIME STEPPING FOR FINITE VOLUME */
#ifndef METHOD_STREAMING
#ifdef LB_FLUID	  
	 time_stepping(p,rhs_p,old_rhs_p,old_old_rhs_p,property.tau_u,p_eq,'p');
#endif
#ifdef LB_TEMPERATURE
	 time_stepping(g,rhs_g,old_rhs_g,old_old_rhs_g,property.tau_t,g_eq,'g');
#endif
#ifdef LB_SCALAR
	 time_stepping(h,rhs_h,old_rhs_h,old_old_rhs_h,property.tau_s,h_eq,'h');
#endif	
#endif

	 /* COMPUTE MACROSCOPIC FIELDS */
#ifdef LB_FLUID
	 hydro_fields('p');
#endif
#ifdef LB_TEMPERATURE
	 hydro_fields('g');
#endif
#ifdef LB_SCALAR
	 hydro_fields('h');
#endif
	 dump_averages();

	 /* Now the lagrangian part */
#ifdef LAGRANGE       
	boundary_conditions_hydro();
        interpolate_vector_at_particles(u,'u');
#ifdef LB_TEMPERATURE
        interpolate_scalar_at_particles(t,'t');
#endif
#ifdef LB_SCALAR
        interpolate_scalar_at_particles(s,'s');
#endif
	move_particles();
        if(itime%10==0)output_particles();
#endif
	}/* loop on time: time_now */

#ifdef TIMING
	t2 = MPI_Wtime();
	measure_time();
#endif

#ifdef OUTPUT_H5
       	write_pop_h5();
#ifdef LAGRANGE
	write_point_particle_h5();
#endif
#endif

	/* Shut down */
       free_fields();
       if(ROOT)fprintf(stdout,"The End\n");
	MPI_Finalize();
	return 0;
}/* end main */
