#include "common_main.h"

int main(int argc, char **argv)
{

	/*        initialization_MPI(&argc, &argv); */
	initialization_MPI(argc, argv);
	assign_parameters();
	processor_splitting();
	allocate_fields();
#ifdef LAGRANGE
	allocate_particles();
#endif
#ifdef LB
	design_lb();
#endif
	read_mesh(); /* the mesh is read */
	compute_volumes();
	compute_interpolation_coefficients();
#ifdef LB_FLUID_FORCING_LANDSCAPE
	read_landscape(); /* is there any topography ? */
#endif
	initial_conditions(resume);
	initialization_forcing();

	/* first recontruction of hydrodynamic fields */
#ifdef LB_FLUID
	hydro_fields('p');
#endif
#ifdef LB_TEMPERATURE
	hydro_fields('g');
#endif
#ifdef LB_SCALAR
	hydro_fields('h');
#endif

	/* second (more accurate ) reconstruction of hydrodynamic fields , taking into account the forcing terms */
#if (defined LB_FLUID_FORCING || defined LB_TEMPERATURE_FORCING || defined LB_SCALAR_FORCING)
	build_forcing();
#ifdef LB_FLUID
	hydro_fields('p');
#endif
#ifdef LB_TEMPERATURE
	hydro_fields('g');
#endif
#ifdef LB_SCALAR
	hydro_fields('h');
#endif
#endif

#ifdef LAGRANGE
	initial_conditions_particles(resume);
#endif

#ifdef SYSTEM_TIMING
	t1 = MPI_Wtime();
#endif

	itime = 0;
	time_now = property.time_start;

#ifdef OUTPUT_AT_START
	dump_averages();
#ifdef LAGRANGE
	output_particles();		  /* binary or h5 */
	dump_particle_averages(); /* ascii diagnostic averages */
#endif
#endif
	for (time_now = property.time_start + property.time_dt; time_now <= property.time_end; time_now += property.time_dt)
	{
		itime++;
		if (itime % 10 == 0 && ROOT)
			fprintf(stderr, "time step %d\n", itime);

#ifndef METHOD_STREAMING
#ifndef METHOD_REDEFINED_POP
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
		/* if REDEFINED_POP is defined the BC are computed for f_aux in advection */
#endif
#endif
#endif /* end of NOT DEFINED method STREAMING */

#ifdef LB_FLUID
		compute_advection(p, rhs_p, property.tau_u, p_eq, 'p');
#ifdef METHOD_HEUN
		copy_pop(rhs_p, old_rhs_p);
#endif
		add_collision(p, rhs_p, property.tau_u, p_eq, 'p');
#endif

#ifdef LB_TEMPERATURE
		compute_advection(g, rhs_g, property.tau_t, g_eq, 'g');
#ifdef METHOD_HEUN
		copy_pop(rhs_g, old_rhs_g);
#endif
		add_collision(g, rhs_g, property.tau_t, g_eq, 'g');

#ifdef LB_TEMPERATURE_CHT_ZIQI						  /* correction for Conjugate Heat Transfer (CHT) */
		temperature_cht(rhs_g, liquid_frac, u, dens); /* Added after collision */
#endif
#endif

#ifdef LB_SCALAR
		compute_advection(h, rhs_h, property.tau_s, h_eq, 'h');
#ifdef METHOD_HEUN
		copy_pop(rhs_h, old_rhs_h);
#endif
		add_collision(h, rhs_h, property.tau_s, h_eq, 'h');
#endif

#if (defined LB_FLUID_FORCING || defined LB_TEMPERATURE_FORCING || defined LB_SCALAR_FORCING)
		build_forcing();
#ifdef LB_TEMPERATURE_MELTING
		melting();
#endif
#ifdef LAGRANGE_POLYMER_FEEDBACK
		/* Add the divergence of the extra stress tensor */
		add_lagrangian_polymer_feedback_on_the_flow();
#endif
		add_forcing();
#endif

#ifdef METHOD_STREAMING
		boundary_and_pbc_conditions_for_streaming();
#endif

#ifdef LB_FLUID
		time_stepping(p, rhs_p, old_rhs_p, old_old_rhs_p, property.tau_u, p_eq, 'p');
#endif
#ifdef LB_TEMPERATURE
		time_stepping(g, rhs_g, old_rhs_g, old_old_rhs_g, property.tau_t, g_eq, 'g');
#endif
#ifdef LB_SCALAR
		time_stepping(h, rhs_h, old_rhs_h, old_old_rhs_h, property.tau_s, h_eq, 'h');
#endif

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
#ifdef LAGRANGE_NUCLEATE
		shake_ghosts();
#endif
		boundary_conditions_hydro();
		interpolate_vector_at_particles(u, 'u');
#ifdef LB_TEMPERATURE
		interpolate_scalar_at_particles(t, 't');
#endif
#ifdef LB_SCALAR
		interpolate_scalar_at_particles(s, 's');
#endif
#ifdef LAGRANGE_NUCLEATE
		nucleate_particles(t);
#endif
		move_particles();
		output_particles();		  /* binary or h5 */
		dump_particle_averages(); /* ascii diagnostic averages */
#endif

	  /* now the extra EULERIAN equations (implemented by Finite-volumes method) */
#ifdef EULER_PARTICLE
compute_rhs_conc();
time_stepping_scalar(conc,rhs_conc,old_rhs_conc);
#endif 

	} /* loop on time: time_now */

#ifdef SYSTEM_TIMING
	t2 = MPI_Wtime();
	measure_time();
#endif

#ifdef OUTPUT_H5
		write_pop_h5();
	#ifdef LB_TEMPERATURE_MELTING
		write_scalar_h5(liquid_frac, "liquid_frac");
	#endif
	#ifdef LAGRANGE
		write_point_particle_h5();
	#endif
	#ifdef EULER_PARTICLE
		write_scalar_h5(conc, "euler_particle_concentration");
	#endif 
#endif

	/* Shut down */
	free_fields();
	if (ROOT)
		fprintf(stdout, "The End\n");
	MPI_Finalize();
	return 0;
} /* end main */
