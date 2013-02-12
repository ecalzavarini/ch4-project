#include "common_main.h"

int 
main(int argc, char **argv)
{

	/* Initialize MPI */
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &me);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	if (ROOT)
		fprintf(stderr, "Hello world! I am ROOT: processor %d of %d\n", me, nprocs);

	assign_parameters();
	processor_splitting();
	allocate_fields();
	read_mesh();
	compute_volumes();
#ifdef LB
	design_lb();
#endif

	/*
	initial_conditions(); 
	boundary_conditions(); 
	*/

#ifdef FLUID
	/*
	hydro_fields();
	compute_advection();
	compute_collision();
	add_forcing();
	boundary_conditions();
	time_stepping();
	 */
#endif

	/* Shut down MPI */
	MPI_Finalize();
	return 0;
}
