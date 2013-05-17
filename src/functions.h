#include "define.h"

/* parallel.c */
void initialization_MPI();

/* parameters.c */
my_double read_parameter();
void assign_parameters();
void processor_splitting();
void allocate_fields();

/* grid.c */
void read_mesh();
void compute_volumes();
void compute_interpolation_coefficients();

/* lb.c */
void design_lb();
void sendrecv_borders_pop(pop *f);
pop equilibrium(pop * f, int i, int j, int k);
pop equilibrium_given_velocity(vector v , my_double rho);
void time_stepping(pop *f, pop *rhs_f,pop *old_rhs_f,my_double tau);

/* boundary_conditions.c */
void boundary_conditions();

/* initial_conditions.c */
void initial_conditions(int restart);

/* fluid.c */
void add_collision(pop * f, pop *rhs_f,my_double tau);
void compute_advection(pop * f, pop *rhs_f);
void hydro_fields();
tensor strain_tensor(pop *f,int i, int j, int k);
vector gradient_scalar(my_double *t, int i, int j, int k);
tensor gradient_vector(vector *t, int i, int j, int k);
void add_forcing();
void build_forcing();

/* output.c */
void dump_averages();

/* output_h5.c */
void write_pop_h5();
void read_pop_h5();