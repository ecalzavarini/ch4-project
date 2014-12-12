#ifdef NO_XMLHEADERS
 #include "define.h"
#endif

/* parallel.c */
void initialization_MPI();
void measure_time();
void sum_output(output *a, output *b, int *length, MPI_Datatype *dtype);
void sum_vector(vector *a, vector *b, int *length, MPI_Datatype *dtype);

/* parameters.c */
my_double read_parameter();
void assign_parameters();
void processor_splitting();
void allocate_fields();
void free_fields();
void set_to_zero_output(output  *f,int size);

/* grid.c */
void read_mesh();
void compute_volumes();
void compute_interpolation_coefficients();
void read_landscape();
void sendrecv_borders_vector(vector *f);
void sendrecv_borders_scalar(my_double *f);

/* lb.c */
void design_lb();
void sendrecv_borders_pop(pop *f);
pop equilibrium(pop * f, int i, int j, int k);
pop equilibrium_given_velocity(vector v , my_double rho);
void time_stepping(pop *f, pop *rhs_f,pop *old_rhs_f,pop *old_old_rhs_f,my_double tau,pop *f_eq,char which_pop);
void copy_pop(pop *f, pop *f_copy);
void streaming(pop *f, pop *rhs_f,int i,int j,int k);

/* boundary_conditions.c */
void boundary_conditions();
void boundary_and_pbc_conditions_for_streaming();
void boundary_conditions_for_equilibrium(char which_pop);
void boundary_conditions_for_hydro(char which_pop);
void boundary_conditions_for_advection(pop * f, char which_pop);

/* initial_conditions.c */
void initial_conditions(int restart);
void initialization_forcing();

/* fluid.c */
void add_collision(pop * f, pop *rhs_f,my_double tau,pop *f_eq,char which_pop);
void compute_advection(pop * f, pop *rhs_f, my_double tau, pop *f_eq,char which_pop);
void hydro_fields(char which_pop);
tensor strain_tensor(pop *f,int i, int j, int k);
vector gradient_scalar(my_double *t, int i, int j, int k);
tensor gradient_vector(vector *t, int i, int j, int k);
void add_forcing();
void build_forcing();
#ifdef METHOD_CENTERED
my_double compute_flux_with_central_difference(pop * f,int i, int j, int k, int pp);
#endif
#ifdef METHOD_UPWIND
my_double compute_flux_with_upwind_first(pop * f, int i, int j, int k, int pp);
#endif
#ifdef METHOD_UPWIND_LINEAR
my_double compute_flux_with_upwind_linear(pop * f, int i, int j, int k, int pp);
#endif
#ifdef METHOD_MYQUICK
my_double compute_flux_with_quick(pop * f, int i, int j, int k, int pp);
#ifdef METHOD_MYQUICK_LIMITER
my_double compute_flux_with_limiters(pop * f, int i, int j, int k, int pp);
my_double limiter(my_double r);
my_double gradient_ratio(pop * f, int i, int j, int k, int pp,int dir);
#endif
#endif

/* output.c */
void dump_averages();

/* output_h5.c */
void write_pop_h5();
void read_pop_h5();
void output_h5();

/* lagrange */
#ifdef LAGRANGE 
void allocate_particles();
void initial_conditions_particles(int restart);
void interpolate_vector_at_particles(vector *f, char which_vector);
void interpolate_scalar_at_particles(my_double *f, char which_scalar); 
void output_particles();
void move_particles();
void sendrecv_particles();
void boundary_conditions_hydro();
void write_point_particle_h5();
void read_point_particle_h5();
#endif

/* random.c */
double myrand();
#ifndef RANDOM48
double ran1();
#endif

/* LES */
#ifdef LB_FLUID_LES
my_double tau_u_les(int i , int j , int k);
#endif
#ifdef LB_TEMPERATURE_LES
my_double tau_t_les(int i , int j , int k);
#endif
#ifdef LB_SCALAR_LES
my_double tau_s_les(int i , int j , int k);
#endif
