/* DEBUG */
//#define DEBUG
//#define DEBUG_HARD


/* OUTPUT */
//#define OUTPUT_ASCII
#define OUTPUT_H5

/* GRID */

#define GRID_REGULAR
//#define GRID_REFINED
//#define GRID_REFINED_ONESIDED
//#define GRID_REGULAR_CARTESIAN
//#define GRID_REGULAR_CARTESIAN_EQUISPACED
//#define GRID_REGULAR_CARTESIAN_EQUISPACED_X
//#define GRID_REGULAR_CARTESIAN_EQUISPACED_Y
//#define GRID_REGULAR_CARTESIAN_EQUISPACED_Z
//#define GRID_REGULAR_CARTESIAN_REFINED_X
//#define GRID_REGULAR_CARTESIAN_REFINED_Y
//#define GRID_REGULAR_CARTESIAN_REFINED_Z

//#define GRID_IRREGULAR_RANDOM

/* LB */
#define LB

/* FLUID */
#define LB_FLUID
//#define LB_FLUID_INITIAL_KOLMOGOROV 
//#define LB_FLUID_INITIAL_POISEUILLE
#define LB_FLUID_FORCING
//#define METHOD_FORCING_GUO
//#define LB_FLUID_FORCING_POISEUILLE
//#define LB_FLUID_FORCING_KOLMOGOROV
#define LB_FLUID_BC
#define LB_FLUID_BC_Y  /* default fluid bc is no-slip */
//#define LB_FLUID_BC_YP_SLIP                                                                        
//#define LB_FLUID_BC_YM_SLIP                         
//#define LB_FLUID_BC_X                                                                              
//#define LB_FLUID_BC_XP_SLIP                                                      
//#define LB_FLUID_BC_XM_SLIP                                                                             
//#define LB_FLUID_BC_XM_INLET                                         
//#define LB_FLUID_BC_XP_OUTLET                                               
//#define LB_FLUID_BC_Z  

/* METHODS for time stepping or convective term */
//#define METHOD_EXPONENTIAL
//#define METHOD_STEPPING_EULER
#define METHOD_STEPPING_AB2
//#define METHOD_MIXED
//#define METHOD_CENTERED
//#define METHOD_UPWIND
#define METHOD_MYQUICK

/* TEMPERATURE */
#define LB_TEMPERATURE
#define LB_TEMPERATURE_INITIAL_LINEAR
//#define LB_TEMPERATURE_INITIAL_ADD_PERTURBATION
//#define LB_TEMPERATURE_INITIAL_CONSTANT
//#define LB_TEMPERATURE_BUOYANCY
#define LB_TEMPERATURE_BC
#define LB_TEMPERATURE_BC_Y

/* EXTRA SCALAR FIELD e.g. CH4 */
//#define LB_SCALAR
