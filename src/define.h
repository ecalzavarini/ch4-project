/* DEBUG */
//#define DEBUG
//#define DEBUG_HARD

/* PARALLEL*/
#define NEW_SENDRECV  /* does not use itermediate transfer buffers */

/* OUTPUT */
#define OUTPUT_ASCII
#define OUTPUT_H5
#define TIMING

/* GRID */
#define GRID_REGULAR
//#define GRID_REFINED /* the default refinement is TANH */ 
//#define GRID_REFINED_ONESIDED
//#define GRID_REFINED_BULK
//#define GRID_REFINED_CHEBYSHEV
//#define GRID_REFINED_SINH
//#define GRID_REFINED_X
//#define GRID_REFINED_Y
//#define GRID_REFINED_Z

//#define GRID_IRREGULAR_RANDOM

/* LB */
#define LB

/* FLUID */
#define LB_FLUID
//#define LB_FLUID_INITIAL_KOLMOGOROV 
#define LB_FLUID_INITIAL_POISEUILLE
//#define LB_FLUID_INITIAL_CHANNEL
//#define LB_FLUID_INITIAL_VORTICES 
//#define LB_FLUID_INITIAL_PERTURBATION
//#define LB_INITIAL_BAROMETRIC
#define LB_FLUID_FORCING
//#define LB_FLUID_FORCING_CHANNEL
#define LB_FLUID_FORCING_POISEUILLE
//#define LB_FLUID_FORCING_CONSTANT_POWER
//#define LB_FLUID_FORCING_KOLMOGOROV
//#define LB_FLUID_FORCING_LANDSCAPE
//#define LB_FLUID_FORCING_PENALIZATION
//#define LB_FLUID_FORCING_DIRECT
#define LB_FLUID_BC
#define LB_FLUID_BC_Y  /* default fluid bc is no-slip */
//#define LB_FLUID_BC_YP_SLIP                                           
//#define LB_FLUID_BC_YM_SLIP                         
//#define LB_FLUID_BC_YP_OUTLET
//#define LB_FLUID_BC_X                                                     
//#define LB_FLUID_BC_XP_SLIP                                             
//#define LB_FLUID_BC_XM_SLIP                                               
//#define LB_FLUID_BC_XM_INLET                                         
//#define LB_FLUID_BC_XP_OUTLET                                               
//#define LB_FLUID_BC_Z  

/* METHODS for time stepping or convective term */
/* For smooth simulations 
Activate either METHOD_FINITE_VOLUME or METHOD_STREAMING */
#define METHOD_FINITE_VOLUME
//#define METHOD_STREAMING

//#define METHOD_REDEFINED_POP
//#define METHOD_LOG
//#define METHOD_EXPONENTIAL
//#define METHOD_STEPPING_EULER
//#define METHOD_STEPPING_AB2
//#define METHOD_STEPPING_AB3
//#define METHOD_COLLISION_IMPLICIT

//#define METHOD_CENTERED
//#define METHOD_UPWIND
//#define METHOD_UPWIND_SKEW
#define METHOD_UPWIND_LINEAR
#define METHOD_UPWIND_LINEAR_IMPROVED  
//#define METHOD_MYQUICK
//#define METHOD_MYQUICK_CARTESIAN
//#define METHOD_MYQUICK_LIMITER
//#define METHOD_MIXED

//#define METHOD_TRAPEZOID

//#define METHOD_EDGES_AND_CORNERS
//#define METHOD_FORCING_GUO

/* TEMPERATURE */
//#define LB_TEMPERATURE
//#define LB_TEMPERATURE_FLUCTUATION
//#define LB_TEMPERATURE_INITIAL_LINEAR
//#define LB_TEMPERATURE_INITIAL_ADD_PERTURBATION
//#define LB_TEMPERATURE_INITIAL_CONSTANT
//#define LB_TEMPERATURE_INITIAL_SPOT
//#define LB_TEMPERATURE_INITIAL_BL
//#define LB_TEMPERATURE_BUOYANCY
//#define LB_TEMPERATURE_BC
//#define LB_TEMPERATURE_BC_Y
//#define LB_TEMPERATURE_BC_Y_VARIABLE
//#define LB_TEMPERATURE_BC_X
//#define LB_TEMPERATURE_BC_X_NOFLUX
//#define LB_TEMPERATURE_BC_KEEP_WITHIN
//#define LB_TEMPERATURE_FORCING
//#define LB_TEMPERATURE_FORCING_BULK 
//#define LB_TEMPERATURE_FORCING_BULK_VARIABLE
//#define LB_TEMPERATURE_FORCING_RADIATION 
//#define LB_TEMPERATURE_FORCING_SOURCE
//#define LB_TEMPERATURE_FORCING_PROFILE
//#define LB_TEMPERATURE_FORCING_REACTION
//#define LB_TEMPERATURE_MELTING
//#define LB_TEMPERATURE_FORCING_DIRECT

/* EXTRA SCALAR FIELD e.g. CH4 */
//#define LB_SCALAR


/***************************************************/
/* Define dependencies not to be modified          */

#ifdef METHOD_FINITE_VOLUME                                            
#define METHOD_STEPPING_EULER                                           
#define METHOD_MYQUICK                                                       
#define METHOD_MYQUICK_CARTESIAN
#endif 

#ifdef METHOD_STREAMING                                                  
#define METHOD_EDGES_AND_CORNERS                                        
#define METHOD_FORCING_GUO                                              
#endif
