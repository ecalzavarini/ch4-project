/* DEBUG */
//#define DEBUG
//#define DEBUG_MEMORY
//#define DEBUG_HARD
//#define DEBUG_MESH
//#define DEBUG_CENTERS
//#define DEBUG_LAGRANGE_INTERPOLATE

/* SYSTEM RELATED */
#define SYSTEM_RANDOM48  /* does enable drand48 , default is our random number generator */
#define SYSTEM_TIMING /* measure computetion time */

/* PARALLEL (MPI related) */
#define NEW_SENDRECV  /* improve speed: does not use itermediate transfer buffers */
//#define PARALLEL_MPICART /* defines an optimized cartesian grid of processor (DO NOT USE YET)*/

/* OUTPUT */
//#define OUTPUT_DIAGN_APPEND
//#define OUTPUT_ASCII
#define OUTPUT_H5
#define OUTPUT_H5_GRID /* write the grid in every h5 file, not essential but useful for VisIt visualizations */
#define OUTPUT_H5_TIMESTAMP_REAL /* enable "time_now" as time stamp on the name of .h5 files instead of itime */
#define OUTPUT_AT_START /* write down the complete intial state */
//#define VERBOSE /* if not defined, the amount of output communication is limited */

/* GRID */
#define GRID_REGULAR
//#define GRID_REFINED /* the default refinement is TANH */ 
//#define GRID_REFINED_ONESIDED
//#define GRID_REFINED_BULK
//#define GRID_REFINED_CHEBYSHEV
//#define GRID_REFINED_SINH
//#define GRID_REFINED_X  /* activate direction to be refined */
//#define GRID_REFINED_Y
//#define GRID_REFINED_Z

//#define GRID_IRREGULAR_RANDOM

/* Our default choice is D3Q19 and is defined at the end of this file */
//#define GRID_POP_D1Q3
//#define GRID_POP_D2Q9 
//#define GRID_POP_D3Q15
#define GRID_POP_D3Q19
//#define GRID_POP_D3Q27

/* LB */
#define LB

/* FLUID */
#define LB_FLUID 
//#define LB_FLUID_PAST
//#define LB_FLUID_INITIAL_KOLMOGOROV 
//#define LB_FLUID_INITIAL_POISEUILLE
//#define LB_FLUID_INITIAL_POISEUILLE_HALF
//#define LB_FLUID_INITIAL_CHANNEL
//#define LB_FLUID_INITIAL_VORTICES 
//#define LB_FLUID_INITIAL_CELLULAR
//#define LB_FLUID_INITIAL_STIRRER
#define LB_FLUID_INITIAL_PERTURBATION
//#define LB_FLUID_INITIAL_LANDSCAPE
//#define LB_INITIAL_BAROMETRIC
//#define LB_INITIAL_BULK
//#define LB_INITIAL_CONSTANT_T_TOP
//#define LB_FLUID_INITIAL_ADD_NOISE
#define LB_FLUID_INITIAL_UNIT_DENSITY
//#define LB_FLUID_AFTER_INIT_PERTURBATION
#define LB_FLUID_FORCING  /* activate force on the fluid */
#define LB_FLUID_FORCING_GRAVITY /* request the gravity parameter vector (without actually using it) */
//#define LB_FLUID_FORCING_POISEUILLE  /* Constant forcing for flow along x direction */
//#define LB_FLUID_FORCING_CHANNEL     /* Constant forcing for turbulent Channnel flow along x */
//#define LB_FLUID_FORCING_CHANNEL_CONSTANT_POWER   /* Constant power forcing for turbulent Channnel flow along x */  
//#define LB_FLUID_FORCING_KOLMOGOROV
//#define LB_FLUID_FORCING_SHEAR_LINEAR
//#define LB_FLUID_FORCING_LANDSCAPE
//#define LB_FLUID_FORCING_PENALIZATION
//#define LB_FLUID_FORCING_PENALIZATION_DIRECTION_X
//#define LB_FLUID_FORCING_DIRECT
//#define LB_FLUID_FORCING_CELLULAR
//#define LB_FLUID_FORCING_CELLULAR_UNSTEADY
//#define LB_FLUID_FORCING_STIRRER
#define LB_FLUID_FORCING_HIT
//#define LB_FLUID_FORCING_HIT_LINEAR
//#define LB_FLUID_FORCING_HIT_ZEROMODE
//#define LB_FLUID_FORCING_HIT_RANDOM
//#define LB_FLUID_FORCING_HIT_TYPE2
//#define LB_FLUID_FORCING_HIT_RECTANGULAR
//#define LB_FLUID_FORCING_HIT_2D
//#define LB_FLUID_FORCING_ABSORB
#define LB_FLUID_FORCING_NOZEROMODE
#define LB_FLUID_FORCING_CONSTANT_POWER
//#define LB_FLUID_FORCING_LAPLACIAN
//#define LB_FLUID_FORCING_CORIOLIS
 
/* Landscape definitions in the domain */
//#define LB_FLUID_FORCING_LANDSCAPE
//#define LB_FLUID_FORCING_LANDSCAPE_CYLINDER
//#define LB_FLUID_FORCING_LANDSCAPE_CUBE
//#define LB_FLUID_FORCING_LANDSCAPE_BUILDINGS

//#define LB_FLUID_BC
//#define LB_FLUID_BC_Y  /* default fluid bc is no-slip */
//#define LB_FLUID_BC_Y_M_SLIP                                           
//#define LB_FLUID_BC_Y_P_SLIP  
//#define LB_FLUID_BC_Y_M_JET                       
//#define LB_FLUID_BC_Y_P_OUTLET
//#define LB_FLUID_BC_Y_M_VELOCITY
//#define LB_FLUID_BC_Y_P_VELOCITY
//#define LB_FLUID_BC_Y_P_GRADIENT

//#define LB_FLUID_BC_X                                                     
//#define LB_FLUID_BC_X_M_SLIP                                             
//#define LB_FLUID_BC_X_P_SLIP                                               
//#define LB_FLUID_BC_X_M_INLET                                         
//#define LB_FLUID_BC_X_M_INLET_POISEUILLE
//#define LB_FLUID_BC_X_M_INLET_POISEUILLE_HALF
//#define LB_FLUID_BC_X_M_INLET_CONSTANT
//#define LB_FLUID_BC_X_M_OUTLET 
//#define LB_FLUID_BC_X_P_OUTLET  
                                             
//#define LB_FLUID_BC_Z
//#define LB_FLUID_BC_Z_M_SLIP                                             
//#define LB_FLUID_BC_Z_P_SLIP  
//#define LB_FLUID_BC_Z_M_OUTLET 
//#define LB_FLUID_BC_Z_P_OUTLET  

/* METHODS for time stepping or convective term */
/* For smooth simulations 
Activate either METHOD_FINITE_VOLUME or METHOD_STREAMING */
//#define METHOD_FINITE_VOLUME
#define METHOD_STREAMING
//#define METHOD_TRT
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
//#define METHOD_UPWIND_LINEAR
//#define METHOD_UPWIND_LINEAR_IMPROVED
//#define METHOD_MYQUICK
//#define METHOD_MYQUICK_CARTESIAN
//#define METHOD_MYQUICK_LIMITER
//#define METHOD_MIXED
//#define METHOD_TRAPEZOID
//#define METHOD_EDGES_AND_CORNERS
//#define METHOD_FORCING_GUO
//#define METHOD_FORCING_MALASPINAS
//#define METHOD_CORRECTION_LATT

/* TEMPERATURE */
//#define LB_TEMPERATURE
//#define LB_TEMPERATURE_PAST
//#define LB_TEMPERATURE_FLUCTUATION
//#define LB_TEMPERATURE_INITIAL_LINEAR
//#define LB_TEMPERATURE_INITIAL_ADD_PERTURBATION
//#define LB_TEMPERATURE_INITIAL_ADD_NOISE
//#define LB_TEMPERATURE_INITIAL_CONSTANT
//#define LB_TEMPERATURE_INITIAL_CONSTANT_MEAN
//#define LB_TEMPERATURE_INITIAL_CONSTANT_BOT
//#define LB_TEMPERATURE_INITIAL_CONSTANT_TOP
//#define LB_TEMPERATURE_INITIAL_SPOT
//#define LB_TEMPERATURE_INITIAL_BL
//#define LB_TEMPERATURE_INITIAL_BULK
//#define LB_TEMPERATURE_AFTER_INIT_PERTURBATION 
//#define LB_TEMPERATURE_BUOYANCY
//#define LB_TEMPERATURE_BUOYANCY_T0_REF
//#define LB_TEMPERATURE_BUOYANCY_T0_REF2  
//#define LB_TEMPERATURE_BUOYANCY_T0_BOT
//#define LB_TEMPERATURE_BUOYANCY_T0_TOP
//#define LB_TEMPERATURE_BUOYANCY_T0_GRAD
//#define LB_TEMPERATURE_BUOYANCY_WATER
//#define LB_TEMPERATURE_BC
//#define LB_TEMPERATURE_BC_Y  /* the default bc is fixed temperature value */
//#define LB_TEMPERATURE_BC_Y_P_OUTLET
//#define LB_TEMPERATURE_BC_Y_P_NOFLUX
//#define LB_TEMPERATURE_BC_Y_P_FLUX
//#define LB_TEMPERATURE_BC_Y_P_FLUX_CONSTANTGLOBALTEMPERATURE
//#define LB_TEMPERATURE_BC_Y_P_VARIABLE
//#define LB_TEMPERATURE_BC_Y_M_OUTLET
//#define LB_TEMPERATURE_BC_Y_M_NOFLUX
//#define LB_TEMPERATURE_BC_Y_M_FLUX 
//#define LB_TEMPERATURE_BC_Y_M_VARIABLE

//#define LB_TEMPERATURE_BC_X
//#define LB_TEMPERATURE_BC_X_P_OUTLET
//#define LB_TEMPERATURE_BC_X_M_OUTLET
//#define LB_TEMPERATURE_BC_X_P_NOFLUX
//#define LB_TEMPERATURE_BC_X_M_NOFLUX

//#define LB_TEMPERATURE_BC_Z
//#define LB_TEMPERATURE_BC_Z_P_OUTLET
//#define LB_TEMPERATURE_BC_Z_M_OUTLET
//#define LB_TEMPERATURE_BC_Z_P_NOFLUX
//#define LB_TEMPERATURE_BC_Z_M_NOFLUX
//#define LB_TEMPERATURE_BC_X_NOFLUX /* obsolete : only for finite-volume */
//#define LB_TEMPERATURE_BC_KEEP_WITHIN
//#define LB_TEMPERATURE_KEEP_CONSTANT /* use with care */
//#define LB_TEMPERATURE_FORCING
//#define LB_TEMPERATURE_FORCING_PAST
//#define LB_TEMPERATURE_FORCING_BULK 
//#define LB_TEMPERATURE_FORCING_BULK_VARIABLE
//#define LB_TEMPERATURE_FORCING_RADIATION 
//#define LB_TEMPERATURE_FORCING_RADIATION_SOLAR
//#define LB_TEMPERATURE_FORCING_RADIATION_REFLECTION
//#define LB_TEMPERATURE_FORCING_SOURCE
//#define LB_TEMPERATURE_FORCING_ABSORB
//#define LB_TEMPERATURE_FORCING_PROFILE
//#define LB_TEMPERATURE_FORCING_REACTION
//#define LB_TEMPERATURE_FORCING_REACTION_FKPP
//#define LB_TEMPERATURE_FORCING_REACTION_FKPP_FLUCTUATION
//#define LB_TEMPERATURE_FORCING_REACTION_ORDER1
//#define LB_TEMPERATURE_FORCING_MONOD 
//#define LB_TEMPERATURE_FORCING_HIT
//#define LB_TEMPERATURE_FORCING_HIT_RANDOM
//#define LB_TEMPERATURE_FORCING_HIT_ZEROMODE
//#define LB_TEMPERATURE_FORCING_HIT_LINEAR
//#define LB_TEMPERATURE_FORCING_VISCOUS
//#define LB_TEMPERATURE_FORCING_VISCOUS_FLUCTUATION
//#define LB_FLUID_TEMPERATURE_NOZEROMODE
//#define LB_TEMPERATURE_FORCING_CONSTANT_POWER
//#define LB_TEMPERATURE_FORCING_GRAD
//#define LB_TEMPERATURE_FORCING_LAPLACIAN

//#define LB_TEMPERATURE_MELTING
//#define LB_TEMPERATURE_MELTING_BOUNCEBACK
//#define LB_TEMPERATURE_MELTING_INITIAL_LIQUID
//#define LB_TEMPERATURE_MELTING_INITIAL_LIQUID_HALF
//#define LB_TEMPERATURE_MELTING_INITIAL_LIQUID_SEMISPHERE
//#define LB_TEMPERATURE_MELTING_INITIAL_LIQUID_SEMISPHERE_TEMPERATURE 
//#define LB_TEMPERATURE_MELTING_INITIAL_LIQUID_CAVITY
//#define LB_TEMPERATURE_MELTING_INITIAL_LIQUID_CAVITY_TEMPERATURE
//#define LB_TEMPERATURE_MELTING_INITIAL_LIQUID_LAYER
//#define LB_TEMPERATURE_MELTING_SOLID_DIFFUSIVITY
//#define LB_TEMPERATURE_MELTING_SOLIDUS
//#define LB_TEMPERATURE_MELTING_SOLIDUS_LINEAR
//#define LB_TEMPERATURE_MELTING_UNDEFORMABLE
//#define LB_TEMPERATURE_MELTING_CHECK_REACH_YP
//#define LB_TEMPERATURE_FORCING_DIRECT
//#define LB_FLUID_TEMPERATURE_NOZEROMODE
//#define LB_TEMPERATURE_FORCING_CONSTANT_POWER
//#define LB_TEMPERATURE_OUTPUT_BABAK

/* ADDED BY ZIQI for WATER-ICE MELTING*/
//#define LB_TEMPERATURE_ZIQI
//#define LB_TEMPERATURE_CHT_ZIQI 
//#define LB_TEMPERATURE_BUOYANCY_WATER_ZIQI 
//#define LB_TEMPERATURE_BUOYANCY_WATER_NOMEANBUOYANCY_ZIQI
//#define LB_TEMPERATURE_NEW_ENTHALPY_ZIQI
//#define LB_TEMPERATURE_INITIAL_CONSTANT_HALFBOT_ZIQI //initial condition
//#define LB_TEMPERATURE_INITIAL_LINEAR_ZERO_TO_TBOT_ZIQI //initial condition


/* EXTRA SCALAR FIELD e.g. CH4 */
//#define LB_SCALAR
//#define LB_SCALAR_INITIAL_LINEAR
//#define LB_SCALAR_INITIAL_ADD_PERTURBATION
//#define LB_SCALAR_INITIAL_BULK
//#define LB_SCALAR_INITIAL_CONSTANT
//#define LB_SCALAR_INITIAL_CONSTANT_MEAN
//#define LB_SCALAR_INITIAL_CONSTANT_BOT
//#define LB_SCALAR_INITIAL_CONSTANT_TOP
//#define LB_SCALAR_INITIAL_SPOT
//#define LB_SCALAR_BUOYANCY
//#define LB_SCALAR_BC
//#define LB_SCALAR_BC_Y
//#define LB_SCALAR_BC_Y_P_OUTLET
//#define LB_SCALAR_BC_Y_M_OUTLET
//#define LB_SCALAR_BC_Y_P_NOFLUX
//#define LB_SCALAR_BC_Y_M_NOFLUX
//#define LB_SCALAR_BC_X
//#define LB_SCALAR_BC_X_P_OUTLET
//#define LB_SCALAR_BC_X_M_OUTLET
//#define LB_SCALAR_BC_X_P_NOFLUX
//#define LB_SCALAR_BC_X_M_NOFLUX
//#define LB_SCALAR_BC_Z
//#define LB_SCALAR_BC_Z_P_OUTLET
//#define LB_SCALAR_BC_Z_M_OUTLET
//#define LB_SCALAR_BC_Z_P_NOFLUX
//#define LB_SCALAR_BC_Z_M_NOFLUX
//#define LB_SCALAR_FLUCTUATION  
//#define LB_SCALAR_FORCING
//#define LB_SCALAR_FORCING_REACTION
//#define LB_SCALAR_FORCING_REACTION_FKPP
//#define LB_SCALAR_FORCING_REACTION_FKPP_FLUCTUATION
//#define LB_SCALAR_FORCING_REACTION_ORDER1
//#define LB_SCALAR_FORCING_MONOD
//#define LB_SCALAR_FORCING_HIT   
//#define LB_SCALAR_FORCING_HIT_RANDOM
//#define LB_SCALAR_FORCING_HIT_ZEROMODE
//#define LB_SCALAR_FORCING_HIT_LINEAR  
//#define LB_SCALAR_FORCING_GRAD
/* Huisman model */
//#define LB_SCALAR_HUISMAN
//#define LB_SCALAR_SETTLING_HUISMAN
//#define LB_SCALAR_FORCING_HUISMAN
//#define LB_SCALAR_BC_Y_P_HUISMAN
//#define LB_SCALAR_BC_Y_M_HUISMAN


/* Lagragian parts */
#define LAGRANGE
//#define LAGRANGE_DEBUG
//#define LAGRANGE_WRAP
//#define LAGRANGE_OUTPUT_DEBUG
//#define LB_LAGRANGE_INITIAL_VELOCITY_FLUID
//#define LAGRANGE_NUCLEATE
//#define LAGRANGE_NUCLEATE_REMELT
//#define LAGRANGE_INITIAL_PAIRS
#define LAGRANGE_ADD_TRACER
//#define LAGRANGE_ADD_TRACER_STATIC
#define LAGRANGE_RADIUSandDENSITY
//#define LAGRANGE_RADIUSandDENSITY_INCREMENT_LOG_RADIUS
//#define LAGRANGE_RADIUSandDENSITY_INCREMENT_LOG_DENSITY
//#define LAGRANGE_TAUDRAG_INCREMENT_LOG
#define LAGRANGE_GRADIENT 
//#define LAGRANGE_GRADIENT_BRD2
#define LAGRANGE_ADDEDMASS
//#define LAGRANGE_ADDEDMASS_INCREMENT_LOG
//#define LAGRANGE_ADDEDMASS_LIFT
//#define LAGRANGE_ADDEDMASS_WAKEDRAG
#define LAGRANGE_GRAVITY
#define LAGRANGE_GRAVITY_VARIABLE
//#define LAGRANGE_GRAVITY_VARIABLE_INCREMENT_LOG
//#define LB_LAGRANGE_BC_INELASTIC   /* if not enabled the particle-wall collisions are elastic */
//#define LAGRANGE_SMALLTAUD /* approximate solution for tau_drag < 10.0 */
//#define LAGRANGE_ORIENTATION
//#define LAGRANGE_ORIENTATION_SECONDORIENTATION
//#define LAGRANGE_ORIENTATION_JEFFREY
//#define LAGRANGE_ORIENTATION_JEFFREY_INCREMENT_LOG
//#define LAGRANGE_ORIENTATION_BC
//#define LAGRANGE_ORIENTATION_JEFFREY_GYROTAXIS
//#define LAGRANGE_ORIENTATION_JEFFREY_GYROTAXIS_LINFENG /* Alternative and correct implementation used in our TAML 2021 paper */
//#define LAGRANGE_ORIENTATION_DIFFUSION
//#define LAGRANGE_ORIENTATION_RANDOM
//#define LAGRANGE_ORIENTATION_ACTIVE
//#define LAGRANGE_ORIENTATION_ACTIVE_JUMP
//#define LAGRANGE_ORIENTATION_ACTIVE_JUMP_HIGHPASS  // the default is LOWPASS
//#define LAGRANGE_ORIENTATION_ACTIVE_JUMP_TEMPERATURE
//#define LAGRANGE_ORIENTATION_ACTIVE_BALLISTIC
//#define LAGRANGE_ORIENTATION_ACTIVE_TEMPERATURE
//#define LAGRANGE_ORIENTATION_ACTIVE_SECONDORIENTATION
//#define LAGRANGE_ORIENTATION_ACTIVE_SELECTIVEORIENTATION
//#define LAGRANGE_ORIENTATION_DRAG
//#define LAGRANGE_POLYMER
//#define LAGRANGE_POLYMER_FEEDBACK

/* LES */
//#define LB_FLUID_LES
//#define LB_FLUID_LES_SISM
//#define LB_FLUID_LES_SISM_KALMAN
//#define LB_FLUID_LES_VANDRIEST
//#define LB_TEMPERATURE_LES
//#define LB_SCALAR_LES


/***************************************************/
/* Define dependencies not to be modified          */

#ifdef METHOD_FINITE_VOLUME                                            
#define METHOD_STEPPING_EULER
#define METHOD_REDEFINED_POP
#define METHOD_REDEFINED_POP_GUO
#define METHOD_FORCING_GUO
#define METHOD_HEUN
#define METHOD_MYQUICK                                                       
#define METHOD_MYQUICK_CARTESIAN
#define OUTPUT_H5_GRID
#endif 

#ifdef METHOD_STREAMING
#define METHOD_EDGES_AND_CORNERS                                        
#define METHOD_FORCING_GUO                                              
#endif

#ifdef LB_TEMPERATURE_BUOYANCY
#ifndef LB_FLUID_FORCING_GRAVITY
#define LB_FLUID_FORCING_GRAVITY
#endif
#endif