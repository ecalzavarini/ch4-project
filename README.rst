=================
General purpose Lattice-Boltzmann code
=================

Introduction
============

A general purpose *Lattice-Boltzmann* code for fluid-dynamics simulations. It includes : 

- **fluid dynamics**  (with several volume forcing terms for Channel flow, Homogeneous Isotropic Turbulence, buoyancy)
- **temperature dynamics** (advection, diffusion , sink/source or reaction terms)
- **phase change** (enthalpy formulation for solid/liquid systems)
- **scalar transport** (same functionalities as temperature)
- **lagrangian dynamics** (tracers, heavy/light & active  point-like particles; non-shperical Jeffrey rotation, gyrotaxis)
- **large eddy simulation** (Smagorinsky, Shear Improved Samgorinsky with Kalman Filter)

Requirements:

- MPI 
- HDF5 
- CMake (optional)

History
=======

This project is a continuation and extension https://github.com/ecalzavarini/ice-project

Contact
=======
This project is based at Mechanics Laboratory of Lille (LML) France. 

Info: 

For more information please contact:

Enrico Calzavarini <enrico.calzavarini@polytech-lille.fr> , www.ecalzavarini.info

*Contributors*: Kalyan Shrestha (Lille University) , Babak Rabbanipour (Lille University)


How to: 
======
To be written

Aknowledgments:
======
Partial support from French National Agency for Research (ANR) by the grant (SEAS: ANR-13-JS09-0010)

Bibliography:
======
This code has been employed in the following studies:

1) *Finite volume versus streaming-based lattice Boltzmann algorithm for fluid-dynamics simulations: A one-to-one accuracy and performance study*, K.Shrestha, G.Mompean and E.Calzavarini, Phys. Rev. E **93**, 023306 (2016).

2) *Micro-bubbles and micro-particles are not faithful tracers of turbulent acceleration*, Varghese Mathai, Enrico Calzavarini,  Jon Brons, Chao Sun and Detlef Lohse, Phys. Rev. Lett. **117**, 024501 (2016).

3) *Basal melting driven by turbulent thermal convection*, Babak Rabbanipour Esfahani, Silvia C. Hirata, Stefano Berti and Enrico Calzavarini, arXiv:1801.03694 (2018).

4) *Propelled micro-probes in turbulence*, E. Calzavarini, Y.  X. Huang, F. G. Schmitt and L. Wang, arXiv:1802.00189 (2018).
