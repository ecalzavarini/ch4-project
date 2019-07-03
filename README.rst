=================
Eulerian-Lagrangian fluid dynamics platform based on the Lattice-Boltzmann method
=================

Introduction
============

A general purpose *Lattice-Boltzmann* code for fluid-dynamics simulations. It includes : 

- **fluid dynamics**  (with several volume forcing terms for Channel flow, Homogeneous Isotropic Turbulence, buoyancy)
- **temperature dynamics** (advection, diffusion , sink/source or reaction terms)
- **phase change** (enthalpy formulation for solid/liquid systems)
- **scalar transport** (same functionalities as temperature)
- **lagrangian dynamics** (tracers, heavy/light & active  point-like particles; non-shperical Jeffery rotation, gyrotaxis)
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
This project is based at Unité de Mécanique de Lille (UML EA 7512, http://uml.univ-lille.fr ) France. 

Info: 

For more information please contact:

Enrico Calzavarini <enrico.calzavarini@polytech-lille.fr> , www.ecalzavarini.info

*Contributors*: Kalyan Shrestha (Lille University) , Babak Rabbanipour Esfahani (Lille University)


How to: 
======
See wiki pages https://github.com/ecalzavarini/ch4-project/wiki (very incomplete)

Aknowledgments:
======
This project received support from the INNOCOLD consortium (innocold.fr) and by the French National Agency for Research (ANR) by the grant (SEAS: ANR-13-JS09-0010).

Bibliography:
======
This code can be cited as:

0) *Eulerian-Lagrangian fluid dynamics platform: The ch4-project* E.Calzavarini, Software Impacts (2019).
   https://doi.org/10.1016/j.simpa.2019.100002

This code has been employed in the following published studies:

1) *Finite volume versus streaming-based lattice Boltzmann algorithm for fluid-dynamics simulations: A one-to-one accuracy and performance study*, K.Shrestha, G.Mompean and E.Calzavarini, Phys. Rev. E **93**, 023306 (2016).
   https://link.aps.org/pdf/10.1103/PhysRevE.93.023306

2) *Micro-bubbles and micro-particles are not faithful tracers of turbulent acceleration*, Varghese Mathai, Enrico Calzavarini,  Jon Brons, Chao Sun and Detlef Lohse, Phys. Rev. Lett. **117**, 024501 (2016).
   https://link.aps.org/doi/10.1103/PhysRevLett.117.024501
   
3) *Propelled microprobes in turbulence*, Enrico Calzavarini, Yongxiang X. Huang, Francois G. Schmitt and Lipo Wang, Phys. Rev. Fluids **3**, 054604 (2018).
   https://link.aps.org/doi/10.1103/PhysRevFluids.3.054604

4) *Basal melting driven by turbulent thermal convection*, Babak Rabbanipour Esfahani, Silvia C. Hirata, Stefano Berti and Enrico Calzavarini, Phys. Rev. Fluids **3**, 053501 (2018).
   https://link.aps.org/doi/10.1103/PhysRevFluids.3.053501

5) *Robustness of heat-transfer in confined inclined convection at high-Prandtl number*, Linfeng Jiang, Chao Sun and Enrico Calzavarini, Phys. Rev. E **99**, 013108 (2019). 
   https://link.aps.org/doi/10.1103/PhysRevE.99.013108
