Key Differences Between IAMR and IAMeX
======================================

IAMR (Base)
-----------
- ``NUM_STATE_MAX = AMREX_SPACEDIM+4`` (velocity + density + 2 tracers + temperature)
- Standard incompressible Navier-Stokes

IAMReX (Extended)
-----------------
- ``NUM_STATE_MAX = AMREX_SPACEDIM+5`` (adds level set field)
- Additional advance methods:
  - ``advance_semistaggered_twophase_ls()``: Two-phase flow with level sets
  - ``advance_semistaggered_fsi_diffusedib()``: Fluid-structure interaction with particle collision

Comparison with Base IAMR
-------------------------

.. list-table:: Difference between IAMR and IAMReX
    :widths: 50 50 50
    :header-rows: 1

    * - Feature
      - IAMR
      - IAMReX
    * - State Variables
      - 4(velocity, density, tracers, temperature)
      - 5(adds level set)
    * - Interface Methods
      - None
      - Level set, immersed boundary, phase field (TBD)
    * - Particle Support
      - None
      - Full 6-DOF dynamics with collisions
    * - Multi-Physics
      - Single-phase only
      - Multi-phase, fluid-solid interaction
    * - Collision Detection
      - None
      - Spatial hashing with DKT model
    * - GPU Support
      - Basic AMReX
      - Enhanced particle operations

Software Architecture
=====================

Core Class Hierarchy
--------------------
- ``NavierStokesBase``: Base class providing AMR infrastructure and core CFD functionality
- ``NavierStokes``: Derived class implementing specific solver methods

Key Components
--------------
- **State Management**: Uses AMReX state descriptors for velocity, pressure, and scalars
- **Boundary Conditions**: Defined in NS_BC.H with physical-to-mathematical BC mapping
- **Time Integration**: Semi-implicit schemes with explicit advection and implicit diffusion
- **Projection Methods**: MAC projection and nodal projection for incompressible flow
- **Level Set Support**: Additional scalar fields for interface tracking of the multiphase flow
- **Fluid-structure interaction**: Implements diffused immersed boundary method for FSI
- **Multiple resolved particles**: Particle collision detection using spatial hashing

Key Files
---------
- ``NavierStokesBase.H/cpp``: Core AMR and CFD infrastructure
- ``NavierStokes.H/cpp``: Problem-specific solver implementation
- ``NS_setup.cpp``: Variable registration and boundary condition setup
- ``DiffusedIB.H/cpp``: Immersed boundary method implementation
- ``Collision.H/cpp``: Particle collision handling

Understanding the Flow
----------------------
1. Variable setup happens in ``NS_setup.cpp``
2. Time advancement occurs in ``NavierStokes::advance()`` methods
3. Boundary conditions are applied via functions in NS_BC.H
