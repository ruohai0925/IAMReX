.. IAMReX documentation master file, created by
   sphinx-quickstart on Sun Apr 27 13:57:52 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

IAMReX documentation
====================

This IAMReX repo extends the capability of original `IAMR <https://amrex-fluids.github.io/IAMR/>`_ codes, aiming at simulating the multiphase incompressible flows and fluid structure interaction problems on both CPUs and GPUs with/without subcycling. The Navier-Stokes equations are solved on an adaptive semi-staggered grid using the projection method. The gas-liquid interface is captured using the level set (LS) method. The fluid-solid interface is resolved using the multidirect forcing immersed boundary method (IBM). The particle-wall as well as the particle-particle collisions are also captured by the adaptive collision time model (ACTM).

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   Introduction_Chapter
   Getting_Started
   Results
   Acknowledgements

.. toctree::
   :caption: References:

   references
