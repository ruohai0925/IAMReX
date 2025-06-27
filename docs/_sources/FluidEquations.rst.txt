Fluid Variables
===============

The fluid variables are defined as follows.

   +-----------------------+--------------------------------------------------+
   | Variable              | Definition                                       |
   +=======================+==================================================+
   | :math:`\rho`          | Fluid density                                    |
   +-----------------------+--------------------------------------------------+
   | :math:`U`             | Fluid velocity                                   |
   +-----------------------+--------------------------------------------------+
   | :math:`\tau`          | Viscous stress tensor                            |
   +-----------------------+--------------------------------------------------+
   | :math:`{\bf H}_U`     | :math:`= (H_x , H_y , H_z )`, External Forces    |
   +-----------------------+--------------------------------------------------+
   | :math:`H_s`           | External sources                                 |
   +-----------------------+--------------------------------------------------+
   | :math:`\phi`          | Level set field                                  |
   +-----------------------+--------------------------------------------------+
   | :math:`\alpha_p`      | Particle volume fraction                         |
   +-----------------------+--------------------------------------------------+
   | :math:`\mathbf{F}_p`  | Eulerian force                                   |
   +-----------------------+--------------------------------------------------+

Compared with `IAMR <https://amrex-fluids.github.io/IAMR/>`_, IAMReX adds a level set field :math:`\phi` for the two-phase flow, and 
the particle volume fraction :math:`\alpha_p`, and the Eulerian force :math:`\mathbf{F}_p` for the particle-fluid interaction.

.. _FluidEquationsPart:

Fluid Equations
===============

Conservation of fluid mass:

.. math::
   :label: eq:div

   \frac{\partial \rho}{\partial t} + \nabla \cdot (\rho U)  = 0

Conservation of fluid momentum:

.. math::
   :label: eq:ns

   \frac{ \partial (\rho U)}{\partial t}
   + \nabla \cdot (\rho U U) + \nabla p = \nabla \cdot \tau + {\bf H}_U

Velocity constraint:

.. math:: \nabla \cdot U = S

where :math:`S` is zero by default, to model incompressible flow.
The :math:`S \ne 0` case is discussed below.

Tracer(s):

.. math:: \frac{\partial \rho s}{\partial t} + \nabla \cdot (\rho s U)  = \nabla \cdot \beta \nabla s + \rho H_s

for conservatively advected scalars and

.. math:: \frac{\partial s}{\partial t} + U \cdot \nabla s  = \nabla \cdot \beta \nabla s + H_s

for passively advected scalars. In general, one could advect an arbitrary number of scalars.

Level set equation:

.. math::
   :label: eq:phi

   \frac{\partial \phi}{\partial t} + U \cdot \nabla \phi = 0

IAMReX has the ability to incorporate general, user-defined external forcing and source terms. The default behaviour is that
:math:`H_s=0`, and :math:`{\bf H}_U` represents gravitational forces, with :math:`{\bf H}_U= (0 , 0 , -\rho g )` in 3d and
:math:`{\bf H}_U= (0 , -\rho g )` in 2d, where :math:`g` is the magnitude of the gravitational acceleration. However, since
by default, :math:`g=0`, :math:`{\bf H}_U = 0` unless ``ns.gravity`` is set.

By default, IAMReX solves the momentum equation in convective form. The inputs parameter ``ns.do_mom_diff = 1`` is used to
switch to conservation form. Tracers are passively advected by default. The inputs parameter ``ns.do_cons_trac = 1``
switches the first tracer to conservative. A second tracer can be included with ``ns.do_trac2 = 1``, and it can be
conservatively advected with ``ns.do_cons_trac2 = 1``.

IAMReX also has the option to solve for temperature, along with a modified divergence constraint on the velocity field:

.. math:: \rho c_p \left( \frac{\partial T}{\partial t} + U \cdot \nabla T \right)  = \nabla \cdot \lambda \nabla T + H_T

      \nabla \cdot U = \frac{1}{\rho c_p T} \nabla \cdot \lambda \nabla T

Here, the divergence constraint captures compressibily effects due to thermal diffusion.
To enable the temperature solve, use ``ns.do_temp = 1`` and set ``ns.temp_cond_coef`` to represent :math:`\lambda / c_p`,
which is taken to be constant. More sophiticated treatments are possible; if interested, please open an issue on github:
https://github.com/ruohai0925/IAMReX/issues


Time Step - Godunov
===================

In IAMReX, the canonical projection is applied to solve the fluid equations. Note that we only use the time-centered Godunov advection :cite:`almgren1998conservative,zeng2022aparallel`, and there no longer needs the predictor and corrector steps.

-  Define the time-centered face-centered (staggered) MAC velocity which is used for advection: :math:`U^{MAC,n+1/2}`

-  Define the new-time density, :math:`\rho^{n+1} = \rho^n - \Delta t (\rho^{n+1/2,pred} U^{MAC,n+1/2})` by setting

-  Define an approximation to the new-time state, :math:`(\rho U)^{\ast}` by setting

   .. math:: (\rho^{n+1} U^{\ast}) &= (\rho^n U^n) -
             \Delta t \nabla \cdot (\rho U^{MAC} U) + \Delta t \nabla {p}^{n-1/2}  \\ &+
             \frac{\Delta t}{2}  (\nabla \cdot \tau^n + \nabla \cdot \tau^\ast) +
             \Delta t \rho g

   (for implicit diffusion, which is the current default)

-  Project :math:`U^{\ast}` by solving

.. math:: \nabla \cdot \frac{1}{\rho} \nabla \phi = \nabla \cdot \left( \frac{1}{\Delta t}
          U^{\ast}+ \frac{1}{\rho} \nabla {p}^{n-1/2} \right)

then defining

.. math:: U^{n+1} = U^{\ast} - \frac{\Delta t}{\rho} \nabla \phi

and

.. math:: {p}^{n+1/2} = \phi

