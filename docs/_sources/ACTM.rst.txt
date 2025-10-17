
.. _ACTM:

Adaptive Collision Time Model
==============================

The Adaptive Collision Time Model (ACTM) is a sophisticated collision model designed for grain-resolving simulations of flows over dense, mobile, polydisperse granular sediment beds. This model combines the immersed boundary method (IBM) to couple fluid and particle dynamics, providing accurate representation of particle-particle and particle-wall interactions in multiphase flow simulations.

Introduction
------------

The ACTM addresses the challenges of simulating dense granular systems where particles are in close proximity and frequently interact through collisions. Traditional collision models often struggle with maintaining numerical stability and physical accuracy in such scenarios, particularly when dealing with:

- Dense particle packing with frequent collisions
- Polydisperse particle size distributions
- Mobile granular sediment beds
- Complex fluid-particle interactions

The ACTM employs an adaptive procedure to obtain the desired restitution coefficient and resolve collisions on the timescale of the fluid solver, ensuring accurate and stable simulations of particle-laden flows.

Model Components
----------------

The ACTM consists of three key components that work together to accurately model particle interactions:

Lubrication Force
~~~~~~~~~~~~~~~~~

When particles approach each other or solid surfaces, the lubrication force accounts for the fluid dynamics in the narrow gap between surfaces. This force prevents direct particle contact and reduces energy dissipation:

.. math::
   :label: eq:lubrication_force

   \mathbf{F}_{lub} = -\frac{6\pi\mu R_{eff}^2}{h} \mathbf{v}_{rel}

where :math:`\mu` is the fluid viscosity, :math:`R_{eff}` is the effective radius, :math:`h` is the gap distance, and :math:`\mathbf{v}_{rel}` is the relative velocity between particles.

Normal Contact Force
~~~~~~~~~~~~~~~~~~~~

The normal contact force prevents particle overlap and models the elastic and damping behavior during collisions using a spring-damper system:

.. math::
   :label: eq:normal_contact_force

   \mathbf{F}_n = k_n \delta_n \mathbf{n} + c_n \mathbf{v}_n

where :math:`k_n` is the normal stiffness, :math:`\delta_n` is the overlap distance, :math:`\mathbf{n}` is the contact normal vector, :math:`c_n` is the normal damping coefficient, and :math:`\mathbf{v}_n` is the normal relative velocity.

Tangential Contact Force
~~~~~~~~~~~~~~~~~~~~~~~~

The tangential contact force accounts for friction and describes the tangential interactions between particles, affecting their rotational and sliding behavior:

.. math::
   :label: eq:tangential_contact_force

   \mathbf{F}_t = \min(\mu_s |\mathbf{F}_n|, k_t \delta_t) \mathbf{t}

where :math:`\mu_s` is the static friction coefficient, :math:`k_t` is the tangential stiffness, :math:`\delta_t` is the tangential displacement, and :math:`\mathbf{t}` is the tangential unit vector.

Adaptive Time Stepping
----------------------

The ACTM employs an adaptive time stepping scheme to ensure numerical stability and physical accuracy. The collision time step is dynamically adjusted based on:

1. **Minimum collision time**: The time until the next collision event
2. **Fluid solver time step**: The time step used by the fluid solver
3. **Contact duration**: The time required to resolve contact interactions

The adaptive time step is calculated as:

.. math::
   :label: eq:adaptive_timestep

   \Delta t_{collision} = \min\left(\frac{t_c}{N_{substeps}}, \Delta t_{fluid}, \Delta t_{contact}\right)

where :math:`t_c` is the collision time, :math:`N_{substeps}` is the number of substeps, and :math:`\Delta t_{contact}` is the contact resolution time.

Enduring Contact Handling
-------------------------

One of the key innovations of the ACTM is its ability to handle enduring contacts, which are common in dense granular systems. The model introduces new criteria for persistent contact situations:

1. **Contact detection**: Identifies when particles remain in contact for extended periods
2. **Momentum conservation**: Maintains complete momentum balance during enduring contacts
3. **Parameter-free approach**: Avoids the use of arbitrary parameter values

The enduring contact force is calculated as:

.. math::
   :label: eq:enduring_contact

   \mathbf{F}_{enduring} = \mathbf{F}_{normal} + \mathbf{F}_{tangential} + \mathbf{F}_{lubrication}

This approach ensures that the model remains physically consistent even in complex contact scenarios.

Mathematical Formulation
------------------------

Collision Detection
~~~~~~~~~~~~~~~~~~~

The ACTM uses a spatial hashing approach to efficiently identify potential collision pairs. For each particle :math:`i` with position :math:`\mathbf{x}_i` and radius :math:`r_i`, the algorithm identifies particles :math:`j` within a search radius:

.. math::
   :label: eq:collision_detection

   |\mathbf{x}_i - \mathbf{x}_j| \leq R_{search} = r_i + r_j + \delta

where :math:`\delta` is a safety margin for robust collision detection.

Collision Time Calculation
~~~~~~~~~~~~~~~~~~~~~~~~~~

For each detected collision pair, the collision time :math:`t_c` is calculated using the relative velocity approach:

.. math::
   :label: eq:collision_time

   t_c = \frac{|\mathbf{x}_i - \mathbf{x}_j| - (r_i + r_j)}{|\mathbf{v}_i - \mathbf{v}_j| \cdot \mathbf{n}_{ij}}

where :math:`\mathbf{v}_i` and :math:`\mathbf{v}_j` are the particle velocities, and :math:`\mathbf{n}_{ij}` is the unit normal vector.

Restitution Coefficient
~~~~~~~~~~~~~~~~~~~~~~~

The restitution coefficient :math:`e` is determined based on the material properties and collision conditions:

.. math::
   :label: eq:restitution_coefficient

   e = \exp\left(-\frac{\pi \zeta}{\sqrt{1-\zeta^2}}\right)

where :math:`\zeta` is the damping ratio:

.. math::
   :label: eq:damping_ratio

   \zeta = \frac{c_n}{2\sqrt{k_n m_{eff}}}

Here, :math:`m_{eff}` is the effective mass:

.. math::
   :label: eq:effective_mass

   m_{eff} = \frac{m_i m_j}{m_i + m_j}

Algorithm Implementation
------------------------

The ACTM algorithm follows these key steps:

1. **Particle Insertion**: Particles are inserted into the collision detection system with their current properties.

2. **Spatial Hashing**: A background mesh is generated for efficient collision pair identification.

3. **Collision Pair Generation**: Potential collision pairs are identified using spatial hashing.

4. **Collision Time Calculation**: Collision times are calculated for each pair using equation :eq:`eq:collision_time`.

5. **Adaptive Time Stepping**: The collision time step is adaptively adjusted using equation :eq:`eq:adaptive_timestep`.

6. **Force Calculation**: Contact forces are calculated using equations :eq:`eq:normal_contact_force` and :eq:`eq:tangential_contact_force`.

7. **Enduring Contact Resolution**: Persistent contacts are handled using equation :eq:`eq:enduring_contact`.

8. **Integration with Fluid Solver**: Collision forces are integrated with the fluid solver for consistent timescales.


This integration ensures that particle collisions are resolved consistently with the fluid dynamics and interface evolution, providing a comprehensive simulation framework for multiphase flows with dense granular systems.

