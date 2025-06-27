.. _LevelSetMethod:

Level Set Method
================

The level set method is a powerful tool for capturing the interface between two immiscible fluids. It is based on the idea of a signed distance function, which is a function that assigns a signed distance to each point in the domain. The level set function is defined as follows:

.. math::

    \phi(\mathbf{x}, t) = 0

where :math:`\phi` is the level set function, :math:`\mathbf{x}` is the position vector, and :math:`t` is the time.

Material Properties
-------------------

The density and viscosity are defined as functions of the level set field:

.. math::
   :label: eq:rho

   \rho(\phi) = \rho_1 H(\phi) + \rho_2 (1 - H(\phi))

.. math::
   :label: eq:mu

   \mu(\phi) = \mu_1 H(\phi) + \mu_2 (1 - H(\phi))

where :math:`H(\phi)` is the Heaviside function:

.. math::

   H(\phi) = \begin{cases}
   1, & \phi > 0 \\
   0, & \phi \leq 0
   \end{cases}

Time Discretization
-------------------

For a single level, the momentum equation :eq:`eq:ns` is advanced by a fractional step method with the approximate projection :cite:`almgren1998conservative, zeng2022aparallel` to enforce the incompressibility condition (equation :eq:`eq:div`). The LS advection equation :eq:`eq:phi` is updated using the Godunov scheme.

At the beginning of each time advancement of level :math:`l`, the velocity :math:`\boldsymbol{u}^{n,l}` and the LS function :math:`\phi^{n,l}` at time :math:`t^{n,l}` are given. Owing to the fractional step method we used, the pressure is staggered in time and thus the pressure :math:`p^{n-1/2,l}` at time :math:`t^{n-1/2,l}` is known. Because the single-level advancement concerns only level :math:`l`, we omit the superscript :math:`{}^l` in this section. To obtain the updated velocity :math:`\boldsymbol{u}^{n+1}`, pressure :math:`p^{n+1/2}`, and LS function :math:`\phi^{n+1}` on level :math:`l`, the solver performs the following steps:

1. **Advance the LS function** as

   .. math::
      :label: eq:s0phin1

      \phi^{n+1} = \phi^{n} - \Delta t \left[\nabla \cdot (\boldsymbol{u}\phi) \right]^{n+1/2}

   The advection term in the above equation is calculated using the Godunov scheme.

2. **Solve the intermediate velocity** :math:`\boldsymbol{u}^{*}` semi-implicitly

   .. math::
      :label: eq:viscsolve

      \begin{aligned}
      &\boldsymbol{u^{*,n+1}}-\frac{\Delta t}{2\rho(\phi^{n+1/2})Re}{\nabla} \cdot {\mu(\phi^{n+1})}{\nabla}\boldsymbol{u^{*,n+1}} =
      \boldsymbol{u^n}-\Delta t \left[\nabla \cdot (\boldsymbol{uu})\right]^{n+1/2} + \\ &\frac{\Delta t}{\rho(\phi^{n+1/2})}\bigg[-\nabla p^{n-1/2}+
      \frac{1}{2Re}{\nabla}\cdot{\mu(\phi^{n})}{\nabla}\boldsymbol{u^n} + \rho(\phi^{n+1/2}) \frac{z}{Fr^2} -  \frac{1}{We}\kappa(\phi^{n+1/2})\delta(x^{n+1/2})\boldsymbol{n}\bigg].
      \end{aligned}

   In equation :eq:`eq:viscsolve`, the detailed discretization of the advection term :math:`\nabla \cdot (\boldsymbol{uu})`, viscous term :math:`{\nabla}\cdot({\mu(\phi)}{\nabla}\boldsymbol{u})`, and surface tension term :math:`\kappa(\phi)\delta(x)\boldsymbol{n}/We` are implemented using standard finite difference schemes. The LS function at :math:`t^{n+1/2}` is calculated by

   .. math::
      :label: eq:ns_half_phi

      \phi^{n+1/2} = \frac{1}{2}(\phi^{n}+\phi^{n+1})

   where :math:`\phi^{n+1}` is obtained from step 1 (equation :eq:`eq:s0phin1`). The :math:`\rho(\phi^{n+1/2})`, :math:`\mu(\phi^{n})`, and :math:`\mu(\phi^{n+1})` are then obtained from equations :eq:`eq:rho` and :eq:`eq:mu`.

3. **Apply the projection method** to obtain the pressure and a solenoidal velocity field. To conduct the level projection, a temporary variable :math:`\boldsymbol{V}` is defined as

   .. math::
      :label: eq:ns_lp_ls1

      \boldsymbol{V} =  \frac{\boldsymbol{{u^{*,n+1}}}}{\Delta t} + \frac{1}{\rho(\phi^{n+1/2})} \nabla p^{n-1/2}

   Then the updated pressure :math:`p^{n+1/2}` is calculated by

   .. math::
      :label: eq:ns_lp_ls2

      L^{cc,\mathrm{level}}_{\rho^{n+1/2}} p^{n+1/2} =  \nabla \cdot \boldsymbol{V}

   where :math:`L^{cc,\mathrm{level}}_{\rho^{n+1/2}}p^{n+1/2}` is a density-weighted approximation to :math:`\nabla \cdot (1/\rho^{n+1/2} \nabla p^{n+1/2})`. Finally, the velocity can be calculated as

   .. math::
      :label: eq:ns_lp_ls3

      \boldsymbol{{u^{n+1}}} = \Delta t \left(\boldsymbol{V} - \frac{1}{\rho^{n+1/2}} \nabla p^{n+1/2}\right)

   As defined in the AMReX framework, :math:`\nabla \cdot` and :math:`\nabla` are the cell-centered level divergence operator :math:`D^{cc,\mathrm{level}}` and level gradient operator :math:`G^{cc,\mathrm{level}}`, respectively. The level gradient operator :math:`G^{cc,\mathrm{level}}` is not the minus transpose of the level divergence operator :math:`D^{cc,\mathrm{level}}`, i.e., :math:`G^{cc,\mathrm{level}} \neq -(D^{cc,\mathrm{level}})^T`. As a result, the idempotency of the approximate projection :math:`\boldsymbol{P} = I - G^{cc,\mathrm{level}}(L^{cc,\mathrm{level}})^{-1}D^{cc,\mathrm{level}}` is not ensured, i.e., :math:`\boldsymbol{P}^{2} \neq \boldsymbol{P}`. Yet, this nonidempotent approximate projection is stable and appears to be well-behaved in various numerical tests and practical applications :cite:`martin2000cell`. Notably, for a uniform single grid with periodic boundary conditions, Lai theoretically proved that this approximate projection method is stable, in that :math:`\|\boldsymbol{P}\| \leq 1`. It should be noted that the approximate projection is applied to the intermediate velocity :math:`\boldsymbol{{u^{*,n+1}}}` (equation :eq:`eq:ns_lp_ls1`). Compared with the form that projects the increment velocity :math:`\boldsymbol{u^{*,n+1}}-\boldsymbol{u^n}`, e.g. as that used in other methods, the projection method used here can reduce the accumulation of pressure errors and lead to a more stable algorithm. The effectiveness and stability of this approximate projection has been validated through various numerical tests.

4. **Reinitialize the LS function** :math:`\phi` to maintain :math:`\phi` as a signed distance function of the interface and guarantee the conservation of the mass of the two phases. In this step, a temporary LS function :math:`d(\boldsymbol{x},\tau)` is updated iteratively using the following pseudo evolution equation:

   .. math::
      :label: eq:ns_reinit1

      \frac{\partial d}{\partial \tau}=S(\phi)(1-|\nabla d|)

   with the initial condition

   .. math::
      :label: eq:ns_reinit2

      d(\boldsymbol{x},\tau = 0)=\phi^{n+1}(\boldsymbol{x})

   where

   .. math::
      :label: eq:ns_reinit3

      S(\phi) = 2\left(H(\phi)-1 / 2\right)

   Here, :math:`\tau` is the pseudo time for iterations. A second-order essentially non-oscillatory (ENO) scheme is used to discretize the distance function and a second-order Runge--Kutta (RK) method is applied for the pseudo time advancing. To ensure the mass conservation, :math:`d(\boldsymbol{x},\tau)` is further corrected by minimizing the differences of the volume of each fluid between :math:`\tau=0` and the final iteration. Finally, the LS function :math:`\phi` is re-initialized by the volume corrected :math:`d` :cite:`zeng2022subcycling,zeng2023consistent,sussman1999adaptive`.

At last, we give a summary of the single-level advancement algorithm as follows :cite:`zeng2022subcycling`.

1. Advance the LS function using equation :eq:`eq:s0phin1`
2. Solve the intermediate velocity using equation :eq:`eq:viscsolve`
3. Apply the projection method to update the pressure and velocity field following equations :eq:`eq:ns_lp_ls1`--:eq:`eq:ns_lp_ls3`
4. Re-initialize the LS function on the single level using equations :eq:`eq:ns_reinit1`--:eq:`eq:ns_reinit3`

