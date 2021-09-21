##################################
The Pencil Code documentation
##################################

.. admonition:: Welcome!

   This is the new homepage of The Pencil Code documentation!
   
   Explore the page hierarchy below (or in the sidebar), and get started with
   :ref:`contributing your own documentation<Contributing to the documentation>`!

The Pencil Code is primarily designed to deal with weakly compressible turbulent flows,  which is why we use high-order first and second derivatives. To achieve good parallelization, we use explicit (as opposed to compact) finite differences. Typical scientific targets include driven MHD turbulence in a periodic box, convection in a slab with non-periodic upper and lower boundaries, a convective star embedded in a fully nonperiodic box, accretion disc turbulence in the shearing sheet approximation, self-gravity, non-local radiation transfer, dust particle evolution with feedback on the gas, etc. A range of artificial viscosity and diffusion schemes can be invoked to deal with supersonic flows. For direct simulations regular viscosity and diffusion is being used.

Please find `more details on our website <http://pencil-code.nordita.org/>`_.


.. toctree::
   :caption: Introduction
   :maxdepth: 2

   intro/getting_started
   intro/usingrst
   intro/links
   intro/discussion

.. toctree::
   :caption: User manuals
   :maxdepth: 2

   Quick Guide <manuals/quick-guide>
   

.. toctree::
   :caption: Tutorials
   :maxdepth: 2

   tutorials/pencil/tutpencil
   tutorials/python/tutpython
   tutorials/mathematica/tutmathematica


.. toctree::
   :caption: Code documentation
   :maxdepth: 2

   toc/modpython
   toc/modidl
   toc/modfortran

  
 


