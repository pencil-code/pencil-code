
.. #########################################
.. Common Substitutions for Pencil Docs
.. #########################################


..   This file defines reusable text and math snippets (substitutions) used throughout
   the Pencil Code documentation.  
   Include it in pages using:

   .. code-block:: rst

      .. include:: /_substitutions.rst

   Or make it global in ``conf.py`` with:

   .. code-block:: python

      rst_prolog = """
      .. include:: /_substitutions.rst
      """

.. ==============================
.. Project and repository info
.. ==============================

.. |PC| replace:: PENCIL CODE
.. |Pencil| replace:: Pencil Code
.. |pencilweb| replace:: `<https://pencil-code.org/>`__
.. |repo| replace:: `GitHub repository <https://github.com/pencil-code/pencil-code>`__
.. |docs| replace:: `online documentation <https://pencil-code.readthedocs.io>`__
.. |manual| replace:: ``manual.tex``
.. |conf| replace:: ``conf.py``
.. |Docker| replace:: `Docker <https://www.docker.com/>`__
.. |Sphinx| replace:: `Sphinx <https://www.sphinx-doc.org/>`__
.. |RTD| replace:: `Read the Docs <https://readthedocs.org/>`__
.. |Github| replace:: `Github <https://github.com>`__

.. ==============================
.. Math and operators
.. ==============================

.. |grad| replace:: :math:`\boldsymbol{\nabla}`
.. |div| replace:: :math:`\boldsymbol{\nabla}\cdot`
.. |curl| replace:: :math:`\boldsymbol{\nabla}\times`
.. |pderiv| replace:: :math:`\frac{\partial}{\partial}`
.. |de| replace:: :math:`\mathrm{d}`
.. |Re| replace:: :math:`\mathrm{Re}`
.. |cs| replace:: :math:`c_{\mathrm{s}}`



.. ==============================
.. Text and references
.. ==============================

.. |seeconf| replace:: See :file:`conf.py` for more details.
.. |latexmanual| replace:: See the LaTeX source in :file:`manual.tex`.

.. ########################################################
.. Operators
.. ########################################################

.. |De|       replace:: :math:`\De`
.. |artanh|   replace:: :math:`\artanh`
.. |const|    replace:: :math:`\const`
.. |Div|      replace:: :math:`\Div`
.. |rot|      replace:: :math:`\rot`
.. |Laplace|  replace:: :math:`\Laplace`
.. |erfc|     replace:: :math:`\erfc`
.. |erf|      replace:: :math:`\erf` 
.. |pderivn|  replace:: :math:`\pderivn{u}{x}{2}`

.. ########################################################
.. Vectors
.. ########################################################

.. |Av|       replace:: :math:`\Av`
.. |Bv|       replace:: :math:`\Bv`
.. |Jv|       replace:: :math:`\Jv`
.. |Uv|       replace:: :math:`\Uv`
.. |Wv|       replace:: :math:`\Wv`
.. |Ev|       replace:: :math:`\Ev`
.. |Fv|       replace:: :math:`\Fv`
.. |fv|       replace:: :math:`\fv`
.. |gv|       replace:: :math:`\gv`
.. |jv|       replace:: :math:`\jv`
.. |kv|       replace:: :math:`\kv`
.. |ov|       replace:: :math:`\ov`
.. |uv|       replace:: :math:`\uv`
.. |vv|       replace:: :math:`\vv`
.. |bv|       replace:: :math:`\bv`
.. |xv|       replace:: :math:`\xv`
.. |zerovect| replace:: :math:`\zerovect`
.. |omv|      replace:: :math:`\omv`
.. |Bhat|     replace:: :math:`\Bhat`
.. |BBhat|    replace:: :math:`\BBhat`

.. ########################################################
.. Units and numbers
.. ########################################################

.. |ns|       replace:: :math:`\ns`
.. |ps|       replace:: :math:`\ps`
.. |EE|       replace:: :math:`\EE{n}`
.. |ttimes|   replace:: :math:`\ttimes{n}`
.. |xtimes|   replace:: :math:`\xtimes{a}{b}`

.. ########################################################
.. Physical numbers
.. ########################################################

.. |Ra|       replace:: :math:`\Ra`
.. |Reynolds| replace:: :math:`\Reynolds`
.. |Rm|       replace:: :math:`\Rm`
.. |vA|       replace:: :math:`\vA`
.. |csnull|   replace:: :math:`\csnull`
.. |Heat|     replace:: :math:`\Heat`
.. |Cool|     replace:: :math:`\Cool`
.. |Heavi|    replace:: :math:`\Heavi`
.. |Strain|   replace:: :math:`\Strain`

.. ===========================
.. Journals
.. ==================================

.. ArXiv categories
.. (used for linking preprints, leave these symbolic if not replaced by URLs)
.. |arXiv| replace:: arXiv preprint
.. |astroph| replace:: astro-ph preprint 
.. |condmat| replace:: cond-mat preprint
.. |physics| replace:: physics preprint
.. |qbio| replace:: q-bio preprint

.. Journals and series

.. |yannr| replace:: *Ann. Rev. Astron. Astrophys.*
.. |yica| replace:: *Icarus*
.. |ysph| replace:: *Solar Phys.*
.. |ysphs| replace:: *Solar Phys.*
.. |ymn| replace:: *Monthly Notices Roy. Astron. Soc.*
.. |yan| replace:: *Astron. Nachr.*
.. |yana| replace:: *Astron. Astrophys.*
.. |yanaN| replace:: *Astron. Astrophys.*
.. |yanas| replace:: *Astron. Astrophys.*
.. |yass| replace:: *Astrophys. Spa. Sci.*
.. |yapj| replace:: *Astrophys. J.*
.. |yapjl| replace:: *Astrophys. J. Lett.*
.. |yapjlS| replace:: *Astrophys. J. Lett.*
.. |yapjs| replace:: *Astrophys. J. Suppl. Series*
.. |yjfm| replace:: *J. Fluid Mech.*
.. |ypepi| replace:: *Phys. Earth Planet. Int.*
.. |ygafd| replace:: *Geophys. Astrophys. Fluid Dyn.*
.. |ypr| replace:: *Phys. Rev.*
.. |yprN| replace:: *Phys. Rev.*
.. |yjour| replace:: *Journal (generic)*
.. |yjourS| replace:: *Journal (generic)*
.. |ybook| replace:: *Book reference*
.. |ypf| replace:: *Phys. Fluids*
.. |ypp| replace:: *Phys. Plasmas*
.. |yepl| replace:: *Europhys. Lett.*
.. |yprl| replace:: *Phys. Rev. Lett.*
.. |ybif| replace:: *Int. J. Bifurc. Chaos*
.. |ycsf| replace:: *Chaos, Solitons & Fractals*
.. |ycsfS| replace:: *Chaos, Solitons & Fractals*
.. |ynat| replace:: *Nature*

.. Submitted / in press / to be submitted versions

.. Physics journals
.. |spf| replace:: *Phys. Fluids* (submitted)
.. |ppf| replace:: *Phys. Fluids* (in press)
.. |ppp| replace:: *Phys. Plasmas* (in press)
.. |spp| replace:: *Phys. Plasmas* (submitted)
.. |tpp| replace:: *Phys. Plasmas* (to be submitted)
.. |tppS| replace:: *Phys. Plasmas* (to be submitted)
.. |pppp| replace:: *Phys. Plasmas* (in press, scheduled)
.. |ppr| replace:: *Phys. Rev.* (in press)
.. |spr| replace:: *Phys. Rev.* (submitted)
.. |tpr| replace:: *Phys. Rev.* (to be submitted)
.. |sprl| replace:: *Phys. Rev. Lett.* (submitted)
.. |pprl| replace:: *Phys. Rev. Lett.* (in press)
.. |pbif| replace:: *Int. J. Bifurc. Chaos* (in press)
.. |sbif| replace:: *Int. J. Bifurc. Chaos* (submitted)

.. Astrophysics journals
.. |sapj| replace:: *Astrophys. J.* (submitted)
.. |sapjS| replace:: *Astrophys. J.* (submitted)
.. |ppapj| replace:: *Astrophys. J.* (in press)
.. |ppapjS| replace:: *Astrophys. J.* (in press)
.. |tapj| replace:: *Astrophys. J.* (to be submitted)
.. |sapjl| replace:: *Astrophys. J. Lett.* (submitted)
.. |sapjlS| replace:: *Astrophys. J. Lett.* (submitted)
.. |ppapjl| replace:: *Astrophys. J. Lett.* (in press)
.. |ppapjlS| replace:: *Astrophys. J. Lett.* (in press)
.. |ppapjs| replace:: *Astrophys. J. Suppl. Series* (in press)
.. |papj| replace:: *Astrophys. J.* (scheduled)
.. |papjS| replace:: *Astrophys. J.* (scheduled)
.. |papjl| replace:: *Astrophys. J. Lett.* (scheduled)
.. |papjlS| replace:: *Astrophys. J. Lett.* (scheduled)

.. Astronomy journals
.. |sana| replace:: *Astron. Astrophys.* (submitted)
.. |sanas| replace:: *Astron. Astrophys.* (submitted)
.. |pana| replace:: *Astron. Astrophys.* (in press)
.. |panas| replace:: *Astron. Astrophys.* (in press)
.. |tana| replace:: *Astron. Astrophys.* (to be submitted)
.. |san| replace:: *Astron. Nachr.* (submitted)
.. |pan| replace:: *Astron. Nachr.* (in press)

.. Geophysics journals
.. |sgafd| replace:: *Geophys. Astrophys. Fluid Dyn.* (submitted)
.. |pgafd| replace:: *Geophys. Astrophys. Fluid Dyn.* (in press)
.. |ppgafd| replace:: *Geophys. Astrophys. Fluid Dyn.*

.. Solar physics
.. |ssph| replace:: *Solar Phys.* (submitted)
.. |psph| replace:: *Solar Phys.* (in press)

.. Fluid mechanics
.. |sjfm| replace:: *J. Fluid Mech.* (submitted)
.. |pjfm| replace:: *J. Fluid Mech.* (in press)
.. |tjfm| replace:: *J. Fluid Mech.* (to be submitted)

.. Monthly Notices (MNRAS)
.. |smn| replace:: *Monthly Notices Roy. Astron. Soc.* (submitted)
.. |pmn| replace:: *Monthly Notices Roy. Astron. Soc.* (in press)
.. |tmn| replace:: *Monthly Notices Roy. Astron. Soc.* (to be submitted)

.. ==============================
.. Example usage
..  ==============================

.. Use any of these in documentation text:

.. .. code-block:: rst

..    The |PC| source is available in the |repo|.  
..    For derivations, see |manual| and note that |grad| and |div| follow standard notation.

.. Result:

..    The |PC| source is available in the |repo|.  
   For derivations, see |manual| and note that |grad| and |div| follow standard notation.
