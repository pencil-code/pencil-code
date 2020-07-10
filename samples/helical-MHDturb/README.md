
 Helical MHD turbulence
==========================

## Maintainer:

Axel Brandenburg <brandenb/nordita[dot]org>

## Added:

08-Jun-2002

## Status:

succeeds

## Recommended resolution:

32x32x32 for nu=eta=5e-3 is fine. For comparison with higher
resolution runs see Brandenburg (2001, ApJ 550, 824), except that
there a forcing wavenumber of kf=5 was used. The forcing function
works with preselected wavevectors that were computed with:
${PENCIL_HOME}/samples/helical-MHDturb/idl/generate_kvectors.pro

## Comments:

This is a helical MHD turbulence run. After about 600 time units
a large scale magnetic field develops from an initially random
magnetic field (if initaa=0 is set).

## Links:
* https://www.nordita.org/~brandenb/projects/LShelicityspec/
* http://pencil-code.nordita.org/samples/turbulence/helical-MHDturb32-4procs/

## References:

*  Brandenburg, A.: 2001 ``The inverse cascade and nonlinear alpha-effect in
   simulations of isotropic helical hydromagnetic turbulence,''
   *Astrophys. J.* **550**, 824-840 |
   [arXiv](http://arXiv.org/abs/astro-ph/0006186) |
   [ADS](http://esoads.eso.org/cgi-bin/nph-bib_query?bibcode=2001ApJ...550..824B)

*  Candelaresi, S., & Brandenburg, A.: 2013, ``Kinetic helicity needed to drive
   large-scale dynamos'' Phys. Rev. E 87, 043104 |
   [arXiv](https://arxiv.org/abs/1208.4529) |
   [ADS](http://adsabs.harvard.edu/abs/2013PhRvE..87d3104C)
