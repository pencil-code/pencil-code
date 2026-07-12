Schwinger effect with MHD
=========================

## Maintainer:

Axel Brandenburg <brandenb/nordita[dot]org>

## Added:

18-Apr-2026

## Status:

succeeds

## Recommended resolution:

512x1x1, but here 64x1x1 is also ok

## Comments:

This sample run covers time ranges with displacement current included
at early times and without displacement current at late times.
The length of the time step is adjusted automatically.

Initial values when lappy_BD_k1D_factor=T and
amplphi=1e0, amplee_BD_prefactor=1e0, kpeak_phi=50e-6
 k1        erms           brms
4e-6  3.174147605E-16  3.23551E-16 
2e-6  1.604501256E-16  1.41110E-16
1e-6  5.748021227E-17  3.81043E-17

Same, but now instead lappy_BD_k1D_factor=F
 k1        erms           brms
4e-6  3.967684507E-08  4.04438E-08
2e-6  5.672768593E-08  4.98898E-08
1e-6  5.748021227E-08  3.81043E-08

3D, but now instead lappy_BD_k1D_factor=T, linv=T currently (NEW)
 k1        erms           brms
4e-6  7.199604133E-15   5.82799E-15

2e-6  7.184665696E-15   5.34148E-15
1e-6  4.897146140E-15   2.64046E-15

3D, but now instead lappy_BD_k1D_factor=T, linv=F currently (NEW **)
 k1        erms           brms
4e-6  2.746431020E-20   2.22320E-20
2e-6  2.740732459E-20   2.03761E-20
1e-6  1.868112999E-20   1.00726E-20

4e-6  5.364618736E-23   4.34369E-23  with 1/N for N=64^3 !nok
4e-6  6.577652321E-21   5.25253E-21  with 1/N for N=64^2 !nok
4e-6  6.199507042E-19   6.31935E-19  with 1/N for N=64   !nok
4e-6  1.212661095E-21   8.97046E-22  with 1/N for N=32^3 !nok

good scaling
4e-6  1.406172682E-17   1.13828E-17  with N for N=64^3
4e-6  2.694206391E-17   2.15144E-17  with N for N=64^2
      2.694206391E-17   2.15144E-17  F (i.e., use forward transform)
      1.103546938E-13   8.81229E-14  T (4e-6)
      7.603928674E-14   5.53677E-14  T (2e-6)
4e-6  3.967684507E-17   4.04438E-17  with N for N=64
      3.967684507E-17   4.04438E-17
      2.539318084E-15   2.58841E-15  T (4e-6)
      1.283601005E-15   1.12888E-15  T (2e-6)
4e-6  3.973647876E-17   2.93944E-17  with N for N=32^3
      3.973647876E-17   2.93944E-17

Same, but now instead lappy_BD_k1D_factor=F currently (NEW)
 k1        erms           brms
4e-6  8.999505167E-07   7.28499E-07
2e-6  2.540162917E-06   1.88850E-06
1e-6  4.897146140E-06   2.64046E-06

Same, but now instead lapply_BD_1ND_factor=T
 k1        erms           brms
4e-6  7.935369014E-11  8.08877E-11
2e-6  8.022506281E-11  7.05548E-11
1e-6  5.748021227E-11  3.81043E-11

Same, but now instead lapply_BD_1ND_factor=T and 
 k1        erms           brms
4e-6  1.239901408E-12  1.26387E-12
2e-6  1.253516606E-12  1.10242E-12 
1e-6  8.981283167E-13  5.95380E-13

Now 3-D WRONG
Same, but now instead lapply_BD_1ND_factor=T and 
 k1        erms           brms
X 4e-6  2.746431020E-20  2.22320E-20
X 2e-6  2.740732459E-20  2.03761E-20 
X 1e-6  1.868112999E-20  1.00726E-20 

Same, but now instead with lapply_BD_ND_factor=T, linv_BD=T
 k1        erms           brms
X 4e-6  8.999505167E-07  7.28499E-07
X 2e-6  2.540162917E-06  1.88850E-06
X 1e-6  4.897146140E-06  2.64046E-06

Same, but now instead with lapply_BD_ND_factor=T, linv_BD=T
 k1        erms           brms
Y 4e-6  4.607746645E-04   3.72992E-04
Y 2e-6  1.300563414E-03   9.66912E-04
Y 1e-6  2.507338824E-03   1.35192E-03

Same, but now instead with lapply_BD_ND_factor=T, linv_BD=F
 k1        erms           brms
Y 4e-6  1.757715853E-09   1.42285E-09
Y 2e-6  4.961255698E-09   3.68848E-09
Y 1e-6  9.564738555E-09   5.15715E-09  (clipped)

updated to match settings in data/oksana/plasma/2D/Hbdn1024alpf90_rho28_ampl1_Gam9o8
updated to minimize result from: pc_diffruns . ~/data/oksana/plasma/2D/Hbdn1024alpf90_rho28_ampl1_Gam9o8 | m

in preparation

## Links:
* http://pencil-code.nordita.org/samples/1d-tests/Schwinger-MHD

## References:

*  Iarygina, O., Sfakianakis, E. I., & Brandenburg, A.: 2025, "Schwinger effect in axion inflation on a lattice," Phys. Rev. Lett., submitted |
   [arXiv](https://arxiv.org/abs/2506.20538) |
   [ADS](http://adsabs.harvard.edu/abs/2025arXiv250620538I)
