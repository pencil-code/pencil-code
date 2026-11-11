SU(2) gauge fields coupled to a Higgs doublet
=============================================

## Maintainer:

Alberto Roper Pol <alberto[dot]roperpol/unige[dot]ch>

## Added:

10-Nov-2026

## Status:

succeeds

## Recommended resolution:

16x16x16

## Comments:

Real-time evolution of the bosonic SU(2) x Higgs sector (no hypercharge) in
temporal gauge, using special/electroweaksu2 together with special/klein_gordon
(lphi_doublet=T, lphi_weakcharge=T).  The nine W^a_i are initialised with
cosine waves of different wavevectors and phases, the Higgs doublet with
phi_up = A (1 + sin x), phi_down = A sin x, and all time derivatives vanish,
so the Gauss constraint is initially satisfied.  The longitudinal (Gamma)
variables are evolved with the default constraint damping
weight_longitudinalW=2, i.e. d(Gamma^a - div W^a)/dt = -2 G^a.

## Links:

## References:
