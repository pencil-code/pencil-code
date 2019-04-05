
 Gravitational Waves from switching on a Beltrami field
=======================================================

## Maintainer:

Axel Brandenburg <brandenb/nordita[dot]org>

## Added:

11-Jun-2018

## Status:

succeeds

## Recommended resolution:

16x16x16 for nu=eta=0 is fine for short times.
For longer runs, one can put nu=eta=1e-3.
To see several oscillations, one should for 400 steps (0.02*400=8).
5-7 us/step/pt on one processor nl6 (laptop), depending on output.

## Comments:

* The gravitational waves are positively circularly polarized.
  Therefore, hel_GWs=GWs and hel_GWh=GWh.
* For EEM=1/2, the amplitudes are EEGW=2pi/4=pi/2; see (28) of Ref.[1]
  and hrms=8pi/2=4pi; see (26) of Ref.[1].
  The spectra denote spec_GWs = S_hdot(k,t), spec_GWh = S_h(k,t),
  and both have amplitude (4pi)^2.
* The time step is fixed to be 0.02 for accuracy reasons (i.e., c*dt/dx=0.051).
  Based on stability alone, the timestep could be 10 times longer [1].
* Alternatively, can use SPECIAL=special/gravitational_waves_hTXk,
  which solves the GW equatin exactly.

## Links:
* https://www.nordita.org/~brandenb/projects/GW

## Reference:
[1] Roper Pol, A., Brandenburg, A., Kahniashvili, T., Kosowsky, A.,
    Mandal, S.: 2018, ``The timestep constraint in solving the
    gravitational wave equations sourced by hydromagnetic turbulence,''
    Geophys. Astrophys. Fluid Dyn., submitted (arXiv:1807.05479v1)

