Relativistic forced turbulence
=======================================================
## Added:

10-September-2026 Antonino S. Midiri (toninomidiri7@virgilio.it)

## Comments:

* this sample simulates hydrodynamic forced turbulence in the relativistic regime for a radiation dominated fluid (with eos = 1/3)
* max_vel (default value 0.999) controls the maximum allowed value for the velocity field and prevents superluminal velocities from appearing
* lvel_limiter controls whether this maximum value is imposed or not
* this velocity limiter has not yet been extended to the conservation form of relativistic hydrodynamics
