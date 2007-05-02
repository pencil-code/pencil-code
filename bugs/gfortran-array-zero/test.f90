!!!!!!!!!!!!!!!!!!!!
!!!   test.f90   !!!
!!!!!!!!!!!!!!!!!!!!

!!!  Author: wd (Wolfgang.Dobler <at> ucalgary.ca)
!!!  Date:   28-Apr-2007
!!!
!!!  Description:
!!!   Demonstrates a compiler bug in gfortran 4.3.0 20070426 (experimental):
!!!   Program will print
!!!       1.000000       1.000000       1.000000       1.000000    
!!!       0.000000       0.000000       0.000000       0.000000    
!!!   instead of
!!!       1.000000       1.000000       1.000000       1.000000    
!!!       0.000000       1.000000       0.000000       1.000000    
!!!
!!!   Fixed in revision 124343 of gcc.

program Test

  implicit none
  real, dimension (2,2) :: smooth_factor

  smooth_factor = 1.
  write(0,*) smooth_factor

  ! Either of the following lines goes wrong:
  smooth_factor(1,:) = 0.
  !smooth_factor(:,1) = 0.
  write(0,*) smooth_factor

endprogram

!!! End of file test.f90
