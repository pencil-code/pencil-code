!!!!!!!!!!!!!!!!!!!!
!!!   seed.f90   !!!
!!!!!!!!!!!!!!!!!!!!

!!!  Author: wd (Wolfgang.Dobler@kis.uni-freiburg.de)
!!!  Date:   26-Aug-2004
!!!
!!!  Description:
!!!   Get random seed before and after getting a first random number

program rand_seed

  implicit none

  integer, parameter :: mseed=1024
  integer :: nseed
  integer, dimension(mseed) :: seed
  real, dimension(4) :: X

  call random_seed(SIZE=nseed)
  if (nseed>mseed) then
    print*, 'ERROR: nseed=', nseed, ' exceeds mseed=', mseed
    STOP
  endif

  call random_seed(GET=seed(1:nseed))
  print*, '---------------------------------------------'
  print*, 'Seed = '
  print*, seed(1:nseed)
  print*, '---------------------------------------------'

  call random_number(X)
  print*, 'X = ', X

  call random_seed(GET=seed(1:nseed))
  print*, '---------------------------------------------'
  print*, 'Seed = '
  print*, seed(1:nseed)
  print*, '---------------------------------------------'

  call random_number(X)
  print*, 'X = ', X

  call random_seed(GET=seed(1:nseed))
  print*, '---------------------------------------------'
  print*, 'Seed = '
  print*, seed(1:nseed)
  print*, '---------------------------------------------'


endprogram rand_seed

!!! End of file seed.f90
