program test

  integer, parameter        :: mx=<N>, my=2, mz=2
  real, dimension(mx,my,mz) :: f=-10.
  real, dimension(mx)       :: fx=1.0

  write(0,*) 'spread_<N>: Starting'

  f(:,:,:) = f(:,:,:) + spread(spread(fx,2,my),3,mz)

  write(0,*) 'spread_<N>: <f> = ', sum(f)/mx/my/mz 
  write(0,*) 'spread_<N>: Done.' 
  write(0,*)

endprogram test
