program test

  integer, parameter :: mx=255, my=2, mz=2
  real,dimension(mx,my,mz) :: f
  real,dimension(mx) :: fx=1.0

  write(0,*) 'spread_255: Starting'

  f(:,:,:) = f(:,:,:) + spread(spread(fx,2,my),3,mz)

  write(0,*) 'spread_255: <f> = ', sum(f)/mx/my/mz 
  write(0,*) 'spread_255: Done.' 
  write(0,*)


endprogram test

