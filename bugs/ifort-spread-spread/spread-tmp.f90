program test

  integer, parameter :: mx=256, my=2, mz=2
  real, dimension(mx,my,mz) :: f
  real, dimension(mx,my) :: tmp
  real, dimension(mx) :: fx=1.0

  write(0,*) 'spread-tmp: Starting'

  tmp = spread(fx,2,my)
  f(:,:,:) = f(:,:,:) + spread(tmp,3,mz)

  write(0,*) 'spread-tmp: <f> = ', sum(f)/mx/my/mz 
  write(0,*) 'spread-tmp: Done.' 
  write(0,*)

endprogram test
