program test

  integer, parameter :: mx=256,my=2,mz=2
  real,dimension(mx,my,mz) :: f
  real,dimension(mx) :: fx=1.0

  f(:,:,:)=f(:,:,:)+spread(spread(fx,2,my),3,mz)

end program test

