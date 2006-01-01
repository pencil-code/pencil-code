! $Id: nopower_spectrum.f90,v 1.12 2006-01-01 15:42:39 ajohan Exp $
!
module  power_spectrum
  !
  use Cdata
  !
  implicit none

  include 'power_spectrum.h'
  !
  contains

!***********************************************************************
  subroutine power(f,sp)
!
  real, dimension (mx,my,mz,mvar+maux) :: f
  character (len=1) :: sp
  !
  if(ip<=15) print*,'Use POWER=power_spectrum in Makefile.local'
  if(sp=='') print*,f(1,1,1,1)  !(to keep compiler happy)
  endsubroutine power
!***********************************************************************
    subroutine power_2d(f,sp)
!
  real, dimension (mx,my,mz,mvar+maux) :: f
  character (len=1) :: sp
  !
  if(ip<=15) print*,'Use POWER=power_spectrum in Makefile.local'
  if(sp=='') print*,f(1,1,1,1)  !(to keep compiler happy)
  endsubroutine power_2d
!***********************************************************************
  subroutine powerhel(f,sp)
!
  real, dimension (mx,my,mz,mvar+maux) :: f
  character (len=3) :: sp
  !
  if(ip<=15) print*,'Use POWER=power_spectrum in Makefile.local'
  if(sp=='') print*,f(1,1,1,1)  !(to keep compiler happy)
  endsubroutine powerhel
!***********************************************************************
  subroutine powerscl(f,sp)
!
  real, dimension (mx,my,mz,mvar+maux) :: f
  character (len=2) :: sp
  !
  if(ip<=15) print*,'Use POWER=power_spectrum in Makefile.local'
  if(sp=='') print*,f(1,1,1,1)  !(to keep compiler happy)
  endsubroutine powerscl
!***********************************************************************
    subroutine power_1d(f,sp,ivec,ivar)
!
    real, dimension (mx,my,mz,mvar+maux) :: f
    character (len=1) :: sp
    integer :: ivec
    integer, optional :: ivar
!
    if(ip<=15) print*,'Use POWER=power_spectrum in Makefile.local'
    if(sp=='') print*,f(1,1,1,1),ivec  !(to keep compiler happy)
!    
  endsubroutine power_1d
!***********************************************************************
    subroutine pdf(f,variabl,pdf_mean,pdf_rms)
!
  real, dimension (mx,my,mz,mvar+maux) :: f
  real :: pdf_mean,pdf_rms
  character (len=*) :: variabl
!
  if(ip<0) print*,f(1,1,1,1),variabl,pdf_mean,pdf_rms !(keep compiler quiet)
  endsubroutine pdf
!***********************************************************************

endmodule power_spectrum
