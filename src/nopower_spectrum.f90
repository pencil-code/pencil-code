! $Id$
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
  real, dimension (mx,my,mz,mfarray) :: f
  character (len=1) :: sp
  !
  if (ip<=15) print*,'Use POWER=power_spectrum in Makefile.local'
  if (sp=='') print*,f  !(to keep compiler happy)
  endsubroutine power
!***********************************************************************
    subroutine power_2d(f,sp)
!
  real, dimension (mx,my,mz,mfarray) :: f
  character (len=1) :: sp
  !
  if (ip<=15) print*,'Use POWER=power_spectrum in Makefile.local'
  if (sp=='') print*,f  !(to keep compiler happy)
  endsubroutine power_2d
!***********************************************************************
  subroutine powerhel(f,sp)
!
  real, dimension (mx,my,mz,mfarray) :: f
  character (len=3) :: sp
  !
  if (ip<=15) print*,'Use POWER=power_spectrum in Makefile.local'
  if (sp=='') print*,f  !(to keep compiler happy)
  endsubroutine powerhel
!***********************************************************************
  subroutine powerscl(f,sp)
!
  real, dimension (mx,my,mz,mfarray) :: f
  character (len=2) :: sp
  !
  if (ip<=15) print*,'Use POWER=power_spectrum in Makefile.local'
  if (sp=='') print*,f  !(to keep compiler happy)
  endsubroutine powerscl
!***********************************************************************
    subroutine power_1d(f,sp,ivec,ivar)
!
    real, dimension (mx,my,mz,mfarray) :: f
    character (len=1) :: sp
    integer :: ivec
    integer, optional :: ivar
!
    if (ip<=15) print*,'Use POWER=power_spectrum in Makefile.local'
!
    if (NO_WARN) print*, f, sp, ivec, ivar
!
  endsubroutine power_1d
!***********************************************************************
    subroutine pdf(f,variabl,pdf_mean,pdf_rms)
!
  real, dimension (mx,my,mz,mfarray) :: f
  real :: pdf_mean,pdf_rms
  character (len=*) :: variabl
!
  if (ip<0) print*,f,variabl,pdf_mean,pdf_rms !(keep compiler quiet)
  endsubroutine pdf
!***********************************************************************
  subroutine power_phi(f,sp)
!
! dummy 
!
    real, dimension (mx,my,mz,mfarray) :: f
    character (len=*) :: sp
!
    if (NO_WARN) print*,f,sp
!
  endsubroutine power_phi
!***********************************************************************
  subroutine powerhel_phi(f,sp)
!
! dummy 
!
    real, dimension (mx,my,mz,mfarray) :: f
    character (len=*) :: sp
!
    if (NO_WARN) print*,f,sp
!
  endsubroutine powerhel_phi
!***********************************************************************
endmodule power_spectrum
