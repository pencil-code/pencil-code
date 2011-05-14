! $Id$
module power_spectrum
!
  use Cdata
  use Sub, only: keep_compiler_quiet
!
  implicit none
!
  integer :: n_spectra=0

  include 'power_spectrum.h'
!
  contains
!***********************************************************************
    subroutine read_power_spectrum_runpars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      call keep_compiler_quiet(unit)
      call keep_compiler_quiet(present(iostat))
!
    endsubroutine read_power_spectrum_runpars
!***********************************************************************
    subroutine write_power_spectrum_runpars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_power_spectrum_runpars

!***********************************************************************
    subroutine power(f,sp)
!
      real, dimension (mx,my,mz,mfarray) :: f
      character (len=1) :: sp
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(sp)
!
    endsubroutine power
!***********************************************************************
    subroutine power_2d(f,sp)
!
      real, dimension (mx,my,mz,mfarray) :: f
      character (len=1) :: sp
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(sp)
!
    endsubroutine power_2d
!***********************************************************************
    subroutine power_xy(f,sp,sp2)
!
      real, dimension (mx,my,mz,mfarray) :: f
      character (len=1) :: sp
      character (len=1), optional :: sp2
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(sp)
      call keep_compiler_quiet(present(sp2))
!
    endsubroutine power_xy
!***********************************************************************
    subroutine powerhel(f,sp)
!
      real, dimension (mx,my,mz,mfarray) :: f
      character (len=3) :: sp
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(sp)
!
    endsubroutine powerhel
!***********************************************************************
    subroutine powerscl(f,sp)
!
      real, dimension (mx,my,mz,mfarray) :: f
      character (len=2) :: sp
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(sp)
!
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
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(sp)
      call keep_compiler_quiet(ivec)
      call keep_compiler_quiet(present(ivar))
!
    endsubroutine power_1d
!***********************************************************************
    subroutine pdf(f,variabl,pdf_mean,pdf_rms)
!
      real, dimension (mx,my,mz,mfarray) :: f
      real :: pdf_mean,pdf_rms
      character (len=*) :: variabl
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(pdf_mean,pdf_rms)
      call keep_compiler_quiet(variabl)
!
    endsubroutine pdf
!***********************************************************************
    subroutine power_phi(f,sp)
!
      real, dimension (mx,my,mz,mfarray) :: f
      character (len=*) :: sp
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(sp)
!
    endsubroutine power_phi
!***********************************************************************
    subroutine powerhel_phi(f,sp)
!
      real, dimension (mx,my,mz,mfarray) :: f
      character (len=*) :: sp
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(sp)
!
    endsubroutine powerhel_phi
!***********************************************************************
  subroutine power_vec(f,sp)
!
      real, dimension (mx,my,mz,mfarray) :: f
      character (len=*) :: sp
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(sp)
!
  endsubroutine power_vec
!***********************************************************************
endmodule power_spectrum
