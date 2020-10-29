! $Id$
!
! MODULE_DOC: reads in full snapshot and calculates power spetrum of u
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lpower_spectrum = .false.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
!***************************************************************
!
module power_spectrum
!
  use Cdata
  use General, only: keep_compiler_quiet
!
  implicit none
!
  integer :: n_spectra=0
!
  include 'power_spectrum.h'
!
  contains
!***********************************************************************
    subroutine initialize_power_spectrum
!
    endsubroutine initialize_power_spectrum
!***********************************************************************
    subroutine read_power_spectrum_run_pars(iostat)
!
      integer, intent(out) :: iostat
!
      iostat = 0
!
    endsubroutine read_power_spectrum_run_pars
!***********************************************************************
    subroutine write_power_spectrum_run_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_power_spectrum_run_pars
!***********************************************************************
    subroutine power(f,sp,iapn_index)
!
      use General, only: ioptest
!
      real, dimension (mx,my,mz,mfarray) :: f
      character (len=1) :: sp
      integer, optional :: iapn_index
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(sp)
      call keep_compiler_quiet(ioptest(iapn_index))
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
    subroutine powerLor(f,sp)
!
      real, dimension (mx,my,mz,mfarray) :: f
      character (len=3) :: sp
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(sp)
!
    endsubroutine powerLor
!***********************************************************************
    subroutine powerEMF(f,sp)
!
      real, dimension (mx,my,mz,mfarray) :: f
      character (len=3) :: sp
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(sp)
!
    endsubroutine powerEMF
!***********************************************************************
    subroutine powerTra(f,sp)
!
      real, dimension (mx,my,mz,mfarray) :: f
      character (len=3) :: sp
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(sp)
!
    endsubroutine powerTra
!***********************************************************************
    subroutine powerGWs(f,sp,lfirstcall)
!
      real, dimension (mx,my,mz,mfarray) :: f
      character (len=3) :: sp
      logical :: lfirstcall
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(sp)
      call keep_compiler_quiet(lfirstcall)
!
    endsubroutine powerGWs
!***********************************************************************
    subroutine powerscl(f,sp,iapn_index,lsqrt)
!
      use General, only: ioptest

      logical, intent(in), optional :: lsqrt
      real, dimension (mx,my,mz,mfarray) :: f
      character (len=2) :: sp
      integer, optional :: iapn_index
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(sp)
      call keep_compiler_quiet(lsqrt)
      call keep_compiler_quiet(ioptest(iapn_index))
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
  subroutine anisoq_diag(f)
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
  endsubroutine anisoq_diag
!***********************************************************************
  subroutine corrfunc_3d(f)
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
  endsubroutine corrfunc_3d
!***********************************************************************
  subroutine k_omega_spectra(f,sp)
!
      real, dimension (mx,my,mz,mfarray) :: f
      character (len=*) :: sp
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(sp)
!
  endsubroutine k_omega_spectra
!***********************************************************************
endmodule power_spectrum
