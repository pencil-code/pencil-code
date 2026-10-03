! $Id$
!
!  A module that solves a single evolution equation, df/dt = 1. This is used to debug the times at which various quantities are output.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lspecial = .true.
!
! MVAR CONTRIBUTION 1
! MAUX CONTRIBUTION 0
!
!***************************************************************
module Special
!
  use Quiet
  use Cdata
!
  implicit none
!
  include '../special.h'
!
  integer :: ispecial=0
!
  integer :: idiag_specialm=0, idiag_specialmz=0, idiag_specialmxy=0
!
  contains
!***********************************************************************
    subroutine register_special
!
      use FArrayManager, only: farray_register_pde
!
      call farray_register_pde('special',ispecial)
!
    endsubroutine register_special
!***********************************************************************
    subroutine init_special(f)
!
!
      real, contiguous, dimension(:,:,:,:) :: f
!
      f(:,:,:,ispecial) = 0.
!
    endsubroutine init_special
!***********************************************************************
    subroutine dspecial_dt(f,df,p)
!
      real, contiguous, dimension(:,:,:,:) :: f
      real, dimension(mx,my,mz,mvar) :: df
      type(pencil_case) :: p
!
      df(l1:l2,m,n,ispecial) = df(l1:l2,m,n,ispecial) + 1.
!
      call calc_diagnostics_special(f,p)
!
    endsubroutine dspecial_dt
!***********************************************************************
    subroutine rprint_special(lreset,lwrite)
!
      use Diagnostics, only: parse_name
!
      logical :: lreset
      logical, optional :: lwrite
!
     integer :: iname
!
      call keep_compiler_quiet(lwrite)
!
      if (lreset) then
          idiag_specialm=0
          idiag_specialmz=0
          idiag_specialmxy=0
      endif
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),&
            'specialm',idiag_specialm)
      enddo
!
      do iname=1,nnamez
        call parse_name(iname,cnamez(iname),cformz(iname),&
            'specialmz',idiag_specialmz)
      enddo
!
      do iname=1,nnamexy
        call parse_name(iname,cnamexy(iname),cformxy(iname),&
            'specialmxy',idiag_specialmxy)
      enddo
!
      if (lwrite_slices) then
        where(cnamev=='special') cformv='DEFINED'
      endif
!
      !TODO: farray_index_append?
!
    endsubroutine rprint_special
!***********************************************************************
    subroutine get_slices_special(f,slices)
!
      use Slices_methods, only: assign_slices_scal
!
      real, contiguous, dimension(:,:,:,:) :: f
      type(slice_data) :: slices
!
!     NOTE: copied from advective_gauge.f90
      select case (trim(slices%name))
      case ('special')
        call assign_slices_scal(slices,f,ispecial)
      endselect
!
    endsubroutine get_slices_special
!***********************************************************************
    subroutine calc_diagnostics_special(f,p)
!
      use Diagnostics
!
      real, contiguous, dimension(:,:,:,:) :: f
      type(pencil_case) :: p
!
      call keep_compiler_quiet(p)
!
      if (ldiagnos) then
        call sum_mn_name(f(l1:l2,m,n,ispecial), idiag_specialm)
      endif
!
      if (l1davgfirst) then
        call xysum_mn_name_z(f(l1:l2,m,n,ispecial), idiag_specialmz)
      endif
!
      if (l2davgfirst) then
        call zsum_mn_name_xy(f(l1:l2,m,n,ispecial), idiag_specialmxy)
      endif
!
    endsubroutine calc_diagnostics_special
!***********************************************************************
    subroutine pushpars2c(p_par)
      use Syscalls, only: copy_addr
!
      integer, parameter :: n_pars=1
      integer(KIND=ikind8), dimension(n_pars) :: p_par
!
      call copy_addr(ispecial, p_par(1)) ! int
!
    endsubroutine pushpars2c
!***********************************************************************
!TODO: may be useful to also test power spectrum output as well; below are just dummy subroutines.
    subroutine special_calc_spectra(f,spectrum,spectrumhel, &
      spectrum_2d,spectrum_2d_hel,&
      lfirstcall,kind)
!
      real, contiguous, dimension(:,:,:,:) :: f
      real, dimension(:) :: spectrum
      real, dimension(:) :: spectrumhel
      real, dimension(:,:) :: spectrum_2d
      real, dimension(:,:) :: spectrum_2d_hel
      logical :: lfirstcall
      character(len=3) :: kind
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(spectrum,spectrumhel)
      call keep_compiler_quiet(spectrum_2d,spectrum_2d_hel)
      call keep_compiler_quiet(lfirstcall)
      call keep_compiler_quiet(kind)
!
    endsubroutine special_calc_spectra
!***********************************************************************
    subroutine special_calc_spectra_byte(f,spectrum,spectrumhel,lfirstcall,kind,len)

      real, contiguous, dimension(:,:,:,:) :: f
      real, dimension(:) :: spectrum
      real, dimension(:) :: spectrumhel
      logical :: lfirstcall
      character, dimension(3) :: kind
      integer :: len
!
      call keep_compiler_quiet(len)
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(spectrum,spectrumhel)
      call keep_compiler_quiet(lfirstcall)
      call keep_compiler_quiet(kind)
!
    endsubroutine special_calc_spectra_byte
!***********************************************************************
!***********************************************************************
!************        DO NOT DELETE THE FOLLOWING       **************
!********************************************************************
!**  This is an automatically generated include file that creates  **
!**  copies dummy routines from nospecial.f90 for any Special      **
!**  routines not implemented in this file                         **
!**                                                                **
    include '../special_dummies.inc'
!********************************************************************
endmodule Special
