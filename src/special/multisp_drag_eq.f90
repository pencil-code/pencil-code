! $Id$
!
!  This module provide a way for users to specify custom
!  (i.e. not in the standard Pencil Code) physics, diagnostics etc.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of special_dummies.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lspecial = .true.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
!***************************************************************
module Special
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages, only: svn_id, fatal_error
!
  implicit none
!
  include "../special.h"
!
! Module Variables
!
  real, dimension(npar_species) :: vpx0 = 0.0, vpy0 = 0.0
  real :: ux0 = 0.0, uy0 = 0.0
!
! Diagnostic variables (needs to be consistent with reset list below).
!
  integer :: idiag_rduxm = 0
  integer :: idiag_rduym = 0
  integer :: idiag_rdux2m = 0
  integer :: idiag_rduy2m = 0
  integer :: idiag_rduxduym = 0
  integer :: idiag_ruzduxm = 0
  integer :: idiag_ruzduym = 0
!
  contains
!****************************************************************************
    subroutine initialize_special(f)
!
!  Called after reading parameters, but before the time loop.
!
!  04-oct-19/ccyang: coded
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
!
      call keep_compiler_quiet(f)
      if (lstart) return
      print *, "initialize_special: hello world! "
!
!  Read the equilibrium velocities from the initial conditions.
!
      open(10, file="data/multisp_drag_eq.dat", form="unformatted", action="read")
      read(10) ux0, uy0, vpx0, vpy0
      close(10)
!
!  Log the velocities.
!
      info: if (lroot) then
        print *, "initialize_special: ux0, uy0 = ", ux0, uy0
        print *, "initialize_special: vpx0 = ", vpx0
        print *, "initialize_special: vpy0 = ", vpy0
      endif info
!
    endsubroutine initialize_special
!***********************************************************************
    subroutine pencil_criteria_special
!
!  All pencils that this special module depends on are specified here.
!
!  04-oct-19/ccyang: coded
!
      diag: if (idiag_rduxm /= 0 .or. idiag_rduym /= 0 .or. &
                idiag_rdux2m /= 0 .or. idiag_rduy2m /= 0 .or. &
                idiag_rduxduym /= 0 .or. idiag_ruzduxm /= 0 .or. idiag_ruzduym /= 0) then
        lpenc_diagnos(i_rho) = .true.
        lpenc_diagnos(i_uu) = .true.
      endif diag
!
    endsubroutine pencil_criteria_special
!***********************************************************************
    subroutine dspecial_dt(f, df, p)
!
!  calculate right hand side of ONE OR MORE extra coupled PDEs
!  along the 'current' Pencil, i.e. f(l1:l2,m,n) where
!  m,n are global variables looped over in equ.f90
!
!  Due to the multi-step Runge Kutta timestepping used one MUST always
!  add to the present contents of the df array.  NEVER reset it to zero.
!
!  Several precalculated Pencils of information are passed for
!  efficiency.
!
!  04-oct-19/ccyang: coded
!
      use Diagnostics, only: sum_mn_name
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(mx,my,mz,mvar), intent(in) :: df
      type(pencil_case), intent(in) :: p
!
      real, dimension(nx) :: dux, duy
!
      call keep_compiler_quiet(f, df)
      call keep_compiler_quiet(p)
!
!  Identify module and boundary conditions.
!
      if (headtt.or.ldebug) print *, "dspecial_dt: special diagnostics"
!
!  Diagnostics
!
      diag: if (ldiagnos) then
        dux = p%uu(:,1) - ux0
        duy = p%uu(:,2) - uy0
        if (idiag_rduxm /= 0) call sum_mn_name(p%rho * dux, idiag_rduxm)
        if (idiag_rduym /= 0) call sum_mn_name(p%rho * duy, idiag_rduym)
        if (idiag_rdux2m /= 0) call sum_mn_name(p%rho * dux**2, idiag_rdux2m)
        if (idiag_rduy2m /= 0) call sum_mn_name(p%rho * duy**2, idiag_rduy2m)
        if (idiag_rduxduym /= 0) call sum_mn_name(p%rho * dux * duy, idiag_rduxduym)
        if (idiag_ruzduxm /= 0) call sum_mn_name(p%rho * p%uu(:,3) * dux, idiag_ruzduxm)
        if (idiag_ruzduym /= 0) call sum_mn_name(p%rho * p%uu(:,3) * duy, idiag_ruzduym)
      endif diag
!
    endsubroutine dspecial_dt
!***********************************************************************
    subroutine rprint_special(lreset, lwrite)
!
!  Reads and registers print parameters relevant to special.
!
!  04-oct-19/ccyang: coded
!
      use Diagnostics, only: parse_name
      use FArrayManager, only: farray_index_append
!
      logical, intent(in) :: lreset
      logical, intent(in), optional :: lwrite
!
      logical :: lwr
      integer :: iname
!
      lwr = .false.
      if (present(lwrite)) lwr = lwrite
!
!  Reset everything in case of reset
!
      reset: if (lreset) then
        idiag_rduxm = 0
        idiag_rduym = 0
        idiag_rdux2m = 0
        idiag_rduy2m = 0
        idiag_rduxduym = 0
        idiag_ruzduxm = 0
        idiag_ruzduym = 0
      endif reset
!
!  Turn on requested diagnostics.
!
      diag: do iname = 1, nname
        call parse_name(iname, cname(iname), cform(iname), "rduxm", idiag_rduxm)
        call parse_name(iname, cname(iname), cform(iname), "rduym", idiag_rduym)
        call parse_name(iname, cname(iname), cform(iname), "rdux2m", idiag_rdux2m)
        call parse_name(iname, cname(iname), cform(iname), "rduy2m", idiag_rduy2m)
        call parse_name(iname, cname(iname), cform(iname), "rduxduym", idiag_rduxduym)
        call parse_name(iname, cname(iname), cform(iname), "ruzduxm", idiag_ruzduxm)
        call parse_name(iname, cname(iname), cform(iname), "ruzduym", idiag_ruzduym)
      enddo diag
!
!  Write column where which magnetic variable is stored
!
      wr: if (lwr) then
        call farray_index_append("idiag_rduxm", idiag_rduxm)
        call farray_index_append("idiag_rduym", idiag_rduym)
        call farray_index_append("idiag_rdux2m", idiag_rdux2m)
        call farray_index_append("idiag_rduy2m", idiag_rduy2m)
        call farray_index_append("idiag_rduxduym", idiag_rduxduym)
        call farray_index_append("idiag_ruzduxm", idiag_ruzduxm)
        call farray_index_append("idiag_ruzduym", idiag_ruzduym)
      endif wr
!
    endsubroutine rprint_special
!***********************************************************************
!***********************************************************************
!
!***********************************************************************
!***********************************************************************
!********************************************************************
!
!********************************************************************
!************        DO NOT DELETE THE FOLLOWING       **************
!********************************************************************
!**  This is an automatically generated include file that creates  **
!**  copies dummy routines from nospecial.f90 for any Special      **
!**  routines not implemented in this file                         **
!**                                                                **
    include '../special_dummies.inc'
!********************************************************************
endmodule Special
