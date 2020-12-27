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
! Diagnostic variables
!
  integer :: idiag_rduxm = 0
  integer :: idiag_rduym = 0
  integer :: idiag_rdux2m = 0
  integer :: idiag_rduy2m = 0
  integer :: idiag_rduxduym = 0
  integer :: idiag_ruzduxm = 0
  integer :: idiag_ruzduym = 0
!
  integer :: idiag_drhopm = 0
  integer :: idiag_drhop2m = 0
  integer :: idiag_rhopdvpxm = 0
  integer :: idiag_rhopdvpym = 0
  integer :: idiag_rhopdvpx2m = 0
  integer :: idiag_rhopdvpy2m = 0
  integer :: idiag_rhopvpz2m = 0
!
! yz-averages
!
  integer :: idiag_rduxmx = 0
  integer :: idiag_rduymx = 0
  integer :: idiag_rdux2mx = 0
  integer :: idiag_rduy2mx = 0
  integer :: idiag_rduxduymx = 0
  integer :: idiag_ruzduxmx = 0
  integer :: idiag_ruzduymx = 0
!
  integer :: idiag_drhopmx = 0
  integer :: idiag_drhop2mx = 0
!
  contains
!****************************************************************************
    subroutine initialize_special(f)
!
! Called after reading parameters, but before the time loop.
!
! 17-jul-20/ccyang: coded
!
      use Mpicomm, only: mpibcast
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
!
      if (lstart) return
      call keep_compiler_quiet(f)
!
!  Read the equilibrium velocities from the initial conditions.
!
      eqvel: if (lroot) then
        open(10, file="data/multisp_drag_eq.dat", form="unformatted", action="read")
        read(10) ux0, uy0, vpx0, vpy0
        close(10)
!
        print *, "initialize_special: ux0, uy0 = ", ux0, uy0
        print *, "initialize_special: vpx0 = ", vpx0
        print *, "initialize_special: vpy0 = ", vpy0
      endif eqvel
!
      call mpibcast(ux0)
      call mpibcast(uy0)
      call mpibcast(vpx0, npar_species)
      call mpibcast(vpy0, npar_species)
!
    endsubroutine initialize_special
!***********************************************************************
    subroutine pencil_criteria_special
!
!  All pencils that this special module depends on are specified here.
!
!  03-aug-20/ccyang: coded
!
      gas: if (idiag_rduxm /= 0 .or. idiag_rduym /= 0 .or. &
               idiag_rdux2m /= 0 .or. idiag_rduy2m /= 0 .or. &
               idiag_rduxduym /= 0 .or. idiag_ruzduxm /= 0 .or. idiag_ruzduym /= 0 .or. &
               idiag_rduxmx /= 0 .or. idiag_rduymx /= 0 .or. &
               idiag_rdux2mx /= 0 .or. idiag_rduy2mx /= 0 .or. &
               idiag_rduxduymx /= 0 .or. idiag_ruzduxmx /= 0 .or. idiag_ruzduymx /= 0) then
        lpenc_diagnos(i_rho) = .true.
        lpenc_diagnos(i_uu) = .true.
      endif gas
!
      if (idiag_drhopm /= 0 .or. idiag_drhop2m /= 0 .or. &
          idiag_drhopmx /= 0 .or. idiag_drhop2mx /= 0) lpenc_diagnos(i_rhop) = .true.
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
!  03-aug-20/ccyang: coded
!
      use Diagnostics, only: sum_mn_name, yzsum_mn_name_x
      use EquationOfState, only: rho0
      use Particles_cdata, only: eps_dtog
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(mx,my,mz,mvar), intent(in) :: df
      type(pencil_case), intent(in) :: p
!
      real, dimension(nx) :: dux, duy, drhop
!
      call keep_compiler_quiet(f, df)
!
!  Identify module and boundary conditions.
!
      if (headtt.or.ldebug) print *, "dspecial_dt: special diagnostics"
!
!  Pencil precalculations
!
      penc: if (ldiagnos .or. l1davgfirst) then
        dux = p%uu(:,1) - ux0
        duy = p%uu(:,2) - uy0
        drhop = p%rhop - eps_dtog * rho0
      endif penc
!
!  Diagnostics
!
      diag: if (ldiagnos) then
        if (idiag_rduxm /= 0) call sum_mn_name(p%rho * dux, idiag_rduxm)
        if (idiag_rduym /= 0) call sum_mn_name(p%rho * duy, idiag_rduym)
        if (idiag_rdux2m /= 0) call sum_mn_name(p%rho * dux**2, idiag_rdux2m)
        if (idiag_rduy2m /= 0) call sum_mn_name(p%rho * duy**2, idiag_rduy2m)
        if (idiag_rduxduym /= 0) call sum_mn_name(p%rho * dux * duy, idiag_rduxduym)
        if (idiag_ruzduxm /= 0) call sum_mn_name(p%rho * p%uu(:,3) * dux, idiag_ruzduxm)
        if (idiag_ruzduym /= 0) call sum_mn_name(p%rho * p%uu(:,3) * duy, idiag_ruzduym)
!
        if (idiag_drhopm /= 0) call sum_mn_name(drhop, idiag_drhopm)
        if (idiag_drhop2m /= 0) call sum_mn_name(drhop**2, idiag_drhop2m)
      endif diag
!
!  1D averages
!
      avg1d: if (l1davgfirst) then
        if (idiag_rduxmx /= 0) call yzsum_mn_name_x(p%rho * dux, idiag_rduxmx)
        if (idiag_rduymx /= 0) call yzsum_mn_name_x(p%rho * duy, idiag_rduymx)
        if (idiag_rdux2mx /= 0) call yzsum_mn_name_x(p%rho * dux**2, idiag_rdux2mx)
        if (idiag_rduy2mx /= 0) call yzsum_mn_name_x(p%rho * duy**2, idiag_rduy2mx)
        if (idiag_rduxduymx /= 0) call yzsum_mn_name_x(p%rho * dux * duy, idiag_rduxduymx)
        if (idiag_ruzduxmx /= 0) call yzsum_mn_name_x(p%rho * p%uu(:,3) * dux, idiag_ruzduxmx)
        if (idiag_ruzduymx /= 0) call yzsum_mn_name_x(p%rho * p%uu(:,3) * duy, idiag_ruzduymx)
!
        if (idiag_drhopmx /= 0) call yzsum_mn_name_x(drhop, idiag_drhopm)
        if (idiag_drhop2mx /= 0) call yzsum_mn_name_x(drhop**2, idiag_drhop2m)
      endif avg1d
!
    endsubroutine dspecial_dt
!***********************************************************************
    subroutine special_calc_particles(f, df, fp, dfp, ineargrid)
!
!  Called before the loop, in case some particle value is needed
!  for the special density/hydro/magnetic/entropy.
!
!  04-oct-19/ccyang: coded
!
      use Particles_sub, only: assign_species, sum_par_name
      use Particles_cdata, only: ipar, npar_loc, ivpx, ivpy, ivpz, irhopswarm
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(mx,my,mz,mvar), intent(in) :: df
      real, dimension(:,:), intent(in) :: fp, dfp
      integer, dimension(:,:), intent(in) :: ineargrid
!
      real, dimension(npar_loc) :: dvpx, dvpy
      integer :: k, jspec
!
      call keep_compiler_quiet(f, df)
      call keep_compiler_quiet(dfp)
      call keep_compiler_quiet(ineargrid)
!
      diag: if (ldiagnos) then
!
!  Compute the deviation from the equilibrium velocities.
!
        vpeq: do k = 1, npar_loc
          jspec = assign_species(ipar(k))
          dvpx(k) = fp(k,ivpx) - vpx0(jspec)
          dvpy(k) = fp(k,ivpy) - vpy0(jspec)
        enddo vpeq
!
!  Compute the diagnostics
!
        rhop: if (lparticles_density) then
          if (idiag_rhopdvpxm /= 0) call sum_par_name(fp(1:npar_loc,irhopswarm) * dvpx, idiag_rhopdvpxm)
          if (idiag_rhopdvpym /= 0) call sum_par_name(fp(1:npar_loc,irhopswarm) * dvpy, idiag_rhopdvpym)
          if (idiag_rhopdvpx2m /= 0) call sum_par_name(fp(1:npar_loc,irhopswarm) * dvpx**2, idiag_rhopdvpx2m)
          if (idiag_rhopdvpy2m /= 0) call sum_par_name(fp(1:npar_loc,irhopswarm) * dvpy**2, idiag_rhopdvpy2m)
          if (idiag_rhopvpz2m /= 0) call sum_par_name(fp(1:npar_loc,irhopswarm) * fp(1:npar_loc,ivpz)**2, idiag_rhopvpz2m)
        endif rhop
!
      endif diag
!
    endsubroutine special_calc_particles
!***********************************************************************
    subroutine rprint_special(lreset, lwrite)
!
!  Reads and registers print parameters relevant to special.
!
!  03-aug-20/ccyang: coded
!
      use Diagnostics, only: parse_name
      use FArrayManager, only: farray_index_append
!
      logical, intent(in) :: lreset
      logical, intent(in), optional :: lwrite
!
      logical :: lwr
      integer :: iname, inamex
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
!
        idiag_drhopm = 0
        idiag_drhop2m = 0
        idiag_rhopdvpxm = 0
        idiag_rhopdvpym = 0
        idiag_rhopdvpx2m = 0
        idiag_rhopdvpy2m = 0
        idiag_rhopvpz2m = 0
!
        idiag_rduxmx = 0
        idiag_rduymx = 0
        idiag_rdux2mx = 0
        idiag_rduy2mx = 0
        idiag_rduxduymx = 0
        idiag_ruzduxmx = 0
        idiag_ruzduymx = 0
!
        idiag_drhopmx = 0
        idiag_drhop2mx = 0
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
!
        call parse_name(iname, cname(iname), cform(iname), "drhopm", idiag_drhopm)
        call parse_name(iname, cname(iname), cform(iname), "drhop2m", idiag_drhop2m)
        call parse_name(iname, cname(iname), cform(iname), "rhopdvpxm", idiag_rhopdvpxm)
        call parse_name(iname, cname(iname), cform(iname), "rhopdvpym", idiag_rhopdvpym)
        call parse_name(iname, cname(iname), cform(iname), "rhopdvpx2m", idiag_rhopdvpx2m)
        call parse_name(iname, cname(iname), cform(iname), "rhopdvpy2m", idiag_rhopdvpy2m)
        call parse_name(iname, cname(iname), cform(iname), "rhopvpz2m", idiag_rhopvpz2m)
      enddo diag
!
!  Turn on requested yz-averages.
!
      mx: do inamex = 1, nnamex
        call parse_name(inamex, cnamex(inamex), cformx(inamex), "rduxmx", idiag_rduxmx)
        call parse_name(inamex, cnamex(inamex), cformx(inamex), "rduymx", idiag_rduymx)
        call parse_name(inamex, cnamex(inamex), cformx(inamex), "rdux2mx", idiag_rdux2mx)
        call parse_name(inamex, cnamex(inamex), cformx(inamex), "rduy2mx", idiag_rduy2mx)
        call parse_name(inamex, cnamex(inamex), cformx(inamex), "rduxduymx", idiag_rduxduymx)
        call parse_name(inamex, cnamex(inamex), cformx(inamex), "ruzduxmx", idiag_ruzduxmx)
        call parse_name(inamex, cnamex(inamex), cformx(inamex), "ruzduymx", idiag_ruzduymx)
!
        call parse_name(inamex, cnamex(inamex), cformx(inamex), "drhopmx", idiag_drhopmx)
        call parse_name(inamex, cnamex(inamex), cformx(inamex), "drhop2mx", idiag_drhop2mx)
      enddo mx
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
!
        call farray_index_append("idiag_drhopm", idiag_drhopm)
        call farray_index_append("idiag_drhop2m", idiag_drhop2m)
        call farray_index_append("idiag_rhopdvpxm", idiag_rhopdvpxm)
        call farray_index_append("idiag_rhopdvpym", idiag_rhopdvpym)
        call farray_index_append("idiag_rhopdvpx2m", idiag_rhopdvpx2m)
        call farray_index_append("idiag_rhopdvpy2m", idiag_rhopdvpy2m)
        call farray_index_append("idiag_rhopvpz2m", idiag_rhopvpz2m)
!
        call farray_index_append("idiag_rduxmx", idiag_rduxmx)
        call farray_index_append("idiag_rduymx", idiag_rduymx)
        call farray_index_append("idiag_rdux2mx", idiag_rdux2mx)
        call farray_index_append("idiag_rduy2mx", idiag_rduy2mx)
        call farray_index_append("idiag_rduxduymx", idiag_rduxduymx)
        call farray_index_append("idiag_ruzduxmx", idiag_ruzduxmx)
        call farray_index_append("idiag_ruzduymx", idiag_ruzduymx)
!
        call farray_index_append("idiag_drhopmx", idiag_drhopmx)
        call farray_index_append("idiag_drhop2mx", idiag_drhop2mx)
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
