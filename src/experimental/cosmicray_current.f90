! $Id: cosmicrayflux.f90 19193 2012-06-30 12:55:46Z wdobler $
!
!  Module for calculating Cosmic Ray Flux.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lcosmicrayflux = .true.
!
! MVAR CONTRIBUTION 3
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED ucr(3); ucrij(3,3); ucrgucr(3); divucr; del2ucr(3)
!
!***************************************************************
module Cosmicrayflux
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
!
  implicit none
!
  include 'cosmicrayflux.h'
!
  character (len=labellen) :: initfcr='zero'
  real :: amplfcr=0., omegahat=0., fcr_const=0., J_param=0., Ma_param=1.
  real :: kx_fcr=1., ky_fcr=1., kz_fcr=1., cs2cr=1., gamma_cr=1., nu_cr=0.
  real, parameter :: rhocr0=1.
  logical :: lupw_ucr=.false.
!
  namelist /cosmicrayflux_init_pars/ &
       omegahat, initfcr, amplfcr, fcr_const, kx_fcr, ky_fcr, kz_fcr, &
       cs2cr, gamma_cr
!
  namelist /cosmicrayflux_run_pars/ &
       omegahat, lupw_ucr, J_param, Ma_param, cs2cr, nu_cr
!
  integer :: idiag_ekincr=0       ! DIAG_DOC: $\left<{1\over2}\varrho\uv_{\rm cr}^2\right>$
  integer :: idiag_ethmcr=0       ! DIAG_DOC: $\left<\varrho_{\rm cr} e_{\rm cr}\right>$
!
  contains
!***********************************************************************
    subroutine register_cosmicrayflux()
!
!  Initialise variables which should know that we solve for the vector
!  potential: ifcr, etc; increase nvar accordingly
!
!  1-may-02/wolf: coded
!
      use FArrayManager
!
      call farray_register_pde('fcr',ifcr,vector=3)
      ifcrx = ifcr; ifcry = ifcr+1; ifcrz = ifcr+2
!
!  Identify version number.
!
      if (lroot) call svn_id( &
          "$Id: cosmicrayflux.f90 19193 2012-06-30 12:55:46Z wdobler $")
!
!  Writing files for use with IDL
!
      if (lroot) then
        if (maux == 0) then
          if (nvar < mvar) write(4,*) ',fcr $'
          if (nvar == mvar) write(4,*) ',fcr'
        else
          write(4,*) ',fcr $'
        endif
        write(15,*) 'fcr = fltarr(mx,my,mz,3)*one'
      endif
!
    endsubroutine register_cosmicrayflux
!***********************************************************************
    subroutine initialize_cosmicrayflux(f)
!
!  Perform any post-parameter-read initialization
!
!  24-nov-02/tony: dummy routine - nothing to do at present
!  20-may-03/axel: reinitalize_aa added
!
      real, dimension (mx,my,mz,mfarray) :: f

      call keep_compiler_quiet(f)
!
    endsubroutine initialize_cosmicrayflux
!***********************************************************************
    subroutine init_fcr(f)
!
!  initialise magnetic field; called from start.f90
!  AB: maybe we should here call different routines (such as rings)
!  AB: and others, instead of accummulating all this in a huge routine.
!  We have an init parameter (initaa) to stear magnetic i.c. independently.
!
!   7-nov-2001/wolf: coded
!
      use Mpicomm
      use Sub
      use Initcond
      use InitialCondition, only: initial_condition_fcr
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      select case (initfcr)

      case ('zero', '0'); f(:,:,:,ifcrx:ifcrz) = 0.
      case ('const_fcr'); f(:,:,:,ifcrz) = fcr_const
      case ('wave-x'); call wave(amplfcr,f,ifcrx,kx=kx_fcr)
      case ('wave-y'); call wave(amplfcr,f,ifcry,ky=ky_fcr)
      case ('wave-z'); call wave(amplfcr,f,ifcrz,kz=kz_fcr)
!
      case default
!
!  Catch unknown values
!
        if (lroot) print*, 'init_fcr: No such such value for initfcr: ', trim(initfcr)
        call stop_it(" ")

      endselect
!
!  Interface for user's own initial condition
!
      if (linitial_condition) call initial_condition_fcr(f)
!
    endsubroutine init_fcr
!***********************************************************************
    subroutine pencil_criteria_cosmicrayflux()
!
!   All pencils that the Magnetic module depends on are specified here.
!
!  19-nov-04/anders: coded
!
      lpenc_requested(i_ucr)=.true.
      lpenc_requested(i_ecr1)=.true.
      lpenc_requested(i_ucrgucr)=.true.
      lpenc_requested(i_del2ucr)=.true.
      if (lhydro) lpenc_requested(i_rho1)=.true.
      lpenc_requested(i_gecr)=.true.
      lpenc_requested(i_bb)=.true.
!
    endsubroutine pencil_criteria_cosmicrayflux
!***********************************************************************
    subroutine pencil_interdep_cosmicrayflux(lpencil_in)
!
      logical, dimension(npencils) :: lpencil_in
!
      if (lpencil_in(i_divucr)) lpencil_in(i_ucrij)=.true.
!
      if (lpencil_in(i_ucrgucr)) then
        lpencil_in(i_ucr)=.true.
        lpencil_in(i_ucrij)=.true.
      endif
!
    endsubroutine pencil_interdep_cosmicrayflux
!***********************************************************************
    subroutine calc_pencils_cosmicrayflux(f,p)
!
!  Calculate Cosmicray Flux pencils - to be done
!
      use Sub, only: u_dot_grad, gij, div_mn, del2v
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      intent(in)  :: f
      intent(inout) :: p
!
! ucr
      if (lpencil(i_ucr)) p%ucr=f(l1:l2,m,n,ifcrx:ifcrz)
! ucrij
      if (lpencil(i_ucrij)) call gij(f,ifcr,p%ucrij,1)
! divucr
      if (lpencil(i_divucr)) call div_mn(p%ucrij,p%divucr,p%ucr)
! ucrgucr
      if (lpencil(i_ucrgucr)) then
        if (headtt.and.lupw_ucr) print *,'calc_pencils_cosmicray_current: upwinding advection term'
        call u_dot_grad(f,ifcr,p%ucrij,p%ucr,p%ucrgucr,UPWIND=lupw_ucr)
      endif
      if (lpencil(i_del2ucr)) call del2v(f,ifcr,p%del2ucr)
!
! fcr
!      if (lpencil(i_fcr)) p%fcr=f(l1:l2,m,n,ifcrx:ifcrz)
!
    endsubroutine calc_pencils_cosmicrayflux
!***********************************************************************
    subroutine dfcr_dt(f,df,p)
!
!  Cosmicray Flux evolution
!
!  08-mar-05/snod: adapted from daa_dt
!
      use Sub
      use Slices
      use Debug_IO, only: output_pencil
      use Mpicomm, only: stop_it
      use Diagnostics, only: sum_mn_name
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx,3) :: delucrxbb, delucrxbb2, gecr_over_ecr
      real, dimension (nx)   :: b2, b21, ucr2
      real, dimension (nx)   :: tmp, ratio
      real :: fact
      integer :: i,j
      type (pencil_case) :: p
!
      intent(in)     :: f
      intent(inout)  :: df
!
!  Identify module and boundary conditions.
!
      if (headtt.or.ldebug) print*,'dfcr_dt: SOLVE'
      if (headtt) then
        call identify_bcs('Fecx',ifcrx)
        call identify_bcs('Fecy',ifcry)
        call identify_bcs('Fecz',ifcrz)
      endif
!
!  Time step control.
!
      if (lfirst.and.ldt) advec_cs2cr=cs2cr*dxyz_2
!
!  Compute auxiliary terms.
!
      call cross(omegahat*(p%ucr-p%uu),p%bb,delucrxbb)
      call multsv(cs2cr*p%ecr1,p%gecr,gecr_over_ecr)
!
!  Take care of gamma_cr factor
!
      if (gamma_cr/=0.) then
        call multsv((p%ecr/rhocr0)**(gamma_cr-1.),gecr_over_ecr,gecr_over_ecr)
      endif
!
!  Cosmic Ray Flux equation.
!
      df(l1:l2,m,n,ifcrx:ifcrz) = df(l1:l2,m,n,ifcrx:ifcrz) &
        +delucrxbb-p%ucrgucr-gecr_over_ecr+nu_cr*p%del2ucr
!
!  Add Lorentz force
!
      if (lhydro) then
        ratio=p%ecr*p%rho1
        call multsv(-ratio,delucrxbb,delucrxbb2)
        df(l1:l2,m,n,iux:iuz)=df(l1:l2,m,n,iux:iuz)+delucrxbb2
      endif
!
!  Calculate diagnostic quantities.
!
      if (ldiagnos) then
        if (idiag_ekincr/=0) then
          call dot2 (p%ucr,ucr2)
          call sum_mn_name(.5*p%ecr*ucr2,idiag_ekincr)
        endif
!
        if (idiag_ethmcr/=0) then
          fact=(1.-1./gamma_cr)*rhocr0*cs2cr
          call sum_mn_name(fact*(p%ecr/rhocr0)**gamma_cr,idiag_ethmcr)
        endif
!
!  cosmicrayflux components at one point (=pt)
!
!        if (lroot.and.m==mpoint.and.n==npoint) then
!          if (idiag_fcrxpt/=0) call save_name(p%fcr(lpoint-nghost,1),idiag_fcrxpt)
!          if (idiag_fcrypt/=0) call save_name(p%fcr(lpoint-nghost,2),idiag_fcrypt)
!         if (idiag_fcrzpt/=0) call save_name(p%fcr(lpoint-nghost,3),idiag_fcrzpt)
!        endif
      endif ! endif (ldiagnos)
!
    endsubroutine dfcr_dt
!***********************************************************************
    subroutine read_cosmicrayflux_init_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=cosmicrayflux_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=cosmicrayflux_init_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_cosmicrayflux_init_pars
!***********************************************************************
    subroutine write_cosmicrayflux_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit,NML=cosmicrayflux_init_pars)
!
    endsubroutine write_cosmicrayflux_init_pars
!***********************************************************************
    subroutine read_cosmicrayflux_run_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=cosmicrayflux_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=cosmicrayflux_run_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_cosmicrayflux_run_pars
!***********************************************************************
    subroutine write_cosmicrayflux_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit,NML=cosmicrayflux_run_pars)
!
    endsubroutine write_cosmicrayflux_run_pars
!***********************************************************************
    subroutine rprint_cosmicrayflux(lreset,lwrite)
!
!  Reads and registers print parameters relevant for cosmicrayflux.
!
!   3-may-02/axel: coded
!  27-may-02/axel: added possibility to reset list
!
      use Diagnostics, only: parse_name
!
      integer :: iname,inamez,ixy,irz
      logical :: lreset,lwr
      logical, optional :: lwrite
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
!  Reset everything in case of RELOAD.
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        idiag_ekincr=0
        idiag_ethmcr=0
      endif
!
!  Check for those quantities that we want to evaluate online.
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'ekincr',idiag_ekincr)
        call parse_name(iname,cname(iname),cform(iname),'ethmcr',idiag_ethmcr)
      enddo
!
!  Check for those quantities for which we want xy-averages.
!
      do inamez=1,nnamez
!        call parse_name(inamez,cnamez(inamez),cformz(inamez),'bxmz',idiag_bxmz)
      enddo
!
!  Check for those quantities for which we want z-averages.
!
      do ixy=1,nnamexy
!        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'bxmxy',idiag_bxmxy)
      enddo
!
!  Check for those quantities for which we want phi-averages.
!
      do irz=1,nnamerz
!        call parse_name(irz,cnamerz(irz),cformrz(irz),'brmphi'  ,idiag_brmphi)
      enddo
!
!  Write column, idiag_XYZ, where our variable XYZ is stored.
!
      if (lwr) then
        write(3,*) 'ifcr=',ifcr
      endif
!
    endsubroutine rprint_cosmicrayflux
!***********************************************************************
endmodule Cosmicrayflux
