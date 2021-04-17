! $Id$
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
!***************************************************************
module Cosmicrayflux
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages

  implicit none

  include 'cosmicrayflux.h'

  character(len=labellen) :: initfcr='zero'
  real, pointer :: K_para, K_perp
  real :: kpara_t=0., kperp_t=0.
  real :: tau=0., tau1=0.
  real :: subgrid_c1=0., subgrid_c2=0.
  real :: subgrid_brms=1.
  real :: subgrid_bmin=0.
  real :: subgrid_s=5./3.
  real :: subgrid_k=0.
  real :: ratio_kpara_kperp=0.
  real, dimension(nx) :: vKperp, vKpara
  real, dimension(nx) :: b_exp
  logical :: lsubgrid_cr = .false.
  logical :: lcosmicrayflux_diffus_dt = .false.
  logical :: ladvect_fcr=.false., lupw_fcr=.false.

  namelist /cosmicrayflux_init_pars/ &
      tau, lsubgrid_cr, &
      lcosmicrayflux_diffus_dt, &
      subgrid_brms, subgrid_bmin, &
      subgrid_c1,subgrid_c2, &
      subgrid_s, subgrid_k, &
      ratio_kpara_kperp



  namelist /cosmicrayflux_run_pars/ &
      tau, lsubgrid_cr, &
      lcosmicrayflux_diffus_dt, ladvect_fcr, lupw_fcr, &
      subgrid_brms, subgrid_bmin, &
      subgrid_c1,subgrid_c2, &
      subgrid_s, subgrid_k, &
      ratio_kpara_kperp

  contains
!*******************************************************************************
    subroutine register_cosmicrayflux()
      ! Initialise variables which should know that we solve for the vector
      ! potential: ifcr, etc; increase nvar accordingly
      !
      ! 1-may-02/wolf: coded

      use FArrayManager

      call farray_register_pde('fcr',ifcr,vector=3)
      ifcrx = ifcr
      ifcry = ifcr+1
      ifcrz = ifcr+2

      ! Identify version number.
      if (lroot) call svn_id("$Id$")

      ! Writing files for use with IDL
      if (lroot) then
        if (maux == 0) then
          if (nvar < mvar) write (4,*) ',fcr $'
          if (nvar == mvar) write (4,*) ',fcr'
        else
          write (4,*) ',fcr $'
        endif
        write (15,*) 'fcr = fltarr(mx,my,mz,3)*one'
      endif
!
    endsubroutine register_cosmicrayflux
!*******************************************************************************
    subroutine initialize_cosmicrayflux(f)
      ! Perform any post-parameter-read initialization
      use Messages, only: warning
      use SharedVariables, only: get_shared_variable
      real, dimension(mx,my,mz,mfarray) :: f

      if (tau /= 0.)  then
        tau1 = 1./tau
      else
        call warning('initialize_cosmicrayflux', &
            'parameter tau was set to zero!')
      endif

      ! Reads shared diffusivities from the cosmicray module
      call get_shared_variable('K_para',K_para)
      call get_shared_variable('K_perp',K_perp)

      kpara_t = K_para * tau1
      kperp_t = K_perp * tau1

      call keep_compiler_quiet(f)

    endsubroutine initialize_cosmicrayflux
!*******************************************************************************
    subroutine init_fcr(f)
      ! initialise magnetic field; called from start.f90
      ! AB: maybe we should here call different routines (such as rings)
      ! AB: and others, instead of accummulating all this in a huge routine.
      ! We have an init parameter (initaa) to stear magnetic i.c. independently.
      !
      ! 7-nov-2001/wolf: coded

      use Mpicomm
      use Sub
      use Initcond
      use InitialCondition, only: initial_condition_fcr

      real, dimension(mx,my,mz,mfarray) :: f

      select case (initfcr)
        case ('zero', '0')
          f(:,:,:,ifcrx:ifcrz) = 0.
          ! probably no more cases needed for fcr
        case default
          !  Catch unknown values
          if (lroot) print*, 'init_fcr: No such such value for initfcr: ', trim(initfcr)
          call stop_it(" ")
      endselect

      ! Interface for user's own initial condition
      if (linitial_condition) call initial_condition_fcr(f)
    endsubroutine init_fcr

!*******************************************************************************
    subroutine pencil_criteria_cosmicrayflux()
      ! All pencils that the Magnetic module depends on are specified here.
      !
      ! 19-11-04/anders: coded

      lpenc_requested(i_gecr) = .true.
      lpenc_requested(i_bb) = .true.
    endsubroutine pencil_criteria_cosmicrayflux
!*******************************************************************************
    subroutine pencil_interdep_cosmicrayflux(lpencil_in)
      logical, dimension(npencils) :: lpencil_in

      call keep_compiler_quiet(lpencil_in)
    endsubroutine pencil_interdep_cosmicrayflux
!*******************************************************************************
    subroutine calc_pencils_cosmicrayflux(f,p)
      ! Calculate Cosmicray Flux pencils - to be done

      real, dimension(mx,my,mz,mfarray) :: f
      type (pencil_case) :: p

      intent(in) :: f
      intent(inout) :: p

      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)
      ! fcr
      !if (lpencil(i_fcr)) p%fcr=f(l1:l2,m,n,ifcrx:ifcrz)
    endsubroutine calc_pencils_cosmicrayflux
!*******************************************************************************
    subroutine dfcr_dt(f,df,p)
      ! Cosmicray Flux evolution
      ! 08-mar-05/snod: adapted from daa_dt
      ! 12-jan-15/luiz: diffusion/advection contribution to the timestep
      ! 22-jun-17/luiz: subgrid prescription for the cosmic ray diffusion
      use Sub
      use Debug_IO, only: output_pencil
      use Mpicomm, only: stop_it
      use General, only: notanumber

      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mx,my,mz,mvar) :: df
      real, dimension(nx,3) :: BuiBujgecr, bunit
      real, dimension(nx)   :: b2, b21, b_brms
      real, dimension(nx)   :: tmp, diffus_cr,advec_kfcr
      real, dimension (nx,3,3) :: gfcr
      real, dimension (nx,3) :: ugfcr
      integer :: i, j
      type (pencil_case) :: p

      intent(in)     :: f
      intent(inout)  :: df

      !  Identifies module and boundary conditions.
      if (headtt .or. ldebug) print*,'dfcr_dt: SOLVE'
      if (headtt) then
        call identify_bcs('Fecx',ifcrx)
        call identify_bcs('Fecy',ifcry)
        call identify_bcs('Fecz',ifcrz)
      endif

      call dot2_mn(p%bb,b2)
      ! with frequency omega_Bz_ext
      b21 = 1./max(tini,b2)
      call multsv_mn(sqrt(b21),p%bb,bunit)

      do i=1,3
        tmp = 0.
        do j=1, 3
          tmp = tmp + bunit(:,i)*bunit(:,j)*p%gecr(:,j)
        enddo
        BuiBujgecr(:,i) = tmp
      enddo

      !  Cosmic Ray Flux equation.
      if (.not.lsubgrid_cr) then
        ! If not using the subgrid prescription, the diffusivities are
        ! constants
        df(l1:l2,m,n,ifcrx:ifcrz) = df(l1:l2,m,n,ifcrx:ifcrz) &
            - tau1*f(l1:l2,m,n,ifcrx:ifcrz)                   &
            - kperp_t*p%gecr                                  &
            - (kpara_t - kperp_t)*BuiBujgecr
      else
        ! Subgrid prescription for parallel diffusion:
        ! kpara = kpara_0 * [1 + c1*(B/Brms)^s + c2*(Brms/B)]
        ! Brms (subgrid_brms) is the rms value of the assumed background field
        ! and s is its spectral index.
        ! (subgrid_)c1 and c2 depend on properties of background field
        vKpara(:) = kpara_t

        b_brms = sqrt(b2)/subgrid_brms
        if (subgrid_c1 /= 0.0) then
          vKpara(:) = vKpara(:) + kpara_t*subgrid_c1*b_brms**subgrid_s
        endif

        if (subgrid_c2 /= 0.0) then
          vKpara(:) = vKpara(:) + kpara_t*subgrid_c2*subgrid_brms*sqrt(b21)
        endif

          ! Subgrid prescription for perpendicular diffusion:
        if (ratio_kpara_kperp /= 0.0) then
          ! If a ration k_para/k_perp was specified, ignore the value of
          ! k_perp and computes the perpendicular diffusivity accordingly
          ! kperp = kpara/r * (B/Brms)^(-k)
          ! where r = ratio_kpara_kperp) and k = subgrid_k
          vKperp(:) = vKpara(:)/ratio_kpara_kperp
          if (subgrid_k /= 0.0) then
            vKperp(:) = vKperp(:) * (b_brms)**(-subgrid_k)
          endif
        else
          ! Otherwise, use constant perpendicular diffusivity previously
          ! specified
          vKperp(:) = kperp_t
        endif

        do i=1,3
          df(l1:l2,m,n,ifcrx+i-1) = df(l1:l2,m,n,ifcrx+i-1) &
              - tau1*f(l1:l2,m,n,ifcrx+i-1)                 &
              - vKperp*p%gecr(:,i)                          &
              - (vKpara - vKperp)*BuiBujgecr(:,i)
        enddo
      endif

      !  Allows optional use of advection term for fcr.
      if (ladvect_fcr) then
        call gij(f,ifcr,gfcr,1)
        call u_dot_grad(f,ifcr,gfcr,p%uu,ugfcr,UPWIND=lupw_fcr)
        df(l1:l2,m,n,ifcrx:ifcrz) = df(l1:l2,m,n,ifcrx:ifcrz) &
            - ugfcr(1:nx,1:3)
      endif

      ! For the timestep calculation, needs maximum diffusion or advection.
      ! Unless the switch lcosmicrayflux_diffus_dt is used, kpara_t and kperp_t
      ! are treated as an advection contribution. Otherwise,
      ! tau*kperp_t/tau*kpara_t are used as cosmic ray diffusivities.
      if (lfirst .and. ldt) then
        if (lcosmicrayflux_diffus_dt) then
          if (lsubgrid_cr) then
            diffus_cr = max(maxval(vKperp)*tau*dxyz_2, &
                            maxval(vKpara)*tau*dxyz_2)
          else
            diffus_cr = max(kperp_t*tau*dxyz_2,kpara_t*tau*dxyz_2)
          endif
          maxdiffus=max(maxdiffus,diffus_cr)
        else
          if (lsubgrid_cr) then
            advec_kfcr = max(maxval(vKperp)*dxyz_2, &
                             maxval(vKpara)*dxyz_2)
          else
            advec_kfcr = max(kperp_t*dxyz_2, kpara_t*dxyz_2)
          endif
          if (notanumber(advec_kfcr)) print*, 'advec_kfcr =',advec_kfcr
          advec2=advec2+advec_kfcr
        endif
      endif

    endsubroutine dfcr_dt
!*******************************************************************************
    subroutine read_cosmicrayflux_init_pars(iostat)
      use File_io, only: parallel_unit
      integer, intent(out) :: iostat
      read(parallel_unit, NML=cosmicrayflux_init_pars, IOSTAT=iostat)
    endsubroutine read_cosmicrayflux_init_pars
!*******************************************************************************
    subroutine write_cosmicrayflux_init_pars(unit)
      integer, intent(in) :: unit
      write(unit, NML=cosmicrayflux_init_pars)
    endsubroutine write_cosmicrayflux_init_pars
!*******************************************************************************
    subroutine read_cosmicrayflux_run_pars(iostat)
      use File_io, only: parallel_unit
      integer, intent(out) :: iostat
      read(parallel_unit, NML=cosmicrayflux_run_pars, IOSTAT=iostat)
    endsubroutine read_cosmicrayflux_run_pars
!*******************************************************************************
    subroutine write_cosmicrayflux_run_pars(unit)
      integer, intent(in) :: unit
      write(unit, NML=cosmicrayflux_run_pars)
    endsubroutine write_cosmicrayflux_run_pars
!*******************************************************************************
    subroutine rprint_cosmicrayflux(lreset,lwrite)
      ! Reads and registers print parameters relevant for cosmicrayflux.
      ! 3-may-02/axel: coded
      ! 27-may-02/axel: added possibility to reset list
      use Sub

      integer :: iname, inamez, ixy, irz
      logical :: lreset
      logical, optional :: lwrite

      ! Reset everything in case of RELOAD.
      ! (this needs to be consistent with what is defined above!)
      if (lreset) then
        !idiag_b2m=0; idiag_bm2=0; idiag_j2m=0; idiag_jm2=0; idiag_abm=0
      endif

      ! Check for those quantities that we want to evaluate online.
      do iname=1,nname
!        call parse_name(iname,cname(iname),cform(iname),'dteta',idiag_dteta)
      enddo

      ! Check for those quantities for which we want xy-averages.
      do inamez=1,nnamez
!        call parse_name(inamez,cnamez(inamez),cformz(inamez),'bxmz',idiag_bxmz)
      enddo

      ! Check for those quantities for which we want z-averages.
      do ixy=1,nnamexy
!        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'bxmxy',idiag_bxmxy)
      enddo

      ! Check for those quantities for which we want phi-averages.
      do irz=1,nnamerz
!        call parse_name(irz,cnamerz(irz),cformrz(irz),'brmphi',idiag_brmphi)
      enddo
    endsubroutine rprint_cosmicrayflux
!*******************************************************************************
endmodule Cosmicrayflux
