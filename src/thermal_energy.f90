! $Id$
!
!  This module can replace the entropy module by using the thermal energy
!  eth as dependent variable. For a perfect gas we have
!
!    deth/dt + div(u*eth) = -P*div(u)
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lentropy = .false.
! CPARAM logical, parameter :: ltemperature = .false.
! CPARAM logical, parameter :: lthermal_energy = .true.
!
! MVAR CONTRIBUTION 1
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED Ma2; fpres(3); transpeth
!
!***************************************************************
module Entropy
!
  use Cdata
  use Cparam
  use Messages
  use Sub, only: keep_compiler_quiet
!
  implicit none
!
  include 'entropy.h'
!
  real :: eth_left, eth_right, widtheth, eth_const=1.0
  logical :: lviscosity_heat=.true.
  character (len=labellen), dimension(ninit) :: initeth='nothing'
!
!  Input parameters.
!
  namelist /entropy_init_pars/ &
      initeth, eth_left, eth_right, widtheth, eth_const
!
!  Run parameters.
!
  namelist /entropy_run_pars/ &
      lviscosity_heat
!
!  Diagnostic variables (needs to be consistent with reset list below).
!
  integer :: idiag_TTmax=0    ! DIAG_DOC: $\max (T)$
  integer :: idiag_TTmin=0    ! DIAG_DOC: $\min (T)$
  integer :: idiag_eem=0      ! DIAG_DOC: $\left< e \right> =
                              ! DIAG_DOC:  \left< c_v T \right>$
                              ! DIAG_DOC: \quad(mean internal energy)
  integer :: idiag_ppm=0
  integer :: idiag_ppmx=0
  integer :: idiag_ppmy=0
  integer :: idiag_ppmz=0
  integer :: idiag_TTm=0
  integer :: idiag_TTmx=0
  integer :: idiag_TTmy=0
  integer :: idiag_TTmz=0
  integer :: idiag_TTmxy=0
  integer :: idiag_TTmxz=0
  integer :: idiag_pdivum=0
  integer :: idiag_ethm=0     ! DIAG_DOC: $\left< e_{\text{th}}\right> =
                              ! DIAG_DOC:  \left< c_v \rho T \right> $
                              ! DIAG_DOC: \quad(mean thermal energy)
!
  real, dimension(nx,nz) :: pp_xz
  real, dimension(ny,nz) :: pp_yz
  real, dimension(nx,ny) :: pp_xy,pp_xy2,pp_xy3,pp_xy4
!
  contains
!***********************************************************************
    subroutine register_entropy()
!
!  Initialise variables which should know that we solve an energy equation.
!
!  04-nov-10/anders+evghenii: adapted
!
      use FArrayManager, only: farray_register_pde
!
      call farray_register_pde('eth',ieth)
!
!  Identify version number.
!
      if (lroot) call svn_id( &
           "$Id$")
!
    endsubroutine register_entropy
!***********************************************************************
    subroutine initialize_entropy(f,lstarting)
!
!  Called by run.f90 after reading parameters, but before the time loop.
!
!  04-nov-10/anders+evghenii: adapted
!
      use EquationOfState, only: select_eos_variable
      use SharedVariables
!
      real, dimension (mx,my,mz,mfarray), intent (inout) :: f
      logical, intent (in) :: lstarting
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(lstarting)
!
      call select_eos_variable('eth',ieth)
!
      call put_shared_variable('lviscosity_heat',lviscosity_heat)
!
    endsubroutine initialize_entropy
!***********************************************************************
    subroutine init_ss(f)
!
!  Initialise thermal energy.
!
!  04-nov-10/anders+evghenii: adapted
!
      use General, only: chn
      use Initcond, only: jump
!
      real, dimension (mx,my,mz,mfarray), intent (inout) :: f
!
      integer :: j
      logical :: lnothing=.true.
      character (len=5) :: iinit_str
!
      do j=1,ninit
!
        if (initeth(j)/='nothing') then
!
          lnothing=.false.
!
          call chn(j,iinit_str)
!
!  Select between various initial conditions.
!
          select case (initeth(j))
!
          case ('zero', '0'); f(:,:,:,ieth) = 0.0
!
          case ('xjump'); call jump(f,ieth,eth_left,eth_right,widtheth,'x')
!
          case ('const_eth'); f(:,:,:,ieth)=f(:,:,:,ieth)+eth_const
!
          case default
!
!  Catch unknown values.
!
            write(unit=errormsg,fmt=*) 'No such value for initeth(' &
                //trim(iinit_str)//'): ',trim(initeth(j))
            call fatal_error('init_ss',errormsg)
!
          endselect
        endif
      enddo
!
    endsubroutine init_ss
!***********************************************************************
    subroutine pencil_criteria_entropy()
!
!  All pencils that the Entropy module depends on are specified here.
!
!  04-nov-10/anders+evghenii: adapted
!
      if (lweno_transport) then
        lpenc_requested(i_transpeth)=.true.
      endif
!  In case of advection we need divu
      lpenc_requested(i_divu)=.true.
!
      lpenc_requested(i_eth)=.true.
      lpenc_requested(i_fpres)=.true.
      if (ldt) lpenc_requested(i_cs2)=.true.
      if (lviscosity.and.lviscosity_heat) lpenc_requested(i_visc_heat)=.true.
!
      if (idiag_ppm/=0) lpenc_diagnos(i_pp)=.true.
      if (idiag_pdivum/=0) then
        lpenc_diagnos(i_pp)=.true.
        lpenc_diagnos(i_divu)=.true.
      endif
      if (idiag_TTm/=0 .or. idiag_TTmin/=0 .or. idiag_TTmax/=0) &
          lpenc_diagnos(i_TT)=.true.
!
    endsubroutine pencil_criteria_entropy
!***********************************************************************
    subroutine pencil_interdep_entropy(lpencil_in)
!
!  Interdependency among pencils from the Entropy module is specified here.
!
!  04-nov-10/anders+evghenii: adapted
!  14-feb-11/bing: update
!
      logical, dimension(npencils) :: lpencil_in
!
      if (lpencil_in(i_fpres)) then
        lpencil_in(i_rho1)=.true.
        lpencil_in(i_geth)=.true.
      endif
!
      if (lpencil_in(i_cs2).or. &
          lpencil_in(i_TT).or. &
          lpencil_in(i_lnTT)) then
        lpencil_in(i_eth)=.true.
        lpencil_in(i_rho1)=.true.
      endif
!
      if (lpencil_in(i_pp)) lpencil_in(i_eth)=.true.
!
    endsubroutine pencil_interdep_entropy
!***********************************************************************
    subroutine calc_pencils_entropy(f,p)
!
!  Calculate Entropy pencils.
!  This routine is called after  calc_pencils_eos
!  Most basic pencils should come first, as others may depend on them.
!
!  04-nov-10/anders+evghenii: adapted
!  14-feb-11/bing: moved eth dependend pecnils to eos_idealgas
!
      use EquationOfState, only: gamma_m1
      use WENO_transport, only: weno_transp
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      integer :: j
! fpres
      if (lpencil(i_fpres)) then
        do j=1,3
          p%fpres(:,j)=-p%rho1*gamma_m1*p%geth(:,j)
        enddo
      endif
! transpeth
      if (lpencil(i_transpeth)) &
          call weno_transp(f,m,n,ieth,-1,iux,iuy,iuz,p%transpeth,dx_1,dy_1,dz_1)
!
    endsubroutine calc_pencils_entropy
!***********************************************************************
    subroutine dss_dt(f,df,p)
!
!  Calculate right hand side of energy equation.
!
!    deth/dt + div(u*eth) = -P*div(u)
!
!  04-nov-10/anders+evghenii: coded
!
      use Diagnostics
      use EquationOfState, only: gamma_m1
      use Sub, only: identify_bcs
      use Viscosity, only: calc_viscous_heat
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      real, dimension (nx) :: Hmax=0.0
!
      intent(inout) :: f,p
      intent(out) :: df
!
!  Identify module and boundary conditions.
!
      if (headtt.or.ldebug) print*, 'dss_dt: solve deth_dt'
      if (headtt) call identify_bcs('eth',ieth)
!
!  Sound speed squared.
!
      if (headtt) print*, 'dss_dt: cs20=', p%cs2(1)
!
!  ``cs2/dx^2'' for timestep
!
      if (lfirst.and.ldt) advec_cs2=p%cs2*dxyz_2
      if (headtt.or.ldebug) print*,'dss_dt: max(advec_cs2) =',maxval(advec_cs2)
!
!  Add pressure gradient term in momentum equation.
!
      if (lhydro) df(l1:l2,m,n,iux:iuz) = df(l1:l2,m,n,iux:iuz) + p%fpres
!
!  Add energy transport term.
!
      if (lweno_transport) then
        df(l1:l2,m,n,ieth) = df(l1:l2,m,n,ieth) - p%transpeth
      else
        call fatal_error('dss_dt','only implemented for WENO transport')
      endif
!
!  Add P*dV work.
!
      df(l1:l2,m,n,ieth) = df(l1:l2,m,n,ieth) - gamma_m1*p%eth*p%divu
!
!  Calculate viscous contribution to temperature.
!
      if (lviscosity.and.lviscosity_heat) call calc_viscous_heat(f,df,p,Hmax)
!
!  Diagnostics.
!
      if (ldiagnos) then
        if (idiag_TTm/=0)    call sum_mn_name(p%TT,idiag_TTm)
        if (idiag_TTmax/=0)  call max_mn_name(p%TT,idiag_TTmax)
        if (idiag_TTmin/=0)  call max_mn_name(-p%TT,idiag_TTmin,lneg=.true.)
        if (idiag_ethm/=0)   call sum_mn_name(p%eth,idiag_ethm)
        if (idiag_pdivum/=0) call sum_mn_name(p%pp*p%divu,idiag_pdivum)
      endif
!
      if (l1davgfirst) then
        if (idiag_ppmx/=0) call yzsum_mn_name_x(p%pp,idiag_ppmx)
        if (idiag_ppmy/=0) call xzsum_mn_name_y(p%pp,idiag_ppmy)
        if (idiag_ppmz/=0) call xysum_mn_name_z(p%pp,idiag_ppmz)
        if (idiag_TTmx/=0) call yzsum_mn_name_x(p%TT,idiag_TTmx)
        if (idiag_TTmy/=0) call xzsum_mn_name_y(p%TT,idiag_TTmy)
        if (idiag_TTmz/=0) call xysum_mn_name_z(p%TT,idiag_TTmz)
      endif
!
!  2-D averages.
!
      if (l2davgfirst) then
        if (idiag_TTmxy/=0) call zsum_mn_name_xy(p%TT,idiag_TTmxy)
        if (idiag_TTmxz/=0) call ysum_mn_name_xz(p%TT,idiag_TTmxz)
      endif
!
      if (lvideo.and.lfirst) then
        pp_yz(m-m1+1,n-n1+1)=p%pp(ix_loc-l1+1)
        if (m==iy_loc)  pp_xz(:,n-n1+1)=p%pp
        if (n==iz_loc)  pp_xy(:,m-m1+1)=p%pp
        if (n==iz2_loc) pp_xy2(:,m-m1+1)=p%pp
        if (n==iz3_loc) pp_xy3(:,m-m1+1)=p%pp
        if (n==iz4_loc) pp_xy4(:,m-m1+1)=p%pp
      endif
!
    endsubroutine dss_dt
!***********************************************************************
    subroutine calc_lentropy_pars(f)
!
!  Dummy routine.
!
!  04-nov-10/anders+evghenii: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      intent(in) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine calc_lentropy_pars
!***********************************************************************
    subroutine read_entropy_init_pars(unit,iostat)
!
!  04-nov-10/anders+evghenii: coded
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=entropy_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=entropy_init_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_entropy_init_pars
!***********************************************************************
    subroutine write_entropy_init_pars(unit)
!
!  04-nov-10/anders+evghenii: coded
!
      integer, intent(in) :: unit
!
      write(unit,NML=entropy_init_pars)
!
    endsubroutine write_entropy_init_pars
!***********************************************************************
    subroutine read_entropy_run_pars(unit,iostat)
!
!  04-nov-10/anders+evghenii: coded
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=entropy_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=entropy_run_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_entropy_run_pars
!***********************************************************************
    subroutine write_entropy_run_pars(unit)
!
!  04-nov-10/anders+evghenii: coded
!
      integer, intent(in) :: unit
!
      write(unit,NML=entropy_run_pars)
!
    endsubroutine write_entropy_run_pars
!***********************************************************************
    subroutine rprint_entropy(lreset,lwrite)
!
!  Reads and registers print parameters relevant to entropy.
!
!  04-nov-10/anders+evghenii: adapted from temperature_idealgas.f90
!
      use Diagnostics, only: parse_name
!
      logical :: lreset
      logical, optional :: lwrite
!
      integer :: iname, inamex, inamey, inamez, inamexy, inamexz
      logical :: lwr
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
!  Reset everything in case of reset.
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        idiag_TTm=0; idiag_TTmax=0; idiag_TTmin=0
        idiag_ethm=0; idiag_pdivum=0; idiag_eem=0; idiag_ppm=0
      endif
!
!  iname runs through all possible names that may be listed in print.in
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'TTm',idiag_TTm)
        call parse_name(iname,cname(iname),cform(iname),'TTmax',idiag_TTmax)
        call parse_name(iname,cname(iname),cform(iname),'TTmin',idiag_TTmin)
        call parse_name(iname,cname(iname),cform(iname),'ethm',idiag_ethm)
        call parse_name(iname,cname(iname),cform(iname),'eem',idiag_eem)
        call parse_name(iname,cname(iname),cform(iname),'ppm',idiag_ppm)
        call parse_name(iname,cname(iname),cform(iname),'pdivum',idiag_pdivum)
      enddo
!
!  Check for those quantities for which we want yz-averages.
!
      do inamex=1,nnamex
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'ppmx',idiag_ppmx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'TTmx',idiag_TTmx)
      enddo
!
!  Check for those quantities for which we want xz-averages.
!
      do inamey=1,nnamey
        call parse_name(inamey,cnamey(inamey),cformy(inamey),'ppmy',idiag_ppmy)
        call parse_name(inamey,cnamey(inamey),cformy(inamey),'TTmy',idiag_TTmy)
      enddo
!
!  Check for those quantities for which we want xy-averages.
!
      do inamez=1,nnamez
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'ppmz',idiag_ppmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'TTmz',idiag_TTmz)
      enddo
!
!  Check for those quantities for which we want z-averages.
!
      do inamexy=1,nnamexy
        call parse_name(inamexy,cnamexy(inamexy),cformxy(inamexy),'TTmxy', &
            idiag_TTmxy)
      enddo
!
!  Check for those quantities for which we want y-averages.
!
      do inamexz=1,nnamexz
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'TTmxz', &
            idiag_TTmxz)
      enddo
!
    endsubroutine rprint_entropy
!***********************************************************************
    subroutine get_slices_entropy(f,slices)
!
!  04-nov-10/anders+evghenii: adapted
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
      call keep_compiler_quiet(f)
!
!  Loop over slices
!
      select case (trim(slices%name))
!  Pressure
        case ('pp')
          slices%yz =pp_yz
          slices%xz =pp_xz
          slices%xy =pp_xy
          slices%xy2=pp_xy2
          if (lwrite_slice_xy3) slices%xy3=pp_xy3
          if (lwrite_slice_xy4) slices%xy4=pp_xy4
          slices%ready=.true.
!
      endselect
!
    endsubroutine get_slices_entropy
!***********************************************************************
    subroutine fill_farray_pressure(f)
!
!  04-nov-10/anders+evghenii: adapted
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine fill_farray_pressure
!***********************************************************************
endmodule Entropy
