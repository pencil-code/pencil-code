! $Id$
!
!  This module takes care of self gravity by solving the Poisson equation
!    (d^2/dx^2 + d^2/dy^2 + d^2/dz^2)phi = 4*pi*G*rho
!  for the potential phi.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 1
! COMMUNICATED AUXILIARIES 1
!
! PENCILS PROVIDED potself; gpotself(3)
!
!***************************************************************
module Selfgravity
!
  use Cdata
  use Cparam
  use Messages
  use Sub, only: keep_compiler_quiet
!
  implicit none
!
  include 'selfgravity.h'
!
  real, target :: rhs_poisson_const=1.0
  real, target :: tstart_selfgrav=0.0
  real :: gravitational_const=0.0
!
  logical :: lselfgravity_gas=.true., lselfgravity_dust=.false.
  logical :: lselfgravity_neutrals=.false.
!
  namelist /selfgrav_init_pars/ &
      rhs_poisson_const, lselfgravity_gas, lselfgravity_dust, &
      lselfgravity_neutrals,tstart_selfgrav, gravitational_const

  namelist /selfgrav_run_pars/ &
      rhs_poisson_const, lselfgravity_gas, lselfgravity_dust, &
      lselfgravity_neutrals,tstart_selfgrav, gravitational_const
!
  integer :: idiag_gpoten=0, idiag_gpotenmxy=0
  integer :: idiag_gpotselfxm=0, idiag_gpotselfym=0, idiag_gpotselfzm=0
  integer :: idiag_gpotselfx2m=0, idiag_gpotselfy2m=0, idiag_gpotselfz2m=0
  integer :: idiag_gxgym=0,idiag_gxgzm=0,idiag_gygzm=0
  integer :: idiag_grgpm=0,idiag_grgzm=0,idiag_gpgzm=0
!
  contains
!***********************************************************************
    subroutine register_selfgravity()
!
!  Initialise self gravity variables.
!
!  15-may-06/anders+jeff: adapted
!
      use FArrayManager
!
!  Set indices for auxiliary variables
!
      call farray_register_auxiliary('potself',ipotself,communicated=.true.)
!
!  Identify version number (generated automatically by CVS)
!
      if (lroot) call cvs_id( &
          "$Id$")
!
!  Set lselfgravity.
!
      lselfgravity=.true.
!
    endsubroutine register_selfgravity
!***********************************************************************
    subroutine initialize_selfgravity()
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  15-may-06/anders+jeff: adapted
!
      use SharedVariables
!
      integer :: ierr=0
!
      if (gravitational_const .ne. 0) then 
        !gravitational const was set, re-define rhs_poisson_const
        rhs_poisson_const=4*pi*gravitational_const
      endif
!
      if (.not.lpoisson) then
        if (lroot) print*, 'initialize_selfgravity: must choose a Poisson '// &
            'solver in Makefile.local for self-gravity'
        call fatal_error('initialize_selfgravity','')
      endif
!
      if (lselfgravity_gas.and..not.(lhydro.and.ldensity)) then
        if (lroot) print*, 'initialize_selfgravity: must choose a hydro '// &
            'and a density module in Makefile.local for self-gravity'
        call fatal_error('initialize_selfgravity','')
      endif
!
      if (lselfgravity_dust.and..not.(ldustvelocity.and.ldustdensity)) then
        if (lroot) then
          print*, 'initialize_selfgravity: must choose a dust velocity '// &
              'and a dust density module in Makefile.local for '// &
              'self-gravity on the dust fluid'
        endif
        call fatal_error('initialize_selfgravity','')
      endif
!
      if (lselfgravity_neutrals.and.&
           .not.(lneutralvelocity.and.lneutraldensity)) then
        if (lroot) then
          print*, 'initialize_selfgravity: must choose a neutral velocity '// &
              'and a neutral density module in Makefile.local for '// &
              'self-gravity on the neutral fluid'
        endif
        call fatal_error('initialize_selfgravity','')
      endif
!
!  Share the variable tstart_selfgrav so that it can be used by other
!  self-gravity modules.
!
      call put_shared_variable('tstart_selfgrav',tstart_selfgrav,ierr)
      if (ierr/=0) then
        if (lroot) print*, 'initialize_selfgravity: there was a problem '// &
            'when sharing tstart_selfgrav!'
        call fatal_error('initialize_selfgravity','')
      endif
!
!  Also share rhs_poisson_const
!      
      call put_shared_variable('rhs_poisson_const',rhs_poisson_const,ierr)
      if (ierr/=0) then
        if (lroot) print*, 'initialize_selfgravity: there was a problem '// &
            'when sharing rhs_poisson_const!'
        call fatal_error('initialize_selfgravity','')
      endif
!
!  Check that density and self-potential have consistent boundary conditions.
!
      if (ldensity) then
        if (bcx(ipotself)=='p' .and. .not.(bcx(ilnrho)=='p')) then
          if (lroot) then
            print*, 'initialize_selfgravity: potself has bcx=''p'', but the density is not'
            print*, '                        periodic! (you must set a proper boundary condition'
            print*, '                        for the potential)'
            print*, 'initialize_selfgravity: bcx=', bcx
          endif
          call fatal_error('initialize_selfgravity','')
        endif
        if (bcy(ipotself)=='p' .and. .not.(bcy(ilnrho)=='p')) then
          if (lroot) then
            print*, 'initialize_selfgravity: potself has bcy=''p'', but the density is not'
            print*, '                        periodic! (you must set a proper boundary condition'
            print*, '                        for the potential)'
            print*, 'initialize_selfgravity: bcy=', bcy
          endif
          call fatal_error('initialize_selfgravity','')
        endif
        if (bcz(ipotself)=='p' .and. .not.(bcz(ilnrho)=='p')) then
          if (lroot) then
            print*, 'initialize_selfgravity: potself has bcz=''p'', but the density is not'
            print*, '                        periodic! (you must set a proper boundary condition'
            print*, '                        for the potential)'
            print*, 'initialize_selfgravity: bcz=', bcz
          endif
          call fatal_error('initialize_selfgravity','')
        endif
      endif
!
      if (lneutraldensity) then
        if (bcx(ipotself)=='p' .and. .not.(bcx(ilnrho)=='p')) then
          if (lroot) then
            print*, 'initialize_selfgravity: potself has bcx=''p'', but the density is not'
            print*, '                        periodic! (you must set a proper boundary condition'
            print*, '                        for the potential)'
            print*, 'initialize_selfgravity: bcx=', bcx
          endif
          call fatal_error('initialize_selfgravity','')
        endif
        if (bcy(ipotself)=='p' .and. .not.(bcy(ilnrho)=='p')) then
          if (lroot) then
            print*, 'initialize_selfgravity: potself has bcy=''p'', but the density is not'
            print*, '                        periodic! (you must set a proper boundary condition'
            print*, '                        for the potential)'
            print*, 'initialize_selfgravity: bcy=', bcy
          endif
          call fatal_error('initialize_selfgravity','')
        endif
        if (bcz(ipotself)=='p' .and. .not.(bcz(ilnrho)=='p')) then
          if (lroot) then
            print*, 'initialize_selfgravity: potself has bcz=''p'', but the density is not'
            print*, '                        periodic! (you must set a proper boundary condition'
            print*, '                        for the potential)'
            print*, 'initialize_selfgravity: bcz=', bcz
          endif
          call fatal_error('initialize_selfgravity','')
        endif
      endif
!
    endsubroutine initialize_selfgravity
!***********************************************************************
    subroutine pencil_criteria_selfgravity()
!
!  All pencils that the Selfgravity module depends on are specified here.
!
!  15-may-06/anders+jeff: adapted
!
      lpenc_requested(i_gpotself)=.true.
      if (idiag_gpoten/=0 .or. idiag_gpotenmxy/=0) then
        lpenc_diagnos(i_rho)=.true.
        lpenc_diagnos(i_potself)=.true.
      endif
!
      if (idiag_grgpm/=0 .or. idiag_grgzm/=0 .or. idiag_gpgzm/=0) then 
        lpenc_diagnos(i_pomx)=.true.
        lpenc_diagnos(i_pomy)=.true.
        lpenc_diagnos(i_phix)=.true.
        lpenc_diagnos(i_phiy)=.true.
      endif
!
    endsubroutine pencil_criteria_selfgravity
!***********************************************************************
    subroutine pencil_interdep_selfgravity(lpencil_in)
!
!  Interdependency among pencils from the Selfgravity module is specified here.
!
!  15-may-06/anders+jeff: adapted
!
      logical, dimension(npencils) :: lpencil_in
!
      call keep_compiler_quiet(lpencil_in)
!
    endsubroutine pencil_interdep_selfgravity
!***********************************************************************
    subroutine calc_pencils_selfgravity(f,p)
!
!  Calculate Selfgravity pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  15-may-06/anders+jeff: coded
!
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      intent(inout) :: f, p
!
      if (lpencil(i_potself)) p%potself = f(l1:l2,m,n,ipotself)
      if (lpencil(i_gpotself)) then
        call grad(f,ipotself,p%gpotself)
        if (igpotselfx/=0) f(l1:l2,m,n,igpotselfx:igpotselfz)=p%gpotself
      endif
!
    endsubroutine calc_pencils_selfgravity
!***********************************************************************
    subroutine calc_selfpotential(f)
!
!  Calculate the potential of the self gravity.
!
!  15-may-06/anders+jeff: coded
!
      use Particles_main, only: particles_calc_selfpotential
      use Poisson
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      real, dimension (nx,ny,nz) :: rhs_poisson
!
      if (t>=tstart_selfgrav) then 
!
!  Consider self-gravity from gas and dust density or from either one.
!
        if (lhydro.and.ldensity.and.lselfgravity_gas) then
          if (ldensity_nolog) then
            rhs_poisson=rhs_poisson_const*f(l1:l2,m1:m2,n1:n2,ilnrho)
          else
            rhs_poisson=rhs_poisson_const*exp(f(l1:l2,m1:m2,n1:n2,ilnrho))
          endif
        endif
!  Dust.
        if (ldustdensity.and.ldustvelocity.and.lselfgravity_dust) then
          if (ldustdensity_log) then
            if (lselfgravity_gas) then  ! No need to zero rhs.
              rhs_poisson = rhs_poisson + &
                  rhs_poisson_const*exp(f(l1:l2,m1:m2,n1:n2,ind(1)))
            else                        ! Must zero rhs.
              rhs_poisson = rhs_poisson_const*exp(f(l1:l2,m1:m2,n1:n2,ind(1)))
            endif
          else
            if (lselfgravity_gas) then  ! No need to zero rhs.
              rhs_poisson = rhs_poisson + &
                  rhs_poisson_const*f(l1:l2,m1:m2,n1:n2,ind(1))
            else                        ! Must zero rhs.
              rhs_poisson = rhs_poisson_const*f(l1:l2,m1:m2,n1:n2,ind(1))
            endif
          endif
        endif
!  Neutrals.
        if (lneutraldensity.and.lneutralvelocity.and.lselfgravity_neutrals) then
          if (lneutraldensity_log) then
            if (lselfgravity_gas.or.lselfgravity_dust) then  ! No need to zero rhs.
              rhs_poisson = rhs_poisson + &
                  rhs_poisson_const*exp(f(l1:l2,m1:m2,n1:n2,ilnrhon))
            else                        ! Must zero rhs.
              rhs_poisson = rhs_poisson_const*exp(f(l1:l2,m1:m2,n1:n2,ilnrhon))
            endif
          else
            if (lselfgravity_gas.or.lselfgravity_dust) then  ! No need to zero rhs.
              rhs_poisson = rhs_poisson + &
                  rhs_poisson_const*f(l1:l2,m1:m2,n1:n2,ilnrhon)
            else                        ! Must zero rhs.
              rhs_poisson = rhs_poisson_const*f(l1:l2,m1:m2,n1:n2,ilnrhon)
            endif
          endif
        endif
!
!  Contribution from particles is taken care of by the particle modules.
!
        if (lparticles) &
            call particles_calc_selfpotential(f,rhs_poisson,rhs_poisson_const, &
            lselfgravity_gas.or.lselfgravity_dust)
!
!  Send the right-hand-side of the Poisson equation to the Poisson solver and
!  receive the self-gravity potential back.
!
        call inverse_laplacian(f,rhs_poisson)
!
!  Put potential into f array.
!
        f(l1:l2,m1:m2,n1:n2,ipotself) = rhs_poisson
!
      endif ! if (t>=tstart_selfgrav) then
!
    endsubroutine calc_selfpotential
!***********************************************************************
    subroutine duu_dt_selfgrav(f,df,p)
!
!  Add self gravity acceleration on gas.
!
!  15-may-06/anders+jeff: coded
!
      use Diagnostics
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      intent(in) :: f,p
      intent(out) :: df
!
!  Add self-gravity acceleration on the gas and on the dust.
!
      if (t>=tstart_selfgrav) then
        if (lhydro.and.lselfgravity_gas) &
            df(l1:l2,m,n,iux:iuz) = df(l1:l2,m,n,iux:iuz) - p%gpotself
        if ( ldustvelocity.and.lselfgravity_dust) &
            df(l1:l2,m,n,iudx(1):iudz(1)) = df(l1:l2,m,n,iudx(1):iudz(1)) - p%gpotself
        if (lneutralvelocity.and.lselfgravity_neutrals) &
            df(l1:l2,m,n,iunx:iunz) = df(l1:l2,m,n,iunx:iunz) - p%gpotself
      endif
!
!  Diagnostic averages.
!
      if (ldiagnos) then
        if (idiag_gpoten/=0) call sum_mn_name(p%potself*p%rho,idiag_gpoten)
        if (idiag_gpotselfxm/=0) &
            call sum_mn_name(p%gpotself(:,1),idiag_gpotselfxm)
        if (idiag_gpotselfym/=0) &
            call sum_mn_name(p%gpotself(:,2),idiag_gpotselfym)
        if (idiag_gpotselfzm/=0) &
            call sum_mn_name(p%gpotself(:,3),idiag_gpotselfzm)
        if (idiag_gpotselfx2m/=0) &
            call sum_mn_name(p%gpotself(:,1)**2,idiag_gpotselfx2m)
        if (idiag_gpotselfy2m/=0) &
            call sum_mn_name(p%gpotself(:,2)**2,idiag_gpotselfy2m)
        if (idiag_gpotselfz2m/=0) &
            call sum_mn_name(p%gpotself(:,3)**2,idiag_gpotselfz2m)
        if (idiag_gxgym/=0) & 
             call sum_mn_name(p%gpotself(:,1)*p%gpotself(:,2),idiag_gxgym)
        if (idiag_gxgzm/=0) &
             call sum_mn_name(p%gpotself(:,1)*p%gpotself(:,3),idiag_gxgzm)
        if (idiag_gygzm/=0) &
             call sum_mn_name(p%gpotself(:,2)*p%gpotself(:,3),idiag_gygzm)
        if (idiag_grgpm/=0 .or. idiag_grgzm/=0 .or. idiag_gpgzm/=0) &
             call calc_cylgrav_stresses(p)
      endif
!
!  2-D averages.
!
      if (l2davgfirst) then
        if (idiag_gpotenmxy/=0) call zsum_mn_name_xy(p%potself*p%rho,idiag_gpotenmxy)
      endif
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)
!
    endsubroutine duu_dt_selfgrav
!***********************************************************************
    subroutine calc_cylgrav_stresses(p)
!
!  Calculates cylindrical gravitational stresses in a cartesian box
!
!  01-jul-07/wlad: coded 
! 
      use Diagnostics, only: sum_mn_name
      use Mpicomm, only: stop_it
!
      real, dimension(nx) :: gpotr,gpotp,gpotz
      type (pencil_case) :: p
!
      gpotr=p%gpotself(:,1)*p%pomx+p%gpotself(:,2)*p%pomy
      gpotp=p%gpotself(:,1)*p%phix+p%gpotself(:,2)*p%phiy
      gpotz=p%gpotself(:,3)
!
      if (idiag_grgpm/=0) call sum_mn_name(gpotr*gpotp,idiag_grgpm)
      if (idiag_grgzm/=0) call sum_mn_name(gpotr*gpotz,idiag_grgzm)
      if (idiag_gpgzm/=0) call sum_mn_name(gpotp*gpotz,idiag_gpgzm)
!     
    endsubroutine calc_cylgrav_stresses
!***********************************************************************
    subroutine read_selfgravity_init_pars(unit,iostat)
!
!  Read self gravity init parameters
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=selfgrav_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=selfgrav_init_pars,ERR=99)
      endif
!
99    return
    endsubroutine read_selfgravity_init_pars
!***********************************************************************
    subroutine write_selfgravity_init_pars(unit)
!
!  Write self gravity init parameters
!
      integer, intent(in) :: unit
!
      write(unit,NML=selfgrav_init_pars)
!
    endsubroutine write_selfgravity_init_pars
!***********************************************************************
    subroutine read_selfgravity_run_pars(unit,iostat)
!
!  Read self gravity run parameters
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=selfgrav_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=selfgrav_run_pars,ERR=99)
      endif
!
99    return
    endsubroutine read_selfgravity_run_pars
!***********************************************************************
    subroutine write_selfgravity_run_pars(unit)
!
!  Write self gravity run parameters
!
      integer, intent(in) :: unit
!
      write(unit,NML=selfgrav_run_pars)
!
    endsubroutine write_selfgravity_run_pars
!***********************************************************************
    subroutine rprint_selfgravity(lreset,lwrite)
!
!  reads and registers print parameters relevant for gravity advance
!  dummy routine
!
!  16-may-06/anders+jeff: adapted
!
      use Diagnostics
!
      logical :: lreset,lwr
      logical, optional :: lwrite
!
      integer :: iname, inamexy
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
      if (lreset) then
        idiag_gpoten=0; idiag_gpotenmxy=0
        idiag_gpotselfxm=0; idiag_gpotselfym=0; idiag_gpotselfzm=0
        idiag_gpotselfx2m=0; idiag_gpotselfy2m=0; idiag_gpotselfz2m=0
        idiag_gxgym=0; idiag_gxgzm=0; idiag_gygzm=0
        idiag_grgpm=0; idiag_grgzm=0; idiag_gpgzm=0
      endif
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'gpoten',idiag_gpoten)
        call parse_name(iname,cname(iname),cform(iname),'gpotselfxm', &
            idiag_gpotselfxm)
        call parse_name(iname,cname(iname),cform(iname),'gpotselfym', &
            idiag_gpotselfym)
        call parse_name(iname,cname(iname),cform(iname),'gpotselfzm', &
            idiag_gpotselfzm)
        call parse_name(iname,cname(iname),cform(iname),'gpotselfx2m', &
            idiag_gpotselfx2m)
        call parse_name(iname,cname(iname),cform(iname),'gpotselfy2m', &
            idiag_gpotselfy2m)
        call parse_name(iname,cname(iname),cform(iname),'gpotselfz2m', &
            idiag_gpotselfz2m)
        call parse_name(iname,cname(iname),cform(iname),'gxgym',idiag_gxgym)
        call parse_name(iname,cname(iname),cform(iname),'gxgzm',idiag_gxgzm)
        call parse_name(iname,cname(iname),cform(iname),'gygzm',idiag_gygzm)
        call parse_name(iname,cname(iname),cform(iname),'grgpm',idiag_grgpm)
        call parse_name(iname,cname(iname),cform(iname),'grgzm',idiag_grgzm)
        call parse_name(iname,cname(iname),cform(iname),'gpgzm',idiag_gpgzm)
      enddo
!
!  check for those quantities for which we want z-averages
!
      do inamexy=1,nnamexy
        call parse_name(inamexy,cnamexy(inamexy),cformxy(inamexy),'gpotenmxy', idiag_gpotenmxy)
      enddo
!
!  write column where which variable is stored
!
      if (lwr) then
        write(3,*) 'i_gpoten=',idiag_gpoten
        write(3,*) 'i_gpotenmxy=',idiag_gpotenmxy
        write(3,*) 'ipotself=', ipotself
      endif
!
    endsubroutine rprint_selfgravity
!***********************************************************************


endmodule Selfgravity
