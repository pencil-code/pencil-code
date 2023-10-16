! $Id$
!
!  This module takes care of everything related to particle self-gravity
!  by inverting the Poisson equation of the particle density mapped on the
!  grid.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
!
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lparticles_selfgravity=.true.
!
! MAUX CONTRIBUTION 3
! COMMUNICATED AUXILIARIES 3
!
!***************************************************************
module Particles_selfgravity
!
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
  use Particles_cdata
  use Particles_map
  use Particles_mpicomm
  use Particles_sub
!
  implicit none
!
  include 'particles_selfgravity.h'
!
  real :: cdtpg=0.2
  real, pointer :: tstart_selfgrav
  logical :: lselfgravity_particles=.true.
  logical :: lnopartingrav=.false.
!
  namelist /particles_selfgrav_init_pars/ &
      lselfgravity_particles,lnopartingrav
!
  namelist /particles_selfgrav_run_pars/ &
      lselfgravity_particles,lnopartingrav
!
  integer :: idiag_gpotenp=0, idiag_potselfpm=0
!
  contains
!***********************************************************************
    subroutine register_particles_selfgrav()
!
!  Set up indices for access to the f, fp and dfp arrays.
!
!  14-jun-06/anders: adapted
!
      use FArrayManager
      use Messages, only: svn_id
!
      if (lroot) call svn_id( &
           "$Id$")
!
!  Index for gradient for the self-potential and for the smooth particle
!  density field.
!
      call farray_register_auxiliary('gpotselfx',igpotselfx,communicated=.true.)
      call farray_register_auxiliary('gpotselfy',igpotselfy,communicated=.true.)
      call farray_register_auxiliary('gpotselfz',igpotselfz,communicated=.true.)
!
    endsubroutine register_particles_selfgrav
!***********************************************************************
    subroutine initialize_particles_selfgrav(f)
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  14-jun-06/anders: adapted
!
      use SharedVariables
!
      real, dimension (mx,my,mz,mfarray) :: f
!
!  For particle block domain decomposition we need to fill blocks with
!  information about the gravitational acceleration.
!
      if (lparticles_blocks) lfill_blocks_gpotself=.true.
!
!  Initialize gravitational acceleration to zero.
!
      f(:,:,:,igpotselfx:igpotselfz)=0.0
!
      call get_shared_variable('tstart_selfgrav',tstart_selfgrav,caller='initialize_particles_selfgrav')
!
!  Boundary condition consistency.
!
      if (any(bcx(igpotselfx:igpotselfz)=='p') .and. .not.(bcx(ipotself)=='p')) then
        if (lroot) then
          print*, 'initialize_particles_selfgrav: igpotself has bcx=''p'', but the potential is not'
          print*, '                               periodic! (you must set a proper boundary'
          print*, '                               condition for the gradient of the potential)'
          print*, 'initialize_particles_selfgrav: bcx=', bcx
        endif
        call fatal_error('initialize_particles_selfgrav','')
      endif
      if (any(bcy(igpotselfx:igpotselfz)=='p') .and. .not.(bcy(ipotself)=='p')) then
        if (lroot) then
          print*, 'initialize_particles_selfgrav: igpotself has bcy=''p'', but the potential is not'
          print*, '                               periodic! (you must set a proper boundary'
          print*, '                               condition for the gradient of the potential)'
          print*, 'initialize_particles_selfgrav: bcy=', bcy
        endif
        call fatal_error('initialize_particles_selfgrav','')
      endif
      if (any(bcz(igpotselfx:igpotselfz)=='p') .and. .not.(bcz(ipotself)=='p')) then
        if (lroot) then
          print*, 'initialize_particles_selfgrav: igpotself has bcz=''p'', but the potential is not'
          print*, '                               periodic! (you must set a proper boundary'
          print*, '                               condition for the gradient of the potential)'
          print*, 'initialize_particles_selfgrav: bcz=', bcz
        endif
        call fatal_error('initialize_particles_selfgrav','')
      endif
!
    endsubroutine initialize_particles_selfgrav
!***********************************************************************
    subroutine calc_selfpotential_particles(f,rhs_poisson,lcontinued)
!
!  Calculate the gravitational potential of the dust particles.
!
!  13-jun-06/anders: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx,ny,nz) :: rhs_poisson
      logical :: lcontinued
!
      if (t>=tstart_selfgrav) then
!
! DM what is this ugly hack that uses lselfgravity_particles ? 
        if (lnopartingrav) then
        else
          if (lselfgravity_particles) then
            if (lcontinued) then  ! Potential has already been zeroed by the gas.
              rhs_poisson = rhs_poisson + f(l1:l2,m1:m2,n1:n2,irhop)
            else                  ! Must zero potential from last time-step.
              rhs_poisson = f(l1:l2,m1:m2,n1:n2,irhop)
            endif
          endif
        endif
!
      endif  ! if (t>=tstart_selfgrav) then
!
    endsubroutine calc_selfpotential_particles
!***********************************************************************
    subroutine pencil_criteria_par_selfgrav()
!
!  All pencils that the Particles_selfgrav module depends on are specified here.
!
!  02-jul-06/anders: adapted
!
      lpenc_requested(i_gpotself)=.true.
      if (idiag_gpotenp/=0) then
        lpenc_diagnos(i_rhop)=.true.
        lpenc_diagnos(i_potself)=.true.
      endif
!
    endsubroutine pencil_criteria_par_selfgrav
!***********************************************************************
    subroutine pencil_interdep_par_selfgrav(lpencil_in)
!
!  Interdependency among pencils provided by the Particles_selfgrav module
!  is specified here.
!
!  02-jul-06/anders: adapted
!
      logical, dimension(npencils) :: lpencil_in
!
      call keep_compiler_quiet(lpencil_in)
!
    endsubroutine pencil_interdep_par_selfgrav
!***********************************************************************
    subroutine calc_pencils_par_selfgrav(f,p)
!
!  Calculate particle pencils.
!
!  02-jul-06/anders: adapted
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)
!
    endsubroutine calc_pencils_par_selfgrav
!***********************************************************************
    subroutine dvvp_dt_selfgrav_pencil(f,df,fp,dfp,p,ineargrid)
!
!  Add self-gravity to particle equation of motion.
!
!  14-jun-06/anders: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mparray) :: fp
      real, dimension (mpar_loc,mpvar) :: dfp
      type (pencil_case) :: p
      integer, dimension (mpar_loc,3) :: ineargrid
!
      intent (in) :: f, df, p, fp, dfp, ineargrid
!
      call calc_diagnostics_particles_selg(p)
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(dfp)
      call keep_compiler_quiet(ineargrid)
!
    endsubroutine dvvp_dt_selfgrav_pencil
!***********************************************************************
    subroutine calc_diagnostics_particles_selg(p)
!
      use Diagnostics

      type (pencil_case) :: p
!
      if (ldiagnos) then
        if (idiag_gpotenp/=0) call sum_mn_name(p%potself*p%rhop,idiag_gpotenp)
      endif

    endsubroutine calc_diagnostics_particles_selg
!***********************************************************************
    subroutine dvvp_dt_selfgrav(f,df,fp,dfp,ineargrid)
!
!  Add self-gravity to particle equation of motion.
!
!  14-jun-06/anders: coded
!
      use Messages, only: fatal_error
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mparray) :: fp
      real, dimension (mpar_loc,mpvar) :: dfp
      integer, dimension (mpar_loc,3) :: ineargrid
!
      real, dimension (3) :: gpotself
      real, dimension (1) :: potself
      real :: rhop_swarm_par
      integer :: k
      logical :: lheader, lfirstcall=.true.
!
      intent (in) :: f, fp
      intent (inout) :: dfp
!
      if (t>=tstart_selfgrav) then
!
!  Print out header information in first time step.
!
        lheader=lfirstcall .and. lroot
!
        lfirstcall=.false.
!
        if (lheader) print*, 'dvvp_dt_selfgrav: add self-gravity'
!
!  Interpolate the gradient of the potential to the location of the particle.
!
        if (lselfgravity_particles) then
          do k=1,npar_loc
            
              if (lparticlemesh_cic) then
                if (lparticles_blocks) then
                  call interpolate_linear(f,igpotselfx,igpotselfz, &
                      fp(k,ixp:izp),gpotself,ineargrid(k,:),inearblock(k),ipar(k))
                else
                  call interpolate_linear(f,igpotselfx,igpotselfz, &
                      fp(k,ixp:izp),gpotself,ineargrid(k,:),0,ipar(k))
                endif
              elseif (lparticlemesh_tsc) then
                if (linterpolate_spline) then
                  if (lparticles_blocks) then
                    call interpolate_quadratic_spline(f,igpotselfx,igpotselfz, &
                        fp(k,ixp:izp),gpotself,ineargrid(k,:),inearblock(k),ipar(k))
                  else
                    call interpolate_quadratic_spline(f,igpotselfx,igpotselfz, &
                        fp(k,ixp:izp),gpotself,ineargrid(k,:),0,ipar(k))
                  endif
                else
!
!  Polynomial interpolation can lead to self acceleration of isolated particle,
!  since one must use the same assignment and interpolation function.
!
                  call fatal_error('dvvp_dt_selfgrav','polynomial interpolation not allowed'// &
                                   'for interpolating self-gravity')
                endif
              else
                gpotself=f(ineargrid(k,1),ineargrid(k,2),ineargrid(k,3),igpotselfx:igpotselfz)
              endif
!
!  Add gravitational acceleration to particles
!
            dfp(k,ivpx:ivpz)=dfp(k,ivpx:ivpz)-gpotself
!
!  Time-step constraint from particle self-gravity. We require that the
!  gravitational acceleration must at most move a particle a fraction of
!  a grid cell in one time-step.
!
            if (lfirst.and.ldt) then
              if (lparticles_blocks) then
                dt1_max(ineargrid(k,1)-nghostb)= &
                    max(dt1_max(ineargrid(k,1)-nghostb), &
                    sqrt(0.5*abs(gpotself(1))*dx_1(ineargrid(k,1)))/cdtpg)
                dt1_max(ineargrid(k,1)-nghostb)= &
                    max(dt1_max(ineargrid(k,1)-nghostb), &
                    sqrt(0.5*abs(gpotself(2))*dy_1(ineargrid(k,2)))/cdtpg)
                dt1_max(ineargrid(k,1)-nghostb)= &
                    max(dt1_max(ineargrid(k,1)-nghostb), &
                    sqrt(0.5*abs(gpotself(3))*dz_1(ineargrid(k,3)))/cdtpg)
              else
                dt1_max(ineargrid(k,1)-nghost)= &
                    max(dt1_max(ineargrid(k,1)-nghost), &
                    sqrt(0.5*abs(gpotself(1))*dx_1(ineargrid(k,1)))/cdtpg)
                dt1_max(ineargrid(k,1)-nghost)= &
                    max(dt1_max(ineargrid(k,1)-nghost), &
                    sqrt(0.5*abs(gpotself(2))*dy_1(ineargrid(k,2)))/cdtpg)
                dt1_max(ineargrid(k,1)-nghost)= &
                    max(dt1_max(ineargrid(k,1)-nghost), &
                    sqrt(0.5*abs(gpotself(3))*dz_1(ineargrid(k,3)))/cdtpg)
              endif
            endif
!
            if (ldiagnos) then
              if (idiag_potselfpm/=0) then
                if (lparticlemesh_cic) then
                  call interpolate_linear(f,ipotself,ipotself, &
                      fp(k,ixp:izp),potself,ineargrid(k,:),0,ipar(k))
                elseif (lparticlemesh_tsc) then
                  call interpolate_quadratic_spline(f,ipotself,ipotself, &
                      fp(k,ixp:izp),potself,ineargrid(k,:),0,ipar(k))
                else
                  potself=f(ineargrid(k,1),ineargrid(k,2),ineargrid(k,3),ipotself)
                endif
                if (lparticles_density) then
                  call sum_par_name(potself*fp(k,irhopswarm),idiag_potselfpm)
                else
                  call get_rhopswarm(mp_swarm,fp,k,ineargrid(k,:),rhop_swarm_par)
                  call sum_par_name(potself*rhop_swarm_par,idiag_potselfpm)
                endif
              endif
            endif
!
          enddo
!         
!  Patch to call sum_par_name from every processor
!  (even those that do not have any particle). 
!
          if (ldiagnos .and. idiag_potselfpm /= 0 .and. npar_loc <= 0) &
              call sum_par_name(fp(1:npar_loc,ixp), idiag_potselfpm)
!
        endif
!
      endif ! if (t>=tstart_selfgrav) then
!
      call keep_compiler_quiet(df)
!
    endsubroutine dvvp_dt_selfgrav
!***********************************************************************
    subroutine read_particles_selfg_init_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=particles_selfgrav_init_pars, IOSTAT=iostat)
!
    endsubroutine read_particles_selfg_init_pars
!***********************************************************************
    subroutine write_particles_selfg_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=particles_selfgrav_init_pars)
!
    endsubroutine write_particles_selfg_init_pars
!***********************************************************************
    subroutine read_particles_selfg_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=particles_selfgrav_run_pars, IOSTAT=iostat)
!
    endsubroutine read_particles_selfg_run_pars
!***********************************************************************
    subroutine write_particles_selfg_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=particles_selfgrav_run_pars)
!
    endsubroutine write_particles_selfg_run_pars
!***********************************************************************
    subroutine rprint_particles_selfgrav(lreset,lwrite)
!
!  Read and register print parameters relevant for particle self-gravity.
!
!  14-jun-06/anders: adapted
!
      use Diagnostics
!
      logical :: lreset
      logical, optional :: lwrite
!
      integer :: iname
!
      if (lreset) then
        idiag_gpotenp=0; idiag_potselfpm=0
      endif
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'gpotenp',idiag_gpotenp)
        call parse_name(iname,cname(iname),cform(iname),'potselfpm',idiag_potselfpm)
      enddo
!
    endsubroutine rprint_particles_selfgrav
!***********************************************************************
endmodule Particles_selfgravity
