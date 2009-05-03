! $Id$
!
!  This module takes care of everything related to particle self-gravity.
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
  use Messages
  use Particles_cdata
  use Particles_sub
!
  implicit none
!
  include 'particles_selfgravity.h'
!
  real, pointer :: tstart_selfgrav
  logical :: lselfgravity_particles=.true.
  logical :: lselfgravity_nbodyparticles=.false.
!
  namelist /particles_selfgrav_init_pars/ &
      lselfgravity_particles,lselfgravity_nbodyparticles
!
  namelist /particles_selfgrav_run_pars/ &
      lselfgravity_particles,lselfgravity_nbodyparticles
!
  integer :: idiag_gpotenp=0
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
      use Messages, only: cvs_id
!
      if (lroot) call cvs_id( &
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
    subroutine initialize_particles_selfgrav(lstarting)
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  14-jun-06/anders: adapted
!
      use SharedVariables
!
      integer :: ierr
      logical :: lstarting
!
      if (NO_WARN) print*, lstarting
!
      call get_shared_variable('tstart_selfgrav',tstart_selfgrav,ierr)
      if (ierr/=0) then
        if (lroot) print*, 'initialize_particles_selfgrav: '// &
            'there was a problem when getting tstart_selfgrav!'
        call fatal_error('initialize_particles_selfgrav','')
      endif
!
! Boundary condition consistency
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
    subroutine calc_selfpotential_particles(f,rhs_poisson,rhs_poisson_const,    lcontinued)
!
!  Calculate the gravity potential of the dust particles.
!
!  13-jun-06/anders: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx,ny,nz) :: rhs_poisson
      real :: rhs_poisson_const
      logical :: lcontinued
!
      if (t>=tstart_selfgrav) then
!
        if (lselfgravity_particles) then
          if (lcontinued) then  ! Potential has already been zeroed by the gas.
            rhs_poisson=rhs_poisson+ &
                rhs_poisson_const*f(l1:l2,m1:m2,n1:n2,irhop)
          else                  ! Must zero potential from last time-step.
            rhs_poisson= &
                rhs_poisson_const*f(l1:l2,m1:m2,n1:n2,irhop)
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
      if (NO_WARN) print*, lpencil_in
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
      if (NO_WARN) print*, f, p
!
    endsubroutine calc_pencils_par_selfgrav
!***********************************************************************
    subroutine dvvp_dt_selfgrav_pencil(f,df,fp,dfp,p,ineargrid)
!
!  Add self-gravity to particle equation of motion.
!
!  14-jun-06/anders: coded
!
      use Diagnostics
      use Messages, only: fatal_error
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mpvar) :: fp, dfp
      type (pencil_case) :: p
      integer, dimension (mpar_loc,3) :: ineargrid
!
      intent (in) :: f, df, p, fp, dfp, ineargrid
!
      if (ldiagnos) then
        if (idiag_gpotenp/=0) call sum_mn_name(p%potself*p%rhop,idiag_gpotenp)
      endif
!
      if (NO_WARN) print*, f, df, fp, dfp, p, ineargrid
!
    endsubroutine dvvp_dt_selfgrav_pencil
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
      real, dimension (mpar_loc,mpvar) :: fp, dfp
      integer, dimension (mpar_loc,3) :: ineargrid
!
      real, dimension (3) :: gpotself
      integer :: k
      logical :: lheader, lnbody, lfirstcall=.true.
!
      intent (in) :: f, fp
      intent (out) :: dfp
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
              call interpolate_linear(f, igpotselfx, igpotselfz, &
                  fp(k,ixp:izp), gpotself, ineargrid(k,:), ipar(k))
            elseif (lparticlemesh_tsc) then
              if (linterpolate_spline) then
                call interpolate_quadratic_spline(f, igpotselfx, igpotselfz, &
                    fp(k,ixp:izp), gpotself, ineargrid(k,:), ipar(k) )
              else
!
!  Polynomial interpolation can lead to self acceleration of isolated particle,
!  since one must use the same assignment and interpolation function.
!                
                if (lroot) then
                  print*, 'dvvp_dt_selfgrav: polynomial interpolation not '// &
                          'allowed'
                  print*, '                  for interpolating self-gravity'
                endif
                call fatal_error('dvvp_dt_selfgrav','')
              endif
            else
              gpotself=f(ineargrid(k,1),ineargrid(k,2),ineargrid(k,3),igpotselfx:igpotselfz)
            endif
!            
!  Possibility of switching off self-gravity for the
!  massive n-body particles
!
            if (lparticles_nbody.and.(.not.lselfgravity_nbodyparticles)) then 
              lnbody=(lparticles_nbody.and.any(ipar(k).eq.ipar_nbody))
              if (lnbody) gpotself=0
            endif
!
            dfp(k,ivpx:ivpz)=dfp(k,ivpx:ivpz)-gpotself
!
          enddo
        endif
!
      endif ! if (t>=tstart_selfgrav) then
!
      if (NO_WARN) print*, df
!
    endsubroutine dvvp_dt_selfgrav
!***********************************************************************
    subroutine read_particles_selfg_init_pars(unit,iostat)
!
      integer, intent (in) :: unit
      integer, intent (inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=particles_selfgrav_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=particles_selfgrav_init_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_particles_selfg_init_pars
!***********************************************************************
    subroutine write_particles_selfg_init_pars(unit)
!
      integer, intent (in) :: unit
!
      write(unit,NML=particles_selfgrav_init_pars)
!
    endsubroutine write_particles_selfg_init_pars
!***********************************************************************
    subroutine read_particles_selfg_run_pars(unit,iostat)
!
      integer, intent (in) :: unit
      integer, intent (inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=particles_selfgrav_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=particles_selfgrav_run_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_particles_selfg_run_pars
!***********************************************************************
    subroutine write_particles_selfg_run_pars(unit)
!
      integer, intent (in) :: unit
!
      write(unit,NML=particles_selfgrav_run_pars)
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
      logical :: lwr
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
      if (lreset) then
        idiag_gpotenp = 0
      endif
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'gpotenp',idiag_gpotenp)
      enddo
!
!  Write information to index.pro
!
      if (lwr) then
        write(3,*) 'igpotselfx=', igpotselfx
        write(3,*) 'igpotselfy=', igpotselfy
        write(3,*) 'igpotselfz=', igpotselfz
      endif
!
    endsubroutine rprint_particles_selfgrav
!***********************************************************************

endmodule Particles_selfgravity
