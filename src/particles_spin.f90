! $Id$
!
!  This module takes care of everything related to particle spin
!  including lifting forces. The module maintains a full f-array
!  vorticity field, to be able to interpolate on the flow vorticity.
!
!  The module should be considered experimental as it is virtually
!  untested (as of aug-08).
!
!  NOTE: all code relating to particle spin or the magnus force
!        have been commented out.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
!
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! MAUX CONTRIBUTION 3
! COMMUNICATED AUXILIARIES 3
!   disabled for now: MPVAR_CONTRIBUTION_3
! CPARAM logical, parameter :: lparticles_spin=.true.
!
!***************************************************************
module Particles_spin

  use Cdata
  use Particles_cdata
  use Particles_sub

  implicit none

  include 'particles_spin.h'

  logical :: lsaffman_lift=.false.
  logical :: lmagnus_lift=.false.
  character (len=labellen), dimension(ninit) :: initsp='nothing'

  namelist /particles_spin_init_pars/ &
    lsaffman_lift, lmagnus_lift, initsp

  namelist /particles_spin_run_pars/ &
    lsaffman_lift, lmagnus_lift

  integer :: idiag_psxm=0, idiag_psym=0, idiag_pszm=0

  contains

!***********************************************************************
    subroutine register_particles_spin()
!
!  Set up indices for access to the fp and dfp arrays
!
!  21-jul-08/kapelrud: coded
!
      use Cdata
      use FArrayManager
      use Mpicomm, only: stop_it
      use Messages, only: svn_id
      use Cparam
!
      if (lroot) call svn_id( &
           "$Id$")
!
!  Indices for flow field vorticity. The vorticity is a communicated auxiliary
!  vector.
!  Assuming that this module is the last to use communicated aux variables,
!  then the three last entries in bc{x,y,z} in start.in sets the boundary
!  conditions of the vorticity.
!
      call farray_register_auxiliary('ox',iox,communicated=.true.)
      call farray_register_auxiliary('oy',ioy,communicated=.true.)
      call farray_register_auxiliary('oz',ioz,communicated=.true.)
!
!  Indices for particle spin
!
!     ipsx=npvar+1
!     ipsy=npvar+2
!     ipsz=npvar+3
!
!     npvar=npvar+3
!
!  Check that the fp and dfp arrays are big enough.
!
      if (npvar > mpvar) then
        if (lroot) write(0,*) 'npvar = ', npvar, ', mpvar = ', mpvar
        call stop_it('register_particles_spin: npvar > mpvar')
      endif
!
!  Make sure that the vorticity field is communicated one time extra
!  before the pencil loop in pde is executed.
!
      lparticles_prepencil_calc=.true.
!
    endsubroutine register_particles_spin
!***********************************************************************
    subroutine initialize_particles_spin(f,lstarting)
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  21-jul-08/kapelrud: coded
!
      use Particles_radius
!
      real, dimension (mx,my,mz,mfarray) :: f
      logical :: lstarting
!
!  Initialize vorticity field to zero.
!
      f(:,:,:,iox)=0.0
      f(:,:,:,ioy)=0.0
      f(:,:,:,ioz)=0.0
!
      if (lroot) print*,'initialize_particles_spin: '// &
          'using communicated auxiliary variables.'
!
!  Request interpolation of variables:
!
      interp%luu=interp%luu.or.lsaffman_lift !.or.lmagnus_lift
      interp%loo=interp%loo.or.lsaffman_lift !.or.lmagnus_lift
      interp%lrho=interp%lrho.or.lsaffman_lift !.or.lmagnus_lift
!
    endsubroutine initialize_particles_spin
!***********************************************************************
    subroutine init_particles_spin(f,fp)
!
!  Initial spin of particles.
!
!  21-jul-08/kapelrud: coded
!
      use Cdata, only: m,n
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mpar_loc,mpvar) :: fp
!
      intent(inout) :: f
!     integer :: j
!
      if (lroot) print*,'init_particles_spin: oo will be generated '// &
          'at run-time.'
!
!     do j=1,ninit
!       select case (initsp(j))
!
!       case ('nothing')
!         if (lroot) print*,'init_particles_spin: nothing (setting initial '// &
!           'particle spin to zero for all particles'
!         fp(1:npar_loc,ipsx:ipsz)=0.0
!
!       endselect
!     enddo
!       
    endsubroutine init_particles_spin
!***********************************************************************
    subroutine particles_spin_prepencil_calc(f)
!
!  Prepare the curl(uu) field here so that ghost zones can be communicated between
!  processors before the spin is calculated in dps_dt_pencil.
!
!  22-jul-08/kapelrud: coded
!
      use Cdata
      use Sub, only: curl
!
      real,dimension(mx,my,mz,mfarray) :: f
!
      intent(inout) :: f
!
!  Calculate curl(uu) along pencils in the internal region of this
!  processor's grid. Ghost zones will have to be set by the boundary conditions
!  and mpi communication as usual.
!
      do m=m1,m2;do n=n1,n2
        call curl(f,iux,f(l1:l2,m,n,iox:ioz))
      enddo;enddo
!
    endsubroutine particles_spin_prepencil_calc
!***********************************************************************
    subroutine pencil_criteria_par_spin()
!
!  All pencils that the Particles_spin module depends on are specified here.
!
!  21-nov-06/anders: coded
!
      use Cdata
!
      !lpenc_requested(i_uu)=.true.
!
    endsubroutine pencil_criteria_par_spin
!***********************************************************************
    subroutine dps_dt_pencil(f,df,fp,dfp,p,ineargrid)
!
!  Evolution of particle spin.
!
!  22-aug-05/anders: coded
!
      use Cdata
      use Viscosity, only: getnu
      use Messages, only: fatal_error
      use Particles_cdata
      use Particles_radius
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mpvar) :: fp, dfp
      type (pencil_case) :: p
      integer, dimension (mpar_loc,3) :: ineargrid
!
      logical :: lheader, lfirstcall=.true.
      real,dimension(3) :: tau
      real :: ip_tilde,nu
      integer :: k
!
      intent (in) :: f,df,fp,ineargrid
      intent (inout) :: dfp
!
      call getnu(nu_imput=nu)
!
!  Print out header information in first time step.
!
      lheader=lfirstcall .and. lroot
      lfirstcall=.false.
!
!  Identify module and boundary conditions.
!
      if (lheader) print*,'dps_dt_pencil: Calculate dps_dt'
!
!  Calculate torque on particle due to the shear flow, and
!  update the particles' spin.
!
!     if (lmagnus_lift) then
!       do k=k1_imn(imn),k2_imn(imn)
!
!  Calculate angular momentum
!
!         ip_tilde=0.4*mpmat*fp(k,iap)**2
!
!         tau=8.0*pi*interp_rho(k)*nu*fp(k,iap)**3* &
!             (0.5*interp_oo(k,:)-fp(k,ipsx:ipsz))
!         dfp(k,ipsx:ipsz)=dfp(k,ipsx:ipsz)+tau/ip_tilde
!       enddo
!     endif
!
    endsubroutine dps_dt_pencil
!***********************************************************************
    subroutine dps_dt(f,df,fp,dfp,ineargrid)
!
!  Evolution of particle spin.
!
!  25-jul-08/kapelrud: coded
!
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mpvar) :: fp, dfp
      integer, dimension (mpar_loc,3) :: ineargrid
!
!  Diagnostic output
!
!      if (ldiagnos) then
!        if (idiag_psxm/=0) call sum_par_name(fp(1:npar_loc,ipsx),idiag_psxm)
!        if (idiag_psym/=0) call sum_par_name(fp(1:npar_loc,ipsy),idiag_psym)
!        if (idiag_pszm/=0) call sum_par_name(fp(1:npar_loc,ipsz),idiag_pszm)
!      endif
!
    endsubroutine dps_dt
!***********************************************************************
    subroutine read_particles_spin_init_pars(unit,iostat)
!
      integer, intent (in) :: unit
      integer, intent (inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=particles_spin_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=particles_spin_init_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_particles_spin_init_pars
!***********************************************************************
    subroutine write_particles_spin_init_pars(unit)
!
      integer, intent (in) :: unit
!
      write(unit,NML=particles_spin_init_pars)
!
    endsubroutine write_particles_spin_init_pars
!***********************************************************************
    subroutine read_particles_spin_run_pars(unit,iostat)
!
      integer, intent (in) :: unit
      integer, intent (inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=particles_spin_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=particles_spin_run_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_particles_spin_run_pars
!***********************************************************************
    subroutine write_particles_spin_run_pars(unit)
!
      integer, intent (in) :: unit
!
      write(unit,NML=particles_spin_run_pars)
!
    endsubroutine write_particles_spin_run_pars
!***********************************************************************
    subroutine rprint_particles_spin(lreset,lwrite)
!
!  Read and register print parameters relevant for particles spin.
!
!  21-jul-08/kapelrud: adapted from particles_radius
!
      use Cdata
      use Diagnostics
!
      logical :: lreset
      logical, optional :: lwrite
!
      logical :: lwr
!
!  Write information to index.pro
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
      if (lwr) write(3,*) 'iox=', iox
!
!  Reset everything in case of reset
!
      if (lreset) then
!        idiag_psxm=0; idiag_psym=0; idiag_pszm=0
      endif
!
!  Run through all possible names that may be listed in print.in
!
!      if (lroot.and.ip<14) &
!          print*, 'rprint_particles_spin: run through parse list'
!
    endsubroutine rprint_particles_spin
!***********************************************************************
    subroutine calc_liftforce(fp,k,rep,liftforce)
!
!  Calculate lifting forces for a given particle. It should be possible to make
!  this a routine operating on pencils.
!
!  22-jul-08/kapelrud: coded
!
      real,dimension(mpvar) :: fp
      integer :: k
      real,dimension(3) :: liftforce
      real :: rep
!
      intent(in) :: fp, k, rep
      intent(out) :: liftforce
!
      real,dimension(3) :: dlift
!
      liftforce=0.0
      if (lsaffman_lift) then
        call calc_saffman_liftforce(fp,k,rep,dlift)
        liftforce=liftforce+dlift
      endif
!     if (lmagnus_lift) then
!       call calc_magnus_liftforce(fp,k,rep,dlift)
!       liftforce=liftforce+dlift
!       endif
!     endif
!
    endsubroutine calc_liftforce
!***********************************************************************
    subroutine calc_saffman_liftforce(fp,k,rep,dlift)
!
!  Calculate the Saffman lifting force for a given particles.
!
!  16-jul-08/kapelrud: coded
!
      use Particles_cdata
      use Sub, only: cross
      use Viscosity, only: getnu
!
      real,dimension(mpvar) :: fp
      integer :: k
      real,dimension(3) :: dlift
      real :: rep
!
      intent(in) :: fp, k, rep
      intent(out) :: dlift
!
      real :: csaff,diameter,beta,oo,nu
!
      call getnu(nu_imput=nu)
!
      if (.not.lparticles_radius) then
        if (lroot) print*,'calc_saffman_liftforce: '//&
             'Particle_radius module must be enabled!'
        call fatal_error('calc_saffman_liftforce','')
      endif
!
      diameter=2*fp(iap)
      oo=sqrt(sum(interp_oo(k,:)**2))
!
      beta=diameter**2*oo/(2.0*rep*nu)
      if (beta<0.005) then
        beta=0.005
      elseif (beta>0.4) then
        beta=0.4
      endif
!
      if (rep<=40) then
        csaff=(1-0.3314*beta**0.5)*exp(-rep/10.0)+0.3314*beta**0.5
      else
        csaff=0.0524*(beta*rep)**0.5
      endif
!
      call cross(interp_uu(k,:)-fp(ivpx:ivpz),interp_oo(k,:),dlift)
      dlift=1.61*csaff*diameter**2*nu**0.5*&
                 interp_rho(k)*oo**(-0.5)*dlift/mpmat
!
    endsubroutine calc_saffman_liftforce
!***********************************************************************
    subroutine calc_magnus_liftforce(fp,k,rep,dlift)
!
!  Calculate the Magnus liftforce for a given spinning particle.
!
!  22-jul-08/kapelrud: coded
!
      use Sub, only: cross
      use Particles_cdata
      use Viscosity, only: getnu
!
      real,dimension(mpvar) :: fp
      integer :: k
      real,dimension(3) :: dlift
      real :: rep
!
      intent(in) :: fp, k, rep
      intent(out) :: dlift
!
      real :: const_lr, spin_omega, area, nu
      real, dimension(3) :: ps_rel,uu_rel
!
      if (.not.lparticles_radius) then
        if (lroot) print*,'calc_magnus_liftforce: '//&
             'Particle_radius module must be enabled!'
        call fatal_error('calc_magnus_liftforce','')
      endif
!
      call getnu(nu_imput=nu)
!
!  Projected area of the particle
!
      area=pi*fp(iap)**2
!
!  Calculate the Magnus lift coefficent
!
      uu_rel=interp_uu(k,:)-fp(ivpx:ivpz)
      spin_omega=fp(iap)*sqrt(sum(fp(ipsx:ipsz)**2))/sqrt(sum(uu_rel**2))
      const_lr=min(0.5,0.5*spin_omega)
!
      ps_rel=fp(ipsx:ipsz)-0.5*interp_oo(k,:)
      call cross(uu_rel,ps_rel,dlift)
      dlift=dlift/sqrt(sum(ps_rel**2))
      dlift=0.25*interp_rho(k)*(rep*nu/fp(iap))*const_lr*area/mpmat*dlift
!
    endsubroutine calc_magnus_liftforce
!***********************************************************************

endmodule Particles_spin
