! $Id: particles_tracers.f90,v 1.3 2005-06-27 00:14:19 mee Exp $
!
!  This module takes care of everything related to tracer particles.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
!
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! MPVAR CONTRIBUTION 3
! CPARAM logical, parameter :: lparticles=.true.
!
!***************************************************************

module Particles

  use Cdata
  use Particles_sub

  implicit none

  include 'particles.h'

  real, dimension (npar,mpvar) :: fp, dfp
  integer, dimension (npar) :: ipar,ipar0
  integer :: npar_loc,npvar

  real :: xp0=0.0, yp0=0.0, zp0=0.0, vpx0=0.0, vpy0=0.0, vpz0=0.0
  real :: tausp=1.0
  character (len=labellen) :: initxxp='origin', initvvp='zero'

  namelist /particles_init_pars/ &
      initxxp, initvvp, xp0, yp0, zp0, vpx0, vpy0, vpz0, tausp

  namelist /particles_run_pars/ &
      tausp

  integer :: idiag_xpm=0, idiag_ypm=0, idiag_zpm=0
  integer :: idiag_vpxm=0, idiag_vpym=0, idiag_vpzm=0

  contains

!***********************************************************************
    subroutine particles_register_modules()
!
!  Register particle modules.
!
!  07-jan-05/anders: coded
!
      call register_particles()
!
    endsubroutine particles_register_modules
!***********************************************************************
    subroutine particles_rprint_list(lreset)
!
!  Read names of diagnostic particle variables to print out during run.
!
!  07-jan-05/anders: coded
!
      logical :: lreset
!
      if (lroot) open(3, file=trim(datadir)//'/index.pro', &
          STATUS='old', POSITION='append')
      call rprint_particles(lreset,LWRITE=lroot)
      if (lroot) close(3)
!
    endsubroutine particles_rprint_list
!***********************************************************************
    subroutine particles_initialize_modules(lstarting)
!
!  Initialize particle modules.
!
!  07-jan-05/anders: coded
!
      logical :: lstarting
!
      call initialize_particles(lstarting)
!
    endsubroutine particles_initialize_modules
!***********************************************************************
    subroutine particles_init(f)
!
!  Set up initial conditios for particle modules.
!
!  07-jan-05/anders: coded
!
      real, dimension (mx,my,mz,mvar+maux) :: f
!
      intent (in) :: f
!
      call init_particles(f)
!
    endsubroutine particles_init
!***********************************************************************
    subroutine particles_read_snapshot(filename)
!
!  Read particle snapshot from file.
!
!  07-jan-05/anders: coded
!
      character (len=*) :: filename
!
      call input_particles(filename,fp,ipar0)
!
    endsubroutine particles_read_snapshot
!***********************************************************************
    subroutine particles_write_snapshot(chsnap,msnap,enum,flist)
!
!  Write particle snapshot to file.
!
!  07-jan-05/anders: coded
!
      integer :: msnap
      logical :: enum
      character (len=*) :: chsnap,flist
      optional :: flist
!
      logical :: lsnap
!
      if (present(flist)) then
        call wsnap_particles_mask(chsnap,fp,msnap,enum,lsnap, &
            npar_loc,ipar,ipar0,flist)
      else
        call wsnap_particles_mask(chsnap,fp,msnap,enum,lsnap, &
            npar_loc,ipar,ipar0)
      endif
!
    endsubroutine particles_write_snapshot
!***********************************************************************
    subroutine particles_write_pdim(filename)
!   
!  Write npar and mpvar to file.
!
!  09-jan-05/anders: coded
!
      character (len=*) :: filename
!
      open(1,file=filename)
        write(1,'(2i7)') npar, mpvar
      close(1)
!
    endsubroutine particles_write_pdim
!***********************************************************************
    subroutine particles_timestep_first()
!
!  Setup dfp in the beginning of each itsub.
!
!  07-jan-05/anders: coded
!
      if (itsub==1) then
        dfp(ipar(1:npar_loc),:)=0.
      else
        dfp(ipar(1:npar_loc),:)=alpha(itsub)*dfp(ipar(1:npar_loc),:)
      endif
!
    endsubroutine particles_timestep_first
    subroutine particles_timestep_second()
!
!  Time evolution of particle variables.
!
!  07-jan-05/anders: coded
!
      fp(ipar(1:npar_loc),:) = &
          fp(ipar(1:npar_loc),:) + dt_beta(itsub)*dfp(ipar(1:npar_loc),:)
!
    endsubroutine particles_timestep_second
!***********************************************************************
    subroutine particles_pde(f)
!
!  07-jan-05/anders: coded
!
      real, dimension (mx,my,mz,mvar+maux) :: f
!
      intent (in) :: f
!
      call boundconds_particles_mask(fp,npar_loc,ipar,ipar0,dfp=dfp)
      call dxxp_dt(f,fp,dfp)
!
    endsubroutine particles_pde
!***********************************************************************
    subroutine register_particles()
!
!  Set up indices for access to the fp and dfp arrays
!
!  29-dec-04/anders: coded
!
      use Mpicomm, only: stop_it
      use Sub, only: cvs_id
!
      logical, save :: first=.true.
!
      if (.not. first) call stop_it('register_particles: called twice')
      first = .false.
!
      if (lroot) call cvs_id( &
          "$Id: particles_tracers.f90,v 1.3 2005-06-27 00:14:19 mee Exp $")
!
!  Indices for particle position.
!
      ixxp=1
      ixp=1
      iyp=2
      izp=3
!
!  Increase npvar accordingly.
!
      npvar=npvar+3
!
!  Check that the fp and dfp arrays are big enough.
!
      if (npvar > mpvar) then
        if (lroot) write(0,*) 'npvar = ', npvar, ', mpvar = ', mpvar
        call stop_it('register_particles: npvar > mpvar')
      endif
!
    endsubroutine register_particles
!***********************************************************************
    subroutine initialize_particles(lstarting)
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  29-dec-04/anders: coded
!
      logical :: lstarting
!
!  Distribute particles evenly among processors to begin with.
!
      if (lstarting) call dist_particles_evenly_mask(npar_loc,ipar,ipar0)
!
!  Size of box at local processor is needed for particle boundary conditions.
!
      Lxyz_loc(1)=Lxyz(1)/nprocx
      Lxyz_loc(2)=Lxyz(2)/nprocy
      Lxyz_loc(3)=Lxyz(3)/nprocz
      xyz0_loc(1)=xyz0(1)
      xyz0_loc(2)=xyz0(2)+ipy*Lxyz_loc(2)
      xyz0_loc(3)=xyz0(3)+ipz*Lxyz_loc(3)
      xyz1_loc(1)=xyz0_loc(1)+Lxyz_loc(1)
      xyz1_loc(2)=xyz0_loc(2)+Lxyz_loc(2)
      xyz1_loc(3)=xyz0_loc(3)+Lxyz_loc(3)
!
    endsubroutine initialize_particles
!***********************************************************************
    subroutine init_particles(f)
!
!  Initial conditions for tracer particles
!
!  29-dec-04/anders: coded
!
      use Cdata
      use General, only: random_number_wrapper
      use Mpicomm, only: stop_it
!
      real, dimension(mx,my,mz,mvar+maux) :: f
!
      integer :: k
!
      intent(in) :: f
!
!  Initial particle position.
!
      select case(initxxp)

      case ('origin')
        print*, 'init_particles: All particles at origin'
        fp(ipar(1:npar_loc),ixp:izp)=0.

      case ('constant')
        print*, 'init_particles: All particles at x,y,z=', xp0, yp0, zp0
        fp(ipar(1:npar_loc),ixp)=xp0
        fp(ipar(1:npar_loc),iyp)=yp0
        fp(ipar(1:npar_loc),izp)=zp0

      case ('random')
        if (lroot) print*, 'init_particles: Random particle positions'
        do k=1,npar_loc
          call random_number_wrapper(fp(ipar(k),ixp))
          call random_number_wrapper(fp(ipar(k),iyp))
          call random_number_wrapper(fp(ipar(k),izp))
        enddo
        fp(ipar(1:npar_loc),ixp)=xyz0(1)+fp(ipar(1:npar_loc),ixp)*Lxyz(1)
        fp(ipar(1:npar_loc),iyp)=xyz0(2)+fp(ipar(1:npar_loc),iyp)*Lxyz(2)
        fp(ipar(1:npar_loc),izp)=xyz0(3)+fp(ipar(1:npar_loc),izp)*Lxyz(3)

      case default
        print*, 'init_particles: No such such value for initxxp: ', &
            trim(initxxp)
        call stop_it("")

      endselect
!
    endsubroutine init_particles
!***********************************************************************
    subroutine dxxp_dt(f,fp,dfp)
!
!  Evolution of tracer particle position.
!
!  29-dec-04/anders: coded
!
      use Cdata
      use Mpicomm, only: stop_it
      use Sub
!
      real, dimension(mx,my,mz,mvar+maux) :: f
      real, dimension(npar,mpvar) :: fp,dfp
!
      real, dimension(3) :: uup
      integer :: k
      logical :: lheader
!
      intent(in) :: f,fp
      intent(out) :: dfp
!
!  Print out header information in first time step.
!
      lheader=(headt .and. lfirst .and. lroot) .or. ldebug
!
!  Identify module and boundary conditions
!
      if (lheader) print*,'dxxp_dt: Calculate dxxp_dt'
      if (lheader) then
        print*, 'dxxp_dt: Particles boundary condition bcpx=', bcpx
        print*, 'dxxp_dt: Particles boundary condition bcpy=', bcpy
        print*, 'dxxp_dt: Particles boundary condition bcpz=', bcpz
      endif
!
!  Interpolate gas velocity to position of particles. 
!  Then set particle velocity equal to the local gas velocity.
!      
      do k=1,npar_loc
        call interpolate_3d_1st(f(:,:,:,iux:iuz),fp(ipar(k),ixp:izp),uup)
        dfp(ipar(k),ixp:iyp) = dfp(ipar(k),ixp:iyp) + uup
      enddo
!
!  Diagnostic output
!
      if (ldiagnos) then
        if (idiag_xpm/=0) &
            call sum_mn_name_par(fp(ipar(1:npar_loc),ixp),idiag_xpm)
        if (idiag_ypm/=0) &
            call sum_mn_name_par(fp(ipar(1:npar_loc),iyp),idiag_ypm)
        if (idiag_zpm/=0) &
            call sum_mn_name_par(fp(ipar(1:npar_loc),izp),idiag_zpm)
      endif
!
    endsubroutine dxxp_dt
!***********************************************************************
    subroutine read_particles_init_pars(unit,iostat)
!    
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=particles_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=particles_init_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_particles_init_pars
!***********************************************************************
    subroutine write_particles_init_pars(unit)
!    
      integer, intent(in) :: unit
!
      write(unit,NML=particles_init_pars)
!
    endsubroutine write_particles_init_pars
!***********************************************************************
    subroutine read_particles_run_pars(unit,iostat)
!    
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=particles_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=particles_run_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_particles_run_pars
!***********************************************************************
    subroutine write_particles_run_pars(unit)
!    
      integer, intent(in) :: unit
!
      write(unit,NML=particles_run_pars)
!
    endsubroutine write_particles_run_pars
!***********************************************************************
    subroutine rprint_particles(lreset,lwrite)
!   
!  Read and register print parameters relevant for particles
!
!  29-dec-04/anders: coded
!
      use Sub, only: parse_name
!
      logical :: lreset
      logical, optional :: lwrite
!
      integer :: iname
      logical :: lwr
! 
!  Write information to index.pro
! 
      lwr = .false.
      if (present(lwrite)) lwr=lwrite

      if (lwr) then
        write(3,*) 'ixxp=',ixxp
        write(3,*) 'ivvp=',ivvp
      endif
!
      if (lreset) then
        idiag_xpm=0; idiag_ypm=0; idiag_zpm=0
      endif
!
!  Run through all possible names that may be listed in print.in
!
      if (lroot .and. ip<14) print*,'rprint_particles: run through parse list'
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'xpm',idiag_xpm)
        call parse_name(iname,cname(iname),cform(iname),'ypm',idiag_ypm)
        call parse_name(iname,cname(iname),cform(iname),'zpm',idiag_zpm)
      enddo
!
    endsubroutine rprint_particles
!***********************************************************************

endmodule Particles
