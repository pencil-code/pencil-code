! $Id: particles_tracers.f90,v 1.4 2005-06-29 12:43:48 ajohan Exp $
!
!  This module takes care of everything related to tracer particles
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

  real, dimension (mpar_loc,mpvar) :: fp=0.0, dfp=0.0
  integer, dimension (mpar_loc) :: ipar
  integer :: npar_loc, npvar

  real :: xp0=0.0, yp0=0.0, zp0=0.0
  real :: dsnap_par_minor=0.0
  character (len=labellen) :: initxxp='origin'

  namelist /particles_init_pars/ &
      initxxp, xp0, yp0, zp0, bcpx, bcpy, bcpz

  namelist /particles_run_pars/ &
      bcpx, bcpy, bcpz

  integer :: idiag_xpm=0, idiag_ypm=0, idiag_zpm=0
  integer :: idiag_xp2m=0, idiag_yp2m=0, idiag_zp2m=0

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
!  Set up initial conditions for particle modules.
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
      call input_particles(filename,fp,npar_loc,ipar)
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
        call wsnap_particles(chsnap,fp,msnap,enum,lsnap,dsnap_par_minor, &
            npar_loc,ipar,flist)
      else
        call wsnap_particles(chsnap,fp,msnap,enum,lsnap,dsnap_par_minor, &
            npar_loc,ipar)
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
      integer :: k
!
      if (itsub==1) then
        do k=1,npar_loc
          dfp(k,:)=0.
        enddo
      else
        do k=1,npar_loc
          dfp(k,:)=alpha(itsub)*dfp(k,:)
        enddo
      endif
!
    endsubroutine particles_timestep_first
!***********************************************************************
    subroutine particles_timestep_second()
!
!  Time evolution of particle variables.
!
!  07-jan-05/anders: coded
!
      integer :: k
!
      do k=1,npar_loc
        fp(k,:) = fp(k,:) + dt_beta(itsub)*dfp(k,:)
      enddo
!
    endsubroutine particles_timestep_second
!***********************************************************************
    subroutine particles_pde(f,df)
!
!  07-jan-05/anders: coded
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
!
      intent (in) :: f
!
      call boundconds_particles(fp,npar_loc,ipar,dfp=dfp)
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
           "$Id: particles_tracers.f90,v 1.4 2005-06-29 12:43:48 ajohan Exp $")
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
      use EquationOfState, only: cs0
!
      logical :: lstarting
!
      real :: rhom
      integer, dimension (0:ncpus-1) :: ipar1, ipar2
!
!  Distribute particles evenly among processors to begin with.
!
      if (lstarting) call dist_particles_evenly_procs(npar_loc,ipar)
!
!  Size of box at local processor is needed for particle boundary conditions.
!
      Lxyz_loc(1)=Lxyz(1)/nprocx
      Lxyz_loc(2)=Lxyz(2)/nprocy
      Lxyz_loc(3)=Lxyz(3)/nprocz
      xyz0_loc(1)=xyz0(1)
      xyz0_loc(2)=xyz0(2)+ipy*Lxyz_loc(2)
      xyz0_loc(3)=xyz0(3)+ipz*Lxyz_loc(3)
      xyz1_loc(1)=xyz1(1)
      xyz1_loc(2)=xyz0(2)+(ipy+1)*Lxyz_loc(2)
      xyz1_loc(3)=xyz0(3)+(ipz+1)*Lxyz_loc(3)
!
    endsubroutine initialize_particles
!***********************************************************************
    subroutine init_particles(f)
!
!  Initial positions and velocities of tracer particles.
!
!  29-dec-04/anders: coded
!
      use Boundcond
      use General, only: random_number_wrapper
      use Mpicomm, only: stop_it
      use Sub
!
      real, dimension (mx,my,mz,mvar+maux) :: f
!
      real, dimension (3) :: uup
      real :: r, p
      integer :: k
!
      intent (in) :: f
!
!  Initial particle position.
!
      select case(initxxp)

      case ('origin')
        if (lroot) print*, 'init_particles: All particles at origin'
        fp(1:npar_loc,ixp:izp)=0.

      case ('constant')
        if (lroot) &
            print*, 'init_particles: All particles at x,y,z=', xp0, yp0, zp0
        fp(1:npar_loc,ixp)=xp0
        fp(1:npar_loc,iyp)=yp0
        fp(1:npar_loc,izp)=zp0

      case ('random')
        if (lroot) print*, 'init_particles: Random particle positions'
        do k=1,npar_loc
          if (nxgrid/=1) call random_number_wrapper(fp(k,ixp))
          if (nygrid/=1) call random_number_wrapper(fp(k,iyp))
          if (nzgrid/=1) call random_number_wrapper(fp(k,izp))
        enddo
        if (nxgrid/=1) fp(1:npar_loc,ixp)=xyz0_loc(1)+fp(1:npar_loc,ixp)*Lxyz_loc(1)
        if (nygrid/=1) fp(1:npar_loc,iyp)=xyz0_loc(2)+fp(1:npar_loc,iyp)*Lxyz_loc(2)
        if (nzgrid/=1) fp(1:npar_loc,izp)=xyz0_loc(3)+fp(1:npar_loc,izp)*Lxyz_loc(3)

      case ('gaussian-z')
        if (lroot) print*, 'init_particles: Gaussian particle positions'
        do k=1,npar_loc
          if (nxgrid/=1) call random_number_wrapper(fp(k,ixp))
          if (nygrid/=1) call random_number_wrapper(fp(k,iyp))
          call random_number_wrapper(r)
          call random_number_wrapper(p)
          fp(k,izp)=zp0*sqrt(-2*alog(r))*cos(2*pi*p)
        enddo
        if (nxgrid/=1) fp(1:npar_loc,ixp)=xyz0_loc(1)+fp(1:npar_loc,ixp)*Lxyz_loc(1)
        if (nygrid/=1) fp(1:npar_loc,iyp)=xyz0_loc(2)+fp(1:npar_loc,iyp)*Lxyz_loc(2)

      case default
        if (lroot) print*, 'init_particles: No such such value for initxxp: ', &
            trim(initxxp)
        call stop_it("")

      endselect
!
!  Redistribute particles among processors (now that positions are determined).
!
      call boundconds_particles(fp,npar_loc,ipar)
!
    endsubroutine init_particles
!***********************************************************************
    subroutine dxxp_dt(f,fp,dfp)
!
!  Evolution of tracer particle position.
!
!  02-jan-05/anders: coded
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mpar_loc,mpvar) :: fp, dfp
!
      real, dimension (3) :: uu
      integer :: k
      logical :: lheader, lfirstcall=.true.
!
      intent (in) :: f, fp
      intent (out) :: dfp
!
!  Print out header information in first time step.
!
      lheader=(headt .and. lfirst .and. lroot) .or. ldebug
!
!  Identify module and boundary conditions.
!
      if (lheader) print*,'dxxp_dt: Calculate dxxp_dt'
      if (lheader) then
        print*, 'dxxp_dt: Particles boundary condition bcpx=', bcpx
        print*, 'dxxp_dt: Particles boundary condition bcpy=', bcpy
        print*, 'dxxp_dt: Particles boundary condition bcpz=', bcpz
      endif
!
      if (lheader) print*, 'dxxp_dt: Set rate of change of particle '// &
          'position equal to gas velocity.'
!
!  Interpolate gas velocity to position of particles. 
!  Then set particle velocity equal to the local gas velocity.
!       
      do k=1,npar_loc
        call interpolate_3d_1st(f,iux,fp(k,ixp:izp),uu)
        dfp(k,ixp:iyp) = dfp(k,ixp:iyp) + uu
      enddo
!
!  Diagnostic output
!
      if (ldiagnos) then
        if (idiag_xpm/=0) call sum_par_name(fp(1:npar_loc,ixp),idiag_xpm)
        if (idiag_ypm/=0) call sum_par_name(fp(1:npar_loc,iyp),idiag_ypm)
        if (idiag_zpm/=0) call sum_par_name(fp(1:npar_loc,izp),idiag_zpm)
        if (idiag_xp2m/=0) call sum_par_name(fp(1:npar_loc,ixp)**2,idiag_xp2m)
        if (idiag_yp2m/=0) call sum_par_name(fp(1:npar_loc,iyp)**2,idiag_yp2m)
        if (idiag_zp2m/=0) call sum_par_name(fp(1:npar_loc,izp)**2,idiag_zp2m)
      endif
!
    endsubroutine dxxp_dt
!***********************************************************************
    subroutine read_particles_init_pars(unit,iostat)
!    
      integer, intent (in) :: unit
      integer, intent (inout), optional :: iostat
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
      integer, intent (in) :: unit
!
      write(unit,NML=particles_init_pars)
!
    endsubroutine write_particles_init_pars
!***********************************************************************
    subroutine read_particles_run_pars(unit,iostat)
!    
      integer, intent (in) :: unit
      integer, intent (inout), optional :: iostat
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
      integer, intent (in) :: unit
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
      use Cdata
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
        write(3,*) 'ivvp=0'
      endif
!
!  Reset everything in case of reset
!
      if (lreset) then
        idiag_xpm=0; idiag_ypm=0; idiag_zpm=0
        idiag_xp2m=0; idiag_yp2m=0; idiag_zp2m=0
      endif
!
!  Run through all possible names that may be listed in print.in
!
      if (lroot .and. ip<14) print*,'rprint_particles: run through parse list'
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'xpm',idiag_xpm)
        call parse_name(iname,cname(iname),cform(iname),'ypm',idiag_ypm)
        call parse_name(iname,cname(iname),cform(iname),'zpm',idiag_zpm)
        call parse_name(iname,cname(iname),cform(iname),'xp2m',idiag_xp2m)
        call parse_name(iname,cname(iname),cform(iname),'yp2m',idiag_yp2m)
        call parse_name(iname,cname(iname),cform(iname),'zp2m',idiag_zp2m)
      enddo
!
    endsubroutine rprint_particles
!***********************************************************************

endmodule Particles
