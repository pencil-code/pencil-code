! $Id: start.f90,v 1.138 2004-07-24 07:47:14 brandenb Exp $
!
!***********************************************************************
      program start
!
!  start, to setup initial condition and file structure
!
!-----------------------------------------------------------------------
!  01-apr-01/axel+wolf: coded
!  01-sep-01/axel: adapted to MPI
!
        use Cdata
        use Grid
        use General
        use Initcond
        use Mpicomm
        use Sub
        use IO
        use Radiation
        use Register
        use Global
        use Param_IO
        use Wsnaps
        use Filter
        use Special, only: init_special
!
        implicit none
!
!  define parameters
!  The f-array includes auxiliary variables
!  Although they are not necessary for start.f90, idl may want to read them,
!  so we define therefore the full array and write it out.
!
        integer :: i,ifilter
        logical :: lnoerase=.false.
!       logical :: lock=.false.
!       logical :: exist
        real, dimension (mx,my,mz,mvar+maux) :: f
        real, dimension (mx,my,mz,mvar) :: df
        real, dimension (mx,my,mz) :: xx,yy,zz
        real :: x00,y00,z00
        real :: c_grid1
!
        lstart = .true.
!
        call register_modules()         ! register modules, etc.
!
!  identify version
!
        if (lroot) call cvs_id( &
             "$Id: start.f90,v 1.138 2004-07-24 07:47:14 brandenb Exp $")
!
!  set default values: box of size (2pi)^3
!
        xyz0 = (/       -pi,        -pi,       -pi /) ! first corner
        xyz1 = (/impossible, impossible, impossible/) ! last corner
        Lxyz = (/impossible, impossible, impossible/) ! box lengths
        lperi =(/.true.,.true.,.true. /) ! all directions periodic
        lequidist=(/.true.,.true.,.true. /) ! all directions equidistant grid
! dgm
        lshift_origin=(/.false.,.false.,.false./) ! don't shift origin
!
!  read parameters from start.in
!  call also rprint_list, because it writes iuu, ilnrho, iss, and iaa to disk.

        call read_startpars(FILE=.true.)
        call rprint_list(.false.)
!
!  Will we write all slots of f?
!
        if (lwrite_aux) then
          mvar_io = mvar+maux
        else
          mvar_io = mvar
        endif
!
!  print resolution
!
        if (lroot) print*, 'nxgrid,nygrid,nzgrid=',nxgrid,nygrid,nzgrid
!
!  postprocess input parameters
!
        gamma1 = gamma-1.
!
!  I don't think there was a good reason to write param.nml twice (but
!  leave this around for some time [wd; rev. 1.71, 5-nov-2002]
!        call wparam()
!
!  Set up directory names and check whether the directories exist
!
        call directory_names()
!
!  Unfortunately the following test for existence of directory fails under
!  OSF1:
!        inquire(FILE=trim(directory_snap), EXIST=exist)
!        if (.not. exist) &
!             call stop_it('Need directory <' // trim(directory_snap) // '>')
!

!
!  Set box dimensions, make sure Lxyz and xyz1 are in sync.
!  Box defaults to [-pi,pi] for all directions if none of xyz1 or Lxyz are set
!
        do i=1,3
          if (Lxyz(i) == impossible) then
            if (xyz1(i) == impossible) then
              Lxyz(i) = 2*pi    ! default value
            else
              Lxyz(i) = xyz1(i)-xyz0(i)
            endif
          else                  ! Lxyz was set
            if (xyz1(i) /= impossible) then ! both Lxyz and xyz1 are set
              call stop_it('Cannot set Lxyz and xyz1 at the same time')
            endif
          endif
        enddo
        xyz1=xyz0+Lxyz
        !
        !  abbreviations
        !
        x0 = xyz0(1) ; y0 = xyz0(2) ; z0 = xyz0(3)
        Lx = Lxyz(1) ; Ly = Lxyz(2) ; Lz = Lxyz(3)

!  check consistency
!
        if (.not.lperi(1).and.nxgrid<2) stop 'for nonperiodic: must have nxgrid>1'
        if (.not.lperi(2).and.nygrid<2) stop 'for nonperiodic: must have nygrid>1'
        if (.not.lperi(3).and.nzgrid<2) stop 'for nonperiodic: must have nzgrid>1'
        if (lperi(1).and.lshift_origin(1)) stop 'for periodic: must have lshift_origin=F'
        if (lperi(2).and.lshift_origin(2)) stop 'for periodic: must have lshift_origin=F'
        if (lperi(3).and.lshift_origin(3)) stop 'for periodic: must have lshift_origin=F'
!
!  Initialise random number generator in processor-dependent fashion for
!  random initial data.
!  Slightly tricky, since setting seed=(/iproc,0,0,0,0,0,0,0,.../)
!  would produce very short-period random numbers with the Intel compiler;
!  so we need to retain most of the initial entropy of the generator. 
!
        call get_nseed(nseed)   ! get state length of random number generator
        call random_seed_wrapper(get=seed(1:nseed))
        seed(1) = -(10+iproc)    ! different random numbers on different CPUs
        call random_seed_wrapper(put=seed(1:nseed))
!
!  generate grid
!
        call construct_grid(x,y,z,dx,dy,dz)
!
!  write grid.dat file
!
        call wgrid(trim(directory)//'/grid.dat')
!
!  write .general file for data explorer
!
          if (lroot) call write_dx_general( &
                            trim(datadir)//'/var.general', &
                            x00-nghost*dx, y00-nghost*dy, z00-nghost*dz)
!
!  as full arrays
!
        xx=spread(spread(x,2,my),3,mz)
        yy=spread(spread(y,1,mx),3,mz)
        zz=spread(spread(z,1,mx),2,my)
!
!  Parameter dependent initialization of module variables and final
!  pre-timestepping setup (must be done before need_XXXX can be used, for
!  example).
!
        call initialize_modules(f,lstarting=.true.)
!
!  Initial conditions: by default, we put f=0 (ss=lnrho=uu=0, etc)
!  alternatively: read existing snapshot and overwrite only some fields
!  by the various procedures below.
!
        if(lread_oldsnap) then
          print*,'read old snapshot file'
          call input(trim(directory_snap)//'/var.dat',f,mvar,1)
!
!  read non-magnetic snapshot into an MHD run
!
        elseif(lread_oldsnap_nomag) then
          print*,'read old snapshot file (but without magnetic field)'
          call input(trim(directory_snap)//'/var.dat',f,mvar-3,1)
          ! shift the rest of the data
          if (iaz<mvar) then
            do i=iaz+1,mvar
              f(:,:,:,i)=f(:,:,:,i-3)
            enddo
            f(:,:,:,iax:iaz)=0.
          endif
!
!  read data without passive scalar into new run with passive scalar
!
        elseif(lread_oldsnap_nopscalar) then
          print*,'read old snapshot file (but without passive scalar)'
          call input(trim(directory_snap)//'/var.dat',f,mvar-1,1)
          ! shift the rest of the data
          if (iaz<mvar) then
            do i=ilncc+1,mvar
              f(:,:,:,i)=f(:,:,:,i-1)
            enddo
            f(:,:,:,ilncc)=0.
          endif
!
!  default: everything zero
!
        else
          f = 0.
        endif
!
!  the following init routines do then only need to add to f.
!  wd: also in the case where we have read in an existing snapshot??
!
        if (lroot) print* !(empty line)
        call init_gg      (f,xx,yy,zz)
        call init_uu      (f,xx,yy,zz)
        call init_lnrho   (f,xx,yy,zz)
        call init_ss      (f,xx,yy,zz)
        call init_aa      (f,xx,yy,zz)
        call init_rad     (f,xx,yy,zz)
        call init_lncc    (f,xx,yy,zz)
        call init_chiral  (f,xx,yy,zz)
        call init_uud     (f)
        call init_nd      (f)
        call init_ecr     (f,xx,yy,zz)
        call init_special (f,xx,yy,zz)
!
!  check whether we want ionization
!
        if(lionization) call ioninit(f)
        if(lradiation_ray) call radtransfer(f)
!
!  filter initial velocity
!  NOTE: this procedure is currently not very efficient,
!  because for all variables boundary conditions and communication
!  are done again and again for each variable.
!
        if(nfilter/=0) then
          do ifilter=1,nfilter
            call rmwig(f,df,iux,iuz,awig)
          enddo
          if(lroot) print*,'DONE: filter initial velocity, nfilter=',nfilter
        endif
!
!  write to disk
!
        if (lwrite_ic) then
          call wsnap( &
               trim(directory_snap)//'/VAR0',f,mvar_io,ENUM=.false.)
        endif
!
!  The option lnowrite writes everything except the actual var.dat file
!  This can be useful if auxiliary files are outdated, and don't want
!  to overwrite an existing var.dat
!
        inquire(FILE="NOERASE", EXIST=lnoerase)
        if (.not.lnowrite .and. .not.lnoerase) then
          if (ip<12) print*,'START: writing to ' // trim(directory_snap)//'/var.dat'
          call wsnap(trim(directory_snap)//'/var.dat',f,mvar_io,ENUM=.false.)
          call wtime(trim(directory)//'/time.dat',t)
        endif
        call wdim(trim(directory)//'/dim.dat')

!
!  also write full dimensions to data/ :
!
        if (lroot) call wdim(trim(datadir)//'/dim.dat', &
             nxgrid+2*nghost,nygrid+2*nghost,nzgrid+2*nghost)
!
!  write global variables:
!
        call wglobal()
!
!  Write input parameters to a parameter file (for run.x and IDL).
!  Do this late enough, so init_entropy etc. can adjust them
!
        call wparam()
!
!  Seed for random number generator to be used in forcing.f90. Have to
!  have the same on each  processor as forcing is applied in (global)
!  Beltrami modes.
!
        seed(1) = 1812
        call outpui(trim(directory)//'/seed.dat',seed,nseed)
!
        call mpifinalize
        if (lroot) print*
        if (lroot) print*,'start.x has completed successfully'
        if (lroot) print* ! (finish with an empty line)
!
      endprogram
