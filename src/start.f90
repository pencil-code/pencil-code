! $Id: start.f90,v 1.47 2002-06-25 14:58:47 dobler Exp $
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
        use General
        use Mpicomm
        use Sub
        use Register
        use Global
        use Param_IO
!
        implicit none
!
!  define parameters
!
        integer :: i
        real, dimension (mx,my,mz,mvar) :: f
        real, dimension (mx,my,mz) :: xx,yy,zz
        real :: x00,y00,z00
!
        call initialize         ! register modules, etc.
!
!  identify version
!
        if (lroot) call cvs_id( &
             "$Id: start.f90,v 1.47 2002-06-25 14:58:47 dobler Exp $")
!
!  set default values: box of size (2pi)^3
!
        xyz0 = (/  -pi,  -pi,  -pi /) ! first corner
        Lxyz = (/ 2*pi, 2*pi, 2*pi /) ! box lengths
        lperi =(/.true.,.true.,.true. /) ! all directions periodic
!
        z1=0.; z2=1.; zref=0.
!
!  read parameters from start.in
!  call also rprint_list, because it writes iuu, ilnrho, ient, and iaa to disk.
!
        call read_inipars()
        call rprint_list(.false.)
!
!  postprocess input parameters
!
        gamma1 = gamma-1.
        call wparam()
        x0 = xyz0(1) ; y0 = xyz0(2) ; z0 = xyz0(3)
        Lx = Lxyz(1) ; Ly = Lxyz(2) ; Lz = Lxyz(3)
!
!  check consistency
!
        if (.not.lperi(1).and.nxgrid<2) stop 'for nonperiodic: must have nxgrid>1'
        if (.not.lperi(2).and.nygrid<2) stop 'for nonperiodic: must have nygrid>1'
        if (.not.lperi(3).and.nzgrid<2) stop 'for nonperiodic: must have nzgrid>1'
!
!  generate mesh, |x| < Lx, and similar for y and z.
!  lperi indicate periodicity of given direction
!
        if (lperi(1)) then; dx = Lx/nxgrid; x00=x0+.5*dx; else; dx = Lx/(nxgrid-1); x00=x0; endif
        if (lperi(2)) then; dy = Ly/nygrid; y00=y0+.5*dy; else; dy = Ly/(nygrid-1); y00=y0; endif
        if (lperi(3)) then; dz = Lz/nzgrid; z00=z0+.5*dz; else; dz = Lz/(nzgrid-1); z00=z0; endif
!
!  set x,y,z arrays
!
        do i=1,mx; x(i)=x00+(i-nghost-1+ipx*nx)*dx; enddo
        do i=1,my; y(i)=y00+(i-nghost-1+ipy*ny)*dy; enddo
        do i=1,mz; z(i)=z00+(i-nghost-1+ipz*nz)*dz; enddo
!
!  write grid.dat file
!
        call wgrid(trim(directory)//'/grid.dat')
!
!  as full arrays
!
        xx=spread(spread(x,2,my),3,mz)
        yy=spread(spread(y,1,mx),3,mz)
        zz=spread(spread(z,1,mx),2,my)
!
!        rr=sqrt(xx**2+yy**2+zz**2)
!        m_pot=(1.+rr**2)/(1.+rr**2+rr**3) ! negative potential
!AB: what's the plan here? move to grav_r?
!
!  different initial conditions
!
        f = 0.
        if (lroot) print* !(empty line)
        call init_hydro(f,xx,yy,zz)
        call init_lnrho(f,xx,yy,zz)
        call init_ent  (f,xx,yy,zz)
        call init_aa   (f,xx,yy,zz)
        call init_grav (f,xx,yy,zz)
!
!  write to disk
!  The option lnowrite writes everything except the actual var.dat file
!  This can be useful if auxiliary files are outdated, and don't want
!  to overwrite an existing var.dat
!
        if (lwrite_ic) call output(trim(directory)//'/VAR0',f,mvar)
        if (.not.lnowrite) call output(trim(directory)//'/var.dat',f,mvar)
        call wdim(trim(directory)//'/dim.dat')
!
!  also write full dimensions to tmp/ :
!
        if (lroot) call wdim('tmp/dim.dat', &
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
!  seed for random number generator, have to have the same on each
!  processor as forcing is applied in (global) Beltrami modes
!
        seed(1)=1000
        call outpui(trim(directory)//'/seed.dat',seed,2)
!
        call mpifinalize
        if (lroot) print* !(finish with an empty line)
!
      endprogram
