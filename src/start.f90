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
        use Gravity
        use Entropy
        use Magnetic
!
        implicit none
!
!  define parameters
!
        integer :: init,i
        real :: ampl
        real, dimension (mx,my,mz,mvar) :: f
        real, dimension (mx,my,mz) :: xx,yy,zz,rr
!
!  identify version and register modules
!
        print*,'$Id: start.f90,v 1.8 2002-01-17 10:02:52 dobler Exp $'
        call initialize         ! register modules, etc.
!
!  read input parameter (by each processor)
!
        pi=2*asin(1.)
        open(1,FILE='start.in',FORM='formatted')
        read(1,*) ip
        read(1,*) Lx,Ly,Lz
        read(1,*) ampl,init
        read(1,*) cs0,rho0,gravz
        close(1)
!
!  output on the console, but only when root processor
!
        if (lroot)then
          print*, 'ip=', ip
          print*, 'Lx,Ly,Lz=', Lx,Ly,Lz
          print*, 'ampl,init=', ampl,init
          print*, 'cs0,gravz=', cs0,gravz
        endif
!
!  generate mesh, |x| < Lx, and similar for y and z.
!
        dx=Lx/nx
        dy=Ly/(ny*nprocy)
        dz=Lz/(nz*nprocz)
        dz=Lz/((nz-1)*nprocz)
        do i=1,mx; x(i)=-Lx/2.+(i-0.5-nghost       )*dx; enddo
        do i=1,my; y(i)=-Ly/2.+(i-0.5-nghost+ipy*ny)*dy; enddo
        do i=1,mz; z(i)=-Lz/2.+(i-0.5-nghost+ipz*nz)*dz; enddo
        do i=1,mz; z(i)=-Lz/2.+(i-1.0-nghost+ipz*nz)*dz; enddo
!
!  as full arrays
!
        xx=spread(spread(x,2,my),3,mz)
        yy=spread(spread(y,1,mx),3,mz)
        zz=spread(spread(z,1,mx),2,my)
!
        rr=sqrt(xx**2+yy**2+zz**2)
        m_pot=(1.+rr**2)/(1.+rr**2+rr**3) ! negative potential
!
        cs20=cs0**2 ! (goes into cdata module)
!
!  different initial conditions
!
        f = 0.
        call init_hydro(f,init,ampl,xx,yy,zz)
        call init_ent  (f,init,ampl,xx,yy,zz)
        call init_aa   (f,init,ampl,xx,yy,zz)
        call init_grav (f,init,ampl,xx,yy,zz)
!
!  write to disk
!
        call output(trim(directory)//'/var.dat',f,mvar)
        call wdim(trim(directory)//'/dim.dat')
!  also write full dimensions to tmp/ :
        if (lroot) call wdim('tmp/dim.dat',nygrid,nzgrid)
!  write global variables:
        call write_global()

!
!  seed for random number generator, have to have the same on each
!  processor as forcing is applied in (global) Beltrami modes
!
        seed(1)=1000
        call outpui(trim(directory)//'/seed.dat',seed,1)
        if (iproc < 10) print*,'iproc,seed=',iproc,seed
!
        call mpifinalize
!
      endprogram
