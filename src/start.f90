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
        real, dimension (mx,my,mz) :: xx,yy,zz
!
        call initialize         ! register modules, etc.
!
!  read input parameter (by each processor)
!
        pi=2*asin(1.)
        open(1,FILE='start.in',FORM='formatted')
        read(1,*) Lx,Ly,Lz
        read(1,*) ampl,init
        read(1,*) cs,gravz
        close(1)
!
!  output on the console, but only when root processor
!
        if (lroot)then
          print*, 'Lx,Ly,Lz=', Lx,Ly,Lz
          print*, 'ampl,init=', ampl,init
          print*, 'cs,gravz=', cs,gravz
        endif
!
!  generate mesh, |x| < Lx, and similar for y and z.
!
        dx=Lx/nx
        dy=Ly/(ny*nprocy)
        dz=Lz/(nz*nprocz)
        do i=1,mx; x(i)=-Lx/2.+(i-0.5-nghost       )*dx; enddo
        do i=1,my; y(i)=-Ly/2.+(i-0.5-nghost+ipy*ny)*dy; enddo
        do i=1,mz; z(i)=-Lz/2.+(i-0.5-nghost+ipz*nz)*dz; enddo
!
!  as full arrays
!
        xx=spread(spread(x,2,my),3,mz)
        yy=spread(spread(y,1,mx),3,mz)
        zz=spread(spread(z,1,mx),2,my)
!
        cs20=cs**2 ! (goes into cdata module)
!
!  different initial conditions
!
        f = 0.
        call init_hydro(f,init,ampl,xx,yy,zz)
        call init_ent  (f,init,ampl,xx,yy,zz)
        call init_aa   (f,init,ampl,xx,yy,zz)
!
!  write to disk
!
        call output(trim(directory)//'/var.dat',f,mvar)
        call wdim(trim(directory)//'/dim.dat')
!  also write full dimensions to tmp/ :
        if (lroot) call wdim('tmp/dim.dat',nygrid,nzgrid)

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
