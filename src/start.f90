! $Id: start.f90,v 1.26 2002-05-13 18:52:54 dobler Exp $
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
        real, dimension (mx,my,mz) :: xx,yy,zz
!
        call initialize         ! register modules, etc.
!
!  identify version
!
        if (lroot) call cvs_id( &
             "$RCSfile: start.f90,v $", &
             "$Revision: 1.26 $", &
             "$Date: 2002-05-13 18:52:54 $")
!
!  read input parameter (by each processor)
!
        open(1,FILE='start.in',FORM='formatted')
        read(1,*) ip
        read(1,*) x0,y0,z0
        read(1,*) Lx,Ly,Lz
        read(1,*) iperx,ipery,iperz
        read(1,*) z1,z2,ztop
        read(1,*) ampl,init,urand
        read(1,*) cs0,gamma,rho0,gravz,grads0
        !
        ! forcing needs no init parameters
        if (lentropy)  read(1,NML=entropy_init_pars )
        if (lmagnetic) read(1,NML=magnetic_init_pars)
        close(1)
!
!  output on the console, but only when root processor
!
        if (lroot) then
          print*, 'ip=', ip
          print*, 'x0,y0,z0=', x0,y0,z0
          print*, 'Lx,Ly,Lz=', Lx,Ly,Lz
          print*, 'iperx,ipery,iperz=', iperx,ipery,iperz 
          print*, 'z1,z2,ztop=', z1,z2,ztop
          print*, 'ampl,init,urand=', ampl,init,urand
          print*, 'cs0,gamma,rho0,gravz,grads0=', cs0,gamma,rho0,gravz,grads0
          !
          ! forcing needs no init parameters
          if (lentropy ) write(*,NML=entropy_init_pars )
          if (lmagnetic) write(*,NML=magnetic_init_pars)

        endif
!
!  write input parameters to a parameter file (for run.x and IDL)
!
        gamma1 = gamma-1.
        call wparam()
!
!  generate mesh, |x| < Lx, and similar for y and z.
!  iper{x,y,z} indicate periodicity of given direction
!
        if (iperx /= 0) then; dx = Lx/nxgrid; else; dx = Lx/(nxgrid-1); endif
        if (ipery /= 0) then; dy = Ly/nygrid; else; dy = Ly/(nygrid-1); endif
        if (iperz /= 0) then; dz = Lz/nzgrid; else; dz = Lz/(nzgrid-1); endif

        do i=1,mx; x(i)=x0+(i-nghost-1+ipx*nx)*dx; enddo
        do i=1,my; y(i)=y0+(i-nghost-1+ipy*ny)*dy; enddo
        do i=1,mz; z(i)=z0+(i-nghost-1+ipz*nz)*dz; enddo
!
        call wgrid(trim(directory)//'/grid.dat')
!
!  as full arrays
!
        xx=spread(spread(x,2,my),3,mz)
        yy=spread(spread(y,1,mx),3,mz)
        zz=spread(spread(z,1,mx),2,my)
!
        rr=sqrt(xx**2+yy**2+zz**2)
!        m_pot=(1.+rr**2)/(1.+rr**2+rr**3) ! negative potential
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
        if (lroot) call wdim('tmp/dim.dat', &
             nxgrid+2*nghost,nygrid+2*nghost,nzgrid+2*nghost)
!  write global variables:
        call wglobal()

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
