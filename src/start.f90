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
!
        implicit none
!
!  define parameters
!
        integer :: init,i
        real :: ampl
        real, dimension (mx,my,mz,mvar) :: f
        real, dimension (mx,my,mz) :: tmp,r,p,xx,yy,zz
!
        call initialize         ! register modules, etc.
!
!  read input parameter (by each processor)
!
        pi=2*asin(1.)
        open(1,FILE='start.in',FORM='formatted')
        read(1,*) Lx,Ly,Lz
        read(1,*) ampl,init
        close(1)
!
!  output on the console, but only when root processor
!
        if (lroot)then
          print*,Lx,Ly,Lz
          print*,ampl,init
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
!  different initial conditions
!
        select case(init)
          case(0)
            call random_number(r)
            call random_number(p)
            tmp=sqrt(-2*alog(r))*sin(2*pi*p)
            f(:,:,:,1)=ampl*tmp
          case(1)
            if (lroot) print*,'stratification'
            f(:,:,:,4)=-zz
          case(2)
            if (lroot) print*,'oblique sound wave'
            f(:,:,:,4)=ampl*cos(xx+2*yy)*sqrt(5.)
            f(:,:,:,1)=ampl*cos(xx+2*yy)
            f(:,:,:,2)=ampl*cos(xx+2*yy)*2.
          case(3)
            f(:,:,:,1)=spread(spread(sin(2*x),2,my),3,mz)* &
                       spread(spread(sin(3*y),1,mx),3,mz)* &
                       spread(spread(cos(1*z),1,mx),2,my)
          case default
            f=0.
        endselect
!
!  write to disk
!
        call output(trim(directory)//'/var.dat',f,mvar)
        call wdim(trim(directory)//'/dim.dat')
!
!  seed for random number generator, have to have the same on each processor
!
        seed(1)=1000
        call outpui(trim(directory)//'/seed.dat',seed,1)
        if (iproc < 10) print*,'iproc,seed=',iproc,seed
!
        call mpifinalize
!
      endprogram
