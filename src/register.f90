
!!!  A module for setting up f and related variables (`register' the
!!!  entropy, magnetic, etc modules). Didn't know where else to put this:
!!!  Entropy uses Sub and Init must be used by both, Start and Run.


module Register

  implicit none 

  contains

!***********************************************************************
    subroutine initialize()
!
!  Call all initialisation hooks, i.e. initialise MPI and register
!  physics modules.
!
!  6-nov-01/wolf: coded
!
      use Cdata
      use Mpicomm
      use Sub
      use Gravity
      use Entropy
      use Magnetic
!
      call mpicomm_init
!
      nvar = 0                  ! to start with
      call register_hydro
      call register_ent
      call register_aa
      call register_grav
!
      if (nvar /= mvar) then
        if (lroot) write(0,*) 'nvar = ', nvar, ', mvar = ', mvar
        call stop_it('Initialize: nvar /= mvar')
      endif

!
    endsubroutine initialize
!***********************************************************************
    subroutine register_hydro
!
!  Initialise variables which should know that we solve the hydro
!  equations: iuu, etc; increase nvar accordingly.
!  Should be located in the Hydro module, if there was one.
!
!  6-nov-01/wolf: coded
!
      use Cdata
      use Mpicomm, only: lroot,stop_it
      use Sub
!
      logical, save :: first=.true.
!
      if (.not. first) call stop_it('init_hydro called twice')
      first = .false.
!
      iuu = nvar+1             ! indices to access uu and lam
      iux = iuu
      iuy = iuu+1
      iuz = iuu+2
      ilnrho = iuu+3
      nvar = nvar+4             ! added 4 variables
!
      if ((ip<=8) .and. lroot) then
        print*, 'Register_hydro:  nvar = ', nvar
        print*, 'iuu,ilnrho = ', iuu,ilnrho
        print*, 'iux,iuy,iuz = ', iux,iuy,iuz
      endif
!
!  identify version number
!
      if (lroot) call cvs_id( &
           "$RCSfile: register.f90,v $", &
           "$Revision: 1.8 $", &
           "$Date: 2002-01-17 11:42:43 $")
!
!
      if (nvar > mvar) then
        if (lroot) write(0,*) 'nvar = ', nvar, ', mvar = ', mvar
        call stop_it('Register_hydro: nvar > mvar')
      endif
!
    endsubroutine register_hydro
!***********************************************************************
    subroutine init_hydro(f,init,ampl,xx,yy,zz)
!
!  initialise uu and lnrho; called from start.f90
!  Should be located in the Hydro module, if there was one.
!
!  7-nov-2001/wolf: coded
!
      use Cdata
      use Sub
!
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (mx,my,mz) :: tmp,r,p,xx,yy,zz
      real :: ampl
      real :: dsdz0 = -0.2      ! (1/c_p)ds/dz
      real :: zmax
      integer :: init
!
      select case(init)
      case(0)               ! random ux (Gaussian distribution)
        if (lroot) print*,'Gaussian-distributed ux'
        call random_number(r)
        call random_number(p)
        tmp=sqrt(-2*alog(r))*sin(2*pi*p)
        f(:,:,:,iux)=ampl*tmp
      case(1)               ! density stratification
        if (lgravz) then
          if (lroot) print*,'vertical density stratification'
!        f(:,:,:,ilnrho)=-zz
! isentropic case:
!        zmax = -cs20/gamma1/gravz
!        print*, 'zmax = ', zmax
!        f(:,:,:,ilnrho) = 1./gamma1 * alog(abs(1-zz/zmax))
! linear entropy gradient;
          f(:,:,:,ilnrho) = -dsdz0*zz &
                            + 1./gamma1*alog( 1 + gamma1*gravz/dsdz0/cs20 &
                                                  *(1-exp(-dsdz0*zz)) )
        endif
        !
        if (lgravr) then
          if (lroot) print*,'radial density stratification (so far trivial)'
          f(:,:,:,ilnrho) = 0.
        endif
        !
!  The following needs better implementation: add random component 
!  *after* major initialisation and make amplitude input parameter
call random_number(tmp)
f(:,:,:,iux) = 0.2*(tmp-0.5)   ! velocity perturbation
call random_number(tmp)
f(:,:,:,iuy) = 0.2*(tmp-0.5)   ! velocity perturbation
call random_number(tmp)
f(:,:,:,iuz) = 0.2*(tmp-0.5)   ! velocity perturbation
!
      case(2)               ! oblique sound wave
        if (lroot) print*,'oblique sound wave'
        f(:,:,:,ilnrho)=ampl*cos(xx+2*yy)*sqrt(5.)
        f(:,:,:,iux)=ampl*cos(xx+2*yy)
        f(:,:,:,iuy)=ampl*cos(xx+2*yy)*2.
      case(3)               ! uu = (sin 2x, sin 3y , cos z)
        if (lroot) print*,'uu harmonic (what is this good for?)'
        f(:,:,:,iux)=spread(spread(sin(2*x),2,my),3,mz)* &
                     spread(spread(sin(3*y),1,mx),3,mz)* &
                     spread(spread(cos(1*z),1,mx),2,my)
      case default
        if (lroot) print*,'Initialising everything to zero'
      endselect
!
!
    endsubroutine init_hydro
!***********************************************************************

endmodule Register

!!! End of file register.f90
