
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
      use Entropy
      use Magnetic
!
      call mpicomm_init
!
      nvar = 0                  ! to start with
      call register_hydro
      call register_ent
      call register_aa
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
      real :: gamma1,zmax
      integer :: init
!
      gamma1 = gamma-1
      select case(init)
      case(0)               ! random ux (Gaussian distribution)
        if (lroot) print*,'Gaussian-distributed ux'
        call random_number(r)
        call random_number(p)
        tmp=sqrt(-2*alog(r))*sin(2*pi*p)
        f(:,:,:,iux)=ampl*tmp
      case(1)               ! density stratification
        if (lroot) print*,'density stratification'
!        f(:,:,:,ilnrho)=-zz
        zmax = -cs20/gamma1/gravz
        print*, 'zmax = ', zmax
        f(:,:,:,ilnrho) = 1./gamma1 * alog(abs(1-zz/zmax))
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
    endsubroutine init_hydro
!***********************************************************************

endmodule Register

!!! End of file init.f90
