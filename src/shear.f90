! $Id: shear.f90,v 1.5 2002-08-16 21:23:48 brandenb Exp $

!  This modules deals with all aspects of shear; if no
!  shear is invoked, a corresponding replacement dummy
!  routine is used instead which absorbs all the calls to the
!  shear relevant subroutines listed in here.
!  Shear can either be given relative to Omega (using qshear),
!  or in absolute fashion via the parameters Sshear.

module Shear

  use Sub
  use Cdata

  implicit none

  namelist /shear_init_pars/ &
       qshear,Sshear

  namelist /shear_run_pars/ &
       qshear,Sshear

  contains

!***********************************************************************
    subroutine register_shear()
!
!  Initialise variables
!
!  2-july-02/nils: coded
!
      use Mpicomm
!
      logical, save :: first=.true.
!
      if (.not. first) call stop_it('register_shear called twice')
      first = .false.
!
      lshear = .true.
!
!  identify version number
!
      if (lroot) call cvs_id( &
           "$Id: shear.f90,v 1.5 2002-08-16 21:23:48 brandenb Exp $")
!
    endsubroutine register_shear
!***********************************************************************
    subroutine shearing(f,df)
!
!  Calculates the shear terms, -uy0*df/dy (shearing sheat approximation)
!
!  2-jul-02/nils: coded
!  6-jul-02/axel: runs through all nvar variables; added timestep check
! 16-aug-02/axel: use now Sshear which is calculated in param_io.f90
!
      use Cparam
      use Deriv
      use Hydro, only:theta
!
      integer :: j
      real, dimension (mx,my,mz,mvar) :: f,df
      real, dimension(nx) :: uy0,dfdy
!
!  print identifier
!
      if (headtt.or.ldebug) print*,'shearing: Sshear,qshear=',Sshear,qshear
!
!  add shear term, -uy0*df/dy, for all variables
!
      uy0=Sshear*x(l1:l2)
      do j=1,nvar
        call der(f,j,dfdy,2)
        df(l1:l2,m,n,j)=df(l1:l2,m,n,j)-uy0*dfdy
      enddo
!
! Taking care of the fact that the Coriolis force changes when 
! we have got shear. The rest of the Coriolis force is calculated 
! in hydro.
!
      if (theta==0) then
        df(l1:l2,m,n,iuy)=df(l1:l2,m,n,iuy)-Sshear*f(l1:l2,m,n,iux)
      else
        if (headtt) print*,'Sure you want Sshear with finite theta??'
        df(l1:l2,m,n,iuy)=df(l1:l2,m,n,iuy)&
             -Sshear*cos(theta*pi/180.)*f(l1:l2,m,n,iux)
      endif
!
!  take shear into account for calculating time step
!
      if (lfirst.and.ldt) maxadvec2=amax1(maxadvec2,uy0**2)
!
    end subroutine shearing
!***********************************************************************
  end module Shear
