! $Id: shear.f90,v 1.4 2002-08-05 08:02:02 nilshau Exp $

!  This modules deals with all aspects of shear; if no
!  shear is invoked, a corresponding replacement dummy
!  routine is used instead which absorbs all the calls to the
!  shear relevant subroutines listed in here.

module Shear

  use Sub
  use Cdata

  implicit none

  namelist /shear_init_pars/ &
       qshear

  namelist /shear_run_pars/ &
       qshear

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
           "$Id: shear.f90,v 1.4 2002-08-05 08:02:02 nilshau Exp $")
!
    endsubroutine register_shear
!***********************************************************************
    subroutine shearing(f,df)
!
!  Calculates the shear terms, -uy0*df/dy (shearing sheat approximation)
!
!  2-july-02/nils: coded
!  6-july-02/axel: runs through all nvar variables; added timestep check
!
      use Cparam
      use Deriv
      use Hydro, only:Omega,theta
!
      integer :: j
      real, dimension (mx,my,mz,mvar) :: f,df
      real, dimension(nx) :: uy0,dfdy
!
!  print identifier
!
      if (headtt.or.ldebug) print*,'shearing: qshear,Omega=',qshear,Omega
!
!  add shear term, -uy0*df/dy, for all variables
!
      if (Omega==0) then 
         uy0=-qshear*x(l1:l2)
         if (headtt.or.ldebug) print*,&
              'shearing: shear without rotation, max(uy0)=',maxval(uy0)
      else 
         uy0=-qshear*Omega*x(l1:l2)
      endif
      do j=1,nvar
        call der(f,j,dfdy,2)
        df(l1:l2,m,n,j)=df(l1:l2,m,n,j)-uy0*dfdy
      enddo
!
! Taking care of the fact that the Coriolis force changes when 
! we have got shear. The rest of the Coriolis force is calculated 
! in hydro.
!
      if (Omega/=0.) then
        if (theta==0) then
          df(l1:l2,m,n,iuy)=df(l1:l2,m,n,iuy)+qshear*Omega*f(l1:l2,m,n,iux)
        else
          df(l1:l2,m,n,iuy)=df(l1:l2,m,n,iuy)&
               +qshear*Omega*cos(theta*pi/180.)*f(l1:l2,m,n,iux)
        endif
      endif
!
!  take shear into account for calculating time step
!
      if (lfirst.and.ldt) maxadvec2=amax1(maxadvec2,uy0**2)
!
    end subroutine shearing
!***********************************************************************
  end module Shear
