! $Id: shear.f90,v 1.1 2002-07-04 10:10:55 nilshau Exp $

!  This modules deals with all aspects of shear; if no
!  shear is invoked, a corresponding replacement dummy
!  routine is used instead which absorbs all the calls to the
!  shear relevant subroutines listed in here.

module Shear

  use Sub
  use Cdata

  implicit none


  namelist /shear_init_pars/ &
       Omega,qshear

  namelist /shear_run_pars/ &
       Omega,qshear

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
           "$Id: shear.f90,v 1.1 2002-07-04 10:10:55 nilshau Exp $")
!
    endsubroutine register_shear
!***********************************************************************
    subroutine shearing(f,df)
!
!  Calculates the shear terms, -uy0*df/dy (shearing sheat approximation)
!
!  2-july-02/nils: coded
!
      use Cparam
      use Deriv
      use Hydro, only:theta
!
      real, dimension (mx,my,mz,mvar) :: f,df
      real, dimension(nx) :: uy0,dlnrhody,shearlnrho,dssdy,shearss
      real, dimension (nx,3) :: dudy, shearu, A_vec, dAdy, shearA
!
      uy0=-qshear*Omega*x(l1:l2)
!
      if (lhydro) then
         call der(f,iux,dudy(:,1),2)
         call der(f,iuy,dudy(:,2),2)
         call der(f,iuz,dudy(:,3),2)
         call multvs_mn(dudy,uy0,shearu)
         df(l1:l2,m,n,iux:iuz)=df(l1:l2,m,n,iux:iuz)-shearu
      end if
!
      if (ldensity) then
         call der(f,ilnrho,dlnrhody,2)
         shearlnrho=dlnrhody*uy0
         df(l1:l2,m,n,ilnrho)=df(l1:l2,m,n,ilnrho)-shearlnrho
      end if
!
      if (lentropy) then
         call der(f,ient,dssdy,2)
         shearss=dssdy*uy0
         df(l1:l2,m,n,ient)=df(l1:l2,m,n,ient)-shearss
      end if
!
      if (lmagnetic) then
         call der(f,iaa+0,dAdy(:,1),2)
         call der(f,iaa+1,dAdy(:,2),2)
         call der(f,iaa+2,dAdy(:,3),2)
         call multvs_mn(dAdy,uy0,shearA)
         df(l1:l2,m,n,iaa:iaa+2)=df(l1:l2,m,n,iaa:iaa+2)-shearA
      end if
!
! Taking care of the fact that the Corriolis force changes when 
! we have got shear. The rest of the Corrilis force is calculated 
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
    end subroutine shearing
!***********************************************************************
  end module Shear
