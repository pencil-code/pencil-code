! $Id: rotation.f90,v 1.1 2002-07-02 17:14:41 nilshau Exp $

!  This modules deals with all aspects of rotation; if no
!  rotation are invoked, a corresponding replacement dummy
!  routine is used instead which absorbs all the calls to the
!  magnetically relevant subroutines listed in here.

module Rotation

  use Sub
  use Cdata

  implicit none


  namelist /rotation_init_pars/ &
       Omega,qshear

  contains

!***********************************************************************
    subroutine register_rot()
!
!  Initialise variables
!
!  2-july-02/nils: coded
!
      use Mpicomm
!
      logical, save :: first=.true.
!
      if (.not. first) call stop_it('register_rot called twice')
      first = .false.
!
      lrotation = .true.
!
!  identify version number
!
      if (lroot) call cvs_id( &
           "$Id: rotation.f90,v 1.1 2002-07-02 17:14:41 nilshau Exp $")
!
    endsubroutine register_rot
!***********************************************************************
    subroutine shearing(f,df)
!
!  Calculates the actuall shear terms
!
!  2-july-02/nils: coded
!
      use Cparam
      use Deriv
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
    end subroutine shearing
!***********************************************************************
  end module Rotation
