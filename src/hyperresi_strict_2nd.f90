! $Id$
!
!  This module applies a sixth order hyperresistivity to the induction
!  equation (following Brandenburg & Sarson 2002). This hyperresistivity
!  ensures that the energy dissipation rate is positive define everywhere.
!
!  Spatial derivatives are accurate to second order.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lhyperresistivity_strict=.true.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 3
!
!***************************************************************
module Hyperresi_strict
!
  use Cparam
  use Cdata
  use Messages
!  use Density
!
  implicit none
!
  include 'hyperresi_strict.h'
!
  contains
!***********************************************************************
    subroutine register_hyperresi_strict()
!
!  Set up indices for hyperresistivity auxiliary slots.
!
!  23-aug-07/anders: adapted from register_hypervisc_strict
!
      use FArrayManager
!
      if (lroot) call svn_id( &
          "$Id$")
!
!  Set indices for auxiliary variables
!
      call farray_register_auxiliary('hypres',ihypres,vector=3)
!
    endsubroutine register_hyperresi_strict
!***********************************************************************
    subroutine hyperresistivity_strict(f)
!
!  Apply sixth order hyperresistivity with positive definite heating rate
!  (see Brandenburg & Sarson 2002).
!
!  To avoid communicating ghost zones after each operator, we use
!  derivatives that are second order in space.
!
!  23-aug-07/anders: adapted from hyperviscosity_strict
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      real, dimension (mx,my,mz,3) :: tmp
!
!  Calculate del2(del2(del2(A))), accurate to second order.
!
      call del2v_2nd(f,tmp,iaa)
      f(:,:,:,ihypres:ihypres+2)=tmp
      call del2v_2nd(f,tmp,ihypres)
      f(:,:,:,ihypres:ihypres+2)=tmp
      call del2v_2nd(f,tmp,ihypres)
      f(:,:,:,ihypres:ihypres+2)=tmp
!
!     [insert resistive heating (eta/mu0)*curl(curl(B))^2 here]
!
    endsubroutine hyperresistivity_strict
!***********************************************************************
    subroutine del2v_2nd(f,del2f,k)
!
!  Calculate Laplacian of a vector, accurate to second order.
!
!  24-nov-03/nils: adapted from del2v
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,3) :: del2f
      real, dimension (mx,my,mz) :: tmp
      integer :: i,k,k1
!
      intent (in) :: f, k
      intent (out) :: del2f
!
      del2f=0.
!
!  Apply Laplacian to each vector component individually.
!
      k1=k-1
      do i=1,3
        call del2_2nd(f,tmp,k1+i)
        del2f(:,:,:,i)=tmp
      enddo
!
    endsubroutine del2v_2nd
!***********************************************************************
    subroutine del2_2nd(f,del2f,k)
!
!  Calculate Laplacian of a scalar.
!  Accurate to second order.
!
!  24-nov-03/nils: adapted from del2
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz) :: del2f,d2fd
      integer :: k,k1
!
      intent (in) :: f, k
      intent (out) :: del2f
!
      k1=k-1
      call der2_2nd(f,d2fd,k,1)
      del2f=d2fd
      call der2_2nd(f,d2fd,k,2)
      del2f=del2f+d2fd
      call der2_2nd(f,d2fd,k,3)
      del2f=del2f+d2fd
!
    endsubroutine del2_2nd
!***********************************************************************
    subroutine der2_2nd(f,der2f,i,j)
!
!  Calculate the second derivative of f.
!  Accurate to second order.
!
!  24-nov-03/nils: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz) :: der2f
      integer :: i,j
!
      intent (in) :: f,i,j
      intent (out) :: der2f
!
      der2f=0.
!
      if (j==1 .and. nxgrid/=1) then
        der2f(2:mx-1,:,:) = (+1.*f(1:mx-2,:,:,i) &
                             -2.*f(2:mx-1,:,:,i) &
                             +1.*f(3:mx  ,:,:,i) ) / (dx**2)
      endif
!
     if (j==2 .and. nygrid/=1) then
        der2f(:,2:my-1,:) = (+1.*f(:,1:my-2,:,i) &
                             -2.*f(:,2:my-1,:,i) &
                             +1.*f(:,3:my  ,:,i) ) / (dy**2)
      endif
!
     if (j==3 .and. nzgrid/=1) then
        der2f(:,:,2:mz-1) = (+1.*f(:,:,1:mz-2,i) &
                             -2.*f(:,:,2:mz-1,i) &
                             +1.*f(:,:,3:mz  ,i) ) / (dz**2)
      endif
!
    endsubroutine der2_2nd
!***********************************************************************
endmodule Hyperresi_strict
