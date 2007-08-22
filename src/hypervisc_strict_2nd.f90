! $Id: hypervisc_strict_2nd.f90,v 1.4 2007-08-22 17:09:22 ajohan Exp $

!
!  This module applies a sixth order hyperviscosity to the equation
!  of motion (following Haugen & Brandenburg 2004). This hyperviscosity
!  ensures that the energy dissipation rate is positive define everywhere.
!
!  The rate of strain tensor
!    S^(3) = (-nab^2)^2*S
!  is a high order generalisation of the first order operator
!    2*S_ij = u_i,j + u_j,i - 2/3*delta_ij*div(u)
!
!  Spatial derivatives are accurate to second order.
!

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lhyperviscosity_strict=.true.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 3
!
!***************************************************************

module Hypervisc_strict

  use Cparam
  use Cdata
  use Messages
  use Density

  implicit none

  include 'hypervisc_strict.h'

  contains

!***********************************************************************
    subroutine register_hypervisc_strict()
!
!  Set up indices for hyperviscosity auxiliary slots.
!
!  19-nov-02/tony: coded
!  24-nov-03/nils: adapted from visc_shock
!
      use Mpicomm, only: stop_it
!
      logical, save :: first=.true.
!
      if (.not. first) call stop_it('register_particles: called twice')
      first = .false.
!
      if (lroot) call cvs_id( &
           "$Id: hypervisc_strict_2nd.f90,v 1.4 2007-08-22 17:09:22 ajohan Exp $")
!
!  Set indices for auxiliary variables
! 
      ihyper = mvar + naux + 1 + (maux_com - naux_com); naux = naux + 3
!
!  Check that we aren't registering too many auxilary variables
!
      if (naux > maux) then
        if (lroot) write(0,*) 'naux = ', naux, ', maux = ', maux
            call stop_it('register_particles: naux > maux')
      endif
! 
    endsubroutine register_hypervisc_strict
!***********************************************************************
    subroutine hyperviscosity_strict(f,k)
!
!  Apply momentum-conserving, symmetric, sixth order hyperviscosity with
!  positive define heating rate (see Haugen & Brandenburg 2004).
!
!  To avoid communicating ghost zones after each operator, we use
!  derivatives that are second order in space.
!
!  24-nov-03/nils: coded
!
      use Cdata, only: lfirst
      use Io
      use Mpicomm
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: k
!
      real, dimension (mx,my,mz,3) :: tmp, tmp2
      real, dimension (nx,ny,nz) :: sij2
      real, dimension (nx,3) :: del6u
      integer :: i,j
!
!  Calculate del2(del2(del2(u))), accurate to second order.
!
      call del2v_2nd(f,tmp,iux)
      f(:,:,:,ihyper:ihyper+2)=tmp
      call del2v_2nd(f,tmp,ihyper)
      f(:,:,:,ihyper:ihyper+2)=tmp
      call del2v_2nd(f,tmp,ihyper)
!
!  Calculate del2(del2(grad(div(u)))), accurate to second order.
!  Probably gives zero derivative at the Nyquist scale, but the del2^3
!  term above gives dissipation at this scale.
!
      call graddivu_2nd(f,tmp2,iux)
      f(:,:,:,ihyper:ihyper+2)=tmp2
      call del2v_2nd(f,tmp2,ihyper)
      f(:,:,:,ihyper:ihyper+2)=tmp2
      call del2v_2nd(f,tmp2,ihyper)
!
!  Add the two terms.
!
      f(:,:,:,ihyper:ihyper+2)=tmp+tmp2/3.
!
! find heating term (yet it only works for ivisc='hyper3')
! the heating term is d/dt(0.5*rho*u^2) = -2*mu3*( S^(2) )^2
!
!      if (ldiagnos) then
!        if (ivisc == 'hyper3') then
!          sij2=0
!          do i=1,2
!            do j=i+1,3
!
! find uij and uji (stored in tmp in order to save space)
!
!              call der_2nd_nof(f(:,:,:,iux+i-1),tmp(:,:,:,1),j) ! uij
!              call der_2nd_nof(f(:,:,:,iux+j-1),tmp(:,:,:,2),i) ! uji
!
! find (nabla2 uij) and (nable2 uji)
!          
!              call del2_2nd_nof(tmp(:,:,:,1),tmp(:,:,:,3))
!              tmp(:,:,:,1)=tmp(:,:,:,3) ! nabla2 uij
!              call del2_2nd_nof(tmp(:,:,:,2),tmp(:,:,:,3))
!              tmp(:,:,:,2)=tmp(:,:,:,3) ! nabla2 uji
!
! find sij2 for full array
!
!              sij2=sij2+tmp(l1:l2,m1:m2,n1:n2,1)*tmp(l1:l2,m1:m2,n1:n2,2) &
!                   +0.5*tmp(l1:l2,m1:m2,n1:n2,1)**2 &
!                   +0.5*tmp(l1:l2,m1:m2,n1:n2,2)**2 
!            enddo
!          enddo
!
!  add diagonal terms
!
!          tmp(:,:,:,1)=0
!          call divu_2nd(f,tmp(:,:,:,1))                     ! divu
!          call del2_2nd_nof(tmp(:,:,:,1),tmp(:,:,:,3))        ! del2(divu)
!          do i=1,3
!            call der_2nd_nof(f(:,:,:,iux+i-1),tmp(:,:,:,1),i) ! uii
!            call del2_2nd_nof(tmp(:,:,:,1),tmp(:,:,:,2))      ! nabla2 uii
!            sij2=sij2+(tmp(l1:l2,m1:m2,n1:n2,2)-tmp(l1:l2,m1:m2,n1:n2,3)/3.)**2
!          enddo
!
!          epsK_hyper=2*nu_hyper3*sum(sij2)
!        endif
!      endif
!
    endsubroutine hyperviscosity_strict
!***********************************************************************
    subroutine divu_2nd(f,df)
!
!  Calculate divergence of a vector u, get scalar.
!  Accurate to 2nd order.
!
!  23-nov-02/tony: coded
!
      use Cdata
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz) :: df
      real :: fac
!
      df=0.
!
      if (nxgrid/=1) then
         fac=1./(2.*dx)
         df(2:mx-1,:,:) =     df(2:mx-1,:,:) &
                           + ( f(3:mx,:,:,iux)-f(1:mx-2,:,:,iux) ) *fac
      endif
!
      if (nygrid/=1) then
         fac=1./(2.*dy)
         df(:,2:my-1,:) =    df(:,2:my-1,:) &  
                          + ( f(:,3:my,:,iuy)-f(:,1:my-2,:,iuy) )*fac
      endif
!
      if (nzgrid/=1) then
         fac=1./(2.*dz)
         df(:,:,2:mz-1) = df(:,:,2:mz-1) &
                          + (f(:,:,3:mz,iuz)-f(:,:,1:mz-2,iuz))*fac
      endif
!      
    endsubroutine divu_2nd
!***********************************************************************
    subroutine der_2nd_nof(var,tmp,j)
!
!  24-nov-03/nils: coded
!
      use Cdata
!
      real, dimension (mx,my,mz) :: var
      real, dimension (mx,my,mz) :: tmp
      integer :: j
!
      intent (in) :: var,j
      intent (out) :: tmp
!
      tmp=0.

      if (j==1 .and. nxgrid/=1) then
          tmp(     1,:,:) = (-3.*var(1,:,:) &
                             +4.*var(2,:,:) &
                             -1.*var(3,:,:))/(2.*dx) 
          tmp(2:mx-1,:,:) = (-1.*var(1:mx-2,:,:) &
                             +1.*var(3:mx  ,:,:))/(2.*dx) 
          tmp(    mx,:,:) = (+1.*var(mx-2,:,:) &
                             -4.*var(mx-1,:,:) &
                             +3.*var(mx  ,:,:))/(2.*dx) 
      endif
!
      if (j==2 .and. nygrid/=1) then
          tmp(:,     1,:) = (-3.*var(:,1,:) &
                             +4.*var(:,2,:) &
                             -1.*var(:,3,:))/(2.*dy) 
          tmp(:,2:my-1,:) = (-1.*var(:,1:my-2,:) &
                             +1.*var(:,3:my  ,:))/(2.*dy) 
          tmp(:,    my,:) = (+1.*var(:,my-2,:) &
                             -4.*var(:,my-1,:) &
                             +3.*var(:,my  ,:))/(2.*dy) 
      endif
!
      if (j==3 .and. nzgrid/=1) then
          tmp(:,:,     1) = (-3.*var(:,:,1) &
                             +4.*var(:,:,2) &
                             -1.*var(:,:,3))/(2.*dz) 
          tmp(:,:,2:mz-1) = (-1.*var(:,:,1:mz-2) &
                             +1.*var(:,:,3:mz  ))/(2.*dz) 
          tmp(:,:,    mz) = (+1.*var(:,:,mz-2) &
                             -4.*var(:,:,mz-1) &
                             +3.*var(:,:,mz  ))/(2.*dz) 
      endif
!
    end subroutine der_2nd_nof
!***********************************************************************
    subroutine del2v_2nd(f,del2f,k)
!
!  24-nov-03/nils: adapted from del2v
!
      use Cdata
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
!  do the del2 diffusion operator
!
      k1=k-1
      do i=1,3
        call del2_2nd(f,tmp,k1+i)
        del2f(:,:,:,i)=tmp
      enddo
!
    end subroutine del2v_2nd
!***********************************************************************
    subroutine del2_2nd(f,del2f,k)
!
!  Calculate del2 of a scalar, get scalar.
!  Accurate to second order.
!
!  24-nov-03/nils: adapted from del2
!
      use Cdata
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz) :: del2f,d2fd
      integer :: i,k,k1
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
    subroutine del2_2nd_nof(f,del2f)
!
!  Calculate del2 of a scalar, get scalar.
!  Same as del2_2nd but for the case where f is a scalar
!
!  24-nov-03/nils: adapted from del2_2nd
!
      use Cdata
!
      real, dimension (mx,my,mz) :: f
      real, dimension (mx,my,mz) :: del2f,d2fd
      integer :: i
!
      intent (in) :: f
      intent (out) :: del2f
!
      call der2_2nd_nof(f,d2fd,1)
      del2f=d2fd
      call der2_2nd_nof(f,d2fd,2)
      del2f=del2f+d2fd
      call der2_2nd_nof(f,d2fd,3)
      del2f=del2f+d2fd
!
    endsubroutine del2_2nd_nof
!***********************************************************************
    subroutine der2_2nd(f,der2f,i,j)
!
!  Calculate the second derivative of f.
!  Accurate to second order.
!
!  24-nov-03/nils: coded
!
      use Cdata
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
    subroutine der2_2nd_nof(f,der2f,j)
!
!  Same as der2_2nd but for the case where f is a scalar.
!
!  07-jan-04/nils: adapted from der2_2nd
!
      use Cdata
!
      real, dimension (mx,my,mz) :: f
      real, dimension (mx,my,mz) :: der2f
      integer :: j
!
      intent (in) :: f,j
      intent (out) :: der2f
!
      der2f=0.
!
      if (j==1 .and. nxgrid/=1) then
        der2f(2:mx-1,:,:) = (+1.*f(1:mx-2,:,:) &
                             -2.*f(2:mx-1,:,:) &
                             +1.*f(3:mx  ,:,:) ) / (dx**2) 
      endif
!
     if (j==2 .and. nygrid/=1) then
        der2f(:,2:my-1,:) = (+1.*f(:,1:my-2,:) &
                             -2.*f(:,2:my-1,:) &
                             +1.*f(:,3:my  ,:) ) / (dy**2) 
      endif
!
     if (j==3 .and. nzgrid/=1) then
        der2f(:,:,2:mz-1) = (+1.*f(:,:,1:mz-2) &
                             -2.*f(:,:,2:mz-1) &
                             +1.*f(:,:,3:mz  ) ) / (dz**2) 
      endif
!
    endsubroutine der2_2nd_nof
!***********************************************************************
    subroutine graddivu_2nd(f,graddivuf,k)
!
!  Calculate the gradient of the divergence of a vector.
!  Accurate to second order.
!
!  24-nov-03/nils: coded
!
      use Cdata
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,3) :: graddivuf
      real, dimension (mx,my,mz) :: tmp,tmp2
      integer :: k,k1
!
      intent (in) :: f, k
      intent (out) :: graddivuf
!
      graddivuf=0.
!
      k1=k-1
!
!  d/dx(div(u))
!
      call der2_2nd(f,tmp,k1+1,1)
      graddivuf(:,:,:,1)=tmp
      call der_2nd_nof(f(:,:,:,k1+2),tmp,1)
      call der_2nd_nof(tmp,tmp2,2)
      graddivuf(:,:,:,1)=graddivuf(:,:,:,1)+tmp2
      call der_2nd_nof(f(:,:,:,k1+3),tmp,1)
      call der_2nd_nof(tmp,tmp2,3)
      graddivuf(:,:,:,1)=graddivuf(:,:,:,1)+tmp2
!
!  d/dy(div(u))
!
      call der2_2nd(f,tmp,k1+2,2)
      graddivuf(:,:,:,2)=tmp
      call der_2nd_nof(f(:,:,:,k1+1),tmp,1)
      call der_2nd_nof(tmp,tmp2,2)
      graddivuf(:,:,:,2)=graddivuf(:,:,:,2)+tmp2
      call der_2nd_nof(f(:,:,:,k1+3),tmp,2)
      call der_2nd_nof(tmp,tmp2,3)
      graddivuf(:,:,:,2)=graddivuf(:,:,:,2)+tmp2
!
!  d/dz(div(u))
!
      call der2_2nd(f,tmp,k1+3,3)
      graddivuf(:,:,:,3)=tmp
      call der_2nd_nof(f(:,:,:,k1+1),tmp,1)
      call der_2nd_nof(tmp,tmp2,3)
      graddivuf(:,:,:,3)=graddivuf(:,:,:,3)+tmp2
      call der_2nd_nof(f(:,:,:,k1+2),tmp,2)
      call der_2nd_nof(tmp,tmp2,3)
      graddivuf(:,:,:,3)=graddivuf(:,:,:,3)+tmp2
!
    endsubroutine graddivu_2nd
!***********************************************************************

endmodule Hypervisc_strict
