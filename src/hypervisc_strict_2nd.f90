! $Id$
!
!  This module applies a sixth order hyperviscosity to the equation
!  of motion (following Haugen & Brandenburg 2004). This hyperviscosity
!  ensures that the energy dissipation rate is positive define everywhere.
!
!  The rate of strain tensor
!
!    S^(3) = (-nab^2)^2*S
!
!  is a high order generalisation of the first order rate of strain tensor
!
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
!
  use Cparam
  use Cdata
  use Messages
!
  implicit none
!
  include 'hypervisc_strict.h'
!
  contains
!***********************************************************************
    subroutine register_hypervisc_strict()
!
!  Set up indices for hyperviscosity auxiliary slots.
!
      use FArrayManager
!
      if (lroot) call svn_id( &
           "$Id$")
!
!  Set indices for auxiliary variables
!
      call farray_register_auxiliary('hypvis',ihypvis,vector=3)
!
    endsubroutine register_hypervisc_strict
!***********************************************************************
    subroutine hyperviscosity_strict(f)
!
!  Apply momentum-conserving, symmetric, sixth order hyperviscosity with
!  positive definite heating rate (see Haugen & Brandenburg 2004).
!
!  The rate of strain tensor of order 3 is defined as:
!
!               /-1/3*div(u)+dux/dx  1/2*(dux/dy+duy/dx) 1/2*(dux/dz+duz/dx) \
!               |                                                            |
!  S^(3)=(d2)^2 |1/2*(dux/dy+duy/dx) -1/3*div(u)+duy/dy  1/2*(duy/dz+duz/dy) |
!               |                                                            |
!               \1/2*(dux/dz+duz/dx) 1/2*(duy/dz+duz/dy) -1/3*div(u)+2*duz/dz/
!
!  where d2 is the Laplacian operator d2=(d^2/dx^2+d2/dy^2+d2/dz^2).
!
!  To avoid communicating ghost zones after each operator, we use
!  derivatives that are second order in space.
!
!  24-nov-03/nils: coded
!
      use SharedVariables, only: get_shared_variable
      use Sub
      use Mpicomm, only: stop_it
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      real, dimension (mx,my,mz,3) :: tmp, tmp2
      real, dimension (mx,my,mz) :: tmp3, tmp4
      integer :: i, j, ierr
      logical, save :: lfirstcall=.true.
      logical, pointer, save :: lvisc_hyper3_nu_const_strict
!
!  Calculate del2(del2(del2(u))), accurate to second order.
!
      call del2v_2nd(f,tmp,iux)
      f(:,:,:,ihypvis:ihypvis+2)=tmp
      call del2v_2nd(f,tmp,ihypvis)
      f(:,:,:,ihypvis:ihypvis+2)=tmp
      call del2v_2nd(f,tmp,ihypvis)
!
!  Calculate del2(del2(grad(div(u)))), accurate to second order.
!  Probably gives zero derivative at the Nyquist scale, but the del2^3
!  term above gives dissipation at this scale.
!
      call graddivu_2nd(f,tmp2,iux)
      f(:,:,:,ihypvis:ihypvis+2)=tmp2
      call del2v_2nd(f,tmp2,ihypvis)
      f(:,:,:,ihypvis:ihypvis+2)=tmp2
      call del2v_2nd(f,tmp2,ihypvis)
!
!  Add the two terms.
!
      f(:,:,:,ihypvis:ihypvis+2)=tmp+tmp2/3.
!
!  For constant nu we also need the term [2*nu*S^3].grad(lnrho).
!
      if (lfirstcall) then
        call get_shared_variable('lvisc_hyper3_nu_const_strict', &
            lvisc_hyper3_nu_const_strict,ierr)
        if (ierr/=0) call stop_it("hyperviscosity_strict: "//&
            "problem getting shared varable lvisc_hyper3_nu_const_strict")
        if (ip<10) &
            print*, 'hyperviscosity_strict: lvisc_hyper3_nu_const_strict=', &
            lvisc_hyper3_nu_const_strict
        lfirstcall=.false.
      endif
!
!  First we calculate the gradient of the logarithmic density (second order
!  accuracy).
!
      if (lvisc_hyper3_nu_const_strict) then
        call grad_2nd(f,ilnrho,tmp)
        if (ldensity_nolog) then
          tmp3=exp(f(:,:,:,irho))
          do i=1,3; tmp(:,:,:,i)=tmp(:,:,:,i)/tmp3; enddo
        endif
!
!  Add [(del2)^2(u_i,j+u_j,i)].grad(lnrho) to hyperviscosity.
!
        do i=1,3; do j=1,3
          call der_2nd(f,iux-1+i,tmp3,j)
          call del2_2nd_nof(tmp3,tmp4)
          call del2_2nd_nof(tmp4,tmp3)
          f(:,:,:,ihypvis-1+i)=f(:,:,:,ihypvis-1+i)+tmp3*tmp(:,:,:,j)
          f(:,:,:,ihypvis-1+j)=f(:,:,:,ihypvis-1+j)+tmp3*tmp(:,:,:,i)
        enddo; enddo
!
!  Add -2/3*[(del2)^2delta_ij*div(u)].grad(lnrho) term.
!
        call div_2nd(f,iux,tmp3)
        call del2_2nd_nof(tmp3,tmp4)
        call del2_2nd_nof(tmp4,tmp3)
        f(:,:,:,ihypvis  )=f(:,:,:,ihypvis  )-2/3.*tmp3*tmp(:,:,:,1)
        f(:,:,:,ihypvis+1)=f(:,:,:,ihypvis+1)-2/3.*tmp3*tmp(:,:,:,2)
        f(:,:,:,ihypvis+2)=f(:,:,:,ihypvis+2)-2/3.*tmp3*tmp(:,:,:,3)
      endif
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
    subroutine der_2nd(f,k,df,j)
!
!  Calculate derivative of scalar, accurate to 2nd order.
!
!  13-sep-07/anders: adapted from div_2nd
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz) :: df
      integer :: j,k
!
      real :: fac
!
      df=0.0
!
      if (j==1 .and. nxgrid/=1) then
        fac=1./(2.*dx)
        df(2:mx-1,:,:) = df(2:mx-1,:,:) &
                         + ( f(3:mx,:,:,k)-f(1:mx-2,:,:,k) ) *fac
      endif
!
      if (j==2 .and. nygrid/=1) then
        fac=1./(2.*dy)
        df(:,2:my-1,:) = df(:,2:my-1,:) &
                         + ( f(:,3:my,:,k)-f(:,1:my-2,:,k) )*fac
      endif
!
      if (j==3 .and. nzgrid/=1) then
        fac=1./(2.*dz)
        df(:,:,2:mz-1) = df(:,:,2:mz-1) &
                         + ( f(:,:,3:mz,k)-f(:,:,1:mz-2,k) )*fac
      endif
!
    endsubroutine der_2nd
!***********************************************************************
    subroutine div_2nd(f,k,df)
!
!  Calculate divergence of a vector, accurate to 2nd order.
!
!  23-nov-02/tony: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz) :: df
      integer :: k
!
      real :: fac
!
      df=0.0
!
      if (nxgrid/=1) then
        fac=1./(2.*dx)
        df(2:mx-1,:,:) = df(2:mx-1,:,:) &
                         + ( f(3:mx,:,:,k  )-f(1:mx-2,:,:,k  ) ) *fac
      endif
!
      if (nygrid/=1) then
        fac=1./(2.*dy)
        df(:,2:my-1,:) = df(:,2:my-1,:) &
                         + ( f(:,3:my,:,k+1)-f(:,1:my-2,:,k+1) )*fac
      endif
!
      if (nzgrid/=1) then
        fac=1./(2.*dz)
        df(:,:,2:mz-1) = df(:,:,2:mz-1) &
                         + ( f(:,:,3:mz,k+2)-f(:,:,1:mz-2,k+2) )*fac
      endif
!
    endsubroutine div_2nd
!***********************************************************************
    subroutine grad_2nd(f,k,df)
!
!  Calculate gradient of scalar, accurate to second order.
!
!  13-sep-07/anders: adapted from div_2nd
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,3) :: df
      integer :: k
!
      real :: fac
!
      df=0.0
!
      if (nxgrid/=1) then
        fac=1./(2.*dx)
        df(2:mx-1,:,:,1) = (f(3:mx,:,:,k)-f(1:mx-2,:,:,k))*fac
      endif
!
      if (nygrid/=1) then
        fac=1./(2.*dy)
        df(:,2:my-1,:,2) = (f(:,3:my,:,k)-f(:,1:my-2,:,k))*fac
      endif
!
      if (nzgrid/=1) then
        fac=1./(2.*dz)
        df(:,:,2:mz-1,3) = (f(:,:,3:mz,k)-f(:,:,1:mz-2,k))*fac
      endif
!
    endsubroutine grad_2nd
!***********************************************************************
    subroutine der_2nd_nof(var,tmp,j)
!
!  24-nov-03/nils: coded
!
      real, dimension (mx,my,mz) :: var
      real, dimension (mx,my,mz) :: tmp
      integer :: j
!
      intent (in) :: var,j
      intent (out) :: tmp
!
      tmp=0.
!
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
    endsubroutine der_2nd_nof
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
!  Calculate del2 of a scalar, get scalar.
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
    subroutine del2_2nd_nof(f,del2f)
!
!  Calculate Laplacian of a scalar, get scalar.
!
!  24-nov-03/nils: adapted from del2_2nd
!
      real, dimension (mx,my,mz) :: f
      real, dimension (mx,my,mz) :: del2f,d2fd
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
!  Calculate the second derivative of f.
!  Accurate to second order.
!
!  07-jan-04/nils: adapted from der2_2nd
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
