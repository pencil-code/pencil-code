! $Id$
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
!  Derivatives are taken in Fourier space.
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
  use Cdata
  use Cparam
  use Fourier
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
!  20-aug-07/anders: coded
!
      use FArrayManager
!
      if (lroot) call svn_id( &
           "$Id$")
!
!  Set indices for auxiliary variables.
!
      call farray_register_auxiliary('hypvis',ihypvis,vector=3)
!
    endsubroutine register_hypervisc_strict
!***********************************************************************
    subroutine hyperviscosity_strict(f)
!
!  Apply momentum-conserving, symmetric, sixth order hyperviscosity with
!  positive define heating rate (see Haugen & Brandenburg 2004).
!
!  20-aug-2007/anders: coded
!
      use Fourier
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      real, dimension (nx,ny,nz) :: uxhat_re, uxhat_im
      real, dimension (nx,ny,nz) :: uyhat_re, uyhat_im
      real, dimension (nx,ny,nz) :: uzhat_re, uzhat_im
      real, dimension (nx,ny,nz) :: duxhat_re, duxhat_im
      real, dimension (nx,ny,nz) :: duyhat_re, duyhat_im
      real, dimension (nx,ny,nz) :: duzhat_re, duzhat_im
      real, dimension (nx,ny,nz) :: tmp
      real :: kx, ky, kz
      complex :: del2del2del2, del2del2divu
      integer :: ikx, iky, ikz
      real, dimension (nx,3) :: del6u
!
!  Identify version
!
      if (lroot .and. ip<10) call svn_id( &
        "$Id$")
!
!  Derivatives are taken in k-space due to the complicated cross terms.
!
      uxhat_re=f(l1:l2,m1:m2,n1:n2,iuu  )
      uxhat_im=0.0
      uyhat_re=f(l1:l2,m1:m2,n1:n2,iuu+1)
      uyhat_im=0.0
      uzhat_re=f(l1:l2,m1:m2,n1:n2,iuu+2)
      uzhat_im=0.0
!
      if (lshear) then
        call fourier_transform_shear(uxhat_re,uxhat_im)
        call fourier_transform_shear(uyhat_re,uyhat_im)
        call fourier_transform_shear(uzhat_re,uzhat_im)
      else
        call fourier_transform(uxhat_re,uxhat_im)
        call fourier_transform(uyhat_re,uyhat_im)
        call fourier_transform(uzhat_re,uzhat_im)
      endif
!
!  Construct hyperviscous acceleration
!
!    f_visc = mu3/rho*(del2(del2(del2(u)))+1/3*del2(del2(grad(div(u))))
!
      do ikz=1,nz; do iky=1,ny; do ikx=1,nx
!
        if (lshear) then
          kx=kx_fft(ikx)+deltay/Lx*ky_fft(iky+ipy*ny)
          ky=ky_fft(iky+ipy*ny)
          kz=kz_fft(ikz)
        else
          kx=kx_fft(ikx)
          ky=ky_fft(iky+ipy*ny)
          kz=kz_fft(ikz+ipz*nz)
        endif
!
        if (lshear) then
!
          if (nzgrid/=1) then  ! 3-D with shear, FFT order (kz,ky',kx)
!
!  del2(del2(del2(u)))
!
            del2del2del2=-(kx**2+ky**2+kz**2)**3
!
            duxhat_re(ikz,iky,ikx)=del2del2del2*uxhat_re(ikz,iky,ikx)
            duxhat_im(ikz,iky,ikx)=del2del2del2*uxhat_im(ikz,iky,ikx)
            duyhat_re(ikz,iky,ikx)=del2del2del2*uyhat_re(ikz,iky,ikx)
            duyhat_im(ikz,iky,ikx)=del2del2del2*uyhat_im(ikz,iky,ikx)
            duzhat_re(ikz,iky,ikx)=del2del2del2*uzhat_re(ikz,iky,ikx)
            duzhat_im(ikz,iky,ikx)=del2del2del2*uzhat_im(ikz,iky,ikx)
!
!  del2(del2(grad(div(u))))
!
            del2del2divu= -1/3.* & ! i*1/3*del2(del2(divu)) operator
                (kx**2+ky**2+kz**2)**2* &
                (kx*cmplx(uxhat_re(ikz,iky,ikx),uxhat_im(ikz,iky,ikx)) + &
                 ky*cmplx(uyhat_re(ikz,iky,ikx),uyhat_im(ikz,iky,ikx)) + &
                 kz*cmplx(uzhat_re(ikz,iky,ikx),uzhat_im(ikz,iky,ikx)) )
!
            duxhat_re(ikz,iky,ikx)=duxhat_re(ikz,iky,ikx)+ real(kx*del2del2divu)
            duxhat_im(ikz,iky,ikx)=duxhat_im(ikz,iky,ikx)+aimag(kx*del2del2divu)
            duyhat_re(ikz,iky,ikx)=duyhat_re(ikz,iky,ikx)+ real(ky*del2del2divu)
            duyhat_im(ikz,iky,ikx)=duyhat_im(ikz,iky,ikx)+aimag(ky*del2del2divu)
            duzhat_re(ikz,iky,ikx)=duzhat_re(ikz,iky,ikx)+ real(kz*del2del2divu)
            duzhat_im(ikz,iky,ikx)=duzhat_im(ikz,iky,ikx)+aimag(kz*del2del2divu)
!
          else  ! 2-D with shear, FFT order (kx,ky',kz)
!
!  del2(del2(del2(u)))
!
            del2del2del2=-(kx**2+ky**2+kz**2)**3
!
            duxhat_re(ikx,iky,ikz)=del2del2del2*uxhat_re(ikx,iky,ikz)
            duxhat_im(ikx,iky,ikz)=del2del2del2*uxhat_im(ikx,iky,ikz)
            duyhat_re(ikx,iky,ikz)=del2del2del2*uyhat_re(ikx,iky,ikz)
            duyhat_im(ikx,iky,ikz)=del2del2del2*uyhat_im(ikx,iky,ikz)
            duzhat_re(ikx,iky,ikz)=del2del2del2*uzhat_re(ikx,iky,ikz)
            duzhat_im(ikx,iky,ikz)=del2del2del2*uzhat_im(ikx,iky,ikz)
!
!  del2(del2(grad(div(u))))
!
            del2del2divu= -1/3.* & ! i*1/3*del2(del2(divu)) operator
                (kx**2+ky**2+kz**2)**2* &
                (kx*cmplx(uxhat_re(ikx,iky,ikz),uxhat_im(ikx,iky,ikz)) + &
                 ky*cmplx(uyhat_re(ikx,iky,ikz),uyhat_im(ikx,iky,ikz)) + &
                 kz*cmplx(uzhat_re(ikx,iky,ikz),uzhat_im(ikx,iky,ikz)) )
!
            duxhat_re(ikx,iky,ikz)=duxhat_re(ikx,iky,ikz)+ real(kx*del2del2divu)
            duxhat_im(ikx,iky,ikz)=duxhat_im(ikx,iky,ikz)+aimag(kx*del2del2divu)
            duyhat_re(ikx,iky,ikz)=duyhat_re(ikx,iky,ikz)+ real(ky*del2del2divu)
            duyhat_im(ikx,iky,ikz)=duyhat_im(ikx,iky,ikz)+aimag(ky*del2del2divu)
            duzhat_re(ikx,iky,ikz)=duzhat_re(ikx,iky,ikz)+ real(kz*del2del2divu)
            duzhat_im(ikx,iky,ikz)=duzhat_im(ikx,iky,ikz)+aimag(kz*del2del2divu)
!
          endif
!
        else ! No shear
!
          if (nzgrid/=1) then !3D - FFT has put in (kz,kx,ky) order.
!
!  del2(del2(del2(u)))
!
            del2del2del2=-(kx**2+ky**2+kz**2)**3
!
            duxhat_re(ikz,ikx,iky)=del2del2del2*uxhat_re(ikz,ikx,iky)
            duxhat_im(ikz,ikx,iky)=del2del2del2*uxhat_im(ikz,ikx,iky)
            duyhat_re(ikz,ikx,iky)=del2del2del2*uyhat_re(ikz,ikx,iky)
            duyhat_im(ikz,ikx,iky)=del2del2del2*uyhat_im(ikz,ikx,iky)
            duzhat_re(ikz,ikx,iky)=del2del2del2*uzhat_re(ikz,ikx,iky)
            duzhat_im(ikz,ikx,iky)=del2del2del2*uzhat_im(ikz,ikx,iky)
!
!  del2(del2(grad(div(u))))
!
            del2del2divu= -1/3.* & ! i*1/3*del2(del2(divu)) operator
                 (kx**2+ky**2+kz**2)**2* &
                 (kx*cmplx(uxhat_re(ikz,ikx,iky),uxhat_im(ikz,ikx,iky)) + &
                  ky*cmplx(uyhat_re(ikz,ikx,iky),uyhat_im(ikz,ikx,iky)) + &
                  kz*cmplx(uzhat_re(ikz,ikx,iky),uzhat_im(ikz,ikx,iky)) )
!
            duxhat_re(ikz,ikx,iky)=duxhat_re(ikz,ikx,iky)+ real(kx*del2del2divu)
            duxhat_im(ikz,ikx,iky)=duxhat_im(ikz,ikx,iky)+aimag(kx*del2del2divu)
            duyhat_re(ikz,ikx,iky)=duyhat_re(ikz,ikx,iky)+ real(ky*del2del2divu)
            duyhat_im(ikz,ikx,iky)=duyhat_im(ikz,ikx,iky)+aimag(ky*del2del2divu)
            duzhat_re(ikz,ikx,iky)=duzhat_re(ikz,ikx,iky)+ real(kz*del2del2divu)
            duzhat_im(ikz,ikx,iky)=duzhat_im(ikz,ikx,iky)+aimag(kz*del2del2divu)
!
            else !2D - order (kx,ky,kz) order
!
!  del2(del2(del2(u)))
!
            del2del2del2=-(kx**2+ky**2+kz**2)**3
!
            duxhat_re(ikx,iky,ikz)=del2del2del2*uxhat_re(ikx,iky,ikz)
            duxhat_im(ikx,iky,ikz)=del2del2del2*uxhat_im(ikx,iky,ikz)
            duyhat_re(ikx,iky,ikz)=del2del2del2*uyhat_re(ikx,iky,ikz)
            duyhat_im(ikx,iky,ikz)=del2del2del2*uyhat_im(ikx,iky,ikz)
            duzhat_re(ikx,iky,ikz)=del2del2del2*uzhat_re(ikx,iky,ikz)
            duzhat_im(ikx,iky,ikz)=del2del2del2*uzhat_im(ikx,iky,ikz)
!
!  del2(del2(grad(div(u))))
!
            del2del2divu= -1/3.* & ! i*1/3*del2(del2(divu)) operator
                 (kx**2+ky**2+kz**2)**2* &
                 (kx*cmplx(uxhat_re(ikx,iky,ikz),uxhat_im(ikx,iky,ikz)) + &
                  ky*cmplx(uyhat_re(ikx,iky,ikz),uyhat_im(ikx,iky,ikz)) + &
                  kz*cmplx(uzhat_re(ikx,iky,ikz),uzhat_im(ikx,iky,ikz)) )
!
            duxhat_re(ikx,iky,ikz)=duxhat_re(ikx,iky,ikz)+ real(kx*del2del2divu)
            duxhat_im(ikx,iky,ikz)=duxhat_im(ikx,iky,ikz)+aimag(kx*del2del2divu)
            duyhat_re(ikx,iky,ikz)=duyhat_re(ikx,iky,ikz)+ real(ky*del2del2divu)
            duyhat_im(ikx,iky,ikz)=duyhat_im(ikx,iky,ikz)+aimag(ky*del2del2divu)
            duzhat_re(ikx,iky,ikz)=duzhat_re(ikx,iky,ikz)+ real(kz*del2del2divu)
            duzhat_im(ikx,iky,ikz)=duzhat_im(ikx,iky,ikz)+aimag(kz*del2del2divu)
!
            endif
!
        endif
!
      enddo; enddo; enddo
!
!  Inverse transform (to real space). The Viscosity module multiplies the
!  result with mu3/rho.
!
      if (lshear) then
        call fourier_transform_shear(duxhat_re,duxhat_im,linv=.true.)
        f(l1:l2,m1:m2,n1:n2,ihypvis)=duxhat_re
        call fourier_transform_shear(duyhat_re,duyhat_im,linv=.true.)
        f(l1:l2,m1:m2,n1:n2,ihypvis+1)=duyhat_re
        call fourier_transform_shear(duzhat_re,duzhat_im,linv=.true.)
        f(l1:l2,m1:m2,n1:n2,ihypvis+2)=duzhat_re
      else
        call fourier_transform(duxhat_re,duxhat_im,linv=.true.)
        f(l1:l2,m1:m2,n1:n2,ihypvis)=duxhat_re
        call fourier_transform(duyhat_re,duyhat_im,linv=.true.)
        f(l1:l2,m1:m2,n1:n2,ihypvis+1)=duyhat_re
        call fourier_transform(duzhat_re,duzhat_im,linv=.true.)
        f(l1:l2,m1:m2,n1:n2,ihypvis+2)=duzhat_re
      endif
!
!  The imaginary part of the functions should be zero once transformed
!  back to real space, since the differential operators preserve
!  complex conjugarity of the scales (kx,ky,kz) and (-kx,-ky,-kz)
!
      if (ip<=10) then
        print'(A,4e15.5)', &
            'hypervisc_strict: min(duxhat_re), max(duxhat_re), '// &
            'min(duxhat_im), max(duxhat_im)', &
            minval(duxhat_re), maxval(duxhat_re), &
            minval(duxhat_im), maxval(duxhat_im)
        print'(A,4e15.5)', &
            'hypervisc_strict: min(duyhat_re), max(duyhat_re), '// &
            'min(duyhat_im), max(duyhat_im)', &
            minval(duyhat_re), maxval(duyhat_re), &
            minval(duyhat_im), maxval(duyhat_im)
        print'(A,4e15.5)', &
            'hypervisc_strict: min(duzhat_re), max(duzhat_re), '// &
            'min(duzhat_im), max(duzhat_im)', &
            minval(duzhat_re), maxval(duzhat_re), &
            minval(duzhat_im), maxval(duzhat_im)
      endif
!
    endsubroutine hyperviscosity_strict
!***********************************************************************
endmodule Hypervisc_strict
