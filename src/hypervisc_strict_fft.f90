! $Id: hypervisc_strict_fft.f90,v 1.2 2007-08-21 12:41:19 ajohan Exp $

!
!  This module applies a sixth order hyperviscosity to the equation
!  of motion (following Haugen & Brandenburg 2004).
!
!  The rate of strain tensor
!    S^(3) = (-nab^2)^2*S
!  is a high order generalisation of the first order operator
!    2*S_ij = u_i,j + u_j,i - 2/3*delta_ij*div(u)
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

  use Cdata
  use Cparam
  use Fourier
  use Messages

  implicit none

  contains

!***********************************************************************
    subroutine register_hypervisc_strict()
!
!  Set up indices for hyperviscosity auxiliary slots.
!
!  20-aug-07/anders: coded
!
      use Mpicomm, only: stop_it
!
      logical, save :: first=.true.
!
      if (.not. first) call stop_it('register_particles: called twice')
      first = .false.
!
      if (lroot) call cvs_id( &
           "$Id: hypervisc_strict_fft.f90,v 1.2 2007-08-21 12:41:19 ajohan Exp $")
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
    subroutine hyperviscosity_strict(f,j)
!
!  Apply momentum-conserving, symmetric, sixth order hyperviscosity with
!  positive define heating rate (see Haugen & Brandenburg 2004).
!
!  20-aug-2007/anders: coded
!
      use Fourier
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: j
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
!
!  Identify version
!
      if (lroot .and. ip<10) call cvs_id( &
        "$Id: hypervisc_strict_fft.f90,v 1.2 2007-08-21 12:41:19 ajohan Exp $")
!
!  Derivatives are taken in k-space due to the complicated cross terms.
!
      uxhat_re=f(l1:l2,m1:m2,n1:n2,j  )
      uxhat_im=0.0
      uyhat_re=f(l1:l2,m1:m2,n1:n2,j+1)
      uyhat_im=0.0
      uzhat_re=f(l1:l2,m1:m2,n1:n2,j+2)
      uzhat_im=0.0
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
!  Construct hyperviscious acceleration
!
!    f_visc = mu3/rho*(del2(del2(del2(u)))+1/3*del2(del2(grad(div(u))))
!
      do ikz=1,nz; do iky=1,ny; do ikx=1,nx
!
        if (lshear) then
          if (nzgrid/=1) then
            kx=kx_fft(ikz+ipz*nz)+deltay/Lx*ky_fft(iky+ipy*ny)
            ky=ky_fft(iky+ipy*ny)
            kz=kz_fft(ikx)
          else
            kx=kx_fft(ikx)+deltay/Lx*ky_fft(iky+ipy*ny)
            ky=ky_fft(iky+ipy*ny)
            kz=kz_fft(ikz)
          endif
        else
          kx=kx_fft(ikx)
          ky=ky_fft(iky+ipy*ny)
          kz=kz_fft(ikz+ipz*nz)        
        endif
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
        del2del2divu=-1/3.* & ! i*1/3*del2(del2(divu)) operator
            (kx**2+ky**2+kz**2)**2* &
            (kx*cmplx(uxhat_re(ikz,iky,ikx),uxhat_im(ikz,iky,ikx)) + &
             ky*cmplx(uyhat_re(ikz,iky,ikx),uyhat_im(ikz,iky,ikx)) + &
             kz*cmplx(uzhat_re(ikz,iky,ikx),uzhat_im(ikz,iky,ikx)) )

        duxhat_re(ikz,iky,ikx)=duxhat_re(ikz,iky,ikx) +  real(kx*del2del2divu)
        duxhat_im(ikz,iky,ikx)=duxhat_im(ikz,iky,ikx) + aimag(kx*del2del2divu)
        duyhat_re(ikz,iky,ikx)=duyhat_re(ikz,iky,ikx) +  real(ky*del2del2divu)
        duyhat_im(ikz,iky,ikx)=duyhat_im(ikz,iky,ikx) + aimag(ky*del2del2divu)
        duzhat_re(ikz,iky,ikx)=duzhat_re(ikz,iky,ikx) +  real(kz*del2del2divu)
        duzhat_im(ikz,iky,ikx)=duzhat_im(ikz,iky,ikx) + aimag(kz*del2del2divu)
!        print'(3f8.1,2e15.5)', kx, ky, kz, duxhat_re(ikx,iky,ikz), duxhat_im(ikx,iky,ikz)
!
      enddo; enddo; enddo
!
!  Inverse transform (to real space). The Viscosity module multiplies the
!  result with mu3/rho.
!
      if (lshear) then
        call fourier_transform_shear(duxhat_re,duxhat_im,linv=.true.)
        f(l1:l2,m1:m2,n1:n2,ihyper)=duxhat_re
        call fourier_transform_shear(duyhat_re,duyhat_im,linv=.true.)
        f(l1:l2,m1:m2,n1:n2,ihyper+1)=duyhat_re
        call fourier_transform_shear(duzhat_re,duzhat_im,linv=.true.)
        f(l1:l2,m1:m2,n1:n2,ihyper+2)=duzhat_re
      else
        call fourier_transform(duxhat_re,duxhat_im,linv=.true.)
        f(l1:l2,m1:m2,n1:n2,ihyper)=duxhat_re
        call fourier_transform(duyhat_re,duyhat_im,linv=.true.)
        f(l1:l2,m1:m2,n1:n2,ihyper+1)=duyhat_re
        call fourier_transform(duzhat_re,duzhat_im,linv=.true.)
        f(l1:l2,m1:m2,n1:n2,ihyper+2)=duzhat_re
      endif
!
!  The imaginary part of the functions should be zero once transformed
!  back to real space, sine the differential operators preserve
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

