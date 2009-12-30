! $Id$
!
!  This module take care of WENO (weighted essentially non oscillatory)
!  transport.
!
!  The key idea of ENO schemes is to use the ‘‘smoothest’’ stencil among
!  several candidates to approximate the fluxes at cell boundaries to a high
!  order accuracy and at the same time to avoid spurious oscillations.
!
!  See e.g. ``Efficient Implementation of Weighted ENO Schemes'' by Jiang†& Shu
!  1996.
!
module WENO_transport
!
  implicit none
!
  private 
!
  public :: weno_transp
!
  integer, parameter :: nw=3
  real, allocatable, dimension(:,:) :: f, df
!
  contains
!***********************************************************************
    subroutine weno_transp(fq,m,n,iq,iq1,iux,iuy,iuz,dq,dx_1,dy_1,dz_1)
!
!  Solve the equation dq/dt + div(u*q) = 0 using the WENO method.
!
!  29-dec-09/evghenii+anders: coded
!
      real, dimension(:,:,:,:), intent(in ) :: fq
      integer, intent(in) :: m, n
      integer, intent(in) :: iq, iq1, iux, iuy, iuz
      real, dimension(:), intent(out) :: dq
      real, dimension(:), intent(in)  :: dx_1, dy_1, dz_1
!
      call weno5(fq,m,n,iq,iq1,iux,iuy,iuz,dq,dx_1,dy_1,dz_1)
!
    endsubroutine weno_transp
!***********************************************************************
    subroutine weno5(fq,m,n,iq,iq1,iux,iuy,iuz,dq_out,dx_1,dy_1,dz_1)
!
!  Fifth order implementation of WENO scheme. This wrapper takes care of all
!  three directions.
!
!  29-dec-09/evghenii: coded
!
      real, dimension(:,:,:,:), intent(in ) :: fq
      integer, intent(in) :: m, n
      integer, intent(in) :: iq, iq1, iux, iuy, iuz
      real, dimension(:), intent(out) :: dq_out
      real, dimension(:), intent(in)  :: dx_1, dy_1, dz_1
!
      real, allocatable, dimension(:) :: vsig, dq, dl_1
      integer :: i, mx, my, mz, nghost
!
!  Possible to multiply transported variable by another variables, e.g. to
!  transport the momentum rho*u.
!
      if (iq1==0) then
        print*, 'weno5: iq1 is zero - are you using ldensity_nolog=T ?'
        STOP
      endif
!
!  Educated guess at grid dimensions.
!
      mx=size(fq(:,1,1,iq))
      my=size(fq(1,:,1,iq))
      mz=size(fq(1,1,:,iq))
!
!  Allocate arrays.
!
      if (.not. allocated(vsig)) allocate(vsig(mx))
      if (.not. allocated(dq))   allocate(dq  (mx))
      if (.not. allocated(dl_1)) allocate(dl_1(mx))
      if (.not. allocated(f))    allocate( f  (-nw:nw,mx))
      if (.not. allocated(df))   allocate(df  (-nw:nw,mx))
!
      dq(:)=0.0
!
!  WENO transport in x-direction.
!
      vsig=abs(fq(:,m,n, iux))
      vsig=max( &
          cshift(vsig,-3), cshift(vsig,-2), cshift(vsig,-1), cshift(vsig, 0), &
          cshift(vsig,+1), cshift(vsig,+2), cshift(vsig,+3))
!      
      do i=-nw-1+1,nw-1
        if (iq1<0) then
          df(i+1,:)=vsig                   *cshift(fq(:,m,n,iq),i)
          f (i+1,:)=cshift(fq(:,m,n,iux),i)*cshift(fq(:,m,n,iq),i)
        else
          df(i+1,:)=vsig                   *cshift(fq(:,m,n,iq)*fq(:,m,n,iq1),i)
          f (i+1,:)=cshift(fq(:,m,n,iux),i)*cshift(fq(:,m,n,iq)*fq(:,m,n,iq1),i)
        endif 
      enddo
!      
      call weno5_1d(dq,dx_1)
!
!  WENO transport in y-direction.
!
      vsig=max( &
          abs(fq(:,m-3,n,iuy)), abs(fq(:,m-2,n,iuy)), abs(fq(:,m-1,n,iuy)), &
          abs(fq(:,m  ,n,iuy)), &
          abs(fq(:,m+1,n,iuy)), abs(fq(:,m+2,n,iuy)), abs(fq(:,m+3,n,iuy)) )
!
      do i=-nw-1+1,nw-1
        if (iq1<0) then
          df(i+1,:)=vsig           *fq(:,m+i,n,iq)
          f (i+1,:)=fq(:,m+i,n,iuy)*fq(:,m+i,n,iq)
        else
          df(i+1,:)=vsig           *fq(:,m+i,n,iq)*fq(:,m+i,n,iq1)
          f (i+1,:)=fq(:,m+i,n,iuy)*fq(:,m+i,n,iq)*fq(:,m+i,n,iq1)
        endif
      enddo
!      
      dl_1(:)=dy_1(m)
!      call weno5_1d(dq,dl_1)
!
!  WENO transport in z-direction.
!
      vsig=max( &
          abs(fq(:,m,n-3,iuz)), abs(fq(:,m,n-2,iuz)), abs(fq(:,m,n-1,iuz)), &
          abs(fq(:,m,n  ,iuz)), &
          abs(fq(:,m,n+1,iuz)), abs(fq(:,m,n+2,iuz)), abs(fq(:,m,n+3,iuz)) )
!
      do i=-nw-1+1,nw-1
        if (iq1<0) then
          df(i+1,:)=vsig           *fq(:,m,n+i,iq)
          f (i+1,:)=fq(:,m,n+i,iuz)*fq(:,m,n+i,iq)
        else
          df(i+1,:)=vsig           *fq(:,m,n+i,iq)*fq(:,m,n+i,iq1)
          f (i+1,:)=fq(:,m,n+i,iuz)*fq(:,m,n+i,iq)*fq(:,m,n+i,iq1)
        endif
      enddo
!      
      dl_1=dz_1(n)
!      call weno5_1d(dq,dl_1)
!
!  Return transport pencil without ghost cells.
!
      nghost=(mx-size(dq_out))/2
      dq_out(:)=dq(nghost+1:mx-nghost)
!
!  Deallocate arrays.
!
      if (allocated(vsig)) deallocate(vsig)
      if (allocated(dq))   deallocate(dq)
      if (allocated(dl_1)) deallocate(dl_1)
      if (allocated(f))    deallocate(f)
      if (allocated(df))   deallocate(df)
!      
    endsubroutine weno5
!***********************************************************************
    subroutine weno5_1d(dq,dx_1)
!
!  Fifth order implementation of WENO scheme (1-D version).
!
!  29-dec-09/evghenii: coded
!
      real, dimension(:), intent(inout) :: dq
      real, dimension(:), intent(in) :: dx_1
!
      real, parameter :: WENO_POW = 2
      real, parameter :: WENO_EPS = 1.0e-6
      real, dimension(size(dq)) :: b1, b2, b3
      real, dimension(size(dq)) :: wh1, wh2, wh3, wh
      real, dimension(size(dq)) :: fh1, fh2, fh3
      real, dimension(size(dq))         :: flux
      real :: g1, g2, g3
!
      f = 0.5 * (f + df)
!      
      b1(:) = &
          13.0/12.0 * (f(-2,:) - 2.0*f(-1,:) +     f(0,:))**2 + &
          1.0 / 4.0 * (f(-2,:) - 4.0*f(-1,:) + 3.0*f(0,:))**2
!         
      b2(:) = &
          13.0/12.0 * (f(-1,:) - 2.0*f(0,:) + f(+1,:))**2 + &
          1.0 / 4.0 * (f(-1,:)              - f(+1,:))**2
!
      b3(:) = &
          13.0/12.0 * (    f(0,:) - 2.0*f(+1,:) + f(+2,:))**2 + &
          1.0 / 4.0 * (3.0*f(0,:) - 4.0*f(+1,:) + f(+2,:))**2
!      
      g1 = 0.1
      g2 = 0.6
      g3 = 0.3
!
      wh1 = g1/(WENO_EPS + b1)**WENO_POW
      wh2 = g2/(WENO_EPS + b2)**WENO_POW
      wh3 = g3/(WENO_EPS + b3)**WENO_POW
!
      wh = 1.0/(wh1 + wh2 + wh3)
      wh1 = wh1 * wh
      wh2 = wh2 * wh
      wh3 = wh3 * wh
!
      fh1(:) =  1.0/3.0*f(-2,:) - 7.0/6.0*f(-1,:) + 11.0/6.0*f( 0,:)
      fh2(:) = -1.0/6.0*f(-1,:) + 5.0/6.0*f( 0,:) + 1.0 /3.0*f(+1,:)
      fh3(:) =  1.0/3.0*f( 0,:) + 5.0/6.0*f(+1,:) - 1.0 /6.0*f(+2,:)
!      
      flux = wh1*fh1 + wh2*fh2 + wh3*fh3
!      
      f = f - df
!      
      b1(:) = &
          13.0/12.0 * (f(+3,:) - 2.0*f(+2,:) +     f(+1,:))**2 + &
          1.0 / 4.0 * (f(+3,:) - 4.0*f(+2,:) + 3.0*f(+1,:))**2
!      
      b2(:) = &
          13.0/12.0 * (f(+2,:) - 2.0*f(+1,:) + f(0,:))**2 + &
          1.0 / 4.0 * (f(+2,:)               - f(0,:))**2
!
      b3(:) = &
          13.0/12.0 * (    f(+1,:) - 2.0*f(0,:) + f(-1,:))**2 +  &
          1.0 / 4.0 * (3.0*f(+1,:) - 4.0*f(0,:) + f(-1,:))**2
!      
      wh1 = g1/(WENO_EPS + b1)**WENO_POW
      wh2 = g2/(WENO_EPS + b2)**WENO_POW
      wh3 = g3/(WENO_EPS + b3)**WENO_POW
!    
      wh  = 1.0/(wh1 + wh2 + wh3)
      wh1 = wh1*wh
      wh2 = wh2*wh
      wh3 = wh3*wh
!      
      fh1(:) =  1.0/3.0*f(+3,:) - 7.0/6.0*f(+2,:) + 11.0/6.0*f(+1,:)
      fh2(:) = -1.0/6.0*f(+2,:) + 5.0/6.0*f(+1,:) + 1.0/ 3.0*f( 0,:)
      fh3(:) =  1.0/3.0*f(+1,:) + 5.0/6.0*f( 0,:) - 1.0/ 6.0*f(-1,:)
!    
      flux = flux + wh1*fh1 + wh2*fh2 + wh3*fh3
!         
      dq = dq + (cshift(flux, 1) - flux)*dx_1
!  
    endsubroutine weno5_1d
!***********************************************************************
endmodule WENO_transport
