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
  real, allocatable, dimension(:,:) :: f, df    ! Why module variables?
!
  contains
!***********************************************************************
    subroutine weno_transp(fq,m,n,iq,iq1,iux,iuy,iuz,dq,dx_1,dy_1,dz_1,ref,ref1)
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
      real, dimension(:), optional, intent(in) :: ref,ref1
!
      call weno5(fq,m,n,iq,iq1,iux,iuy,iuz,dq,dx_1,dy_1,dz_1,ref,ref1)
!
    endsubroutine weno_transp
!***********************************************************************
    subroutine weno5(fq,m,n,iq,iq1,iux,iuy,iuz,dq_out,dx_1,dy_1,dz_1,ref,ref1)
!
!  Fifth order implementation of WENO scheme. This wrapper takes care of all
!  three directions.
!
!  29-dec-09/evghenii: coded
!  19-feb-15/MR: adapted to use of reference state; for x direction ghost
!                values of it still missing
!
      real, dimension(:,:,:,:), intent(in ) :: fq
      integer, intent(in) :: m, n
      integer, intent(in) :: iq, iq1, iux, iuy, iuz
      real, dimension(:), intent(out) :: dq_out
      real, dimension(:), intent(in)  :: dx_1, dy_1, dz_1
      real, dimension(:), optional, intent(in) :: ref,ref1
!
      real, dimension(size(fq,1)) :: vsig, dq, fl, fr, tmp
      integer :: i, mx, nghost
      logical :: lref,lref1
!
!  Possible to multiply transported variable by another variables, e.g. to
!  transport the momentum rho*u.
!
      if (iq1==0) then
        print*, 'weno5: iq1 is zero - are you using ldensity_nolog=T ?'
        STOP
      endif
!
      lref=present(ref); lref1=present(ref1)
!
!  Educated guess at grid dimensions.
!
      mx=size(fq,1)
      nghost=(mx-size(dq_out))/2
!
!  Allocate arrays.
!
      if (.not. allocated(f))  allocate( f(-nw:nw,mx))
      if (.not. allocated(df)) allocate(df(-nw:nw,mx))
!
!  WENO transport in x-direction.
!
      vsig=abs(fq(:,m,n, iux))
      vsig=max( &
          cshift(vsig,-3), cshift(vsig,-2), cshift(vsig,-1), cshift(vsig, 0), &
          cshift(vsig,+1), cshift(vsig,+2), cshift(vsig,+3))
!
      do i=-nw-1+1,nw-1

        tmp = fq(:,m,n,iq)
        if (lref) tmp(nghost+1:mx-nghost)=tmp(nghost+1:mx-nghost)+ref  ! ghost zones not corrected!

        if (iq1>0) then
          fr = fq(:,m,n,iq1)
          if (lref1) fr(nghost+1:mx-nghost) = fr(nghost+1:mx-nghost) + ref1
          tmp=tmp*fr
        endif
 
        df(i+1,:)=vsig                   *cshift(tmp,i)
        f (i+1,:)=cshift(fq(:,m,n,iux),i)*cshift(tmp,i)

      enddo
!
      call weno5_1d(fl)
!
!  Time derivative for x-transport.
!
      fr = cshift(fl,1)
      dq = -(fl - fr) * dx_1
!
!  WENO transport in y-direction.
!
      vsig=max( &
          abs(fq(:,m-3,n,iuy)), abs(fq(:,m-2,n,iuy)), abs(fq(:,m-1,n,iuy)), &
          abs(fq(:,m  ,n,iuy)), &
          abs(fq(:,m+1,n,iuy)), abs(fq(:,m+2,n,iuy)), abs(fq(:,m+3,n,iuy)) )
!
!  Left fluxes.
!
      do i=-nw-1+1,nw-1

        tmp = fq(:,m+i,n,iq)
        if (lref) tmp(nghost+1:mx-nghost)=tmp(nghost+1:mx-nghost)+ref

        if (iq1>0) then
          fr = fq(:,m+i,n,iq1)
          if (lref1) fr(nghost+1:mx-nghost) = fr(nghost+1:mx-nghost) + ref1
          tmp=tmp*fr
        endif
 
        df(i+1,:)=vsig           *tmp
        f (i+1,:)=fq(:,m+i,n,iuy)*tmp

      enddo
!
      call weno5_1d(fl)
!
!  Right fluxes.
!
      do i=-nw-1+1,nw-1

        tmp = fq(:,m+i+1,n,iq)
        if (lref) tmp(nghost+1:mx-nghost)=tmp(nghost+1:mx-nghost)+ref

        if (iq1>0) then
          fr=fq(:,m+i+1,n,iq1)
          if (lref1) fr(nghost+1:mx-nghost) = fr(nghost+1:mx-nghost) + ref1
          tmp=tmp*fr
        endif
 
        df(i+1,:)=vsig             *tmp
        f (i+1,:)=fq(:,m+i+1,n,iuy)*tmp

      enddo
!
      call weno5_1d(fr)
!
!  Time derivative for y-transport.
!
      dq = dq - (fl - fr) * dy_1(m)
!
!  WENO transport in z-direction.
!
      vsig=max( &
          abs(fq(:,m,n-3,iuz)), abs(fq(:,m,n-2,iuz)), abs(fq(:,m,n-1,iuz)), &
          abs(fq(:,m,n  ,iuz)), &
          abs(fq(:,m,n+1,iuz)), abs(fq(:,m,n+2,iuz)), abs(fq(:,m,n+3,iuz)) )
!
      do i=-nw-1+1,nw-1

        tmp = fq(:,m,n+i,iq)
        if (lref) tmp(nghost+1:mx-nghost)=tmp(nghost+1:mx-nghost)+ref

        if (iq1>0) then
          fr=fq(:,m,n+i,iq1)
          if (lref1) fr(nghost+1:mx-nghost) = fr(nghost+1:mx-nghost) + ref1
          tmp=tmp*fr
        endif

        df(i+1,:)=vsig           *tmp
        f (i+1,:)=fq(:,m,n+i,iuz)*tmp

      enddo
!
!      dl_1=dz_1(n)
!      call weno5_1d(dq,dl_1)
!
!  Return transport pencil without ghost cells.
!
      dq_out(:)=dq(nghost+1:mx-nghost)
!
!  Deallocate arrays.
!
      deallocate(f,df)
!
    endsubroutine weno5
!***********************************************************************
    subroutine weno5_1d(flux)
!
!  Fifth order implementation of WENO scheme (1-D version).
!
!  29-dec-09/evghenii: coded
!
      real, dimension(:), intent(inout) :: flux
!
      real, parameter :: WENO_POW = 2
      real, parameter :: WENO_EPS = 1.0e-6
      real, dimension(size(flux)) :: b1, b2, b3
      real, dimension(size(flux)) :: wh1, wh2, wh3, wh
      real, dimension(size(flux)) :: fh1, fh2, fh3
      real :: g1, g2, g3
!
      f(:,:) = 0.5 * (f(:,:) + df(:,:))
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
    endsubroutine weno5_1d
!***********************************************************************
endmodule WENO_transport
