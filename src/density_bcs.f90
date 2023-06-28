  module DensityBcs

    use Cdata
    use EquationOfState, only: cs20, cs2bot, cs2top
    use Messages

    include 'density_bcs.h'

    private

    integer, parameter :: XBOT=1, XTOP=nx

    real, dimension(:,:), pointer :: reference_state

    contains
!**************************************************************************************************
    subroutine initialize_density_bcs

      use SharedVariables, only: get_shared_variable
!
!  Get the shared variables
!
      if (lreference_state) call get_shared_variable('reference_state',reference_state, &
                                                     caller='initialize_density_bcs')
    endsubroutine initialize_density_bcs
!**************************************************************************************************
    subroutine bc_lnrho_cfb_r_iso(f,topbot)
!
!  Boundary condition for radial centrifugal balance
!
!  This sets
!    \partial_{r} \ln\rho
!  such that
!    (\partial_{r} p)/\rho = cs^2 \partial_{r} \ln\rho} = uphi**2/rad - \partial_{r} Phi
!  where Phi is the gravitational potential
!
!  i.e. it enforces centrifugal balance at the boundary.
!
!  As it is, works only for isobaric, isothermal and cylindrical coordinates
!
!  21-aug-2006/wlad: coded
!
      use DensityMethods, only: getrho
      use Gravity, only: potential
      use Sub, only: div
!
      real, dimension (:,:,:,:), intent (inout) :: f
      integer, intent(IN) :: topbot
      real, dimension (size(f,2),size(f,3)) :: centterm,uphi,rho
      real :: potp,potm,rad,step,gravterm
      integer :: i
!
      select case (topbot)
!
!  Bottom boundary
!
      case(BOT)
!
        if (ldensity_nolog) call getrho(f(l1,:,:,irho),XBOT,rho)
!
        do i=1,nghost
!
          call potential(R=x(l1-i),pot=potm)
          call potential(R=x(l1+i),pot=potp)
!
          gravterm= -(potm-potp)/cs20
!
          step=-dx2_bound(-i)
          rad=x(l1-i)
!
          centterm = f(l1-i,:,:,iuy)**2 * step/(rad*cs20)  !???
          if (ldensity_nolog) then
            f(l1-i,:,:,irho) = f(l1+i,:,:,irho) + rho*(gravterm + centterm)
            if (lreference_state) &
              f(l1-i,:,:,irho)=  f(l1-i,:,:,irho) + dx2_bound(-i)*reference_state(XBOT,iref_grho)
          else
            f(l1-i,:,:,ilnrho)=f(l1+i,:,:,ilnrho) + gravterm + centterm
          endif
!
          !print*,'potentials',potm,potp,-(potm-potp)
          !print*,'centrifugal',f(l1-i,mpoint,npoint,iuy)**2 *step/rad
          !stop
!
        enddo
!
!  Top boundary
!
      case(TOP)
!
        if (ldensity_nolog) call getrho(f(l2,:,:,irho),XTOP,rho)
!
        do i=1,nghost
!
          call potential(R=x(l2+i),pot=potp)
          call potential(R=x(l2-i),pot=potm)
!
          gravterm = -(potp-potm)/cs20
!
          step=dx2_bound(i)
          rad=x(l2+i)
!
          centterm = f(l2+i,:,:,iuy)**2 * step/(rad*cs20)
          if (ldensity_nolog) then
            f(l2+i,:,:,irho) = f(l2-i,:,:,irho) + rho*(gravterm + centterm)
            if (lreference_state) &
              f(l2+i,:,:,irho)=  f(l2+i,:,:,irho) - dx2_bound(i)*reference_state(XTOP,iref_grho)
          else
            f(l2+i,:,:,ilnrho) = f(l2-i,:,:,ilnrho) + gravterm + centterm
          endif
!
          !if (i==nghost) then
          !  print*,'potentials',potp,potm,-potp+potm,-(potp-potm)
          !  print*,'centrifugal',f(l2+i,mpoint,npoint,iuy)**2 *step/rad
          !  stop
          !endif
        enddo
!
      case default
!
      endselect
!
    endsubroutine bc_lnrho_cfb_r_iso
!***********************************************************************
    subroutine bc_lnrho_hdss_z_iso(f,topbot)
!
!  Smooth out density perturbations with respect to hydrostatic
!  stratification in Fourier space.
!
!  Note: Since boundconds_x and boundconds_y are called first,
!  this doesn't set the corners properly. However, this is
!  not a problem since cross derivatives of density are never
!  needed.
!
!  05-jul-07/tobi: Adapted from bc_aa_pot3
!
      use Fourier, only: fourier_transform_xy_xy, fourier_transform_other, kx_fft, ky_fft
      use Gravity, only: potential
!
      real, dimension (:,:,:,:), intent (inout) :: f
      integer, intent(IN) :: topbot
!
      real, dimension (nx,ny) :: kx,ky,kappa,exp_fact
      real, dimension (nx,ny) :: tmp_re,tmp_im
      real :: pot
      integer :: i
!
!  Get local wave numbers
!
      kx = spread(kx_fft(ipx*nx+1:ipx*nx+nx),2,ny)
      ky = spread(ky_fft(ipy*ny+1:ipy*ny+ny),1,nx)
!
!  Calculate 1/k^2, zero mean
!
      if (lshear) then
        kappa = sqrt((kx+ky*deltay/Lx)**2+ky**2)
      else
        kappa = sqrt(kx**2 + ky**2)
      endif
!
!  Check whether we want to do top or bottom (this is processor dependent)
!
      select case (topbot)
!
!  Potential field condition at the bottom
!
      case(BOT)
!
        do i=1,nghost
!
! Calculate delta_z based on z(), not on dz to improve behavior for
! non-equidistant grid (still not really correct, but could be OK)
!
          exp_fact = exp(-kappa*(z(n1+i)-z(n1-i)))
!
!  Determine potential field in ghost zones
!
          !  Fourier transforms of x- and y-components on the boundary
          call potential(z=z(n1+i),pot=pot)
          if (ldensity_nolog) then
            tmp_re = f(l1:l2,m1:m2,n1+i,irho)*exp(+pot/cs2bot)
          else
            tmp_re = f(l1:l2,m1:m2,n1+i,ilnrho)   +pot/cs2bot
          endif
          tmp_im = 0.0
          if (nxgrid>1 .and. nygrid>1) then
            call fourier_transform_xy_xy(tmp_re,tmp_im)
          else
            call fourier_transform_other(tmp_re,tmp_im)
          endif
          tmp_re = tmp_re*exp_fact
          tmp_im = tmp_im*exp_fact
          ! Transform back
          if (nxgrid>1 .and. nygrid>1) then
            call fourier_transform_xy_xy(tmp_re,tmp_im,linv=.true.)
          else
            call fourier_transform_other(tmp_re,tmp_im,linv=.true.)
          endif
          call potential(z=z(n1-i),pot=pot)
          if (ldensity_nolog) then
            f(l1:l2,m1:m2,n1-i,irho)   = tmp_re*exp(-pot/cs2bot)
          else
            f(l1:l2,m1:m2,n1-i,ilnrho) = tmp_re     -pot/cs2bot
          endif
!
        enddo
!
!  Potential field condition at the top
!
      case(TOP)
!
        do i=1,nghost
!
! Calculate delta_z based on z(), not on dz to improve behavior for
! non-equidistant grid (still not really correct, but could be OK)
!
          exp_fact = exp(-kappa*(z(n2+i)-z(n2-i)))
!
!  Determine potential field in ghost zones
!
          !  Fourier transforms of x- and y-components on the boundary
          call potential(z=z(n2-i),pot=pot)
          if (ldensity_nolog) then
            tmp_re = f(l1:l2,m1:m2,n2-i,irho)*exp(+pot/cs2top)
          else
            tmp_re = f(l1:l2,m1:m2,n2-i,ilnrho)   +pot/cs2top
          endif
          tmp_im = 0.0
          if (nxgrid>1 .and. nygrid>1) then
            call fourier_transform_xy_xy(tmp_re,tmp_im)
          else
            call fourier_transform_other(tmp_re,tmp_im)
          endif
          tmp_re = tmp_re*exp_fact
          tmp_im = tmp_im*exp_fact
          ! Transform back
          if (nxgrid>1 .and. nygrid>1) then
            call fourier_transform_xy_xy(tmp_re,tmp_im,linv=.true.)
          else
            call fourier_transform_other(tmp_re,tmp_im,linv=.true.)
          endif
          call potential(z=z(n2+i),pot=pot)
          if (ldensity_nolog) then
            f(l1:l2,m1:m2,n2+i,irho)   = tmp_re*exp(-pot/cs2top)
          else
            f(l1:l2,m1:m2,n2+i,ilnrho) = tmp_re     -pot/cs2top
          endif
!
        enddo
!
      case default
!
        if (lroot) call fatal_error("bc_lnrho_hdss_z_iso","invalid argument")
!
      endselect
!
    endsubroutine bc_lnrho_hdss_z_iso
!***********************************************************************
    subroutine bc_lnrho_hds_z_iso_dens(f,topbot)
!
!  Boundary condition for density *and* entropy.
!
!  This sets
!    \partial_{z} \ln\rho
!  such that
!    \partial_{z} p = \rho g_{z},
!  i.e. it enforces hydrostatic equlibrium at the boundary.
!
!  Currently this is only correct if
!    \partial_{z} lnT = 0
!  at the boundary.
!
!  12-Juil-2006/dintrans: coded
!
      use General, only: itoa
      use Gravity, only: potential, gravz
      use Sub, only: div
!
      real, dimension (:,:,:,:), intent (inout) :: f
      integer, intent(IN) :: topbot
!
      real, dimension (size(f,1),size(f,2)) :: cs2
      real, dimension (nx) :: divu
      real :: potp,potm
      integer :: i
!
      if (ltemperature) return    !???

      if (bcz12(ilnTT,topbot)/='s') call fatal_error("bc_lnrho_hds_z_iso", &
         "This boundary condition for density is "// &
         "currently only correct for bcz"//trim(itoa(topbot))//"(i[ln]TT)='s'")

      select case (topbot)
!
!  Bottom boundary
!
      case(BOT)
!
!  Isothermal or polytropic equations of state.
!
        do i=1,nghost

          call potential(z=z(n1-i),pot=potm)
          call potential(z=z(n1+i),pot=potp)
!
          if (.true.) then
            cs2 = cs2bot
          else
            ! Note: Since boundconds_x and boundconds_y are called first,
            ! this doesn't set the corners properly. However, this is
            ! not a problem since cross derivatives of density are never
            ! needed.
            n = n1+i
            do m = m1,m2
              call div(f,iuu,divu)
              cs2(l1:l2,m) = cs2bot - f(l1:l2,m,n,ishock)*divu
            enddo
          endif
!
          if (ldensity_nolog) then
            f(:,:,n1-i,irho)   = f(:,:,n1+i,irho)*exp(-(potm-potp)/cs2)
          else
            f(:,:,n1-i,ilnrho) = f(:,:,n1+i,ilnrho)   -(potm-potp)/cs2
          endif
!
        enddo
!
!  Top boundary
!
      case(TOP)
!
!  Isothermal or polytropic equations of state.
!
        do i=1,nghost

          call potential(z=z(n2+i),pot=potp)
          call potential(z=z(n2-i),pot=potm)
          if (.true.) then
            cs2 = cs2top
          else
            ! Note: Since boundconds_x and boundconds_y are called first,
            ! this doesn't set the corners properly. However, this is
            ! not a problem since cross derivatives of density are never
            ! needed.
            n = n2-i
            do m = m1,m2
              call div(f,iuu,divu)
              cs2(l1:l2,m) = cs2top - f(l1:l2,m,n,ishock)*divu
            enddo
          endif
          if (ldensity_nolog) then
            f(:,:,n2+i,irho)   = f(:,:,n2-i,irho)*exp(-(potp-potm)/cs2)
          else
            f(:,:,n2+i,ilnrho) = f(:,:,n2-i,ilnrho)   -(potp-potm)/cs2
          endif
        enddo
!
      case default
!
      endselect
!
    endsubroutine bc_lnrho_hds_z_iso_dens
!***********************************************************************
    subroutine bc_ism_dens(f,topbot,j)
!
!  30-nov-15/fred: Replaced bc_ctz and bc_cdz.
!  Apply observed scale height locally from Reynolds 1991, Manchester & Taylor
!  1981 for warm ionized gas - dominant scale height above 500 parsecs.
!  Apply constant local temperature across boundary for entropy.
!  Motivation to prevent numerical spikes in shock fronts, which cannot be
!  absorbed in only three ghost cells, but boundary thermodynamics still
!  responsive to interior dynamics.
!
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
      integer :: j,k
      real :: density_scale1
!
      if (j/=irho .and. j/=ilnrho) call fatal_error('bc_ism_dens','only for irho, ilnrho')
!
      density_scale1=1./density_scale_factor
!
      select case (topbot)
!
      case(BOT)               ! bottom boundary
        do k=1,nghost
          if (ldensity_nolog) then
            f(:,:,k,j)=f(:,:,n1,j)*exp(-(z(n1)-z(k))*density_scale1)
          else
            f(:,:,k,j)=f(:,:,n1,j)     -(z(n1)-z(k))*density_scale1
          endif
        enddo
!
      case(TOP)               ! top boundary

        do k=1,nghost
          if (ldensity_nolog) then
            f(:,:,n2+k,j)=f(:,:,n2,j)*exp(-(z(n2+k)-z(n2))*density_scale1)
          else
            f(:,:,n2+k,j)=f(:,:,n2,j)     -(z(n2+k)-z(n2))*density_scale1
          endif
        enddo
!
      case default
        print*, "bc_ism_dens topbot should be BOT or TOP"
      endselect
!
    endsubroutine bc_ism_dens
!***********************************************************************
  endmodule DensityBcs

