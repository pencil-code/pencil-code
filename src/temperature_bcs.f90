  module EnergyBcs

    use Cdata
    use EquationOfState, only: cs2bot, cs2top
    use Messages

    private

    include 'energy_bcs.h'
    include 'eos_params.h'

    real :: gamma, gamma_m1, cp, density_scale, density_scale1

    contains
!**************************************************************************************************
    subroutine initialize_energy_bcs

      use EquationOfState, only: get_gamma_etc

      call get_gamma_etc(gamma,cp=cp)
      gamma_m1=gamma-1.

    endsubroutine initialize_energy_bcs
!***********************************************************************
    subroutine bc_lnrho_hds_z_iso_energ(f,topbot)
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
      use EquationOfState, only: eoscalc
      use Gravity, only: gravz
!
      real, dimension (:,:,:,:), intent (inout) :: f
      integer, intent(IN) :: topbot
!
      real :: dlnrhodz, cs2_point
      integer :: i
!
      select case (topbot)
!
!  Bottom boundary
!
      case(BOT)
!
!  Energy equation formulated in logarithmic temperature.
!
        if (ltemperature_nolog) then
          call eoscalc(ilnrho_TT,f(l1,m1,n1,ilnrho),f(l1,m1,n1,iTT),cs2=cs2_point)
        else
          call eoscalc(ilnrho_lnTT,f(l1,m1,n1,ilnrho),f(l1,m1,n1,ilnTT),cs2=cs2_point)
        endif
        dlnrhodz = gamma * gravz/cs2_point    ! dT/dz?  (see comment above)
!
        do i=1,nghost
          f(:,:,n1-i,ilnrho) = f(:,:,n1+i,ilnrho) - dz2_bound(-i)*dlnrhodz
        enddo
!
!  Top boundary
!
      case(TOP)
!
!  Energy equation formulated in logarithmic temperature.
!
        if (ltemperature_nolog) then
          call eoscalc(ilnrho_TT,f(l2,m2,n2,ilnrho),f(l2,m2,n2,iTT),cs2=cs2_point)
        else
          call eoscalc(ilnrho_lnTT,f(l2,m2,n2,ilnrho),f(l2,m2,n2,ilnTT),cs2=cs2_point)
        endif
        dlnrhodz =  gamma *gravz/cs2_point
!
        do i=1,nghost
          f(:,:,n2+i,ilnrho) = f(:,:,n2-i,ilnrho) + dz2_bound(i)*dlnrhodz
        enddo
!
      case default
!
      endselect
!
    endsubroutine bc_lnrho_hds_z_iso_energ
!***********************************************************************
    subroutine bc_ss_temp_x(f,topbot)
!
!  boundary condition for entropy: constant temperature
!
!   3-aug-2002/wolf: coded
!  26-aug-2003/tony: distributed across ionization modules
!
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
      integer :: i
!
!  Constant temperature/sound speed for entropy, i.e. antisymmetric
!  ln(cs2) relative to cs2top/cs2bot.
!  This assumes that the density is already set (ie density _must_ register
!  first!)
!
!  check whether we want to do top or bottom (this is processor dependent)
!
      select case (topbot)
!
!  bottom boundary
!
      case(BOT)
        if (ldebug) print*, 'bc_ss_temp_x: set x bottom temperature: cs2bot=',cs2bot
        if (cs2bot<=0.) call fatal_error('bc_ss_temp_x','cannot have cs2bot<=0')
!
         f(l1,:,:,ilnTT) = log(cs2bot/(gamma_m1*cp))    !! factor 1/cp corrected
         do i=1,nghost; f(l1-i,:,:,ilnTT)=2*f(l1,:,:,ilnTT)-f(l1+i,:,:,ilnTT); enddo
!
!  top boundary
!
      case(TOP)
        if (ldebug) print*, 'bc_ss_temp_x: set x top temperature: cs2top=',cs2top
        if (cs2top<=0.) call fatal_error('bc_ss_temp_x','cannot have cs2top<=0')
!
        f(l2,:,:,ilnTT) = log(cs2top/(gamma_m1*cp))
        do i=1,nghost; f(l2+i,:,:,ilnTT)=2*f(l2,:,:,ilnTT)-f(l2-i,:,:,ilnTT); enddo
!
      case default
        call fatal_error('bc_ss_temp_x','invalid argument')
      endselect
!
    endsubroutine bc_ss_temp_x
!***********************************************************************
    subroutine bc_ss_temp_z(f,topbot,lone_sided)
!
!  boundary condition for entropy: constant temperature
!
!   3-aug-2002/wolf: coded
!  26-aug-2003/tony: distributed across ionization modules
!  11-oct-2016/MR: changes for use of one-sided BC formulation (chosen by setting new optional switch lone_sided)
!
      use General, only: loptest
      use Deriv, only: set_ghosts_for_onesided_ders
!
      real, dimension (:,:,:,:) :: f
      integer, intent(IN) :: topbot
      logical, optional :: lone_sided
!
      integer :: i
!
!  Constant temperature/sound speed for entropy, i.e. antisymmetric
!  ln(cs2) relative to cs2top/cs2bot.
!  This assumes that the density is already set (ie density _must_ register
!  first!)
!
!  check whether we want to do top or bottom (this is processor dependent)
!
      select case (topbot)
!
!  bottom boundary
!
      case(BOT)
        if (ldebug) print*, 'bc_ss_temp_z: set z bottom temperature: cs2bot=',cs2bot
        if (cs2bot<=0.) call fatal_error('bc_ss_temp_z','cannot have cs2bot<=0')

        if (ltemperature_nolog) then
          f(:,:,n1,iTT)   = cs2bot/(gamma_m1*cp)      ! factor 1/cp corrected
        else
          f(:,:,n1,ilnTT) = log(cs2bot/(gamma_m1*cp))
        endif
        if (loptest(lone_sided)) then
          call set_ghosts_for_onesided_ders(f,topbot,ilnTT,3,.true.)
        else
          do i=1,nghost; f(:,:,n1-i,ilnTT)=2*f(:,:,n1,ilnTT)-f(:,:,n1+i,ilnTT); enddo
        endif
!
!  top boundary
!
      case(TOP)
        if (ldebug) print*,'bc_ss_temp_z: set z top temperature: cs2top=',cs2top
        if (cs2top<=0.) call fatal_error('bc_ss_temp_z','cannot have cs2top<=0')

        if (ltemperature_nolog) then
          f(:,:,n2,iTT)   = cs2top/(gamma_m1*cp)
        else
          f(:,:,n2,ilnTT) = log(cs2top/(gamma_m1*cp))
        endif
        if (loptest(lone_sided)) then
          call set_ghosts_for_onesided_ders(f,topbot,ilnTT,3,.true.)
        else
          do i=1,nghost; f(:,:,n2+i,ilnTT)=2*f(:,:,n2,ilnTT)-f(:,:,n2-i,ilnTT); enddo
        endif
!
      case default
        call fatal_error('bc_ss_temp_z','invalid argument')
      endselect
!
    endsubroutine bc_ss_temp_z
!***********************************************************************
    subroutine bc_ism_energ(f,topbot,j)
!
!  30-nov-15/fred: Replaced bc_ctz and bc_cdz.
!  Apply observed scale height locally from Reynolds 1991, Manchester & Taylor
!  1981 for warm ionized gas - dominant scale height above 500 parsecs.
!  Apply constant local temperature across boundary for entropy.
!  Motivation to prevent numerical spikes in shock fronts, which cannot be
!  absorbed in only three ghost cells, but boundary thermodynamics still
!  responsive to interior dynamics.
!
      real, dimension (:,:,:,:) :: f
      integer, intent(IN) :: topbot
      integer :: j,k

      if (j/=ilnTT) call fatal_error('bc_ism_energ','only for log temperature')

      select case (topbot)
!
      case(BOT)               ! bottom boundary

        do k=1,nghost
          f(:,:,n1-k,ilnTT)=f(:,:,n1,ilnTT)+log((z(n1)-z(n1-k))*density_factor+1.)
        enddo
!
      case(TOP)               ! top boundary

        do k=1,nghost
          f(:,:,n1-k,ilnTT)=f(:,:,n1,ilnTT)+log((z(n1)-z(n1-k))*density_factor+1.)
        enddo
!
      case default
        print*, "bc_ism_energ topbot should be BOT or TOP"
      endselect
!
    endsubroutine bc_ism_energ
!***********************************************************************
!********************************************************************
!********************************************************************
!************        DO NOT DELETE THE FOLLOWING        *************
!********************************************************************
!**  This is an automatically generated include file that creates  **
!**  copies dummy routines from nospecial.f90 for any Special      **
!**  routines not implemented in this file                         **
!**                                                                **
    include 'energy_bcs_dummies.inc'
!***********************************************************************
  endmodule EnergyBcs
