! $Id$
!
!  This module is the interface for the coala module
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lcoala= .true.
!
!***************************************************************
module Coala
   use Cdata

   implicit none

   real :: tabflux_coag_k0(ndustspec,ndustspec,ndustspec)
   real, dimension(ndustspec+1) :: massgrid
   
   include 'coala.h'
   contains
!***********************************************************************
   subroutine initialize_coala(kernel,rhograin,Q,kpol)

   use Dustvelocity, only: mdminus,mdplus
   use functions_tabflux_tabintflux
   use functions_flux_intflux
   use GQLeg_nodes_weights, only: GQLeg_nodes,GQLeg_weights
   use polynomials_legendre,only: compute_mat_coeffs

   character(len=30) :: kernel
   real          :: rhograin
   integer :: Q,kpol

   real, allocatable :: mat_coeffs_leg(:,:),vecnodes(:),vecweights(:)
   integer :: k
   real :: K0

   !generate Legendre polynomial coefficient
   allocate(mat_coeffs_leg(kpol+1,kpol+1))
   call compute_mat_coeffs(kpol,mat_coeffs_leg)

   !for Gauss quadrature
   allocate(vecnodes(Q),vecweights(Q))
   call GQLeg_nodes(Q,vecnodes)
   call GQLeg_weights(Q,vecweights)

   do k=1,ndustspec
     massgrid(k) = mdminus(k)
   enddo
   massgrid(ndustspec+1) = mdplus(ndustspec)

   K0 = pi*(4./3.*pi*rhograin)**(-2./3.)
   if(kpol == 0) then
      call compute_coagtabflux_GQ_k0(kernel,K0,Q,vecnodes,vecweights,ndustspec,kpol,massgrid,mat_coeffs_leg,tabflux_coag_k0)
   endif

   endsubroutine initialize_coala 
!***********************************************************************
   subroutine coala_advance(rhodust,rhodust_floor,dv,new_rhodust)

    use interface_coala_coag,only: coala_coag_k0

    real :: rhodust_floor
    real, dimension(ndustspec) :: rhodust, new_rhodust
    real, dimension(ndustspec,ndustspec) :: dv

    call coala_coag_k0(ndustspec,massgrid,tabflux_coag_k0,rhodust,rhodust_floor,dv,dt,new_rhodust)

   endsubroutine coala_advance
!***********************************************************************
    subroutine pushpars2c(p_par)

    use Syscalls, only: copy_addr
    use General , only: string_to_enum

    integer, parameter :: n_pars=2
    integer(KIND=ikind8), dimension(n_pars) :: p_par

    call copy_addr(tabflux_coag_k0,p_par(1)) ! (ndustspec) (ndustspec) (ndustspec)
    call copy_addr(massgrid,p_par(2)) ! (ndustspec+1)

    endsubroutine pushpars2c
!***********************************************************************
endmodule Coala
