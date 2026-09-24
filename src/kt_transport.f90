! $Id$
!
!  Kurganov-Tadmor (KT) flux-limited transport for the CONSERVATIVE
!  RELATIVISTIC hydro equations
!
!      d_t T^00 + d_i S_i    = 0
!      d_t S_j  + d_i T^ij   = 0,   T^ij = S_i S_j/(T^00+p) + p delta_ij,
!
!  with the bag-model pressure p = w/4 - eps,  w = (4/3)(T^00-eps)(2 sqrt(1-lam)-1),
!  lam = (3/4)|S|^2/(T^00-eps)^2.
!
!  eps(t,x) is an externally prescribed vacuum energy.  With lhiggsless=T it is
!  the bag constant of the Higgsless scheme (Jinno et al.), reconstructed here
!  from the bubble-wall crossing-time auxiliary (hless) exactly as in hydro.f90;
!  see doc/implementation/higgsless.  Higgsless is one application of the
!  scheme, not a requirement: with lhiggsless=F, eps=0 throughout and the
!  closure reduces to the ordinary radiation equation of state, so KT is
!  available to any conservative relativistic run.
!
!  Scheme: MUSCL reconstruction with the generalized minmod limiter
!  (theta-limiter) + Rusanov (local Lax-Friedrichs) flux whose wave speed is
!  evaluated per face from the two reconstructed states (state_signal_speed):
!  a = max over the two sides of (|v|+cs)/(1+|v|cs), the relativistic
!  composition of the flow speed with the sound speed, so cs <= a <= 1,
!  plus projection of the reconstructed face states onto the admissible set
!  (positivity of T^00-eps and causality |S| bound).
!
!  This is a line-by-line port of two cross-validated references:
!    - jax_lattice/dynamics_files/relativistic_perfect_fluid.py
!      (_reconstruct_muscl_kt, _project_admissible, _flux, _kt_rhs), and
!    - the C++ Higgsless code (Higgsless-Simulation-V2, KTsolver::update_H,
!      del_x, minmod, pN_staggered), including its face treatment of eps
!      (al_ph = 0.5*(al + al_p): one interpolated value per face, seen by
!      both sides).
!
!  It is selected at compile time with KT_TRANSPORT = kt_transport in
!  Makefile.local (default: nokt_transport) and activated at run time with
!  lkt_transport=T in &hydro_run_pars; both defaults leave every existing
!  setup untouched.  Needs 2 ghost zones (Pencil has 3).
!
!  02-sep-2026/Isak Stomberg: created for the Higgsless application
!  16-sep-2026/Isak Stomberg: decoupled from lhiggsless; eps=0 without it
!
module KT_transport
!
  use Cdata
!
  implicit none
!
  include "kt_transport.h"
!
!  Configuration, set once by hydro's initialize via kt_init.
!
  integer :: iux_kt=0, ihless_kt=0
  real :: eps_kt=0.0, width_abs_kt=0.0, safety_kt=1e-8, theta_kt=2.0
  real, parameter :: cs_kt=sqrt(1.0/3.0)  ! sound speed of the radiation/bag EOS
  real, parameter :: tiny_kt=1e-30
!
  contains
!***********************************************************************
    subroutine kt_init(irho_in,iux_in,ihless_in,eps_in,width_abs_in,theta_in)
!
!  Store farray indices and the limiter parameter, plus the Higgsless vacuum
!  energy if that source is active.  Called from initialize_hydro when
!  lkt_transport=T.  Without Higgsless, pass ihless_in=0 and eps_in=0: the
!  scheme then runs on the plain conservative relativistic system.
!
!  02-sep-2026/Isak Stomberg: coded
!
      use Quiet, only: keep_compiler_quiet
!
      integer, intent(in) :: irho_in, iux_in, ihless_in
      real, intent(in) :: eps_in, width_abs_in, theta_in
!
!  iux and ihless are already valid when initialize_hydro calls this; irho
!  (Cdata) is not yet populated at that point, so it is read from Cdata at
!  call time in kt_div/kt_div_tensor instead of being stored here.
!
      call keep_compiler_quiet(irho_in)
      iux_kt=iux_in; ihless_kt=ihless_in
      eps_kt=eps_in; width_abs_kt=width_abs_in; theta_kt=theta_in
      ! safety_kt keeps its 1e-8 default (internal face-state admissible-projection floor)
      if (lroot) print*, 'kt_init: KT flux-limited transport enabled; theta=', theta_kt
!
    endsubroutine kt_init
!***********************************************************************
    subroutine kt_div(f,divergence)
!
!  KT counterpart of Sub's div for the conservative relativistic energy
!  equation: the flux divergence of the energy flux,
!    divergence = d_i S_i,
!  summed over all three directions; the caller subtracts it from df.
!  kt_div_tensor is the counterpart of div_tensor for the three momentum
!  components, and the two share the reconstruction described there.
!
!  Like div and kt_div_tensor, this takes no m,n and no time: the mn-loop
!  indices m,n and the clock t are the Cdata globals.
!
!  02-sep-2026/Isak Stomberg: coded
!  16-sep-2026/Isak Stomberg: renamed from kt_transp; m, n, mu and tcur dropped
!
      real, contiguous, dimension(:,:,:,:), intent(in) :: f
      real, dimension(nx), intent(out) :: divergence
!
      integer, parameter :: mu_energy=1
      real, dimension(nx,4,-2:2) :: ucons
      real, dimension(nx,-1:1) :: epsc
      real, dimension(nx,4) :: slope_m, slope_0, slope_p, ul, ur
      real, dimension(nx) :: hph, hmh, aph, amh, pl, pr, dl1, kt_cfl_rate
      integer :: dir
!
      divergence=0.0
      kt_cfl_rate=0.0
!
! Each direction x,y,z below contributes one term: (H_{+1/2} - H_{-1/2}) / dx_dir to divergence.
!
      do dir=1,3
        call gather_stencil(f,dir,ucons,epsc,slope_m,slope_0,slope_p)
!
!  Face +1/2, then face -1/2: MUSCL states either side, one eps per face.
!
        call muscl_faces(ucons,slope_0,slope_p, 0,ul,ur)
        call face_state(ul,ur,0.5*(epsc(:,0)+epsc(:,1)),pl,pr,aph)
        hph=flux_component(ul,ur,pl,pr,aph,dir,mu_energy)
!
        call muscl_faces(ucons,slope_m,slope_0,-1,ul,ur)
        call face_state(ul,ur,0.5*(epsc(:,-1)+epsc(:,0)),pl,pr,amh)
        hmh=flux_component(ul,ur,pl,pr,amh,dir,mu_energy)
!
!  Divergence contribution: (H_{+1/2} - H_{-1/2}) / dx_dir.
!
        select case (dir)
          case (1); dl1=dx_1(l1:l2)
          case (2); dl1=dy_1(m)
          case (3); dl1=dz_1(n)
        endselect
        divergence=divergence+(hph-hmh)*dl1
        kt_cfl_rate=kt_cfl_rate+max(aph,amh)*dl1
      enddo
!
!  Register the multidimensional Rusanov CFL rate.  Using the same local face
!  speeds as the numerical flux keeps the dissipation and timestep estimates
!  consistent.
!
      if (lupdate_courant_dt.and..not.ldt_paronly) &
        dt1_max=max(dt1_max,kt_cfl_rate/cdt)
!
    endsubroutine kt_div
!***********************************************************************
    subroutine kt_div_tensor(f,divergence)
!
!  KT counterpart of Sub's div_tensor: returns the momentum flux divergence
!  d_i T^{ij} for all three components at once, so that the caller can swap
!
!    call kt_div_tensor(f,divTij)
!    call div_tensor(f,divTij,iTij,lyz_first=.true.)
!
!  without a component loop.  Note the difference in input: div_tensor
!  differentiates the STORED tensor at iTij, whereas this reconstructs the
!  fluxes from the conserved state, so the two are different discretizations
!  of the same continuum quantity, not the same operator.
!
!  Like div_tensor, this takes no m,n and no time: the mn-loop indices m,n and
!  the clock t are the Cdata globals (der, and hence div_tensor, use them the
!  same implicit way), so passing them would only shadow them with equal values.
!
!  This is the single-pass form: the stencil gather (gather_stencil), the face
!  states (muscl_faces) and the projection, pressure and signal speed
!  (face_state) are all independent of the component mu, so they are evaluated
!  ONCE per face and shared by the three momentum fluxes.  Only
!  flux_component depends on mu, and kt_div calls exactly the same routines
!  with mu=1 -- so the two are the same reconstruction by construction rather
!  than by keeping two copies in step.
!
!  The energy component has no div_tensor analogue and is obtained instead from
!  kt_div, which density.f90 calls for the continuity equation.
!
!  09-sep-2026/Isak Stomberg: coded
!  14-sep-2026/Isak Stomberg: single-pass (shared reconstruction over mu)
!
      real, contiguous, dimension(:,:,:,:), intent(in) :: f
      real, dimension(nx,3), intent(out) :: divergence
!
      real, dimension(nx,4,-2:2) :: ucons
      real, dimension(nx,-1:1) :: epsc
      real, dimension(nx,4) :: slope_m, slope_0, slope_p, ul, ur
      real, dimension(nx,3) :: hph, hmh
      real, dimension(nx) :: pl, pr, a, dl1
      integer :: dir, j, mu
!
      divergence=0.0
      do dir=1,3
        call gather_stencil(f,dir,ucons,epsc,slope_m,slope_0,slope_p)
!
!  Face +1/2 and face -1/2, each projected and thermodynamically evaluated once.
!
        call muscl_faces(ucons,slope_0,slope_p, 0,ul,ur)
        call face_state(ul,ur,0.5*(epsc(:,0)+epsc(:,1)),pl,pr,a)
        do mu=2,4
          hph(:,mu-1)=flux_component(ul,ur,pl,pr,a,dir,mu)
        enddo
!
        call muscl_faces(ucons,slope_m,slope_0,-1,ul,ur)
        call face_state(ul,ur,0.5*(epsc(:,-1)+epsc(:,0)),pl,pr,a)
        do mu=2,4
          hmh(:,mu-1)=flux_component(ul,ur,pl,pr,a,dir,mu)
        enddo
!
!  Divergence contribution: (H_{+1/2} - H_{-1/2}) / dx_dir.
!
        select case (dir)
          case (1); dl1=dx_1(l1:l2)
          case (2); dl1=dy_1(m)
          case (3); dl1=dz_1(n)
        endselect
        do j=1,3
          divergence(:,j)=divergence(:,j)+(hph(:,j)-hmh(:,j))*dl1
        enddo
      enddo
!
    endsubroutine kt_div_tensor
!***********************************************************************
    subroutine gather_stencil(f,dir,ucons,epsc,slope_m,slope_0,slope_p)
!
!  Gather the stencil along direction dir and limit its slopes.  Shared by
!  kt_div and kt_div_tensor: everything up to the face states is independent
!  of which flux is being formed, so the two cannot drift apart.
!
!  getcell returns a whole x-row (nx values) at a time, so ucons holds nx
!  INDEPENDENT one-dimensional stencils side by side: a 5 x nx band for
!  dir=2,3, and for dir=1 the same row slid by o.  Each x-index is
!  reconstructed on its own -- the neighbour pattern for one cell is always
!  the five points ucons(i,:,-2:2), never anything two-dimensional.
!
!  The conserved variables need 5 points (-2..2) because the limited slopes
!  are wanted at offsets -1,0,+1 and each of those reaches one cell further
!  out.  eps is only averaged to the faces, never reconstructed, so 3 suffice.
!
!  Component order is (K^0, K^x, K^y, K^z), energy first, so that physflux can
!  select the flux-direction momentum as u4(:,1+dir).
!
!  24-sep-2026/Isak Stomberg: factored out of kt_div and kt_div_tensor
!
      real, contiguous, dimension(:,:,:,:), intent(in) :: f
      integer, intent(in) :: dir
      real, dimension(nx,4,-2:2), intent(out) :: ucons
      real, dimension(nx,-1:1), intent(out) :: epsc
      real, dimension(nx,4), intent(out) :: slope_m, slope_0, slope_p
!
      integer :: comp, o
!
      do o=-2,2
        ucons(:,1,o)=getcell(f,o,dir,m,n,irho)
        do comp=2,4
          ucons(:,comp,o)=getcell(f,o,dir,m,n,iux_kt+comp-2)
        enddo
      enddo
      do o=-1,1
        epsc(:,o)=epscell(f,o,dir,m,n,real(t))
      enddo
!
!  Limited slope-differences (slope*dx) at cell offsets -1, 0, +1:
!  minmod( theta*backward, centered, theta*forward )  [jax _limited_slope,
!  C++ del_x with theta=2].
!
      do comp=1,4
        slope_m(:,comp)=minmod3(theta_kt*(ucons(:,comp,-1)-ucons(:,comp,-2)), &
                                     0.5*(ucons(:,comp, 0)-ucons(:,comp,-2)), &
                                theta_kt*(ucons(:,comp, 0)-ucons(:,comp,-1)))
        slope_0(:,comp)=minmod3(theta_kt*(ucons(:,comp, 0)-ucons(:,comp,-1)), &
                                     0.5*(ucons(:,comp, 1)-ucons(:,comp,-1)), &
                                theta_kt*(ucons(:,comp, 1)-ucons(:,comp, 0)))
        slope_p(:,comp)=minmod3(theta_kt*(ucons(:,comp, 1)-ucons(:,comp, 0)), &
                                     0.5*(ucons(:,comp, 2)-ucons(:,comp, 0)), &
                                theta_kt*(ucons(:,comp, 2)-ucons(:,comp, 1)))
      enddo
!
    endsubroutine gather_stencil
!***********************************************************************
    subroutine muscl_faces(ucons,sl,sr,o,ul,ur)
!
!  MUSCL states either side of the face at o+1/2:
!    U^L = U_o + sigma_o/2,   U^R = U_{o+1} - sigma_{o+1}/2,
!  with no factor dx, because minmod is homogeneous of degree one and the
!  slopes are already undivided.  Called as (slope_0,slope_p, 0) for the
!  +1/2 face and (slope_m,slope_0,-1) for the -1/2 face.
!
!  24-sep-2026/Isak Stomberg: factored out of kt_div and kt_div_tensor
!
      real, dimension(nx,4,-2:2), intent(in) :: ucons
      real, dimension(nx,4), intent(in) :: sl, sr
      integer, intent(in) :: o
      real, dimension(nx,4), intent(out) :: ul, ur
!
      integer :: comp
!
      do comp=1,4
        ul(:,comp)=ucons(:,comp,o)  +0.5*sl(:,comp)
        ur(:,comp)=ucons(:,comp,o+1)-0.5*sr(:,comp)
      enddo
!
    endsubroutine muscl_faces
!***********************************************************************
    subroutine bag_pressure(u4,epsf,p)
!
!  Bag-model pressure of a (projected) state; identical to the thermodynamic
!  part of physflux, factored out so it can be shared across components.
!
!  14-sep-2026/Isak Stomberg: coded
!
      real, dimension(nx,4), intent(in) :: u4
      real, dimension(nx), intent(in) :: epsf
      real, dimension(nx), intent(out) :: p
!
      real, dimension(nx) :: k0e, ki2, lam, w
!
      k0e=u4(:,1)-epsf
      ki2=u4(:,2)**2+u4(:,3)**2+u4(:,4)**2
      lam=min(0.75*ki2/max(k0e**2,tiny_kt),0.75)
      w=(4.0/3.0)*k0e*(2.0*sqrt(max(1.0-lam,0.0))-1.0)
      p=0.25*w-epsf
!
    endsubroutine bag_pressure
!***********************************************************************
    function state_signal_speed(u4,p) result(a)
!
!  Conservative upper bound on the coordinate speed of a sound signal from
!  one state.  The fastest signal is obtained by relativistically composing
!  the fluid speed magnitude with the radiation-fluid sound speed.  Using
!  |v| rather than its face-normal component makes this an upper bound in all
!  coordinate directions without requiring the full multidimensional
!  eigensystem.
!
      real, dimension(nx,4), intent(in) :: u4
      real, dimension(nx), intent(in) :: p
      real, dimension(nx) :: a, vmag
!
      vmag=sqrt(u4(:,2)**2+u4(:,3)**2+u4(:,4)**2)/max(u4(:,1)+p,tiny_kt)
      vmag=min(max(vmag,0.0),1.0)
      a=(vmag+cs_kt)/(1.0+vmag*cs_kt)
!
    endfunction state_signal_speed
!***********************************************************************
    function getcell(f,o,dir,m,n,ind) result(q)
!
!  One x-row (nx values) of the conserved variable ind, at cell offset o along
!  direction dir from the current row (m,n).  For dir=1 the row is slid along
!  itself; for dir=2,3 it is a parallel row displaced in y or z.
!
      real, contiguous, dimension(:,:,:,:), intent(in) :: f
      integer, intent(in) :: o, dir, m, n, ind
      real, dimension(nx) :: q
!
      select case (dir)
        case (1); q=f(l1+o:l2+o,m,n,ind)
        case (2); q=f(l1:l2,m+o,n,ind)
        case (3); q=f(l1:l2,m,n+o,ind)
      endselect
!
    endfunction getcell
!***********************************************************************
    function epscell(f,o,dir,m,n,tcur) result(eps)
!
!  One x-row (nx values) of the vacuum energy eps(t,x) at cell offset o along
!  direction dir, indexed as in getcell: the same crossing-time ramp as
!  hydro.f90's hydro_after_boundary_conservative (linear ramp of temporal
!  width width_abs_kt centred on the crossing time stored in the hless aux).
!
      real, contiguous, dimension(:,:,:,:), intent(in) :: f
      integer, intent(in) :: o, dir, m, n
      real, intent(in) :: tcur
      real, dimension(nx) :: eps, tau
!
!  Without the Higgsless source there is no vacuum energy and no hless slot to
!  read; the bag closure then reduces to the ordinary radiation equation of
!  state.  ihless_kt=0 is the sentinel kt_init sets for that case.
!
      if (ihless_kt==0) then
        eps=0.0
        return
      endif
!
      tau=getcell(f,o,dir,m,n,ihless_kt)
      if (width_abs_kt==0.0) then
        where (tcur<tau)
          eps=eps_kt
        elsewhere
          eps=0.0
        endwhere
      else
        eps=eps_kt*max(0.0,min(1.0,(tau+0.5*width_abs_kt-tcur)/width_abs_kt))
      endif
!
    endfunction epscell
!***********************************************************************
    elemental function minmod3(a,b,c) result(s)
!
!  Three-argument minmod: 0 unless all arguments share a sign, else the one
!  of smallest magnitude [C++ KTsolver::minmod, jax _minmod3].
!
      real, intent(in) :: a, b, c
      real :: s
!
      if (a>0.0 .and. b>0.0 .and. c>0.0) then
        s=min(a,b,c)
      elseif (a<0.0 .and. b<0.0 .and. c<0.0) then
        s=max(a,b,c)
      else
        s=0.0
      endif
!
    endfunction minmod3
!***********************************************************************
    subroutine project_admissible(u4,epsf)
!
!  Project a reconstructed face state onto the admissible set
!  [jax _project_admissible]:
!    1. positivity: T^00 - eps >= safety, and
!    2. causality:  |S| <= (T^00-eps) (1-safety), by rescaling S.
!
      real, dimension(nx,4), intent(inout) :: u4
      real, dimension(nx), intent(in) :: epsf
!
      real, dimension(nx) :: k0e, kmag, scal
      integer :: comp
!
      k0e=max(u4(:,1)-epsf,safety_kt)
      u4(:,1)=k0e+epsf
      kmag=sqrt(u4(:,2)**2+u4(:,3)**2+u4(:,4)**2+tiny_kt)
      scal=min(1.0,k0e*(1.0-safety_kt)/kmag)
      do comp=2,4
        u4(:,comp)=u4(:,comp)*scal
      enddo
!
    endsubroutine project_admissible
!***********************************************************************
    subroutine face_state(ul,ur,epsf,pl,pr,a)
!
!  Everything at one face that does NOT depend on which flux component is being
!  formed: both sides projected onto the admissible set, the bag-EOS pressure of
!  each, and the larger of the two local signal speeds.  Evaluated ONCE per face
!  and shared by every component -- which is what makes kt_div_tensor a single
!  pass rather than three single-component passes.
!
!  ul and ur are projected in place.
!
!  24-sep-2026/Isak Stomberg: factored out of face_flux and face_flux3
!
      real, dimension(nx,4), intent(inout) :: ul, ur
      real, dimension(nx), intent(in) :: epsf
      real, dimension(nx), intent(out) :: pl, pr, a
!
      call project_admissible(ul,epsf)
      call project_admissible(ur,epsf)
      call bag_pressure(ul,epsf,pl)
      call bag_pressure(ur,epsf,pr)
      a=max(state_signal_speed(ul,pl),state_signal_speed(ur,pr))
!
    endsubroutine face_state
!***********************************************************************
    function flux_component(ul,ur,pl,pr,a,dir,mu) result(h)
!
!  Rusanov (local Lax-Friedrichs) flux for ONE component mu at one face:
!    H = (F(U^L)+F(U^R))/2 - (a/2) (U^R - U^L)
!  [jax _kt_rhs / C++ H_temp assembly].  This is the only part of the flux that
!  depends on mu; the face must already have been prepared by face_state.
!
!  kt_div calls it once with mu=1 (energy), kt_div_tensor three times with
!  mu=2,3,4 on the same prepared face.
!
!  24-sep-2026/Isak Stomberg: coded
!
      real, dimension(nx,4), intent(in) :: ul, ur
      real, dimension(nx), intent(in) :: pl, pr, a
      integer, intent(in) :: dir, mu
      real, dimension(nx) :: h
!
      h=0.5*(physflux(ul,pl,dir,mu)+physflux(ur,pr,dir,mu)) &
        -0.5*a*(ur(:,mu)-ul(:,mu))
!
    endfunction flux_component
!***********************************************************************
    function physflux(u4,p,dir,mu) result(fl)
!
!  Physical flux of component mu in direction dir [jax _flux, C++ f_i_mu]:
!    energy   (mu=1)   : F = S_dir
!    momentum (mu=2..4): F = S_dir S_{mu-1}/(T^00+p) + p delta_{dir,mu-1}
!
!  The pressure is passed in rather than re-derived here: bag_pressure is the
!  single place the bag EOS is evaluated [jax compute_thermo, C++ pN_staggered],
!  and face_state has already called it for this face.
!
!  24-sep-2026/Isak Stomberg: takes p instead of epsf (was an inlined copy of
!                             bag_pressure)
!
      real, dimension(nx,4), intent(in) :: u4
      real, dimension(nx), intent(in) :: p
      integer, intent(in) :: dir, mu
      real, dimension(nx) :: fl
!
      if (mu==1) then
        fl=u4(:,1+dir)
      else
        fl=u4(:,1+dir)*u4(:,mu)/max(u4(:,1)+p,tiny_kt)
        if (mu==1+dir) fl=fl+p
      endif
!
    endfunction physflux
!***********************************************************************
    subroutine pushpars2c(p_par)
!
!  For review: restored from upstream adea2d9afe (GPU parameter marshalling).
!  cs_kt and tiny_kt are parameters, so they need no entry here.
!
    use Syscalls, only: copy_addr
    use General , only: string_to_enum

    integer, parameter :: n_pars=100
    integer(KIND=ikind8), dimension(n_pars) :: p_par
    call copy_addr(iux_kt,p_par(1)) ! int
    call copy_addr(ihless_kt,p_par(2)) ! int
    call copy_addr(eps_kt,p_par(3))
    call copy_addr(width_abs_kt,p_par(4))
    call copy_addr(safety_kt,p_par(5))
    call copy_addr(theta_kt,p_par(6))

    endsubroutine pushpars2c
!***********************************************************************
endmodule KT_transport
