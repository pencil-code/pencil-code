! $Id$
!
!  Kurganov-Tadmor (KT) flux-limited transport for the CONSERVATIVE
!  RELATIVISTIC (Higgsless) hydro equations
!
!      d_t T^00 + d_i S_i    = 0
!      d_t S_j  + d_i T^ij   = 0,   T^ij = S_i S_j/(T^00+p) + p delta_ij,
!
!  with the bag-model pressure p = w/4 - eps,  w = (4/3)(T^00-eps)(2 sqrt(1-lam)-1),
!  lam = (3/4)|S|^2/(T^00-eps)^2, and the externally prescribed vacuum energy
!  eps(t,x) of the Higgsless scheme (Jinno et al.), reconstructed here from the
!  bubble-wall crossing-time auxiliary (hless) exactly as in hydro.f90.
!
!  Scheme: MUSCL reconstruction with the generalized minmod limiter
!  (theta-limiter) + Rusanov (local Lax-Friedrichs) flux with a = 1
!  (speed of light; relativistic characteristics obey cs <= |v| <= 1),
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
!
module KT_transport
!
  use Cparam
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
  real, parameter :: a_face_kt=1.0        ! Rusanov speed (amax_kind='unity' in jax)
  real, parameter :: tiny_kt=1e-30
!
  contains
!***********************************************************************
    subroutine kt_init(irho_in,iux_in,ihless_in,eps_in,width_abs_in,theta_in)
!
!  Store farray indices and Higgsless/limiter parameters.
!  Called from initialize_hydro when lkt_transport=T.
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
!  call time in kt_transp instead of being stored here.
!
      call keep_compiler_quiet(irho_in)
      iux_kt=iux_in; ihless_kt=ihless_in
      eps_kt=eps_in; width_abs_kt=width_abs_in; theta_kt=theta_in
      ! safety_kt keeps its 1e-8 default (internal face-state admissible-projection floor)
      if (lroot) print*, 'kt_init: KT flux-limited transport enabled; theta=', theta_kt
!
    endsubroutine kt_init
!***********************************************************************
    subroutine kt_transp(f,m,n,mu,tcur,dq)
!
!  KT flux divergence for conserved component mu of one pencil row:
!    mu=1        : energy    dq = d_i S_i          (KT form)
!    mu=2,3,4    : momentum  dq = d_i T^{i,mu-1}   (KT form)
!  i.e. the caller subtracts dq from df.  All three directions are summed.
!
!  02-sep-2026/Isak Stomberg: coded
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      integer, intent(in) :: m, n, mu
      real, intent(in) :: tcur
      real, dimension(nx), intent(out) :: dq
!
      real, dimension(nx,4,-2:2) :: u
      real, dimension(nx,-1:1) :: epsc
      real, dimension(nx,4) :: sdm, sd0, sdp, ul, ur
      real, dimension(nx) :: hph, hmh, dl1
      integer :: dir, comp, o
!
      dq=0.0
      do dir=1,3
!
!  Gather the 5-point stencil of conserved variables and the 3-point
!  stencil of eps along direction dir.
!
        do o=-2,2
          u(:,1,o)=getcell(f,o,dir,m,n,irho)
          do comp=2,4
            u(:,comp,o)=getcell(f,o,dir,m,n,iux_kt+comp-2)
          enddo
        enddo
        do o=-1,1
          epsc(:,o)=epscell(f,o,dir,m,n,tcur)
        enddo
!
!  Limited slope-differences (slope*dx) at cell offsets -1, 0, +1:
!  minmod( theta*backward, centered, theta*forward )  [jax _limited_slope,
!  C++ del_x with theta=2].
!
        do comp=1,4
          sdm(:,comp)=minmod3(theta_kt*(u(:,comp,-1)-u(:,comp,-2)), &
                              0.5*(u(:,comp,0)-u(:,comp,-2)), &
                              theta_kt*(u(:,comp,0)-u(:,comp,-1)))
          sd0(:,comp)=minmod3(theta_kt*(u(:,comp,0)-u(:,comp,-1)), &
                              0.5*(u(:,comp,1)-u(:,comp,-1)), &
                              theta_kt*(u(:,comp,1)-u(:,comp,0)))
          sdp(:,comp)=minmod3(theta_kt*(u(:,comp,1)-u(:,comp,0)), &
                              0.5*(u(:,comp,2)-u(:,comp,0)), &
                              theta_kt*(u(:,comp,2)-u(:,comp,1)))
        enddo
!
!  Face +1/2: left/right MUSCL states and one eps per face.
!
        ul=u(:,:,0)+0.5*sd0
        ur=u(:,:,1)-0.5*sdp
        call face_flux(ul,ur,0.5*(epsc(:,0)+epsc(:,1)),dir,mu,hph)
!
!  Face -1/2.
!
        ul=u(:,:,-1)+0.5*sdm
        ur=u(:,:,0)-0.5*sd0
        call face_flux(ul,ur,0.5*(epsc(:,-1)+epsc(:,0)),dir,mu,hmh)
!
!  Divergence contribution: (H_{+1/2} - H_{-1/2}) / dx_dir.
!
        select case (dir)
          case (1); dl1=dx_1(l1:l2)
          case (2); dl1=dy_1(m)
          case (3); dl1=dz_1(n)
        endselect
        dq=dq+(hph-hmh)*dl1
      enddo
!
    endsubroutine kt_transp
!***********************************************************************
    function getcell(f,o,dir,m,n,ind) result(q)
!
!  Conserved variable at cell offset o along direction dir, for the row (m,n).
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
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
!  Vacuum energy eps(t,x) at cell offset o: the same crossing-time ramp as
!  hydro.f90's hydro_after_boundary_conservative (linear ramp of temporal
!  width width_abs_kt centred on the crossing time stored in the hless aux).
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      integer, intent(in) :: o, dir, m, n
      real, intent(in) :: tcur
      real, dimension(nx) :: eps, tau
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
!    2. causality:  |S| <= sqrt(4/3) (T^00-eps) (1-safety), by rescaling S.
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
      scal=min(1.0,sqrt(4.0/3.0)*k0e*(1.0-safety_kt)/kmag)
      do comp=2,4
        u4(:,comp)=u4(:,comp)*scal
      enddo
!
    endsubroutine project_admissible
!***********************************************************************
    subroutine face_flux(ul,ur,epsf,dir,mu,h)
!
!  Rusanov flux for component mu at one set of faces:
!    H = (F(U_L)+F(U_R))/2 - (a/2) (U_R - U_L),  a = a_face_kt
!  [jax _kt_rhs / C++ H_temp assembly].  Both states are projected first.
!
      real, dimension(nx,4), intent(inout) :: ul, ur
      real, dimension(nx), intent(in) :: epsf
      integer, intent(in) :: dir, mu
      real, dimension(nx), intent(out) :: h
!
      real, dimension(nx) :: fl, fr
!
      call project_admissible(ul,epsf)
      call project_admissible(ur,epsf)
      fl=physflux(ul,epsf,dir,mu)
      fr=physflux(ur,epsf,dir,mu)
      h=0.5*(fl+fr)-0.5*a_face_kt*(ur(:,mu)-ul(:,mu))
!
    endsubroutine face_flux
!***********************************************************************
    function physflux(u4,epsf,dir,mu) result(fl)
!
!  Physical flux of component mu in direction dir [jax _flux, C++ f_i_mu]:
!    energy   (mu=1)   : F = S_dir
!    momentum (mu=2..4): F = S_dir S_{mu-1}/(T^00+p) + p delta_{dir,mu-1}
!  with p from the bag EOS [jax compute_thermo, C++ pN_staggered].
!
      real, dimension(nx,4), intent(in) :: u4
      real, dimension(nx), intent(in) :: epsf
      integer, intent(in) :: dir, mu
      real, dimension(nx) :: fl
!
      real, dimension(nx) :: k0e, ki2, lam, w, p
!
      if (mu==1) then
        fl=u4(:,1+dir)
      else
        k0e=u4(:,1)-epsf
        ki2=u4(:,2)**2+u4(:,3)**2+u4(:,4)**2
        lam=min(0.75*ki2/max(k0e**2,tiny_kt),0.75)
        w=(4.0/3.0)*k0e*(2.0*sqrt(max(1.0-lam,0.0))-1.0)
        p=0.25*w-epsf
        fl=u4(:,1+dir)*u4(:,mu)/max(u4(:,1)+p,tiny_kt)
        if (mu==1+dir) fl=fl+p
      endif
!
    endfunction physflux
!***********************************************************************
endmodule KT_transport
