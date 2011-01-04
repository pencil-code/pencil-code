! $Id$
!
!  Module for NSCBC (Navier-Stokes Characteristic Boundary Conditions).
!  NSCBC is an alternative way of imposing (time-dependent) boundary
!  conditions through solving differential equations on the boundaries.
!
! 2009.09.09 (Nils Erland) : Moved all NSCBC stuff from boundcond.f90 to
!                            this module.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lnscbc = .true.
!
!***************************************************************
module NSCBC
!
  use Cdata
  use Cparam
  use Messages
  use Mpicomm
  use EquationOfState
!
  implicit none
!
include 'NSCBC.h'
!
! Format: nscbc_bc = 'bottom_x:top_x','bottom_y:top_y','bottom_z:top_z'
! for top and bottom boundary treatment at x, y and z-boundaries.
!
! nscbc_bc1(1) refers to boundary treatment at bottom x-boundary, nscbc_bc2(1)
! at top x-boundary etc.
! fbcx1, fbcx2 etc. are still used to impose values of variables at the
! boundaries.
!
! nscbc_sigma is a parameter describing how fast the velocity will reach
! an imposed velocity at an inlet in the absence of outgoing waves. How
! to set it is case dependent.
!
  character(len=2*nscbc_len+1), dimension(3) :: nscbc_bc=''
  character(len=nscbc_len), dimension(3) :: nscbc_bc1,nscbc_bc2
  character (len=labellen), dimension(ninit) :: inlet_profile='nothing'
  character(len=40) :: turb_inlet_dir=''
  character(len=40) :: turb_inlet_file='var.dat'
  real :: nscbc_sigma_out = 1.,nscbc_sigma_in = 1., p_infty=1.
  real :: transversal_damping=0.2, turb_profile_mag=1.0
  logical :: inlet_from_file=.false., jet_inlet=.false.
  logical :: first_NSCBC=.true.,onesided_inlet=.true.
  logical :: notransveral_terms=.true.
!
!  Variables to be used when getting timevarying inlet from file
!
  real, allocatable, dimension(:,:,:,:) :: f_in
  real, allocatable, dimension(:) :: x_in
  real, allocatable, dimension(:) :: y_in
  real, allocatable, dimension(:) :: z_in
  character :: prec_in
  real :: t_in,dx_in,dy_in,dz_in
  integer :: mx_in,my_in,mz_in,nv_in
  integer :: l1_in, nx_in, ny_in, nz_in
  integer :: mvar_in,maux_in,mglobal_in
  integer :: nghost_in,ipx_in, ipy_in, ipz_in
  integer :: m1_in
  integer :: n1_in
  integer :: l2_in
  integer :: m2_in
  integer :: n2_in
  real :: Lx_in
  real :: Ly_in
  real :: Lz_in
  real :: inlet_zz1, inlet_zz2
  real :: smooth_time=0.
!
!  Variables for inlet profiles
!
  real, dimension(2) :: radius_profile=(/0.0182,0.0364/)
  real, dimension(2) :: momentum_thickness=(/0.014,0.0182/)
  real, dimension(2) :: jet_center=(/0.,0./)
  real :: velocity_ratio=3.3, sigma=1.
  character (len=labellen), dimension(ninit) :: velocity_profile='nothing'
  character (len=labellen), dimension(ninit) :: zz_profile='nothing'
  logical :: lfinal_velocity_profile=.false.
!
  namelist /NSCBC_init_pars/  &
      nscbc_bc, nscbc_sigma_in, nscbc_sigma_out, p_infty, inlet_from_file,&
      turb_inlet_dir, jet_inlet,radius_profile,momentum_thickness,jet_center,&
      velocity_ratio,turb_inlet_file,turb_profile_mag,inlet_zz1,inlet_zz2,&
          zz_profile
!
  namelist /NSCBC_run_pars/  &
      nscbc_bc, nscbc_sigma_in, nscbc_sigma_out, p_infty, inlet_from_file,&
      turb_inlet_dir,jet_inlet,inlet_profile,smooth_time,onesided_inlet,&
      notransveral_terms, transversal_damping,radius_profile,momentum_thickness,&
      jet_center,velocity_ratio,turb_inlet_file,sigma,lfinal_velocity_profile,&
      velocity_profile,turb_profile_mag,inlet_zz1,inlet_zz2,zz_profile
!
  contains
!***********************************************************************
    subroutine nscbc_boundtreat(f,df)
!
!  Boundary treatment of the df-array.
!
!  This is a way to impose (time-
!  dependent) boundary conditions by solving a so-called characteristic
!  form of the fluid equations on the boundaries, as opposed to setting
!  actual values of the variables in the f-array. The method is called
!  Navier-Stokes characteristic boundary conditions (NSCBC).
!  The current implementation only solves a simplified version of the
!  equations, namely a set of so-called local one-dimensional inviscid
!  (LODI) relations. This means that transversal and viscous terms are
!  dropped on the boundaries.
!
!  The treatment should be done after the y-z-loop, but before the Runge-
!  Kutta solver adds to the f-array.
!
!   7-jul-08/arne: coded.
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df

      intent(inout) :: f
      intent(inout) :: df
!
      if (nscbc_bc1(1) /= '' .or. nscbc_bc2(1) /= '') &
          call nscbc_boundtreat_xyz(f,df,1)
      if (nscbc_bc1(2) /= '' .or. nscbc_bc2(2) /= '') &
          call nscbc_boundtreat_xyz(f,df,2)
      if (nscbc_bc1(3) /= '' .or. nscbc_bc2(3) /= '') &
          call nscbc_boundtreat_xyz(f,df,3)
!
    endsubroutine nscbc_boundtreat
!***********************************************************************
    subroutine nscbc_boundtreat_xyz(f,df,j)
!
!   NSCBC boundary treatment.
!   j = 1, 2 or 3 for x, y or z-boundaries respectively.
!
!   7-jul-08/arne: coded.

!
      use General, only: safe_character_assign, chn

      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      character (len=nscbc_len), dimension(3) :: bc12
      character (len=3) :: topbot
      character (len=60) :: turbfile
      integer i,j,k,ip_ok,ip_test
      integer :: imin,imax,jmin,jmax,igrid,jgrid
      real, dimension(mcom) :: valx,valy,valz
      logical :: proc_at_inlet
      integer :: ipx_in, ipy_in, ipz_in, iproc_in, nprocx_in, nprocy_in, nprocz_in
      character (len=120) :: directory_in
      character (len=5) :: chproc_in
      real :: T_t,u_t
      real, dimension(nchemspec) :: YYk

      intent(inout) :: f
      intent(inout) :: df
      intent(in)    :: j
!
      proc_at_inlet=.false.
!
!  Create the final velocity profile
!
      if (lfinal_velocity_profile) then
        if (j==1) then
          imin=m1; imax=m2; jmin=n1; jmax=n2
          igrid=ny; jgrid=nz
        elseif (j==2) then
          imin=l1; imax=l2; jmin=n1; jmax=n2
          igrid=nx; jgrid=nz
        elseif (j==3) then
          imin=l1; imax=l2; jmin=m1; jmax=m2
          igrid=nx; jgrid=ny
        endif
      endif
!
      do k=1,2                ! loop over 'bot','top'
        if (k==1) then
          topbot='bot'; bc12(j)=nscbc_bc1(j);!val=bt_val1(j)
          valx=fbcx1; valy=fbcy1; valz=fbcz1; ip_ok=0
        else
          topbot='top'; bc12(j)=nscbc_bc2(j);!val=bt_val2(j)
          valx=fbcx2; valy=fbcy2; valz=fbcz2
          if (j==1) ip_ok=nprocx-1
          if (j==2) ip_ok=nprocy-1
          if (j==3) ip_ok=nprocz-1
        endif
        if (j==1) ip_test=ipx
        if (j==2) ip_test=ipy
        if (j==3) ip_test=ipz
!
! Read data to be used at the inlet from file
!
        if (inlet_from_file) then
!
! Read data from file only initially
! At later times this is stored in processor memory.
! NILS: Is it really a good idea to store all this in memory?
!
          if (first_NSCBC) then
!
! Check which processor we want to read from.
! In the current implementation it is required that:
!   1) The number of mesh points and processors at the interface between
!      the two computational domains are equal. The two comp. domains I
!      am refering to here is the domain of the current simulation and the
!      domain of the pre-run isotropic turbulence simulation defining the
!      turbulence at the inlet.
!   2) The pre-run simulaion can not have multiple processors in the flow
!      direction of the current simulation.
!
            if (lprocz_slowest) then
              ipx_in=ipx
              ipy_in=ipy
              ipz_in=ipz
              nprocx_in=nprocx
              nprocy_in=nprocy
              nprocz_in=nprocz
              if (j==1) then
                if ((topbot=='bot'.and.lfirst_proc_x).or. &
                    (topbot=='top'.and.llast_proc_x)) then
                  proc_at_inlet=.true.
                  ipx_in=0
                  nprocx_in=1
                endif
              elseif (j==2) then
                if ((topbot=='bot'.and.lfirst_proc_y).or. &
                    (topbot=='top'.and.llast_proc_y)) then
                  proc_at_inlet=.true.
                  ipy_in=0
                  nprocy_in=1
                endif
              elseif (j==3) then
                if ((topbot=='bot'.and.lfirst_proc_z).or. &
                    (topbot=='top'.and.llast_proc_z)) then
                  proc_at_inlet=.true.
                  ipz_in=0
                  nprocz_in=1
                endif
              else
                call fatal_error("nscbc_boundtreat_xyz",'No such direction!')
              endif
              iproc_in=ipz_in*nprocy_in*nprocx_in+ipy_in*nprocx_in+ipx_in
            else
              call fatal_error("nscbc_boundtreat_xyz",&
                  'lprocz_slowest=F not implemeted for inlet from file!')
            endif
!
!  Read data only if required, i.e. if we are at a processor handling inlets
!
            if (proc_at_inlet) then
              call chn(iproc_in,chproc_in)
              call safe_character_assign(directory_in,&
                  trim(turb_inlet_dir)//'/data/proc'//chproc_in)
              call safe_character_assign(turbfile,&
                  trim(directory_in)//'/'//trim(turb_inlet_file))
              open(1,FILE=turbfile,FORM='unformatted')
              if (ip<=8) print*,'input: open, mx_in,my_in,mz_in,nv_in=',&
                  mx_in,my_in,mz_in,nv_in
              read(1) f_in
              read(1) t_in,x_in,y_in,z_in,dx_in,dy_in,dz_in
              nx_in=mx_in-2*nghost
              ny_in=my_in-2*nghost
              nz_in=mz_in-2*nghost
              l1_in=nghost+1
              m1_in=nghost+1
              n1_in=nghost+1
              l2_in=mx_in-nghost
              m2_in=my_in-nghost
              n2_in=mz_in-nghost
              Lx_in=x_in(l2_in+1)-x_in(l1_in)
              Ly_in=y_in(m2_in+1)-y_in(m1_in)
              Lz_in=z_in(n2_in+1)-z_in(n1_in)
              first_NSCBC=.false.
              close(1)
            endif
          endif
        endif
!
!  Check if this is a physical boundary
!
        if (ip_test==ip_ok) then
!
!  Set the values of T_t and u_t dependent on direction
!
          T_t=0
          if (j==1) then
            if (ilnTT > 0) T_t=valx(ilnTT)
            u_t=valx(j)
          elseif (j==2) then
            if (ilnTT > 0) T_t=valy(ilnTT)
            u_t=valy(j)
          elseif (j==3) then
            if (ilnTT > 0) T_t=valz(ilnTT)
            u_t=valz(j)
          endif
!
! Set the values of species
!
       if (ichemspec(1)>0) then
         YYk=0.
         do i=1,nchemspec
           YYk(i)=valx(ichemspec(i))
         enddo
       endif
!
!  Do the NSCBC boundary
!
          select case (bc12(j))
!
          case ('part_ref_outlet')
            call bc_nscbc_prf(f,df,j,topbot,.false.)
!
          case ('part_ref_inlet')
            call bc_nscbc_prf(f,df,j,topbot,.true.,linlet=.true.,u_t=u_t,T_t=T_t,YYk=YYk)
!
          case ('ref_inlet')
            call bc_nscbc_prf(f,df,j,topbot,.false.,linlet=.true.,u_t=u_t,T_t=T_t,YYk=YYk)
!
          case ('subsonic_inflow')
            if (j==1) then
              call bc_nscbc_subin_x(f,df,topbot,valx)
            elseif (j==2) then
            endif
          case ('subsonic_inflow_new')
            if (j==1) then
              call bc_nscbc_subin_x_new(f,df,topbot,valx)
            elseif (j==2) then
            endif
!
          case ('subson_nref_outflow')
            if (j==1) then
              call bc_nscbc_nref_subout_x(f,df,topbot,nscbc_sigma_out)
            elseif (j==2) then
              call bc_nscbc_nref_subout_y(f,df,topbot,nscbc_sigma_out)
            elseif (j==3) then
              call bc_nscbc_nref_subout_z(f,df,topbot,nscbc_sigma_out)
            endif
          case ('no_nscbc')
              call no_nscbc(f,topbot,valx)
          case ('')
!   Do nothing.
          case ('none')
            print*,'nscbc_boundtreat_xyz: doing nothing!'
!   Do nothing.
          case default
            call fatal_error("nscbc_boundtreat_xyz",&
                'You must specify nscbc bouncond!')
          endselect
          if (lfinal_velocity_profile) then
            call final_velocity_profile(f,j,topbot, &
             imin,imax,jmin,jmax,igrid,jgrid)
          endif
        endif
      enddo


!
    endsubroutine nscbc_boundtreat_xyz
!***********************************************************************
    subroutine bc_nscbc_prf(f,df,dir,topbot,non_reflecting_inlet,linlet,u_t,T_t,YYk)
!
!   Calculate du and dlnrho at a partially reflecting outlet/inlet normal to
!   x-direction acc. to LODI relations. Uses a one-sided finite diff. stencil.
!
!   7-jul-08/arne: coded.
!  25-nov-08/nils: extended to work in multiple dimensions and with cross terms
!                  i.e. not just the LODI equations.
!  22-jan-10/nils: made general with respect to direction such that we do not
!                  need three different routines for the three different
!                  directions.
!  26-jui-10/julien: added the full composition array as an optional argument to
!                    enable a relaxation of the species mass fractions at the inlet
!                    in the case of a multi-species flow
!
      use Deriv, only: der_onesided_4_slice, der_pencil, der2_pencil
      use Chemistry
      use General, only: random_number_wrapper

      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      character (len=3) :: topbot
      logical, optional :: linlet
      logical :: llinlet, non_reflecting_inlet
      real, optional :: u_t, T_t
      real, optional, dimension(nchemspec) :: YYk
      real :: Mach,KK,nu, cs0_average
      integer, dimension(33) :: stat
      integer lll,k,ngridpoints,imin,imax,jmin,jmax,i
      integer sgn,dir,iused,dir1,dir2,dir3,igrid,jgrid
      integer :: imass=1
      logical :: non_zero_transveral_velo
      real, allocatable, dimension(:,:,:,:) :: dui_dxj
      real, allocatable, dimension(:,:,:) :: &
          fslice, dfslice,grad_rho,u_in,grad_T,grad_P
      real, allocatable, dimension(:,:) :: &
          TT,mu1,grad_mu1,rho0,P0,L_1,L_2,L_3,L_4,L_5,&
          prefac1,prefac2,T_1,T_2,T_3,T_4,T_5,cs,&
          cs2,gamma,dY_dx,scalar_profile
      real, allocatable, dimension(:,:,:) :: L_k, dYk_dx
      real, allocatable, dimension(:,:,:) :: YYk_full
      real, allocatable, dimension(:,:)   :: sum_Lk
!
      intent(inout) :: f
      intent(inout) :: df
!
!  Define the direction of this boundary
!
      llinlet = .false.
      if (present(linlet)) llinlet = linlet
      if (llinlet.and..not.present(u_t)) call stop_it(&
           'bc_nscbc_prf: when using linlet=T, you must also specify u_t)')
      if (llinlet.and.ilnTT>0.and..not.present(T_t)) call stop_it(&
           'bc_nscbc_prf: when using linlet=T, you must also specify T_t)')
      select case (topbot)
      case ('bot')
        if (dir == 1) lll = l1
        if (dir == 2) lll = m1
        if (dir == 3) lll = n1
        sgn = 1
      case ('top')
        if (dir == 1) lll = l2
        if (dir == 2) lll = m2
        if (dir == 3) lll = n2
        sgn = -1
      case default
        print*, "bc_nscbc_prf: ", topbot, " should be `top' or `bot'"
      end select
!
!  Set some auxillary variables
!
      if (dir==1) then
        dir1=1; dir2=2; dir3=3
        igrid=ny; jgrid=nz
        imin=m1; imax=m2
        jmin=n1; jmax=n2
      elseif (dir==2) then
        dir1=2; dir2=1; dir3=3
        igrid=nx; jgrid=nz
        imin=l1; imax=l2
        jmin=n1; jmax=n2
     elseif (dir==3) then
        dir1=3; dir2=1; dir3=2
        igrid=nx; jgrid=ny
        imin=l1; imax=l2
        jmin=m1; jmax=m2
      else
        call fatal_error('bc_nscbc_prf:','No such dir!')
      endif
!
!  Allocate all required arrays
!
      stat=0
      allocate(  fslice(igrid,jgrid,mfarray),STAT=stat(1))
      allocate( dfslice(igrid,jgrid,mvar   ),STAT=stat(2))
      allocate(      TT(igrid,jgrid),        STAT=stat(3))
      allocate(     mu1(igrid,jgrid),        STAT=stat(4))
      allocate(grad_mu1(igrid,jgrid),        STAT=stat(5))
      allocate(    rho0(igrid,jgrid),        STAT=stat(6))
      allocate(      P0(igrid,jgrid),        STAT=stat(7))
      allocate(     L_1(igrid,jgrid),        STAT=stat(8))
      allocate(     L_2(igrid,jgrid),        STAT=stat(9))
      allocate(     L_3(igrid,jgrid),        STAT=stat(10))
      allocate(     L_4(igrid,jgrid),        STAT=stat(11))
      allocate(     L_5(igrid,jgrid),        STAT=stat(12))
      allocate( prefac1(igrid,jgrid),        STAT=stat(13))
      allocate( prefac2(igrid,jgrid),        STAT=stat(14))
      allocate(     T_1(igrid,jgrid),        STAT=stat(15))
      allocate(     T_2(igrid,jgrid),        STAT=stat(16))
      allocate(     T_3(igrid,jgrid),        STAT=stat(17))
      allocate(     T_4(igrid,jgrid),        STAT=stat(18))
      allocate(     T_5(igrid,jgrid),        STAT=stat(19))
      allocate(      cs(igrid,jgrid),        STAT=stat(20))
      allocate(     cs2(igrid,jgrid),        STAT=stat(21))
      allocate(   gamma(igrid,jgrid),        STAT=stat(22))
      allocate(dYk_dx(igrid,jgrid,nchemspec),  STAT=stat(23))
      allocate(grad_rho(igrid,jgrid,3),      STAT=stat(24))
      allocate(    u_in(igrid,jgrid,3),      STAT=stat(25))
      allocate(  grad_T(igrid,jgrid,3),      STAT=stat(26))
      allocate(  grad_P(igrid,jgrid,3),      STAT=stat(27))
      allocate( dui_dxj(igrid,jgrid,3,3),    STAT=stat(28))
      if (lpscalar_nolog) allocate(scalar_profile(igrid,jgrid),STAT=stat(29))
      allocate(L_k(igrid,jgrid,nchemspec),   STAT=stat(30))
      allocate(YYk_full(igrid,jgrid,nchemspec), STAT=stat(31))
      allocate(dY_dx(igrid,jgrid),           STAT=stat(32))
      allocate(sum_Lk(igrid,jgrid),          STAT=stat(33))
!
      if (maxval(stat) > 0) &
          call stop_it("Couldn't allocate memory for all vars in bc_nscbc_prf")
!
!  Initialize fslice and dfslice
!
      if (dir == 1) then
        fslice=f(lll,m1:m2,n1:n2,:)
        dfslice=df(lll,m1:m2,n1:n2,:)
      elseif (dir == 2) then
        fslice=f(l1:l2,lll,n1:n2,:)
        dfslice=df(l1:l2,lll,n1:n2,:)
      elseif (dir == 3) then
        fslice=f(l1:l2,m1:m2,lll,:)
        dfslice=df(l1:l2,m1:m2,lll,:)
      else
        call fatal_error('bc_nscbc_prf','No such dir!')
      endif
!
      if (ichemspec(1)>0) then
        i = 1
        do k=ichemspec(1), ichemspec(nchemspec)
          call der_onesided_4_slice(f,sgn,k,dY_dx,lll,dir1)
          dYk_dx(:,:,i) = dY_dx(:,:)
          i=i+1
        enddo
      endif
!
!  Find derivatives at boundary
!
      call derivate_boundary(f,sgn,dir1,lll,dui_dxj,grad_T,grad_rho)
!
!  Get some thermodynamical variables
!  This includes values and gradients of; density, temperature, pressure
!  In addition also speed of sound, gamma and mu1 is found.
!
      call get_thermodynamics(mu1,grad_mu1,gamma,cs2,cs,rho0,TT,P0,nu,&
          grad_rho,grad_P,grad_T,lll,sgn,dir1,fslice,T_t,llinlet,&
          imin,imax,jmin,jmax)
!
!  Define some prefactors to be used later
!
      prefac1 = -1./(2.*cs2)
      prefac2 = -1./(2.*rho0*cs)
!
!  Find Mach number
!  (NILS: I do not think this is a good way to determine the Mach
!  number since this is a local Mach number for this processor. Furthermore
!  by determining the Mach number like this we will see the Mach number varying
!  with the phase of an acoustic wave as the wave pass through the boundary.
!  I think that what we really want is a Mach number averaged over the
!  timescale of several acoustic waves. How could this be done????)
!
      ngridpoints=igrid*jgrid
      Mach=sum(fslice(:,:,dir1)/cs)/ngridpoints
!
!  We will need the transversal terms of the waves entering the domain
!
      call transversal_terms(T_1,T_2,T_3,T_4,T_5,rho0,cs,&
          fslice,grad_rho,grad_P,dui_dxj,dir1,dir2,dir3)
!
      if (llinlet) then
!
!  Find the velocity and composition to be used at the inlet.
!
        if (dir==1) then
          call find_velocity_at_inlet(u_in,non_zero_transveral_velo,&
              Lx_in,nx_in,u_t,dir,m1_in,m2_in,n1_in,n2_in,imin,imax,jmin,jmax,&
              igrid,jgrid,scalar_profile)
        elseif (dir==2) then
          call find_velocity_at_inlet(u_in,non_zero_transveral_velo,&
              Ly_in,ny_in,u_t,dir,l1_in,l2_in,n1_in,n2_in,imin,imax,jmin,jmax,&
              igrid,jgrid,scalar_profile)
        elseif (dir==3) then
          call find_velocity_at_inlet(u_in,non_zero_transveral_velo,&
              Lz_in,nz_in,u_t,dir,l1_in,l2_in,m1_in,m2_in,imin,imax,jmin,jmax,&
              igrid,jgrid,scalar_profile)
        endif
!
        call find_composition_at_inlet(nchemspec,YYk,YYk_full,&
             dir,imin,imax,jmin,jmax,igrid,jgrid)
!
        if (lroot .and. ip<5) then
          print*,'bc_nscbc_prf: Finalized reading velocity and composition profiles at the inlet.'
        endif

!
!  Having found the velocity at the inlet we are now ready to start
!  defining the L's, which are really the Lodi equations.
!
        if (non_reflecting_inlet) then
          L_1 = (fslice(:,:,dir1) - sgn*cs)&
              *(grad_P(:,:,dir1) - sgn*rho0*cs*dui_dxj(:,:,dir1,dir1))
          if (ilnTT>0) then
!            L_2=nscbc_sigma_in*(TT-T_t)&
!              *cs*rho0*Rgas/Lxyz(dir1)-(cs2*T_1-T_5)
!  Julien: The above formula seems erroneous which could explain the following
!          remark from Nils. The following correction is proposed but needs
!          some tests
!  Corrected according to Yoo et al. (Combustion Theory and Modelling, 2005)
!           L_2 = nscbc_sigma_in*(TT-T_t)&
!                *rho0*Rgas/(cs*Lxyz(dir1))-(cs2*T_1-T_5)
!
!  NILS: There seems to be something wrong with the calculation of L_2, as the
!  NILS: coded crashes with the version above. For now L_2 is therefore
!  NILS: set to zero. This must be fixed!!!!!!!
!
            L_2=0
          else
            L_2=0
          endif
!
!  The inlet is non-reflecting only when nscbc_sigma_in is set to 0, this
!  might however lead to problems as the inlet velocity will tend to drift
!  away from the target velocity u_t. This problem should be overcome by
!  setting a small but non-zero nscbc_sigma_in.
!
          L_3=nscbc_sigma_in*(fslice(:,:,dir2)-u_in(:,:,dir2))&
              *cs/Lxyz(dir1)-T_3
          L_4=nscbc_sigma_in*(fslice(:,:,dir3)-u_in(:,:,dir3))&
              *cs/Lxyz(dir1)-T_4
          L_5 = nscbc_sigma_in*cs2*rho0&
              *sgn*(fslice(:,:,dir1)-u_in(:,:,dir1))*(1-Mach**2)/Lxyz(dir1)&
              -(T_5+sgn*rho0*cs*T_2)
          if (ichemspec(1)>0) then
            do k = 1, nchemspec
              L_k(:,:,k)=nscbc_sigma_in*(fslice(:,:,ichemspec(k))-YYk_full(:,:,k))&
                  *cs/Lxyz(dir1)
            enddo
          else
            L_k=0.
          endif
        else
          L_3=0
          L_4=0
          L_5 = L_1
          L_2=0.5*(gamma-1)*(L_5+L_1)
          L_k=0.
        endif
      else
!
!  Find the L_i's.
!
        cs0_average=sum(cs)/ngridpoints
        KK=nscbc_sigma_out*(1.-Mach**2)*cs0_average/Lxyz(dir1)
        L_1 = KK*(P0-p_infty)-(T_5-sgn*rho0*cs*T_2)*(1-transversal_damping)
        if (ilnTT > 0) then
          L_2=fslice(:,:,dir1)*(cs2*grad_rho(:,:,dir1)-grad_P(:,:,dir1))
        else
          L_2=0
        endif
        L_3 = fslice(:,:,dir1)*dui_dxj(:,:,dir2,dir1)
        L_4 = fslice(:,:,dir1)*dui_dxj(:,:,dir3,dir1)
        L_5 = (fslice(:,:,dir1) - sgn*cs)*(grad_P(:,:,dir1)&
             - sgn*rho0*cs*dui_dxj(:,:,dir1,dir1))
        if (ichemspec(1)>0) then
          do k = 1, nchemspec
            L_k(:,:,k)=fslice(:,:,dir1)*dYk_dx(:,:,k)
          enddo
        else
          L_k=0.
        endif
      endif
!
      if (lroot .and. ip<5) then
        print*,'bc_nscbc_prf: Finalized setting up the Ls.'
      endif
!
!  Find the evolution equation for the normal velocity at the boundary
!  For 'top' L_1 plays the role of L5 and L_5 the role of L1
!
      select case (topbot)
      case ('bot')
        if (llinlet) then
          dfslice(:,:,dir1) = prefac2*( L_5 - L_1)-T_2
        else
          dfslice(:,:,dir1) = prefac2*( L_1 - L_5)+T_2
        endif
      case ('top')
        if (llinlet) then
          dfslice(:,:,dir1) = prefac2*( L_1 - L_5)+T_2
        else
          dfslice(:,:,dir1) = prefac2*(-L_1 + L_5)-T_2
        endif
      endselect
!
!  Find the evolution equation for the other variables at the boundary
!
      dfslice(:,:,ilnrho) = prefac1*(2*L_2 + L_1 + L_5)-T_1
      if (ilnTT>0) then
        sum_Lk=0.
        if (ichemspec(1)>0) then
          do k = 1,nchemspec
            sum_Lk = sum_Lk + (rho0*cs2)/(species_constants(k,imass)*mu1)*L_k(:,:,k)
          enddo
        endif
        dfslice(:,:,ilnTT) = -1./(rho0*cs2)*(-L_2+0.5*(gamma-1.)*(L_5+L_1)&
                             - sum_Lk)*TT + TT*(T_1/rho0-T_5/P0)
      endif
      dfslice(:,:,dir2) = -L_3-T_3
      dfslice(:,:,dir3) = -L_4-T_4
      if (ichemspec(1)>0) then
        do k = 1, nchemspec
          dfslice(:,:,ichemspec(k)) = -L_k(:,:,k)
        enddo
      endif
!
!  Check if we are solving for logrho or rho
!
      if (.not. ldensity_nolog) then
        dfslice(:,:,ilnrho)=dfslice(:,:,ilnrho)/rho0
      endif
!
!  Check if we are solving for logT or T
!
      if (.not. ltemperature_nolog .and. ilnTT>0) then
        dfslice(:,:,ilnTT)=dfslice(:,:,ilnTT)/TT
      endif
!
      if (lroot .and. ip<5) then
        print*,'bc_nscbc_prf: Finalized setting up the dfslice array.'
      endif
!
! Impose required variables at the boundary for reflecting inlets
!
      if (llinlet) then
        if (.not. non_reflecting_inlet) then
          fslice(:,:,dir1) = u_in(:,:,dir1)
          if (non_zero_transveral_velo) then
            fslice(:,:,dir2) = u_in(:,:,dir2)
            fslice(:,:,dir3) = u_in(:,:,dir3)
          else
            fslice(:,:,dir2) = 0.
            fslice(:,:,dir3) = 0.
          endif
          if (iTT>0) then
            fslice(:,:,iTT) = T_t
          elseif (ilnTT>0) then
            fslice(:,:,ilnTT) = log(T_t)
          endif
!
! Set the derivatives of these variables to zero
! NILS: Is this really required?
!
          dfslice(:,:,dir1)=0
          dfslice(:,:,dir2)=0
          dfslice(:,:,dir3)=0
          if (ilnTT>0) dfslice(:,:,ilnTT)=0
          if (ichemspec(1)>0) then
            do k = 1, nchemspec
              fslice(:,:,ichemspec(k)) = YYk_full(:,:,k)
              dfslice(:,:,ichemspec(k)) = 0
            enddo
          endif
        endif
      endif
!
! Treat all other variables as passive scalars
!
      iused=max(ilnTT,ilnrho)
      if (mvar>iused) then
        do k=iused+nchemspec+1,mvar
          if (llinlet) then
            if (k == icc) then
              fslice(:,:,icc)=scalar_profile
              dfslice(:,:,icc)=0
            else
              call der_onesided_4_slice(f,sgn,k,dY_dx,lll,dir1)
              dfslice(:,:,k)=-fslice(:,:,dir1)*dY_dx
            endif
          else
            call der_onesided_4_slice(f,sgn,k,dY_dx,lll,dir1)
            dfslice(:,:,k)=-fslice(:,:,dir1)*dY_dx
          endif
        enddo
      endif
!
!  Put everything that has been temporarily stored in dfslice back into
!  the df array
!
      if (dir == 1) then
        df(lll,m1:m2,n1:n2,:)=dfslice
        f( lll,m1:m2,n1:n2,:)=fslice
      elseif (dir == 2) then
        df(l1:l2,lll,n1:n2,:)=dfslice
        f( l1:l2,lll,n1:n2,:)=fslice
      elseif (dir == 3) then
        df(l1:l2,m1:m2,lll,:)=dfslice
        f( l1:l2,m1:m2,lll,:)=fslice
      else
        call fatal_error('bc_nscbc_prf','No such dir!')
      endif
!
      if (lroot .and. ip<5) then
        print*,'bc_nscbc_prf: Finalized bc_nscbc_prf.'
      endif
!
    endsubroutine bc_nscbc_prf
!***********************************************************************
    subroutine read_NSCBC_init_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=NSCBC_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=NSCBC_init_pars,ERR=99)
      endif
!
      if (lnscbc) call parse_nscbc(nscbc_bc,nscbc_bc1,nscbc_bc2)
!
99    return
!
    endsubroutine read_NSCBC_init_pars
!***********************************************************************
    subroutine read_NSCBC_run_pars(unit,iostat)

      use Sub, only : rdim

      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
      integer :: stat
      logical :: exist
      character (len=500) :: file
!
! Define default for the inlet_profile for backward compatibility
!
      inlet_profile(1)='uniform'
      zz_profile(1)='uniform'
!
      if (present(iostat)) then
        read(unit,NML=NSCBC_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=NSCBC_run_pars,ERR=99)
      endif
!
      if (lnscbc) call parse_nscbc(nscbc_bc,nscbc_bc1,nscbc_bc2)
!
! Check if we will read turbulent inlet data from data file
!
      if (inlet_from_file) then
        print*,'Read inlet data from file!'
!
! Read the size of the data to be found on the file.
!
        file=trim(turb_inlet_dir)//'/data/proc0/dim.dat'
        inquire(FILE=trim(file),EXIST=exist)
        if (exist) then
          call rdim(file,&
              mx_in,my_in,mz_in,mvar_in,maux_in,mglobal_in,prec_in,&
              nghost_in,ipx_in, ipy_in, ipz_in)
          nv_in=mvar_in+maux_in+mglobal_in
        else
          print*,'file=',file
          call stop_it('read_NSCBC_run_pars: Could not find file!')
        endif
!
! Allocate array for data to be used at the inlet.
! For now every processor reads all the data - this is clearly an overkill,
! but we leav it like this during the development of this feature.
!
        allocate( f_in(mx_in,my_in,mz_in,nv_in),STAT=stat)
        if (stat>0) call stop_it("Couldn't allocate memory for f_in ")
        allocate( x_in(mx_in),STAT=stat)
        if (stat>0) call stop_it("Couldn't allocate memory for x_in ")
        allocate( y_in(my_in),STAT=stat)
        if (stat>0) call stop_it("Couldn't allocate memory for y_in ")
        allocate( z_in(mz_in),STAT=stat)
        if (stat>0) call stop_it("Couldn't allocate memory for z_in ")
      endif
!
99    return
!
    endsubroutine read_NSCBC_run_pars
!***********************************************************************
    subroutine write_NSCBC_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit,NML=NSCBC_init_pars)
!
    endsubroutine write_NSCBC_init_pars
!***********************************************************************
    subroutine write_NSCBC_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit,NML=NSCBC_run_pars)
!
    endsubroutine write_NSCBC_run_pars
!***********************************************************************
    subroutine parse_nscbc(bc,bc1,bc2)
!
!  Parse boundary conditions, which may be in the form `a' (applies to
!  both `lower' and `upper' boundary) or `a:s' (use `a' for lower,
!  `s' for upper boundary.
!
!  Be aware of the difference in format when giving NSCBC boundary
!  conditions and usual boundary conditions through the bcx, bcy, bcz
!  variables:
!
!  bc = 'top_x:bottom_x','top_y:bottom_y','top_z:bottom_z'
!  for top and bottom boundary treatment at x, y and z-boundaries.
!
!  bc1(1) then refers to boundary treatment at bottom x-boundary, bc2(1)
!  at top x-boundary etc.
!  Ie. this routine sets the boundary condition variables, bc1 and bc2,
!  for x, y and z-boundaries all at once.
!
!   7-jul-08/arne: adapted from parse_bc
!
      character (len=2*nscbc_len+1), dimension(3) :: bc
      character (len=nscbc_len), dimension(3) :: bc1,bc2
      integer :: j,isep
!
      intent(in) :: bc
      intent(out) :: bc1,bc2
!
      do j=1,3
        isep = index(bc(j),':')
        if (isep > 0) then
          bc1(j) = bc(j)(1:isep-1)
          bc2(j) = bc(j)(isep+1:)
        else
          bc1(j) = bc(j)(1:nscbc_len)
          bc2(j) = bc(j)(1:nscbc_len)
        endif
      enddo
!
    endsubroutine
!***********************************************************************
    subroutine NSCBC_clean_up
!
!  Deallocate all allocatable arrays
!
      if (allocated(f_in)) deallocate(f_in)
      if (allocated(x_in)) deallocate(x_in)
      if (allocated(y_in)) deallocate(y_in)
      if (allocated(z_in)) deallocate(z_in)
!
    endsubroutine NSCBC_clean_up
!***********************************************************************
    subroutine find_velocity_at_inlet(u_in,non_zero_transveral_velo,&
        domain_length,grid_points,u_t,dir,imin_turb,imax_turb,&
        jmin_turb,jmax_turb,imin,imax,jmin,jmax,igrid,jgrid,scalar_profile)
!
!  Find velocity at inlet.
!
!  2010.01.20/Nils Erland L. Haugen: coded
!
      logical, intent(out) :: non_zero_transveral_velo
      real, dimension(:,:,:), intent(out) :: u_in
      real, dimension(:,:), intent(out) :: scalar_profile
      real, intent(in) :: domain_length,u_t
      integer, intent(in) :: grid_points,dir,imin,imax,jmin,jmax,igrid,jgrid
      integer, intent(in) :: imin_turb,imax_turb,jmin_turb,jmax_turb
!
      real :: shift, grid_shift, weight, round
      integer :: iround, lowergrid, uppergrid
      real, allocatable, dimension(:,:)   :: u_profile
      real, allocatable, dimension(:,:,:) :: u_turb
      real, dimension(3) :: velo
      real :: radius_mean, smooth, rad
      integer :: i,j,kkk,jjj
      integer, dimension(10) :: stat
!
!  Allocate allocatables
!
      stat=0
      allocate(u_profile(igrid,jgrid  ),STAT=stat(1))
      allocate(u_turb   (igrid,jgrid,3),STAT=stat(2))
      if (maxval(stat)>0) &
          call stop_it("Couldn't allocate memory for all vars in find_velocity_at_inlet")
!
! Define velocity profile at inlet
!
        u_in=0
        do j=1,ninit
          select case (inlet_profile(j))
!
          case ('nothing')
            if (lroot .and. it==1 .and. j == 1 .and. lfirst) &
                print*,'inlet_profile: nothing'
            non_zero_transveral_velo=.false.
            u_profile=1.
!
          case ('uniform')
            if (lroot .and. it==1 .and. lfirst) &
                print*,'inlet_profile: uniform'
            non_zero_transveral_velo=.false.
            u_in(:,:,dir)=u_in(:,:,dir)+u_t
            u_profile=1.
!
          case ('coaxial_jet')
            if (lroot .and. it==1 .and. lfirst) &
                print*,'inlet_profile: coaxial_jet'
            non_zero_transveral_velo=.true.
            velo(1)=u_t
            velo(2)=velo(1)*velocity_ratio
            velo(3)=0.04*velo(2)
            radius_mean=(radius_profile(1)+radius_profile(2))/2.
            do jjj=imin,imax
              do kkk=jmin,jmax
                if (dir==1) then
                  rad=sqrt((y(jjj)-jet_center(1))**2+(z(kkk)-jet_center(2))**2)
                elseif (dir==2) then
                  rad=sqrt((x(jjj)-jet_center(1))**2+(z(kkk)-jet_center(2))**2)
                elseif (dir==3) then
                  rad=sqrt((x(jjj)-jet_center(1))**2+(y(kkk)-jet_center(2))**2)
                endif
                  ! Add mean velocity profile
                if (rad < radius_mean) then
                  u_in(jjj-imin+1,kkk-jmin+1,dir)&
                      =u_in(jjj-imin+1,kkk-jmin+1,dir)&
                      +(velo(1)+velo(2))/2&
                      +(velo(2)-velo(1))/2*tanh((rad-radius_profile(1))/&
                      (2*momentum_thickness(1)))
                else
                  u_in(jjj-imin+1,kkk-jmin+1,dir)&
                      =u_in(jjj-imin+1,kkk-jmin+1,dir)&
                      +(velo(2)+velo(3))/2&
                      +(velo(3)-velo(2))/2*tanh((rad-radius_profile(2))/&
                      (2*momentum_thickness(2)))
                endif
                ! Define profile for turbulence on inlet
                u_profile(jjj-imin+1,kkk-jmin+1)&
                    =0.5-0.5*&
                    tanh((rad-radius_profile(2)*1.2)/(0.5*momentum_thickness(1)))
                if (rad < 1.2*radius_profile(2)) then 
                  u_profile(jjj-imin+1,kkk-jmin+1)&
                      =u_profile(jjj-imin+1,kkk-jmin+1)&
                      -0.5+0.5*&
                      tanh((rad-radius_profile(2))/(0.5*momentum_thickness(1)))
                endif
                if (rad < 1.2*radius_profile(1)) then
                  u_profile(jjj-imin+1,kkk-jmin+1)&
                      =u_profile(jjj-imin+1,kkk-jmin+1)&
                      +0.5-0.5*&
                      tanh((rad-radius_profile(1))/(0.5*momentum_thickness(1)))
                endif
                if (rad < 0.8*radius_profile(1)) then
                  u_profile(jjj-imin+1,kkk-jmin+1)&
                      =u_profile(jjj-imin+1,kkk-jmin+1)&
                      -0.5+0.5*&
                      tanh((rad-radius_profile(1)*0.6)/(&
                      0.5*momentum_thickness(1)))  
                endif
                ! Define profile for passive scalar
                if (lpscalar_nolog) then
                  scalar_profile(jjj-imin+1,kkk-jmin+1)&
                      =0.5-0.5&
                      *tanh((rad-radius_profile(2)*1.2)/(&
                      0.5*momentum_thickness(2)))
                  if (rad < radius_mean) then
                    scalar_profile(jjj-imin+1,kkk-jmin+1)&
                        =scalar_profile(jjj-imin+1,kkk-jmin+1)&
                        -0.5+0.5*&
                        tanh((rad-radius_profile(1)*0.9)/(&
                        0.5*momentum_thickness(2)))
                  endif
                endif
              enddo
            enddo
!
          case ('single_jet')
            if (lroot .and. it==1 .and. lfirst) &
                print*,'inlet_profile: single_jet'
            non_zero_transveral_velo=.true.
            velo(1)=u_t
            velo(2)=velo(1)/velocity_ratio
            do jjj=imin,imax
              do kkk=jmin,jmax
                if (dir==1) then
                  rad=sqrt((y(jjj)-jet_center(1))**2+(z(kkk)-jet_center(1))**2)
                elseif (dir==2) then
                  rad=sqrt((x(jjj)-jet_center(1))**2+(z(kkk)-jet_center(1))**2)
                elseif (dir==3) then
                  rad=sqrt((x(jjj)-jet_center(1))**2+(y(kkk)-jet_center(1))**2)
                endif
                  ! Add mean velocity profile
                u_in(jjj-imin+1,kkk-jmin+1,dir)&
                    =u_in(jjj-imin+1,kkk-jmin+1,dir)&
                    +velo(1)*(1-tanh((rad-radius_profile(1))/&
                    momentum_thickness(1)))*0.5+velo(2)
                ! Define profile for turbulence on inlet
                u_profile(jjj-imin+1,kkk-jmin+1)&
                    =0.2*((velo(1))/2-(velo(1))/2&
                    *tanh((rad-radius_profile(1)*1.3)/(2*momentum_thickness(1))))
                ! Define profile for passive scalar
                if (lpscalar_nolog) then
                  scalar_profile(jjj-imin+1,kkk-jmin+1)&
                      =0.5-0.5&
                      *tanh((rad-radius_profile(1)*1.2)/(momentum_thickness(1)))
                endif
              enddo
            enddo
!
          end select
!
        enddo
!
!  Check if we are using a pre-run data file as inlet
!  condition, if not the chosen inlet profile will be used.
!
        if (inlet_from_file) then
          non_zero_transveral_velo=.true.
          if (domain_length == 0) call fatal_error('find_velocity_at_inlet',&
              'domain_length=0. Check that the precisions are the same.')
          round=t*u_t/domain_length
          iround=int(round)
          shift=round-iround
          grid_shift=shift*grid_points
          lowergrid=l1_in+int(grid_shift)
          uppergrid=lowergrid+1
          weight=grid_shift-int(grid_shift)
!
!  Do we want a smooth start
!
          if (smooth_time > 0) then
            smooth=min(t/smooth_time,1.)
          else
            smooth=1.
          endif
!
!  Set the turbulent inlet velocity
!
          if (dir==1) then
            call turbulent_vel_x(u_turb,lowergrid,imin_turb,imax_turb,&
                jmin_turb,jmax_turb,weight,smooth)
          elseif (dir==2) then
            call turbulent_vel_y(u_turb,lowergrid,imin_turb,imax_turb,&
                jmin_turb,jmax_turb,weight,smooth)
          elseif (dir==3) then
            call turbulent_vel_z(u_turb,lowergrid,imin_turb,imax_turb,&
                jmin_turb,jmax_turb,weight,smooth)
          endif
!
!  Add the mean inlet velocity to the turbulent one
!
          do i=1,3
            u_in(:,:,i)=u_in(:,:,i)+u_turb(:,:,i)*u_profile*turb_profile_mag
          enddo

        endif
!
!  Deallocate
!
        deallocate(u_profile)
        deallocate(u_turb)
!
      end subroutine find_velocity_at_inlet
!***********************************************************************
    subroutine find_composition_at_inlet(nchemspec,YYi,YYi_full,dir,imin,imax,jmin,jmax,igrid,jgrid)
!
!  Find composition at inlet.
!
!  2010.07.26/Julien Savre: coded
!
  use  EquationOfState
!
      real, dimension(:,:,:), intent(inout) :: YYi_full
      real, dimension(:), intent(in) :: YYi
      real, dimension(igrid,jgrid) :: zz
      real :: zz1, zz2, init_y1, init_y2
      real :: beta
      integer, intent(in) :: imin,imax,jmin,jmax,igrid,jgrid,nchemspec,dir
!
      real, dimension(igrid) :: dim
      integer :: i_CH4=0, i_O2=0, i_N2=0
      integer :: ichem_CH4=0, ichem_O2=0, ichem_N2=0
      integer :: i,jj,k
      logical :: lO2, lN2, lCH4
!
! Define composition profile at inlet
!
      call keep_compiler_quiet(jmin)
      call keep_compiler_quiet(jmax)
!
        do jj=1,ninit
          select case (zz_profile(jj))
!
          case ('nothing')
            if (lroot .and. it==1 .and. jj == 1 .and. lfirst) &
                print*,'inlet_YY_profile: nothing'
!
          case ('uniform')
            if (lroot .and. it==1 .and. lfirst) &
                print*,'inlet_YY_profile: uniform,'
            do k = 1, nchemspec
              YYi_full(:,:,k)=YYi(k)
            enddo
!
          case ('triple_flame')
            if (lroot .and. it==1 .and. lfirst) &
                print*,'inlet_YY_profile: triple flame, constant gradient'
            zz1=inlet_zz1
            zz2=inlet_zz2
!
            if (dir == 1) then
              dim(1:igrid) = y(imin:imax)
              init_y1 = xyz0(2) + Lxyz(2)/3.
              init_y2 = xyz0(2) + 2.*Lxyz(2)/3
            else if (dir == 2) then
              dim(1:igrid) = x(imin:imax)
              init_y1 = xyz0(1) + Lxyz(1)/3.
              init_y2 = xyz0(1) + 2.*Lxyz(1)/3
            else if (dir == 3) then
              dim(1:igrid) = x(imin:imax)
              init_y1 = xyz0(1) + Lxyz(1)/3.
              init_y2 = xyz0(1) + 2.*Lxyz(1)/3
            else
              call fatal_error('find_composition_at_inlet','No such dir!')
              init_y1 = 0.
              init_y2 = 0.
            endif
            do i = 1, igrid
              if (dim(i) >= init_y1 .and. dim(i) <= init_y2) then
                zz(i,:) = zz2 - (zz2-zz1) * (init_y2 - dim(i)) / &
                                          (init_y2 - init_y1)
              else if (dim(i) <= init_y1) then
                zz(i,:) = zz1
              else if (dim(i) >= init_y2) then
                zz(i,:) = zz2
              endif
            enddo
!
            call find_species_index('O2' ,ichem_O2,i_O2 ,lO2)
            call find_species_index('N2' ,ichem_N2,i_N2,lN2)
            call find_species_index('CH4',ichem_CH4,i_CH4,lCH4)
            do k = 1, nchemspec
              YYi_full(:,:,k)=YYi(k)
            enddo
            beta = YYi(i_N2)/YYi(i_O2)
            do i = 1, igrid
              YYi_full(i,:,i_CH4)=zz(i,:)
              YYi_full(i,:,i_O2)=(1.-YYi_full(i,:,i_CH4)) / (1.+beta)
              YYi_full(i,:,i_N2)=1.-YYi_full(i,:,i_O2)-YYi_full(i,:,i_CH4)
            enddo
!
          end select
!
        enddo
!
      end subroutine find_composition_at_inlet
!***********************************************************************
      subroutine turbulent_vel_x(u_turb,lowergrid,imin,imax,jmin,jmax,weight,smooth)
!
!  Set the turbulent inlet velocity
!
!  2010.01.21/Nils Erland: coded
!
        real, dimension(ny,nz,3), intent(out) :: u_turb
        integer, intent(in) :: lowergrid,imin,imax,jmin,jmax
        real, intent(in) :: weight,smooth
!
        u_turb(:,:,:)&
            =(f_in(lowergrid,imin:imax,jmin:jmax,iux:iuz)*(1-weight)&
            +f_in(lowergrid+1,imin:imax,jmin:jmax,iux:iuz)*weight)*smooth
!
      end subroutine turbulent_vel_x
!***********************************************************************
      subroutine turbulent_vel_y(u_turb,lowergrid,imin,imax,jmin,jmax,weight,smooth)
!
!  Set the turbulent inlet velocity
!
!  2010.01.21/Nils Erland: coded
!
        real, dimension(nx,nz,3), intent(out) :: u_turb
        integer, intent(in) :: lowergrid,imin,imax,jmin,jmax
        real, intent(in) :: weight,smooth
!
        u_turb(:,:,:)&
            =(f_in(imin:imax,lowergrid,jmin:jmax,iux:iuz)*(1-weight)&
            +f_in(imin:imax,lowergrid+1,jmin:jmax,iux:iuz)*weight)*smooth
!
      end subroutine turbulent_vel_y
!***********************************************************************
      subroutine turbulent_vel_z(u_turb,lowergrid,imin,imax,jmin,jmax,weight,smooth)
!
!  Set the turbulent inlet velocity
!
!  2010.01.21/Nils Erland: coded
!
        real, dimension(nx,ny,3), intent(out) :: u_turb
        integer, intent(in) :: lowergrid,imin,imax,jmin,jmax
        real, intent(in) :: weight,smooth
!
        u_turb(:,:,:)&
            =(f_in(imin:imax,jmin:jmax,lowergrid,iux:iuz)*(1-weight)&
            +f_in(imin:imax,jmin:jmax,lowergrid+1,iux:iuz)*weight)*smooth
!
      end subroutine turbulent_vel_z
!***********************************************************************
      subroutine get_thermodynamics(mu1,grad_mu1,gamma_,cs2,cs,rho0,TT,P0,nu,&
          grad_rho,grad_P,grad_T,lll,sgn,direction,fslice,T_t,llinlet,&
          imin,imax,jmin,jmax)
!
!  Find thermodynamical quantities, including density and temperature
!
! 2010.01.21/Nils Erland: coded
!
        Use Chemistry
        use Viscosity
        use EquationOfState, only: cs0, cs20, eoscalc, irho_TT, gamma
!
        integer, intent(in) :: direction, sgn,lll,imin,imax,jmin,jmax
        real, dimension(:,:), intent(out)  :: mu1,cs2,cs,gamma_, grad_mu1
        real, dimension(:,:), intent(out)  :: rho0,TT,P0
        real, dimension(:,:,:), intent(inout)  :: grad_rho,grad_T,grad_P
        real, dimension(:,:,:), intent(in) :: fslice
        real, intent(out) :: nu
        real, intent(inout) :: T_t
        logical :: llinlet
!
        integer :: i,ii,jj
!
!  Find density
!
      if (ldensity_nolog) then
        rho0 = fslice(:,:,irho)
      else
        rho0 = exp(fslice(:,:,ilnrho))
      endif
!
!  Find temperature
!
      if (iTT>0) then
        TT = fslice(:,:,iTT)
      elseif (ilnTT>0) then
        TT = exp(fslice(:,:,ilnTT))
        if (llinlet) T_t=exp(T_t)
      endif
!
!  Get viscoity
!
      call getnu(nu)
!
!  Find mu1, grad_mu1 and gamma
!
        if (ilnTT>0 .or. iTT>0) then
          if (leos_chemistry) then
            call get_mu1_slice(mu1,grad_mu1,lll,sgn,direction)
            call get_gamma_slice(gamma_,direction,lll)
          else
            gamma_=gamma
          endif
        else
          gamma_=1.
        endif
!
!  Set arrays for the speed of sound and for the speed of sound squared (is it
!  really necessarry to have both arrays?)
!  Set prefactors to be used later.
!
      if (leos_idealgas) then
        if (ilnTT>0 .or. iTT>0) then
          do ii=1,imax-imin+1
            do jj=1,jmax-jmin+1
              call eoscalc(irho_TT,rho0(ii,jj),TT(ii,jj),pp=P0(ii,jj),&
                  cs2=cs2(ii,jj))
            enddo
          enddo
          cs=sqrt(cs2)
        else
          cs2=cs20
          cs=cs0
        endif
      elseif (leos_chemistry) then
        call get_cs2_slice(cs2,direction,lll)
        cs=sqrt(cs2)
      else
        print*,"bc_nscbc_prf_x: leos_idealgas=",leos_idealgas,"."
        print*,"bc_nscbc_prf_x: leos_chemistry=",leos_chemistry,"."
        print*,"NSCBC boundary treatment only implemented for ideal gas or"
        print*,"chemistry. Boundary treatment skipped."
        return
      endif
!
!  Find gradient of rho and temperature
!
      if (.not. ldensity_nolog) then
        do i=1,3
          grad_rho(:,:,i)=grad_rho(:,:,i)*rho0
        enddo
      endif
      if (.not. ltemperature_nolog .and. ilnTT>0) then
        do i=1,3
          grad_T(:,:,i)=grad_T(:,:,i)*TT
        enddo
      endif
!
!  Find pressure and the gradient of pressure
!
      do i=1,3
        if (ilnTT>0) then
          if (lchemistry) then
            grad_P(:,:,i)&
                =grad_rho(:,:,i)*TT*Rgas*mu1&
                +grad_T(:,:,i)*rho0*Rgas*mu1&
                +Rgas*grad_mu1*TT*rho0
            P0=rho0*Rgas*mu1*TT
          else
            grad_P(:,:,i)=cs2*(grad_rho(:,:,i)+grad_T(:,:,i)*rho0/TT)/gamma_
          endif
        else
          grad_P(:,:,i)=grad_rho(:,:,i)*cs2
          P0=rho0*cs2
        endif
      enddo
!
    end subroutine get_thermodynamics
!***********************************************************************
      subroutine derivate_boundary(f,sgn,dir,lll,dui_dxj,grad_T,grad_rho)
!
!  Find all derivatives at the boundary
!
!  2010.01.21/Nils Erland: coded
!
        use Deriv, only: der_onesided_4_slice, der_pencil, der2_pencil
!
        integer, intent(in) :: sgn,dir,lll
        real, dimension(:,:,:,:), intent(out) :: dui_dxj
        real, dimension(:,:,:)  , intent(out) :: grad_T,grad_rho
        real, dimension (mx,my,mz,mfarray), intent(in) :: f
!
!        real, dimension (my,mz) :: tmp1,tmp2,tmp3,tmp_lnrho
        real, dimension (:,:), allocatable :: tmp1,tmp2,tmp3,tmp_lnrho
        integer :: i,stat
!
!  Allocate arrays
!
      if (dir == 1) then
        allocate(tmp1(my,mz),STAT=stat)
        if (stat>0) call stop_it("Couldn't allocate memory for tmp1 ")
        allocate(tmp2(my,mz),STAT=stat)
        if (stat>0) call stop_it("Couldn't allocate memory for tmp2 ")
        allocate(tmp3(my,mz),STAT=stat)
        if (stat>0) call stop_it("Couldn't allocate memory for tmp3 ")
        allocate(tmp_lnrho(my,mz),STAT=stat)
        if (stat>0) call stop_it("Couldn't allocate memory for tmp_lnrho ")
      elseif (dir == 2) then
        allocate(tmp1(mx,mz),STAT=stat)
        if (stat>0) call stop_it("Couldn't allocate memory for tmp1 ")
        allocate(tmp2(mx,mz),STAT=stat)
        if (stat>0) call stop_it("Couldn't allocate memory for tmp2 ")
        allocate(tmp3(mx,mz),STAT=stat)
        if (stat>0) call stop_it("Couldn't allocate memory for tmp3 ")
        allocate(tmp_lnrho(mx,mz),STAT=stat)
        if (stat>0) call stop_it("Couldn't allocate memory for tmp_lnrho ")
      elseif (dir == 3) then
        allocate(tmp1(mx,my),STAT=stat)
        if (stat>0) call stop_it("Couldn't allocate memory for tmp1 ")
        allocate(tmp2(mx,my),STAT=stat)
        if (stat>0) call stop_it("Couldn't allocate memory for tmp2 ")
        allocate(tmp3(mx,my),STAT=stat)
        if (stat>0) call stop_it("Couldn't allocate memory for tmp3 ")
        allocate(tmp_lnrho(mx,my),STAT=stat)
        if (stat>0) call stop_it("Couldn't allocate memory for tmp_lnrho ")
      else
        call fatal_error('bc_nscbc_prf_x','No such dir!')
      endif
!
!  Initialize arrays
!
        dui_dxj =0
        grad_rho=0
!
!  Find the derivatives in the direction normal to the boudary by
!  one-sided stencils
!
        call der_onesided_4_slice(f,sgn,ilnrho,grad_rho(:,:,dir),lll,dir)
        call der_onesided_4_slice(f,sgn,iux,dui_dxj(:,:,1,dir),lll,dir)
        call der_onesided_4_slice(f,sgn,iuy,dui_dxj(:,:,2,dir),lll,dir)
        call der_onesided_4_slice(f,sgn,iuz,dui_dxj(:,:,3,dir),lll,dir)
        if (ilnTT>0 .or. iTT>0) then
          call der_onesided_4_slice(f,sgn,ilnTT,grad_T(:,:,dir),lll,dir)
        endif
!
!  Do central differencing in the directions parallell to the boundary
!
        if (dir == 1) then
          if (nygrid /= 1) then
            do i=n1,n2
              call der_pencil(2,f(lll,:,i,iux),tmp1(:,i))
              call der_pencil(2,f(lll,:,i,iuy),tmp2(:,i))
              call der_pencil(2,f(lll,:,i,iuz),tmp3(:,i))
              call der_pencil(2,f(lll,:,i,ilnrho),tmp_lnrho(:,i))
            enddo
            dui_dxj(:,:,1,2)=tmp1(m1:m2,n1:n2)
            dui_dxj(:,:,2,2)=tmp2(m1:m2,n1:n2)
            dui_dxj(:,:,3,2)=tmp3(m1:m2,n1:n2)
            grad_rho(:,:,2)=tmp_lnrho(m1:m2,n1:n2)
          endif
          if (nzgrid /= 1) then
            do i=m1,m2
              call der_pencil(3,f(lll,i,:,iux),tmp1(i,:))
              call der_pencil(3,f(lll,i,:,iuy),tmp2(i,:))
              call der_pencil(3,f(lll,i,:,iuz),tmp3(i,:))
              call der_pencil(3,f(lll,i,:,ilnrho),tmp_lnrho(i,:))
            enddo
            dui_dxj(:,:,1,3)=tmp1(m1:m2,n1:n2)
            dui_dxj(:,:,2,3)=tmp2(m1:m2,n1:n2)
            dui_dxj(:,:,3,3)=tmp3(m1:m2,n1:n2)
            grad_rho(:,:,3)=tmp_lnrho(m1:m2,n1:n2)
          endif
        elseif (dir == 2) then
          if (nxgrid /= 1) then
            do i=n1,n2
              call der_pencil(1,f(:,lll,i,iux),tmp1(:,i))
              call der_pencil(1,f(:,lll,i,iuy),tmp2(:,i))
              call der_pencil(1,f(:,lll,i,iuz),tmp3(:,i))
              call der_pencil(1,f(:,lll,i,ilnrho),tmp_lnrho(:,i))
            enddo
            dui_dxj(:,:,1,1)=tmp1(l1:l2,n1:n2)
            dui_dxj(:,:,2,1)=tmp2(l1:l2,n1:n2)
            dui_dxj(:,:,3,1)=tmp3(l1:l2,n1:n2)
            grad_rho(:,:,1)=tmp_lnrho(l1:l2,n1:n2)
          endif
          if (nzgrid /= 1) then
            do i=l1,l2
              call der_pencil(3,f(i,lll,:,iux),tmp1(i,:))
              call der_pencil(3,f(i,lll,:,iuy),tmp2(i,:))
              call der_pencil(3,f(i,lll,:,iuz),tmp3(i,:))
              call der_pencil(3,f(i,lll,:,ilnrho),tmp_lnrho(i,:))
            enddo
            dui_dxj(:,:,1,3)=tmp1(l1:l2,n1:n2)
            dui_dxj(:,:,2,3)=tmp2(l1:l2,n1:n2)
            dui_dxj(:,:,3,3)=tmp3(l1:l2,n1:n2)
            grad_rho(:,:,3)=tmp_lnrho(l1:l2,n1:n2)
          endif
        elseif (dir == 3) then
          if (nxgrid /= 1) then
            do i=m1,m2
              call der_pencil(1,f(:,i,lll,iux),tmp1(:,i))
              call der_pencil(1,f(:,i,lll,iuy),tmp2(:,i))
              call der_pencil(1,f(:,i,lll,iuz),tmp3(:,i))
              call der_pencil(1,f(:,i,lll,ilnrho),tmp_lnrho(:,i))
            enddo
            dui_dxj(:,:,1,1)=tmp1(l1:l2,m1:m2)
            dui_dxj(:,:,2,1)=tmp2(l1:l2,m1:m2)
            dui_dxj(:,:,3,1)=tmp3(l1:l2,m1:m2)
            grad_rho(:,:,1)=tmp_lnrho(l1:l2,m1:m2)
          endif
          if (nygrid /= 1) then
            do i=l1,l2
              call der_pencil(2,f(i,:,lll,iux),tmp1(i,:))
              call der_pencil(2,f(i,:,lll,iuy),tmp2(i,:))
              call der_pencil(2,f(i,:,lll,iuz),tmp3(i,:))
              call der_pencil(2,f(i,:,lll,ilnrho),tmp_lnrho(i,:))
            enddo
            dui_dxj(:,:,1,2)=tmp1(l1:l2,m1:m2)
            dui_dxj(:,:,2,2)=tmp2(l1:l2,m1:m2)
            dui_dxj(:,:,3,2)=tmp3(l1:l2,m1:m2)
            grad_rho(:,:,2)=tmp_lnrho(l1:l2,m1:m2)
          endif
        endif
!
      end subroutine derivate_boundary
!***********************************************************************
      subroutine transversal_terms(T_1,T_2,T_3,T_4,T_5,rho0,cs,&
          fslice,grad_rho,grad_P,dui_dxj,dir1,dir2,dir3)
!
!  Find the transversal terms.
!  This correspond to the T's in Lodato et al. JCP (2008)
!
!  2010.01.21/Nils Erland: coded
!
        integer,                  intent(in) :: dir1,dir2,dir3
        real, dimension(:,:),     intent(out):: T_1, T_2, T_3, T_4, T_5
        real, dimension(:,:),     intent(in) :: rho0,cs
        real, dimension(:,:,:),   intent(in) :: fslice,grad_rho,grad_P
        real, dimension(:,:,:,:), intent(in) :: dui_dxj
!
!  Calculate the T's
!
        if (.not. notransveral_terms) then
          T_1= rho0*dui_dxj(:,:,dir2,dir2)+fslice(:,:,dir2)*grad_rho(:,:,dir2)&
              +rho0*dui_dxj(:,:,dir3,dir3)+fslice(:,:,dir3)*grad_rho(:,:,dir3)
          T_2= fslice(:,:,dir2)*dui_dxj(:,:,dir1,dir2)&
              +fslice(:,:,dir3)*dui_dxj(:,:,dir1,dir3)
          T_3= fslice(:,:,dir2)*dui_dxj(:,:,dir2,dir2)&
              +fslice(:,:,dir3)*dui_dxj(:,:,dir2,dir3)&
              +grad_P(:,:,dir2)/rho0
          T_4= fslice(:,:,dir2)*dui_dxj(:,:,dir3,dir2)&
              +fslice(:,:,dir3)*dui_dxj(:,:,dir3,dir3)&
              +grad_P(:,:,dir3)/rho0
          T_5= fslice(:,:,dir2)*grad_P(:,:,dir2)&
              +fslice(:,:,dir3)*grad_P(:,:,dir3)&
              +rho0*cs*cs*(dui_dxj(:,:,dir2,dir2)+dui_dxj(:,:,dir3,dir3))
        else
          T_1=0
          T_2=0
          T_3=0
          T_4=0
          T_5=0
        endif
!
      endsubroutine transversal_terms
!***********************************************************************
   subroutine bc_nscbc_subin_x(f,df,topbot,val)
!
!   nscbc case
!   subsonic inflow boundary conditions
!   now it is 2D case (should be corrected for 3D case)
!    ux,uy,T are fixed, drho  is calculated
!
!   16-nov-08/natalia: coded.
!
      use EquationOfState, only: cs0, cs20
      use Deriv, only: der_onesided_4_slice, der_pencil
      use Chemistry
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz) :: mom2
      real, dimension (mx,my,mz,mvar) :: df
      character (len=3) :: topbot
      real, dimension(ny,nz) :: dux_dx, L_1, L_2, L_5, dpp_dx
      real, dimension(my,mz) :: rho0, gamma0, dmom2_dy, TT0
      real, dimension(my,mz) :: cs0_ar,cs20_ar
    !  real, dimension (my,mz) :: cs2x
      real, dimension (mx,my,mz) :: cs2_full, gamma_full, rho_full, pp
      real, dimension(nchemspec) :: YYi
      real, dimension (mcom), optional :: val
      integer :: lll,i, sgn,k
      real :: u_t, T_t
!
      intent(inout) :: f
      intent(out) :: df
!

 !    if (.not.present(val)) call stop_it(&
 !          'bc_nscbc_subin_x: you must specify fbcx)')

      u_t=val(iux)
      T_t=val(ilnTT)
      do k=1,nchemspec
       YYi(k)=val(ichemspec(k))
      enddo

      if (leos_chemistry) then
        call get_cs2_full(cs2_full)
        call get_gamma_full(gamma_full)
      endif
!
      select case (topbot)
      case ('bot')
        lll = l1
        sgn = 1
      case ('top')
        lll = l2
        sgn = -1
      case default
        print*, "bc_nscbc_subin_x: ", topbot, " should be `top' or `bot'"
      endselect
!
      if (leos_chemistry) then
         cs20_ar=cs2_full(lll,:,:)
         cs0_ar=cs2_full(lll,:,:)**0.5
         gamma0=gamma_full(lll,:,:)
       !   TT0=TT_full(lll,:,:)
         if (ltemperature_nolog) then
          TT0=f(lll,:,:,iTT)
         else
          TT0=exp(f(lll,:,:,ilnTT))
         endif

         if (ldensity_nolog) then
          rho_full=f(:,:,:,irho)
         else
          rho_full=exp(f(:,:,:,ilnrho))
         endif

         rho0(:,:) = rho_full(lll,:,:)
         mom2(lll,:,:)=rho0(:,:)*f(lll,:,:,iuy)
         do i=1,my
         do k=1,mz
          if (minval(gamma_full(:,i,k))<=0.) then
           pp(:,i,k)=0.
          else
           pp(:,i,k)=cs2_full(:,i,k)*rho_full(:,i,k)/gamma_full(:,i,k)
          endif
         enddo
         enddo
      else
        print*,"bc_nscbc_subin_x: leos_idealgas=",leos_idealgas,"."
        print*,"NSCBC subsonic inflos is only implemented for "//&
            "the chemistry case."
        print*,"Boundary treatment skipped."
        return
      endif
!
!
      !  call der_onesided_4_slice(f,sgn,ilnrho,dlnrho_dx,lll,1)
        call der_onesided_4_slice(f,sgn,iux,dux_dx,lll,1)
        call der_onesided_4_slice(pp,sgn,dpp_dx,lll,1)
!
        do i=1,mz
          call der_pencil(2,mom2(lll,:,i),dmom2_dy(:,i))
        enddo

        select case (topbot)
        case ('bot')
          L_1 = (f(lll,m1:m2,n1:n2,iux) - cs0_ar(m1:m2,n1:n2))*&
              (dpp_dx - rho0(m1:m2,n1:n2)*cs0_ar(m1:m2,n1:n2)*dux_dx)
          L_5 =L_1-2.*rho0(m1:m2,n1:n2)*cs0_ar(m1:m2,n1:n2)*&
              df(lll,m1:m2,n1:n2,iux)
        case ('top')
          L_5 = (f(lll,m1:m2,n1:n2,iux) + cs0_ar(m1:m2,n1:n2))*&
              (dpp_dx + rho0(m1:m2,n1:n2)*cs0_ar(m1:m2,n1:n2)*dux_dx)
          L_1 = L_5+2.*rho0(m1:m2,n1:n2)*cs0_ar(m1:m2,n1:n2)*&
              df(lll,m1:m2,n1:n2,iux)
        endselect

        if (ltemperature_nolog) then
        L_2 = 0.5*(gamma0(m1:m2,n1:n2)-1.)*(L_5+L_1) &
            +rho0(m1:m2,n1:n2)*cs20_ar(m1:m2,n1:n2) &
            *df(lll,m1:m2,n1:n2,iTT)/TT0(m1:m2,n1:n2)
        else
        L_2 = 0.5*(gamma0(m1:m2,n1:n2)-1.)*(L_5+L_1) &
            +rho0(m1:m2,n1:n2)*cs20_ar(m1:m2,n1:n2)*df(lll,m1:m2,n1:n2,ilnTT)
        endif
        if (ldensity_nolog) then
          df(lll,m1:m2,n1:n2,ilnrho) = -1./cs20_ar(m1:m2,n1:n2)*&
              (L_2+0.5*(L_5 + L_1)) ! -dmom2_dy(m1:m2,n1:n2)
        else
          df(lll,m1:m2,n1:n2,ilnrho) = &
              -1./rho0(m1:m2,n1:n2)/cs20_ar(m1:m2,n1:n2) &
              *(L_2+0.5*(L_5 + L_1)) !&
           !-1./rho0(m1:m2,n1:n2)*dmom2_dy(m1:m2,n1:n2)
        endif

!
! natalia: this subroutine is still under construction
! Please, do not remove this commented part !!!!!
!

!
!  this conditions can be important!
!  check without them
!

   !     do k=1,nchemspec
   !      f(lll,m1:m2,n1:n2,ichemspec(k))=YYi(k)
   !     enddo
   !      f(lll,m1:m2,n1:n2,iux) = u_t
   !     f(lll,m1:m2,n1:n2,ilnTT) = T_t
   !      df(lll,m1:m2,n1:n2,ilnTT)=0.
    endsubroutine bc_nscbc_subin_x
!***********************************************************************
  subroutine bc_nscbc_subin_x_new(f,df,topbot,val)
!
!   nscbc case
!   subsonic inflow boundary conditions
!   now it is 2D case (should be corrected for 3D case)
!    ux,uy,T are fixed, drho  is calculated
!
!   16-nov-08/natalia: coded.
!
      use EquationOfState, only: cs0, cs20
      use Deriv, only: der_onesided_4_slice, der_pencil
      use Chemistry
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz) :: mom2
      real, dimension (mx,my,mz,mvar) :: df
      character (len=3) :: topbot
      real, dimension(ny,nz) :: dux_dx, L_1, L_2,  L_3, L_4, L_5, dpp_dx
      real, dimension(ny,nz,nchemspec) :: L_k
      real, dimension(my,mz) :: rho0, gamma0, TT0 !, dmom2_dy
      real, dimension(my,mz) :: cs0_ar,cs20_ar
    !  real, dimension (my,mz) :: cs2x
      real, dimension (mx,my,mz) :: cs2_full, gamma_full, rho_full, pp
      real, dimension(nchemspec) :: YYi
      real, dimension (mcom), optional :: val
      integer :: lll,i, sgn,k
      real :: u_t, T_t, Mach_num
!
      intent(inout) :: f
      intent(out) :: df
!
      logical :: non_zero_transveral_velo
      integer :: dir, dir2, dir3, igrid, jgrid, imin, imax, jmin, jmax
      real, allocatable, dimension(:,:,:) :: u_in
      real, allocatable, dimension(:,:) :: scalar_profile
!
! reading data from file
!
      if (inlet_from_file) then

        dir=1; dir2=2; dir3=3
        igrid=ny; jgrid=nz
        imin=m1; imax=m2
        jmin=n1; jmax=n2
        non_zero_transveral_velo=.false.
!
        allocate(u_in(igrid,jgrid,3))
!
        call find_velocity_at_inlet(u_in,non_zero_transveral_velo,&
              Lx_in,nx_in,u_t,dir,m1_in,m2_in,n1_in,n2_in,imin,imax,jmin,jmax,&
              igrid,jgrid,scalar_profile)
      endif
!
!    if (.not.present(val)) call stop_it(&
!          'bc_nscbc_subin_x: you must specify fbcx)')
!
      u_t=val(iux)
      T_t=val(ilnTT)
      do k=1,nchemspec
       YYi(k)=val(ichemspec(k))
      enddo

      if (leos_chemistry) then
        call get_cs2_full(cs2_full)
        call get_gamma_full(gamma_full)
      endif
!
      select case (topbot)
      case ('bot')
        lll = l1
        sgn = 1
      case ('top')
        lll = l2
        sgn = -1
      case default
         print*, "bc_nscbc_subin_x: ", topbot, " should be `top' or `bot'"
      endselect
!
      if (leos_chemistry) then
         cs20_ar=cs2_full(lll,:,:)
         cs0_ar=cs2_full(lll,:,:)**0.5
         gamma0=gamma_full(lll,:,:)
       !   TT0=TT_full(lll,:,:)
         if (ltemperature_nolog) then
          TT0=f(lll,:,:,iTT)
         else
          TT0=exp(f(lll,:,:,ilnTT))
         endif

         if (ldensity_nolog) then
          rho_full=f(:,:,:,irho)
         else
          rho_full=exp(f(:,:,:,ilnrho))
         endif

         rho0(:,:) = rho_full(lll,:,:)
         mom2(lll,:,:)=rho0(:,:)*f(lll,:,:,iuy)
         do i=1,my
         do k=1,mz
          if (minval(gamma_full(:,i,k))<=0.) then
           pp(:,i,k)=0.
          else
           pp(:,i,k)=cs2_full(:,i,k)*rho_full(:,i,k)/gamma_full(:,i,k)
          endif
         enddo
         enddo
      else
        print*,"bc_nscbc_subin_x: leos_idealgas=",leos_idealgas,"."
        print*,"NSCBC subsonic inflos is only implemented for "//&
            "the chemistry case."
        print*,"Boundary treatment skipped."
        return
      endif
!
      Mach_num=maxval(f(lll,m1:m2,n1:n2,iux)/cs0_ar(m1:m2,n1:n2))
!
      !  call der_onesided_4_slice(f,sgn,ilnrho,dlnrho_dx,lll,1)
        call der_onesided_4_slice(f,sgn,iux,dux_dx,lll,1)
        call der_onesided_4_slice(pp,sgn,dpp_dx,lll,1)
!
!        do i=1,mz
!          call der_pencil(2,mom2(lll,:,i),dmom2_dy(:,i))
!        enddo

        select case (topbot)
        case ('bot')
          L_1 = (f(lll,m1:m2,n1:n2,iux) - cs0_ar(m1:m2,n1:n2))*&
              (dpp_dx - rho0(m1:m2,n1:n2)*cs0_ar(m1:m2,n1:n2)*dux_dx)
          L_5=nscbc_sigma_in*rho0(m1:m2,n1:n2)*cs20_ar(m1:m2,n1:n2) &
              *(1.-Mach_num**2)/Lxyz(1)*(f(lll,m1:m2,n1:n2,iux)-(u_t+u_in(:,:,1)))


        case ('top')
          L_5 = (f(lll,m1:m2,n1:n2,iux) + cs0_ar(m1:m2,n1:n2))*&
              (dpp_dx + rho0(m1:m2,n1:n2)*cs0_ar(m1:m2,n1:n2)*dux_dx)
          L_1=-nscbc_sigma_in*rho0(m1:m2,n1:n2)*cs20_ar(m1:m2,n1:n2) &
              *(1.-Mach_num**2)/Lxyz(1)*(f(lll,m1:m2,n1:n2,iux)-(u_t+u_in(:,:,1)))
        endselect

        if (ltemperature_nolog) then

        else
         L_2 = -nscbc_sigma_in*pp(lll,m1:m2,n1:n2)*cs0_ar(m1:m2,n1:n2)/Lxyz(1) &
               *(1.-exp(T_t)/TT0(m1:m2,n1:n2))
        endif
        if (ldensity_nolog) then
          df(lll,m1:m2,n1:n2,ilnrho) = -1./cs20_ar(m1:m2,n1:n2)*&
              (L_2+0.5*(L_5 + L_1)) ! -dmom2_dy(m1:m2,n1:n2)
        else

          df(lll,m1:m2,n1:n2,ilnrho) = &
              -1./rho0(m1:m2,n1:n2)/cs20_ar(m1:m2,n1:n2) &
              *(L_2+0.5*(L_5 + L_1)) !&
           !-1./rho0(m1:m2,n1:n2)*dmom2_dy(m1:m2,n1:n2)
        endif

        if (ltemperature_nolog) then
            call stop_it('bc_nscbc_subin_x:ltemperature_nolog case does not work now!')
        else
           df(lll,m1:m2,n1:n2,ilnTT) = &
          -1./(rho0(m1:m2,n1:n2)*cs20_ar(m1:m2,n1:n2))*(-L_2 &
          +0.5*(gamma0(m1:m2,n1:n2)-1.)*(L_5+L_1))
        endif

        L_3=nscbc_sigma_in*cs0_ar(m1:m2,n1:n2)/Lxyz(1) &
                *(f(lll,m1:m2,n1:n2,iuy)-u_in(:,:,2))
        L_4=nscbc_sigma_in*cs0_ar(m1:m2,n1:n2)/Lxyz(1) &
                *(f(lll,m1:m2,n1:n2,iuz)-u_in(:,:,3))
       do k=1,nchemspec
        L_k(:,:,k)=nscbc_sigma_in*cs0_ar(m1:m2,n1:n2)/Lxyz(1) &
                *(f(lll,m1:m2,n1:n2,ichemspec(k))-YYi(k))
       enddo
!
        df(lll,m1:m2,n1:n2,iux) =  &
            -0.5/rho0(m1:m2,n1:n2)/cs0_ar(m1:m2,n1:n2)*(L_5 - L_1)

        df(lll,m1:m2,n1:n2,iuy) = -L_3
        df(lll,m1:m2,n1:n2,iuz) = -L_4
!
        do k=1,nchemspec
         df(lll,m1:m2,n1:n2,ichemspec(k))=-L_k(:,:,k)
         f(lll,m1:m2,n1:n2,ichemspec(k))=YYi(k)
        enddo
         f(lll,m1:m2,n1:n2,ilnTT) = T_t
!
        if (inlet_from_file) then
         f(lll,m1:m2,n1:n2,iux) = u_t + u_in(:,:,1)
         f(lll,m1:m2,n1:n2,iuy) = u_in(:,:,2)
         f(lll,m1:m2,n1:n2,iuz) = u_in(:,:,3)
        else
         f(lll,m1:m2,n1:n2,iux) = u_t
         f(lll,m1:m2,n1:n2,iuy) = 0.
         f(lll,m1:m2,n1:n2,iuz) = 0.
        endif
!
        if (allocated(u_in)) deallocate(u_in)
!
 endsubroutine bc_nscbc_subin_x_new
!***********************************************************************
    subroutine bc_nscbc_nref_subout_x(f,df,topbot,nscbc_sigma_out)

!
!   nscbc case
!   subsonic non-reflecting outflow boundary conditions
!   now it is 2D case (should be corrected for 3D case)
!    drho. dT, dux, duy  are calculated, p_inf can be
!   fixed (if nscbc_sigma <>0)
!
!   16-nov-08/natalia: coded.
!
      use EquationOfState, only: cs0, cs20
      use Deriv, only: der_onesided_4_slice,der_pencil, der2_pencil
      use Chemistry
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      character (len=3) :: topbot
      real, dimension (my,mz) :: rho0,gamma0, TT0
      real, dimension (mx,my,mz) :: mom2,mom3!, rho_ux2, rho_uy2
      real, dimension (mx,my,mz) :: rhoE_p
      real, dimension (mx,my,mz,2) ::  rhoE_pU
      real, dimension (my,mz) ::  dYk_dx,dYk_dy,dYk_dz
      real, dimension (ny,nz) :: drho_prefac, KK

      real, dimension (ny,nz) :: L_1, L_2, L_3, L_4, L_5
      real, dimension (ny,nz) :: M_1, M_2, M_3, M_4, M_5
      real, dimension (ny,nz)  :: N_1, N_2, N_3, N_4, N_5
  !
      real, dimension (my,mz) :: cs0_ar,cs20_ar,dmom2_dy,dmom3_dz
      real, dimension (my,mz,2) :: drhoE_pU!,dYk_dy
   !   real, dimension (mx,my,mz,nchemspec) :: bound_rhs_Y
    !  real, dim(:,nnn)ension (ny,nz) :: bound_rhs_T
      real, dimension (mx,my,mz) :: cs2_full, gamma_full,pp, rho_full
      real, dimension (mx,my,mz,nchemspec) ::  RHS_Y
      real, dimension (nx,ny,nz) :: p_inf

      real, dimension (ny,nz,3,3) :: dui_dxj
      real, dimension (my,mz)     :: tmp1,tmp2,tmp3, tmp_rho,tmp_pp
      real, dimension (ny,nz,3)   :: grad_rho, grad_pp
     ! real, dimension(ny,nz) ::     d2u1_dy2,d2u1_dz2
     ! real, dimension(ny,nz) ::     d2u2_dy2,d2u2_dz2,d2u3_dy2,d2u3_dz2
   !   real, dimension (ny,nz,3) :: dlnT_dxj
      real, dimension (ny,nz) :: T_1_y, T_2_y, T_3_y, T_4_y, T_5_y
      real, dimension (ny,nz) :: T_1_z, T_2_z, T_3_z, T_4_z, T_5_z
!
      integer :: lll, sgn,i,j,k,irho_tmp, nn, nnn, mm, mmm
      real :: Mach_num, nscbc_sigma_out
!
      intent(inout) :: f
      intent(out) :: df
      intent(in) :: nscbc_sigma_out
       logical :: lcorner_y=.false.,lcorner_z=.false.

      if (leos_chemistry) then
        call get_cs2_full(cs2_full)
        call get_gamma_full(gamma_full)
      endif

      call get_RHS_Y_full(RHS_Y)
!
      select case (topbot)
      case ('bot')
        lll = l1; sgn = 1
        if (leos_chemistry) then
          p_inf(1,:,:)=p_infty
        endif
      case ('top')
        lll = l2; sgn = -1
        if (leos_chemistry) then
          p_inf(nx,:,:)=p_infty
        endif
      case default
        print*, "bc_nscbc_subin_x: ", topbot, " should be `top' or `bot'"
      endselect

      if (leos_chemistry) then
         cs20_ar=cs2_full(lll,:,:)
         cs0_ar=cs2_full(lll,:,:)**0.5
         gamma0=gamma_full(lll,:,:)
         if (ltemperature_nolog) then
          TT0=f(lll,:,:,iTT)
         else
          TT0=exp(f(lll,:,:,ilnTT))
         endif

        if (ldensity_nolog) then
          rho_full = f(:,:,:,irho)
          rho0(:,:) = f(lll,:,:,irho)
          drho_prefac=-1./cs20_ar(m1:m2,n1:n2)
          irho_tmp=irho
        !  call stop_it('bc_nscbc_nref_subout_x: NSCBC works now only for lnrho')
        else
          rho_full = exp(f(:,:,:,ilnrho))
          rho0(:,:) = rho_full(lll,:,:)
          drho_prefac=-1./rho0(m1:m2,n1:n2)/cs20_ar(m1:m2,n1:n2)
          irho_tmp=ilnrho
        endif

         do i=1,my
         do k=1,mz
          if (minval(gamma_full(:,i,k))<=0.) then
           pp(:,i,k)=0.
          else
           pp(:,i,k)=cs2_full(:,i,k)*rho_full(:,i,k)/gamma_full(:,i,k)
          endif
         enddo
         enddo

         mom2(lll,:,:)=rho0(:,:)*f(lll,:,:,iuy)
         mom3(lll,:,:)=rho0(:,:)*f(lll,:,:,iuz)
         rhoE_p(lll,:,:)=0.5*rho_full(lll,:,:) &
            *(f(lll,:,:,iux)**2+f(lll,:,:,iuy)**2+f(lll,:,:,iuz)**2) &
             +gamma_full(lll,:,:)/(gamma_full(lll,:,:)-1)*pp(lll,:,:)
         rhoE_pU(lll,:,:,1)=rhoE_p(lll,:,:)*f(lll,:,:,iuy)
         rhoE_pU(lll,:,:,2)=rhoE_p(lll,:,:)*f(lll,:,:,iuz)


         Mach_num=maxval(f(lll,m1:m2,n1:n2,iux)/cs0_ar(m1:m2,n1:n2))
         KK=nscbc_sigma_out*(1.-Mach_num*Mach_num)*cs0_ar(m1:m2,n1:n2)/Lxyz(1)
      else
        print*,"bc_nscbc_subin_x: leos_idealgas=",leos_idealgas,"."
        print*,"NSCBC subsonic inflos is only implemented "//&
            "for the chemistry case."
        print*,"Boundary treatment skipped."
        return
      endif

        call der_onesided_4_slice(rho_full,sgn,grad_rho(:,:,1),lll,1)
        call der_onesided_4_slice(pp,sgn,grad_pp(:,:,1),lll,1)
        call der_onesided_4_slice(f,sgn,iux,dui_dxj(:,:,1,1),lll,1)
        call der_onesided_4_slice(f,sgn,iuy,dui_dxj(:,:,2,1),lll,1)
        call der_onesided_4_slice(f,sgn,iuz,dui_dxj(:,:,3,1),lll,1)
    !    call der_onesided_4_slice(f,sgn,ilnTT,dlnT_dxj(:,:,1),lll,1)


       if (nygrid /= 1) then
         do i=n1,n2
          call der_pencil(2,f(lll,:,i,iux),tmp1(:,i))
          call der_pencil(2,f(lll,:,i,iuy),tmp2(:,i))
          call der_pencil(2,f(lll,:,i,iuz),tmp3(:,i))
          call der_pencil(2,rho_full(lll,:,i),tmp_rho(:,i))
          call der_pencil(2,pp(lll,:,i),tmp_pp(:,i))
          call der_pencil(2,mom2(lll,:,i),dmom2_dy(:,i))
          call der_pencil(2,rhoE_pU(lll,:,i,1),drhoE_pU(:,i,1))
         enddo
      !  do i=n1,n2
      !    call der2_pencil(2,f(lll,:,i,iux),tmpy)
      !    d2u1_dy2(:,i-n1+1)=tmpy(:)
      !    call der2_pencil(2,f(lll,:,i,iuy),tmpy)
      !    d2u2_dy2(:,i-n1+1)=tmpy(:)
      !    call der2_pencil(2,f(lll,:,i,iuz),tmpy)
      !    d2u3_dy2(:,i-n1+1)=tmpy(:)
      !    enddo
      else
        tmp3=0
        tmp2=0
        tmp1=0
        tmp_rho=0
        tmp_pp=0
        dmom2_dy=0
        drhoE_pU(:,:,1)=0
     !   d2u1_dy2=0
     !   d2u2_dy2=0
     !   d2u3_dy2=0
      endif
        dui_dxj(:,:,1,2)=tmp1(m1:m2,n1:n2)
        dui_dxj(:,:,2,2)=tmp2(m1:m2,n1:n2)
        dui_dxj(:,:,3,2)=tmp3(m1:m2,n1:n2)
        grad_rho(:,:,2)=tmp_rho(m1:m2,n1:n2)
        grad_pp(:,:,2)=tmp_pp(m1:m2,n1:n2)

      if (nzgrid /= 1) then
         do j=m1,m2
          call der_pencil(3,f(lll,j,:,iux),tmp1(j,:))
          call der_pencil(3,f(lll,j,:,iuy),tmp2(j,:))
          call der_pencil(3,f(lll,j,:,iuz),tmp3(j,:))
          call der_pencil(3,rho_full(lll,j,:),tmp_rho(j,:))
          call der_pencil(3,pp(lll,j,:),tmp_pp(j,:))
          call der_pencil(3,mom3(lll,j,:),dmom3_dz(j,:))
          call der_pencil(3,rhoE_pU(lll,j,:,2),drhoE_pU(j,:,2))
         enddo
      else
        tmp3=0
        tmp2=0
        tmp1=0
        tmp_rho=0
        tmp_pp=0
        dmom3_dz=0
        drhoE_pU(:,:,2)=0
      endif
        dui_dxj(:,:,1,3)=tmp1(m1:m2,n1:n2)
        dui_dxj(:,:,2,3)=tmp2(m1:m2,n1:n2)
        dui_dxj(:,:,3,3)=tmp3(m1:m2,n1:n2)
        grad_rho(:,:,3)=tmp_rho(m1:m2,n1:n2)
        grad_pp(:,:,3)=tmp_pp(m1:m2,n1:n2)


!
      select case (topbot)
      case ('bot')
        L_5=KK*(cs20_ar(m1:m2,n1:n2)/gamma0(m1:m2,n1:n2)*&
            rho0(m1:m2,n1:n2)-p_inf(1,1:ny,1:nz))
        L_1 = (f(lll,m1:m2,n1:n2,iux) - cs0_ar(m1:m2,n1:n2))*&
           (grad_pp(:,:,1)-rho0(m1:m2,n1:n2)*cs0_ar(m1:m2,n1:n2)*dui_dxj(:,:,1,1))
      case ('top')
        L_1=KK*(cs20_ar(m1:m2,n1:n2)/gamma0(m1:m2,n1:n2)*&
            rho0(m1:m2,n1:n2)-p_inf(nx,1:ny,1:nz))
        L_5 = (f(lll,m1:m2,n1:n2,iux) + cs0_ar(m1:m2,n1:n2))*&
            (grad_pp(:,:,1)+ rho0(m1:m2,n1:n2)*cs0_ar(m1:m2,n1:n2)*dui_dxj(:,:,1,1))
     endselect
!
       L_2 = f(lll,m1:m2,n1:n2,iux)*(cs20_ar(m1:m2,n1:n2)*grad_rho(:,:,1)-grad_pp(:,:,1))
       L_3 = f(lll,m1:m2,n1:n2,iux)*dui_dxj(:,:,2,1)
       L_4 = f(lll,m1:m2,n1:n2,iux)*dui_dxj(:,:,3,1)
!
       df(lll,m1:m2,n1:n2,irho_tmp) = &
         drho_prefac*(L_2+0.5*(L_5 + L_1))
       df(lll,m1:m2,n1:n2,iux) =  &
         -1./(2.*rho0(m1:m2,n1:n2)*cs0_ar(m1:m2,n1:n2))*(L_5 - L_1)
       df(lll,m1:m2,n1:n2,iuy) = -L_3
       df(lll,m1:m2,n1:n2,iuz) = -L_4
       if (ltemperature_nolog) then
        df(lll,m1:m2,n1:n2,ilnTT) = &
          -1./(rho0(m1:m2,n1:n2)*cs20_ar(m1:m2,n1:n2))*(-L_2 &
          +0.5*(gamma0(m1:m2,n1:n2)-1.)*(L_5+L_1))*TT0(m1:m2,n1:n2) !&
      ! + RHS_T_full(lll,m1:m2,n1:n2)
       else
        df(lll,m1:m2,n1:n2,ilnTT) = &
          -1./(rho0(m1:m2,n1:n2)*cs20_ar(m1:m2,n1:n2))*(-L_2 &
          +0.5*(gamma0(m1:m2,n1:n2)-1.)*(L_5+L_1))
       endif
!NNNNNNNNNN

        if ((nygrid /= 1) .or.(nzgrid /= 1))  then

          T_1_y(:,:)=-dmom2_dy(m1:m2,n1:n2)/rho0(m1:m2,n1:n2)
          T_1_z(:,:)=-dmom3_dz(m1:m2,n1:n2)/rho0(m1:m2,n1:n2)
          T_2_y(:,:)=-f(lll,m1:m2,n1:n2,iuy)*dui_dxj(:,:,1,2)
          T_2_z(:,:)=-f(lll,m1:m2,n1:n2,iuz)*dui_dxj(:,:,1,3) !&
                !   -nu_full(lll,m1:m2,n1:n2)*(d2u1_dy2+d2u1_dz2)
          T_3_y(:,:)=-f(lll,m1:m2,n1:n2,iuy)*dui_dxj(:,:,2,2) &
                     -grad_pp(:,:,2)/rho0(m1:m2,n1:n2)
         T_3_z(:,:)=-f(lll,m1:m2,n1:n2,iuz)*dui_dxj(:,:,2,3)
                !   -nu_full(lll,m1:m2,n1:n2)*(d2u2_dy2+d2u2_dz2)
          T_4_y(:,:)=-f(lll,m1:m2,n1:n2,iuy)*dui_dxj(:,:,3,2)
         T_4_z(:,:)=-f(lll,m1:m2,n1:n2,iuz)*dui_dxj(:,:,3,3) &
                     -grad_pp(:,:,3)/rho0(m1:m2,n1:n2)
          T_5_y(:,:)= drho_prefac(:,:)*(gamma0(m1:m2,n1:n2)-1.) &
                      *gamma0(m1:m2,n1:n2)*drhoE_pU(m1:m2,n1:n2,1)
          T_5_z(:,:)= drho_prefac(:,:)*(gamma0(m1:m2,n1:n2)-1.) &
                      *gamma0(m1:m2,n1:n2)*(drhoE_pU(m1:m2,n1:n2,2))

         df(lll,m1:m2,n1:n2,irho_tmp) = df(lll,m1:m2,n1:n2,irho_tmp) + T_1_y + T_1_z
         df(lll,m1:m2,n1:n2,iux)      = df(lll,m1:m2,n1:n2,iux)      + T_2_y + T_2_z
         df(lll,m1:m2,n1:n2,iuy)      = df(lll,m1:m2,n1:n2,iuy)      + T_3_y + T_3_z
         df(lll,m1:m2,n1:n2,iuz)      = df(lll,m1:m2,n1:n2,iuz)      + T_4_y + T_4_z
        if (ltemperature_nolog) then
         df(lll,m1:m2,n1:n2,ilnTT)    = df(lll,m1:m2,n1:n2,ilnTT)   &
                                + (T_5_y + T_5_z)*TT0(m1:m2,n1:n2)
        else
         df(lll,m1:m2,n1:n2,ilnTT)    = df(lll,m1:m2,n1:n2,ilnTT)    + T_5_y + T_5_z
        endif

 !!!
!!!  Corner points are described in bc_nscbc_nref_subout_y and bc_nscbc_nref_subout_z
!!!


       endif
!NATALIA

      if (nygrid /= 1) then
           M_1(:,:)=(f(lll,m1:m2,n1:n2,iuy) - cs0_ar(m1:m2,n1:n2))&
             *(grad_pp(:,:,2)-rho0(m1:m2,n1:n2)*cs0_ar(m1:m2,n1:n2)*dui_dxj(:,:,2,2))
           M_2(:,:)=f(lll,m1:m2,n1:n2,iuy)*(cs20_ar(m1:m2,n1:n2) &
                  *grad_rho(:,:,2)-grad_pp(:,:,2))
           M_3(:,:)=f(lll,m1:m2,n1:n2,iuy)*dui_dxj(:,:,1,2)
           M_4(:,:)=f(lll,m1:m2,n1:n2,iuy)*dui_dxj(:,:,3,2)
           M_5(:,:)=(f(lll,m1:m2,n1:n2,iuy) + cs0_ar(m1:m2,n1:n2))&
            *(grad_pp(:,:,2)+ rho0(m1:m2,n1:n2)*cs0_ar(m1:m2,n1:n2)*dui_dxj(:,:,2,2))


           if (y(m2)==Lxyz(2)) then
            M_1(ny,:)=KK(ny,:)*(cs20_ar(m2,n1:n2)/gamma0(m2,n1:n2) &
               *rho0(m2,n1:n2)-p_inf(lll-3,ny,:))
           endif
           if  (y(m2)==xyz0(2)) then
            M_5(1,:)=KK(1,:)*(cs20_ar(m1,n1:n2)/gamma0(m1,n1:n2)*&
             rho0(m1,n1:n2)-p_inf(lll-3,1,:))
           endif

        else
        M_1=0; M_2=0; M_3=0; M_4=0; M_5=0
       endif
!
       if (nzgrid /= 1)  then
         N_1(:,:)=(f(lll,m1:m2,n1:n2,iuz) - cs0_ar(m1:m2,n1:n2))&
          *(grad_pp(:,:,3)-rho0(m1:m2,n1:n2)*cs0_ar(m1:m2,n1:n2)*dui_dxj(:,:,3,3))

         N_2(:,:)=f(lll,m1:m2,n1:n2,iuz)*dui_dxj(:,:,1,3)

         N_3(:,:)=f(lll,m1:m2,n1:n2,iuz)*dui_dxj(:,:,2,3)

         N_4(:,:)=f(lll,m1:m2,n1:n2,iuz) &
            *(cs20_ar(m1:m2,n1:n2)*grad_rho(:,:,3)-grad_pp(:,:,3))

         N_5(:,:)=(f(lll,m1:m2,n1:n2,iuz) + cs0_ar(m1:m2,n1:n2))&
            *(grad_pp(:,:,3)+ rho0(m1:m2,n1:n2)*cs0_ar(m1:m2,n1:n2)*dui_dxj(:,:,3,3))

         if (z(n2)==Lxyz(3)) then
          N_1(:,nz)=KK(:,nz)*(cs20_ar(m1:m2,n2)/gamma0(m1:m2,n2)*&
            rho0(m1:m2,n2)-p_inf(lll-3,:,nz))
         endif
         if (z(n1)==xyz0(3)) then
          N_5(:,1)=KK(:,1)*(cs20_ar(m1:m2,n1)/gamma0(m1:m2,n1)*&
             rho0(m1:m2,n1)-p_inf(lll-3,:,1))
         endif


       else
        N_1=0; N_2=0; N_3=0; N_4=0; N_5=0
       endif


         do i=1,2
         if (i==1) then
            mm=m1; mmm=1
            if (y(m1)==xyz0(2)) lcorner_y=.true.
         elseif (i==2) then
            mm=m2; mmm=ny
            if (y(m2)==Lxyz(2)) lcorner_y=.true.
         endif
         if (lcorner_y)  then

!          df(lll,mm,n1:n2,irho_tmp) = T_1_z(mmm,:)&
!          + drho_prefac(mmm,:)*(M_2(mmm,:)+0.5*(M_5(mmm,:) + M_1(mmm,:))) &
!          + drho_prefac(mmm,:)*(L_2(mmm,:)+0.5*(L_5(mmm,:) + L_1(mmm,:)))

!         df(lll,mm,n1:n2,iux) =  T_2_z(mmm,:) -M_3(mmm,:)- 1./&
!             (2.*rho0(lll,mm)*cs0_ar(lll,mm))*(L_5(mmm,:) - L_1(mmm,:))

!         df(lll,mm,n1:n2,iuy) =  T_3_z(mmm,:) - 1./&
!           (2.*rho0(lll,mm)*cs0_ar(lll,mm))*(M_5(mmm,:) - M_1(mmm,:))-L_3(mmm,:)
!         df(lll,mm,n1:n2,iuz) = T_4_z(mmm,:)- M_4(mmm,:) -L_4(mmm,:)
!         df(lll,mm,n1:n2,ilnTT) =T_5_z(mmm,:) &
!           +drho_prefac(mmm,:)*(-M_2(mmm,:) &
!           +0.5*(gamma0(lll,mm)-1.)*(M_5(mmm,:)+M_1(mmm,:))) &
!           +drho_prefac(mmm,:)*(-L_2(mmm,:) &
!           +0.5*(gamma0(lll,mm)-1.)*(L_5(mmm,:)+L_1(mmm,:)))
!          lcorner_y=.false.


        endif
        enddo

        do i=1,2
         if (i==1) then
          nn=n1; nnn=1
          if (z(n1)==xyz0(3)) lcorner_z=.true.
         elseif (i==2) then
          nn=n2; nnn=nz
          if (z(n2)==Lxyz(3))  lcorner_z=.true.
         endif

         if  (lcorner_z) then
!          df(lll,m1:m2,nn,irho_tmp) = T_1_y(:,nnn) &
!          + drho_prefac(:,nnn)*(L_2(:,nnn)+0.5*(L_5(:,nnn) + L_1(:,nnn))) &
!          + drho_prefac(:,nnn)*(N_4(:,nnn)+0.5*(N_5(:,nnn) + N_1(:,nnn)))
!          df(lll,m1:m2,nn,iuy) = T_3_y(:,nnn) -L_3(:,nnn)-N_3(:,nnn)
!          df(lll,m1:m2,nn,iux) = T_2_y(:,nnn) - 1./&
!           (2.*rho0(m1:m2,nn)*cs0_ar(m1:m2,nn))*(L_5(:,nnn) &
!                - L_1(:,nnn))-N_2(:,nnn)
!          df(lll,m1:m2,nn,iuz) =  T_4_y(:,nnn) - L_4(:,nnn)  - N_4(:,nnn)
!          df(lll,m1:m2,nn,ilnTT) = T_5_y(:,nnn) &
!           +drho_prefac(:,nnn)*(-L_2(:,nnn) &
!           +0.5*(gamma0(lll,nn)-1.)*(L_5(:,nnn)+L_1(:,nnn))) &
!           +drho_prefac(:,nnn)*(-N_4(:,nnn) &
!           +0.5*(gamma0(lll,nn)-1.)*(N_5(:,nnn)+N_1(:,nnn)))
!           lcorner_z=.false.
         endif
        enddo


        do i=1,2
         if (i==1) then
          nn=n1; nnn=1
          if (z(n1)==xyz0(3)) lcorner_z=.true.
         elseif (i==2) then
          nn=n2; nnn=nz
          if (z(n2)==Lxyz(3))  lcorner_z=.true.
         endif

        do j=1,2
         if (j==1) then
           mm=m1; mmm=1
           if (y(m1)==xyz0(2)) lcorner_y=.true.
         elseif (j==2) then
           mm=m2; mmm=ny
           if (y(m2)==Lxyz(2)) lcorner_y=.true.
         endif

       if ((lcorner_y)  .and. (lcorner_z)) then
!        df(lll,mm,nn,irho_tmp) = &
!            drho_prefac(mmm,nnn)*(L_2(mmm,nnn)+0.5*(L_5(mmm,nnn) + L_1(mmm,nnn))) &
!          + drho_prefac(mmm,nnn)*(M_2(mmm,nnn)+0.5*(M_5(mmm,nnn) + M_1(mmm,nnn))) &
!          + drho_prefac(mmm,nnn)*(N_4(mmm,nnn)+0.5*(N_5(mmm,nnn) + N_1(mmm,nnn)))
!        df(lll,mm,nn,iux) =  -1./&
!          (2.*rho0(mm,nn)*cs0_ar(mm,nn))*(L_5(mmm,nnn) - L_1(mmm,nnn)) &
!          -M_2(mmm,nnn)-N_2(mmm,nnn)
!        df(lll,mm,nn,iuy) =  - L_3(mmm,nnn) - 1./&
!          (2.*rho0(mm,nn)*cs0_ar(mm,nn))*(M_5(mmm,nnn) - M_1(mmm,nnn))-N_3(mmm,nnn)
!        df(lll,mm,nn,iuz) =  - L_4(mmm,nnn) - M_4(mmm,nnn)  - 1./&
!          (2.*rho0(mm,nn)*cs0_ar(mm,nn))*(N_5(mmm,nnn) - N_1(mmm,nnn))
!        df(lll,mm,nn,ilnTT) = drho_prefac(mmm,nnn)*(-L_2(mmm,nnn) &
!          +0.5*(gamma0(mm,nn)-1.)*(L_5(mmm,nnn)+L_1(mmm,nnn)))  &
!          +drho_prefac(mmm,nnn)*(-M_2(mmm,nnn) &
!          +0.5*(gamma0(mm,nn)-1.)*(M_5(mmm,nnn)+M_1(mmm,nnn))) &
!           +drho_prefac(mmm,nnn)*(-N_4(mmm,nnn) &
!          +0.5*(gamma0(mm,nn)-1.)*(N_5(mmm,nnn)+N_1(mmm,nnn)))
!           lcorner_y=.false.
!           lcorner_z=.false.
       endif
      enddo
      enddo


!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (nchemspec>1) then
       do k=1,nchemspec
          call der_onesided_4_slice(f,sgn,ichemspec(k),dYk_dx(m1:m2,n1:n2),lll,1)
         do i=n1,n2
          call der_pencil(2,f(lll,:,i,ichemspec(k)),dYk_dy(:,i))
         enddo
         do i=m1,m2
           call der_pencil(3,f(lll,i,:,ichemspec(k)),dYk_dz(i,:))
         enddo
!NMNMNMN
           df(lll,m1:m2,n1:n2,ichemspec(k))=&
              -f(lll,m1:m2,n1:n2,iux)*dYk_dx(m1:m2,n1:n2) &
              -f(lll,m1:m2,n1:n2,iuy)*dYk_dy(m1:m2,n1:n2) &
              -f(lll,m1:m2,n1:n2,iuz)*dYk_dz(m1:m2,n1:n2) &
              + RHS_Y(lll,m1:m2,n1:n2,k)

       ! if (lfilter) then
         do i=m1,m2
         do j=n1,n2
           if ((f(lll,i,j,ichemspec(k))+df(lll,i,j,ichemspec(k))*dt)<-1e-25 ) then
             df(lll,i,j,ichemspec(k))=-1e-25*dt
           endif
           if ((f(lll,i,j,ichemspec(k))+df(lll,i,j,ichemspec(k))*dt)>1.) then
             f(lll,i,j,ichemspec(k))=1.*dt
           endif
         enddo
         enddo
       ! endif
       enddo
      endif
!
    endsubroutine bc_nscbc_nref_subout_x
!***********************************************************************
    subroutine bc_nscbc_nref_subout_y(f,df,topbot,nscbc_sigma_out)
!
!   nscbc case
!   subsonic non-reflecting outflow boundary conditions
!
!   16-jun-09/natalia: coded.
!
      use EquationOfState, only: cs0, cs20
      use Deriv, only: der_onesided_4_slice,der_pencil
      use Chemistry
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      character (len=3) :: topbot
      real, dimension (mx,mz) :: rho0,gamma0
      real, dimension (mx,my,mz) :: mom1, mom3, pp!, rho_ux2, rho_uy2
      real, dimension (mx,mz) ::  dYk_dy,dYk_dx,dYk_dz
      real, dimension (nx,nz) :: drho_prefac, KK
      real, dimension (nx,nz) :: M_1, M_2, M_3, M_4, M_5
      real, dimension (nx,nz)  :: L_1, L_2, L_3, L_4, L_5
      real, dimension (nx,nz)  :: N_1, N_2, N_3, N_4, N_5
      real, dimension (nx,nz) :: T_1_x, T_2_x, T_3_x, T_4_x, T_5_x
      real, dimension (nx,nz) :: T_1_z, T_2_z, T_3_z, T_4_z, T_5_z
      real, dimension (mx,mz) :: cs0_ar,cs20_ar,dmom1_dx,dmom3_dz
      real, dimension (mx,my,mz) :: rhoE_p
      real, dimension (mx,my,mz,2) ::   rhoE_pU
      real, dimension (mx,mz,2) :: drhoE_pU
      real, dimension (mx,my,mz,nchemspec) ::  RHS_Y
      !real, dimension (mx,mz) ::  !dYk_dy
  !    real, dimension (nx,nz,nchemspec) :: bound_rhs_Y
    !  real, dimension (nx,nz) :: bound_rhs_T
      real, dimension (mx,my,mz) :: cs2_full, gamma_full, rho_full
      real, dimension (nx,ny,nz) :: p_inf
      real, dimension (nx,nz,3,3) :: dui_dxj
      real, dimension (mx,mz)     :: tmp11,tmp21,tmp31,tmp13,tmp23,tmp33
      real, dimension (mx,mz)     :: tmp1_rho,tmp3_rho,tmp1_pp,tmp3_pp
      real, dimension (nx,nz,3)   :: grad_rho, grad_pp

!
      integer :: mmm, sgn,i,j,k, irho_tmp, nn, nnn, ll, lll
      real :: Mach_num,nscbc_sigma_out
      logical :: lcorner_x=.false.,lcorner_z=.false.
!
      intent(inout) :: f
      intent(out) :: df
      intent(in) :: nscbc_sigma_out
!
      if (leos_chemistry) then
        call get_cs2_full(cs2_full)
        call get_gamma_full(gamma_full)
      endif

      call get_RHS_Y_full(RHS_Y)
!

      select case (topbot)
      case ('bot')
        mmm = m1; sgn = 1
        if (leos_chemistry) then
          p_inf(:,1,:)=p_infty
        endif
      case ('top')
        mmm = m2; sgn = -1
        if (leos_chemistry) then
          p_inf(:,ny,:)=p_infty
        endif
      case default
        print*, "bc_nscbc_subout_y: ", topbot, " should be `top' or `bot'"
      endselect


      if (leos_chemistry) then
         cs20_ar=cs2_full(:,mmm,:)
         cs0_ar=cs2_full(:,mmm,:)**0.5
         gamma0=gamma_full(:,mmm,:)

        if (ldensity_nolog) then
          rho_full = f(:,:,:,irho)
          rho0(:,:) = f(:,mmm,:,irho)
          drho_prefac=-1./cs20_ar(l1:l2,n1:n2)
          irho_tmp=irho

         call stop_it('bc_nscbc_nref_subout_y: NSCBC works now only for lnrho')

        else
          rho_full = exp(f(:,:,:,ilnrho))
          rho0(:,:) = rho_full(:,mmm,:)
          drho_prefac=-1./rho0(l1:l2,n1:n2)/cs20_ar(l1:l2,n1:n2)
          irho_tmp=ilnrho
        endif

         do i=1,my
         do k=1,mz
          if (minval(gamma_full(:,i,k))<=0.) then
           pp(:,i,k)=0.
          else
           pp(:,i,k)=cs2_full(:,i,k)*rho_full(:,i,k)/gamma_full(:,i,k)
          endif
         enddo
         enddo

         rhoE_p(:,mmm,:)=0.5*rho_full(:,mmm,:) &
            *(f(:,mmm,:,iux)**2+f(:,mmm,:,iuy)**2+f(:,mmm,:,iuz)**2) &
             +gamma_full(:,mmm,:)/(gamma_full(:,mmm,:)-1)*pp(:,mmm,:)
         rhoE_pU(:,mmm,:,1)=rhoE_p(:,mmm,:)*f(:,mmm,:,iux)
         rhoE_pU(:,mmm,:,2)=rhoE_p(:,mmm,:)*f(:,mmm,:,iuz)
         mom1(:,mmm,:)=rho0(:,:)*f(:,mmm,:,iux)
         mom3(:,mmm,:)=rho0(:,:)*f(:,mmm,:,iuz)

         Mach_num=maxval(f(l1:l2,mmm,n1:n2,iuy)/cs0_ar(l1:l2,n1:n2))
         KK=nscbc_sigma_out*(1.-Mach_num*Mach_num)*cs0_ar(l1:l2,n1:n2)/Lxyz(2)


      else
        print*,"bc_nscbc_subin_y: leos_idealgas=",leos_idealgas,"."
        print*,"NSCBC subsonic inflos is only implemented "//&
            "for the chemistry case."
        print*,"Boundary treatment skipped."
        return
      endif

      call der_onesided_4_slice(rho_full,sgn,grad_rho(:,:,2),mmm,2)
      call der_onesided_4_slice(pp,sgn,grad_pp(:,:,2),mmm,2)
      call der_onesided_4_slice(f,sgn,iux,dui_dxj(:,:,1,2),mmm,2)
      call der_onesided_4_slice(f,sgn,iuy,dui_dxj(:,:,2,2),mmm,2)
      call der_onesided_4_slice(f,sgn,iuz,dui_dxj(:,:,3,2),mmm,2)

      if (nxgrid /= 1) then
         do i=n1,n2
          call der_pencil(1,f(:,mmm,i,iux),tmp11(:,i))
          call der_pencil(1,f(:,mmm,i,iuy),tmp21(:,i))
          call der_pencil(1,f(:,mmm,i,iuz),tmp31(:,i))
          call der_pencil(1,rho_full(:,mmm,i),tmp1_rho(:,i))
          call der_pencil(1,pp(:,mmm,i),tmp1_pp(:,i))
          call der_pencil(1,mom1(:,mmm,i),dmom1_dx(:,i))
          call der_pencil(1,rhoE_pU(:,mmm,i,1),drhoE_pU(:,i,1))
         enddo

      else
        tmp31=0
        tmp21=0
        tmp11=0
        tmp1_rho=0
        tmp1_pp=0
        dmom1_dx=0
        drhoE_pU(:,:,1)=0
      endif
        dui_dxj(:,:,1,1)=tmp11(l1:l2,n1:n2)
        dui_dxj(:,:,2,1)=tmp21(l1:l2,n1:n2)
        dui_dxj(:,:,3,1)=tmp31(l1:l2,n1:n2)
        grad_rho(:,:,1)=tmp1_rho(l1:l2,n1:n2)
        grad_pp(:,:,1)=tmp1_pp(l1:l2,n1:n2)

      if (nzgrid /= 1) then
         do j=l1,l2
          call der_pencil(3,f(j,mmm,:,iux),tmp13(j,:))
          call der_pencil(3,f(j,mmm,:,iuy),tmp23(j,:))
          call der_pencil(3,f(j,mmm,:,iuz),tmp33(j,:))
          call der_pencil(3,rho_full(j,mmm,:),tmp3_rho(j,:))
          call der_pencil(3,pp(j,mmm,:),tmp3_pp(j,:))
          call der_pencil(3,mom3(j,mmm,:),dmom3_dz(j,:))
          call der_pencil(3,rhoE_pU(j,mmm,:,2),drhoE_pU(j,:,2))

         enddo

      else
        tmp33=0
        tmp23=0
        tmp13=0
        tmp3_rho=0
        tmp3_pp=0
        dmom3_dz=0
        drhoE_pU(:,:,2)=0
      endif
        dui_dxj(:,:,1,3)=tmp13(l1:l2,n1:n2)
        dui_dxj(:,:,2,3)=tmp23(l1:l2,n1:n2)
        dui_dxj(:,:,3,3)=tmp33(l1:l2,n1:n2)
        grad_rho(:,:,3)=tmp3_rho(l1:l2,n1:n2)
        grad_pp(:,:,3)=tmp3_pp(l1:l2,n1:n2)


!
      select case (topbot)
      case ('bot')
        M_5=KK*(cs20_ar(l1:l2,n1:n2)/gamma0(l1:l2,n1:n2)*&
            rho0(l1:l2,n1:n2)-p_inf(1:nx,1,1:nz))
        M_1 = (f(l1:l2,mmm,n1:n2,iuy) - cs0_ar(l1:l2,n1:n2))*&
            (grad_pp(:,:,2)- rho0(l1:l2,n1:n2)*cs0_ar(l1:l2,n1:n2)*dui_dxj(:,:,2,2))
      case ('top')
        M_1=KK*(cs20_ar(l1:l2,n1:n2)/gamma0(l1:l2,n1:n2)*&
            rho0(l1:l2,n1:n2)-p_inf(1:nx,ny,1:nz))
        M_5 = (f(l1:l2,mmm,n1:n2,iuy) + cs0_ar(l1:l2,n1:n2))*&
            (grad_pp(:,:,2)+ rho0(l1:l2,n1:n2)*cs0_ar(l1:l2,n1:n2)*dui_dxj(:,:,2,2))
      endselect
!
      M_2 = f(l1:l2,mmm,n1:n2,iuy)*(cs20_ar(l1:l2,n1:n2)*grad_rho(:,:,2)-grad_pp(:,:,2))
      M_3 = f(l1:l2,mmm,n1:n2,iuy)*dui_dxj(:,:,1,2)
      M_4 = f(l1:l2,mmm,n1:n2,iuy)*dui_dxj(:,:,3,2)

      df(l1:l2,mmm,n1:n2,irho_tmp) = drho_prefac*(M_2+0.5*(M_5 + M_1))
      df(l1:l2,mmm,n1:n2,iux) = -M_3
      df(l1:l2,mmm,n1:n2,iuy) = -1./&
          (2.*rho0(l1:l2,n1:n2)*cs0_ar(l1:l2,n1:n2))*(M_5 - M_1) !&
      df(l1:l2,mmm,n1:n2,iuz) = -M_4
      df(l1:l2,mmm,n1:n2,ilnTT) = drho_prefac(:,:)*(-M_2 &
          +0.5*(gamma0(l1:l2,n1:n2)-1.)*(M_5+M_1))


!NNNNNNNN
       if ((nxgrid /= 1) .or. (nzgrid /= 1)) then

          T_1_x(:,:)=-dmom1_dx(l1:l2,n1:n2)/rho0(l1:l2,n1:n2)
          T_1_z(:,:)=-dmom3_dz(l1:l2,n1:n2)/rho0(l1:l2,n1:n2)
          T_3_x(:,:)=-f(l1:l2,mmm,n1:n2,iux)*dui_dxj(:,:,2,1)
          T_3_z(:,:)=-f(l1:l2,mmm,n1:n2,iuz)*dui_dxj(:,:,2,3)
          T_2_x(:,:)=-f(l1:l2,mmm,n1:n2,iux)*dui_dxj(:,:,1,1) &
                    -grad_pp(:,:,1)/rho0(l1:l2,n1:n2)
          T_2_z(:,:)=-f(l1:l2,mmm,n1:n2,iuz)*dui_dxj(:,:,1,3)
          T_4_x(:,:)=-f(l1:l2,mmm,n1:n2,iux)*dui_dxj(:,:,3,1)
          T_4_z(:,:)=-f(l1:l2,mmm,n1:n2,iuz)*dui_dxj(:,:,3,3) &
                   -grad_pp(:,:,3)/rho0(l1:l2,n1:n2)
          T_5_x(:,:)=+drho_prefac(:,:)*(gamma0(l1:l2,n1:n2)-1.)*gamma0(l1:l2,n1:n2) &
                   *(drhoE_pU(l1:l2,n1:n2,1))
          T_5_z(:,:)=+drho_prefac(:,:)*(gamma0(l1:l2,n1:n2)-1.)*gamma0(l1:l2,n1:n2) &
                   *(drhoE_pU(l1:l2,n1:n2,2))


        df(l1:l2,mmm,n1:n2,irho_tmp) = df(l1:l2,mmm,n1:n2,irho_tmp) + T_1_x+T_1_z
        df(l1:l2,mmm,n1:n2,iux) =      df(l1:l2,mmm,n1:n2,iux)      + T_2_x+T_2_z
        df(l1:l2,mmm,n1:n2,iuy) =      df(l1:l2,mmm,n1:n2,iuy)      + T_3_x+T_3_z
        df(l1:l2,mmm,n1:n2,iuz) =      df(l1:l2,mmm,n1:n2,iuz)      + T_4_x+T_4_z
        df(l1:l2,mmm,n1:n2,ilnTT) =    df(l1:l2,mmm,n1:n2,ilnTT)    + T_5_x+T_5_z

   !    if ((nxgrid /= 1) .and. (nzgrid /= 1)) then
    ! if ((nxgrid /= 1) .and. (nzgrid == 1)) then
!!!
!!! Corner points
!!!
!!!

        if (nxgrid /= 1) then
           L_1(:,:)=(f(l1:l2,mmm,n1:n2,iux) - cs0_ar(l1:l2,n1:n2))&
             *(grad_pp(:,:,1)-rho0(l1:l2,n1:n2)*cs0_ar(l1:l2,n1:n2)*dui_dxj(:,:,1,1))
           L_2(:,:)=f(l1:l2,mmm,n1:n2,iux)*(cs20_ar(l1:l2,n1:n2) &
                  *grad_rho(:,:,1)-grad_pp(:,:,1))
           L_3(:,:)=f(l1:l2,mmm,n1:n2,iux)*dui_dxj(:,:,2,1)
           L_4(:,:)=f(l1:l2,mmm,n1:n2,iux)*dui_dxj(:,:,3,1)
           L_5(:,:)=(f(l1:l2,mmm,n1:n2,iux) + cs0_ar(l1:l2,n1:n2))&
            *(grad_pp(:,:,1)+ rho0(l1:l2,n1:n2)*cs0_ar(l1:l2,n1:n2)*dui_dxj(:,:,1,1))


           if (x(l2)==Lxyz(1)) then
            L_1(nx,:)=KK(nx,:)*(cs20_ar(l2,n1:n2)/gamma0(l2,n1:n2) &
               *rho0(l2,n1:n2)-p_inf(nx,mmm-3,:))
           endif
           if  (x(l1)==xyz0(1)) then
            L_5(1,:)=KK(1,:)*(cs20_ar(l1,n1:n2)/gamma0(l1,n1:n2)*&
             rho0(l1,n1:n2)-p_inf(1,mmm-3,:))
           endif

        else
        L_1=0; L_2=0; L_3=0; L_4=0; L_5=0
       endif
!
       if (nzgrid /= 1)  then
         N_1(:,:)=(f(l1:l2,mmm,n1:n2,iuz) - cs0_ar(l1:l2,n1:n2))&
          *(grad_pp(:,:,3)-rho0(l1:l2,n1:n2)*cs0_ar(l1:l2,n1:n2)*dui_dxj(:,:,3,3))

         N_2(:,:)=f(l1:l2,mmm,n1:n2,iuz)*dui_dxj(:,:,1,3)

         N_3(:,:)=f(l1:l2,mmm,n1:n2,iuz)*dui_dxj(:,:,2,3)

         N_4(:,:)=f(l1:l2,mmm,n1:n2,iuz) &
            *(cs20_ar(l1:l2,n1:n2)*grad_rho(:,:,3)-grad_pp(:,:,3))

         N_5(:,:)=(f(l1:l2,mmm,n1:n2,iuz) + cs0_ar(l1:l2,n1:n2))&
            *(grad_pp(:,:,3)+ rho0(l1:l2,n1:n2)*cs0_ar(l1:l2,n1:n2)*dui_dxj(:,:,3,3))

         if (z(n2)==Lxyz(3)) then
          N_1(:,nz)=KK(:,nz)*(cs20_ar(l1:l2,n2)/gamma0(l1:l2,n2)*&
            rho0(l1:l2,n2)-p_inf(:,mmm-3,nz))
         endif
         if (z(n1)==xyz0(3)) then
          N_5(:,1)=KK(:,1)*(cs20_ar(l1:l2,n1)/gamma0(l1:l2,n1)*&
             rho0(l1:l2,n1)-p_inf(:,mmm-3,1))
         endif


       else
        N_1=0; N_2=0; N_3=0; N_4=0; N_5=0
       endif


         do i=1,2
         if (i==1) then
            ll=l1; lll=1
            if (x(l1)==xyz0(1)) lcorner_x=.true.
         elseif (i==2) then
            ll=l2; lll=nx
            if (x(l2)==Lxyz(1)) lcorner_x=.true.
         endif
         if (lcorner_x)  then
          df(ll,mmm,n1:n2,irho_tmp) = &
            drho_prefac(lll,:)*(L_2(lll,:)+0.5*(L_5(lll,:) + L_1(lll,:))) &
          + drho_prefac(lll,:)*(M_2(lll,:) &
          +0.5*(M_5(lll,:) + M_1(lll,:))) +T_1_z(lll,:)

          df(ll,mmm,n1:n2,iux) =  -1./&
           (2.*rho0(ll,n1:n2)*cs0_ar(ll,n1:n2))*(L_5(lll,:) - L_1(lll,:)) &
           -M_3(lll,:)+T_2_z(lll,:)
          df(ll,mmm,n1:n2,iuy) =  - L_3(lll,:) - 1./&
           (2.*rho0(ll,n1:n2)*cs0_ar(ll,n1:n2))*(M_5(lll,:) - M_1(lll,:)) &
           + T_3_z(lll,:)
          df(ll,mmm,n1:n2,iuz) =  - L_4(lll,:) - M_4(lll,:)  +T_4_z(lll,:)

          df(ll,mmm,n1:n2,ilnTT) = drho_prefac(lll,:)*(-L_2(lll,:) &
           +0.5*(gamma0(ll,n1:n2)-1.)*(L_5(lll,:)+L_1(lll,:)))  &
           +drho_prefac(lll,:)*(-M_2(lll,:) &
           +0.5*(gamma0(ll,n1:n2)-1.)*(M_5(lll,:)+M_1(lll,:))) +T_5_z(lll,:)
          lcorner_x=.false.
        endif
        enddo

        do i=1,2
         if (i==1) then
          nn=n1; nnn=1
          if (z(n1)==xyz0(3)) lcorner_z=.true.
         elseif (i==2) then
          nn=n2; nnn=nz
          if (z(n2)==Lxyz(3))  lcorner_z=.true.
         endif

         if  (lcorner_z) then
          df(l1:l2,mmm,nn,irho_tmp) = T_1_x(:,nnn) &
          + drho_prefac(:,nnn)*(M_2(:,nnn)+0.5*(M_5(:,nnn) + M_1(:,nnn))) &
          + drho_prefac(:,nnn)*(N_4(:,nnn)+0.5*(N_5(:,nnn) + N_1(:,nnn)))
          df(l1:l2,mmm,nn,iux) = T_2_x(:,nnn) -M_3(:,nnn)-N_2(:,nnn)
          df(l1:l2,mmm,nn,iuy) = T_3_x(:,nnn) - 1./&
           (2.*rho0(l1:l2,nn)*cs0_ar(l1:l2,nn))*(M_5(:,nnn) &
                - M_1(:,nnn))-N_3(:,nnn)
          df(l1:l2,mmm,nn,iuz) =  T_4_x(:,nnn) - M_4(:,nnn)  - 1./&
           (2.*rho0(l1:l2,nn)*cs0_ar(l1:l2,nn))*(N_5(:,nnn) - N_1(:,nnn))
          df(l1:l2,mmm,nn,ilnTT) = T_5_x(:,nnn) &
           +drho_prefac(:,nnn)*(-M_2(:,nnn) &
           +0.5*(gamma0(l1:l2,nn)-1.)*(M_5(:,nnn)+M_1(:,nnn))) &
           +drho_prefac(:,nnn)*(-N_4(:,nnn) &
           +0.5*(gamma0(l1:l2,nn)-1.)*(N_5(:,nnn)+N_1(:,nnn)))
           lcorner_z=.false.
         endif
        enddo


        do i=1,2
         if (i==1) then
          nn=n1; nnn=1
          if (z(n1)==xyz0(3)) lcorner_z=.true.
         elseif (i==2) then
          nn=n2; nnn=nz
          if (z(n2)==Lxyz(3))  lcorner_z=.true.
         endif

        do j=1,2
         if (j==1) then
           ll=l1; lll=1
           if (x(l1)==xyz0(1)) lcorner_x=.true.
         elseif (j==2) then
           ll=l2; lll=nx
           if (x(l2)==Lxyz(1)) lcorner_x=.true.
         endif

!NNNNNNNN

    !   if ((x(l1)==xyz0(1))  .or. (x(l2)==Lxyz(1)))  lcorner_x=.true.
    !   if ((z(n1)==xyz0(3))  .or. (z(n2)==Lxyz(3)))  lcorner_z=.true.


       if ((lcorner_x)  .and. (lcorner_z)) then
        df(ll,mmm,nn,irho_tmp) = &
          drho_prefac(lll,nnn)*(L_2(lll,nnn)+0.5*(L_5(lll,nnn) + L_1(lll,nnn))) &
          + drho_prefac(lll,nnn)*(M_2(lll,nnn)+0.5*(M_5(lll,nnn) + M_1(lll,nnn))) &
          + drho_prefac(lll,nnn)*(N_4(lll,nnn)+0.5*(N_5(lll,nnn) + N_1(lll,nnn)))
        df(ll,mmm,nn,iux) =  -1./&
          (2.*rho0(ll,nn)*cs0_ar(ll,nn))*(L_5(lll,nnn) - L_1(lll,nnn)) &
          -M_3(lll,nnn)-N_2(lll,nnn)
        df(ll,mmm,nn,iuy) =  - L_3(lll,nnn) - 1./&
          (2.*rho0(ll,nn)*cs0_ar(ll,nn))*(M_5(lll,nnn) - M_1(lll,nnn))-N_3(lll,nnn)
        df(ll,mmm,nn,iuz) =  - L_4(lll,nnn) - M_4(lll,nnn)  - 1./&
          (2.*rho0(ll,nn)*cs0_ar(ll,nn))*(N_5(lll,nnn) - N_1(lll,nnn))
        df(ll,mmm,nn,ilnTT) = drho_prefac(lll,nnn)*(-L_2(lll,nnn) &
          +0.5*(gamma0(ll,nn)-1.)*(L_5(lll,nnn)+L_1(lll,nnn)))  &
          +drho_prefac(lll,nnn)*(-M_2(lll,nnn) &
          +0.5*(gamma0(ll,nn)-1.)*(M_5(lll,nnn)+M_1(lll,nnn))) &
           +drho_prefac(lll,nnn)*(-N_4(lll,nnn) &
          +0.5*(gamma0(ll,nn)-1.)*(N_5(lll,nnn)+N_1(lll,nnn)))
           lcorner_x=.false.
           lcorner_z=.false.
       endif
      enddo
      enddo

    endif



      if (nchemspec>1) then
       do k=1,nchemspec
          call der_onesided_4_slice(f,sgn,ichemspec(k),dYk_dy(l1:l2,n1:n2),mmm,2)
         do i=n1,n2
          call der_pencil(1,f(:,mmm,i,ichemspec(k)),dYk_dx(:,i))
         enddo
         do i=l1,l2
          call der_pencil(3,f(i,mmm,:,ichemspec(k)),dYk_dz(i,:))
         enddo


          df(l1:l2,mmm,n1:n2,ichemspec(k))=&
              -f(l1:l2,mmm,n1:n2,iux)*dYk_dx(l1:l2,n1:n2) &
              -f(l1:l2,mmm,n1:n2,iuy)*dYk_dy(l1:l2,n1:n2) &
              -f(l1:l2,mmm,n1:n2,iuz)*dYk_dz(l1:l2,n1:n2) &
              +RHS_Y(l1:l2,mmm,n1:n2,k)

       ! if (lfilter) then
         do i=l1,l2
         do j=n1,n2
           if ((f(i,mmm,j,ichemspec(k))+df(i,mmm,j,ichemspec(k))*dt)<-1e-25 ) then
             df(i,mmm,j,ichemspec(k))=-1e-25*dt
           endif
           if ((f(i,mmm,j,ichemspec(k))+df(i,mmm,j,ichemspec(k))*dt)>1.) then
             f(i,mmm,j,ichemspec(k))=1.*dt
           endif
         enddo
         enddo
       ! endif
       enddo
      endif
!
    endsubroutine bc_nscbc_nref_subout_y
!***********************************************************************
   subroutine bc_nscbc_nref_subout_z(f,df,topbot,nscbc_sigma_out)
!
!   nscbc case
!   subsonic non-reflecting outflow boundary conditions
!
!   16-jun-09/natalia: coded.
!
      use EquationOfState, only: cs0, cs20
      use Deriv, only: der_onesided_4_slice,der_pencil
      use Chemistry
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      character (len=3) :: topbot
      real, dimension (mx,my) :: rho0,gamma0
      real, dimension (mx,my,mz) :: mom1,mom2,pp
      real, dimension (mx,my,mz) :: rhoE_p
      real, dimension (mx,my,mz,nchemspec) ::  RHS_Y
      real, dimension (mx,my,mz,2) ::   rhoE_pU=0.
      real, dimension (mx,my,2) :: drhoE_pU
      real, dimension (mx,my) ::  dYk_dz=0.,dYk_dx=0.,dYk_dy=0.
      real, dimension (nx,ny) :: drho_prefac, KK, N_1, N_2, N_3,N_4, N_5
      real, dimension (nx,ny)  :: M_1, M_2, M_3, M_4, M_5
      real, dimension (nx,ny)  :: L_1, L_2, L_3, L_4, L_5

      real, dimension (mx,my) :: cs0_ar,cs20_ar,dmom1_dx,dmom2_dy
      real, dimension (mx,my,mz) :: cs2_full, gamma_full, rho_full
      real, dimension (nx,ny,nz) :: p_inf
      real, dimension(nx,ny,3,3) :: dui_dxj=0.
       real, dimension (mx,my)     :: tmp11,tmp21,tmp31,tmp12,tmp22,tmp32
      real, dimension (mx,my)     :: tmp1_rho,tmp2_rho,tmp1_pp,tmp2_pp
      real, dimension (nx,ny,3)   :: grad_rho=0., grad_pp=0.
      real, dimension (nx,ny) :: T_1_x, T_2_x, T_3_x, T_4_x, T_5_x
      real, dimension (nx,ny) :: T_1_y, T_2_y, T_3_y, T_4_y, T_5_y
!
      integer :: nnn, sgn,i,j,k, mmm,mm,lll,ll, irho_tmp
      real :: Mach_num,nscbc_sigma_out

      logical :: lcorner_x=.false.,lcorner_y=.false.
!
      intent(inout) :: f
      intent(out) :: df
      intent(in) :: nscbc_sigma_out
!
      if (leos_chemistry) then
        call get_cs2_full(cs2_full)
        call get_gamma_full(gamma_full)
      endif

        call get_RHS_Y_full(RHS_Y)
!
      select case (topbot)
      case ('bot')
        nnn = n1; sgn = 1
        if (leos_chemistry) then
          p_inf(:,:,1)=p_infty
        endif
      case ('top')
        nnn = n2; sgn = -1
        if (leos_chemistry) then
          p_inf(:,:,nz)=p_infty
        endif
      case default
        print*, "bc_nscbc_subout_z: ", topbot, " should be `top' or `bot'"
      endselect

        if (leos_chemistry) then
         cs20_ar=cs2_full(:,:,nnn)
         cs0_ar=cs2_full(:,:,nnn)**0.5
         gamma0=gamma_full(:,:,nnn)

        if (ldensity_nolog) then
          rho_full = f(:,:,:,irho)
          rho0(:,:) = f(:,:,nnn,irho)
          drho_prefac=-1./cs20_ar(l1:l2,m1:m2)
           irho_tmp=irho
        else
          rho_full = exp(f(:,:,:,ilnrho))
          rho0(:,:) = rho_full(:,:,nnn)
          drho_prefac=-1./rho0(l1:l2,m1:m2)/cs20_ar(l1:l2,m1:m2)
           irho_tmp=ilnrho
        endif

         do i=1,my
         do k=1,mz
          if (minval(gamma_full(:,i,k))<=0.) then
           pp(:,i,k)=0.
          else
           pp(:,i,k)=cs2_full(:,i,k)*rho_full(:,i,k)/gamma_full(:,i,k)
          endif
         enddo
         enddo

         mom1(:,:,nnn)=rho0(:,:)*f(:,:,nnn,iux)
         mom2(:,:,nnn)=rho0(:,:)*f(:,:,nnn,iuy)
         rhoE_p(:,:,nnn)=0.5*rho_full(:,:,nnn) &
            *(f(:,:,nnn,iux)**2+f(:,:,nnn,iuy)**2+f(:,:,nnn,iuz)**2) &
             +gamma_full(:,:,nnn)/(gamma_full(:,:,nnn)-1)*pp(:,:,nnn)
         rhoE_pU(:,:,nnn,1)=rhoE_p(:,:,nnn)*f(:,:,nnn,iux)
         rhoE_pU(:,:,nnn,2)=rhoE_p(:,:,nnn)*f(:,:,nnn,iuy)

        Mach_num=maxval(f(l1:l2,m1:m2,nnn,iuz)/cs0_ar(l1:l2,m1:m2))
        KK=nscbc_sigma_out*(1.-Mach_num*Mach_num)*cs0_ar(l1:l2,m1:m2)/Lxyz(3)

      else
        print*,"bc_nscbc_subin_y: leos_idealgas=",leos_idealgas,"."
       print*,"NSCBC subsonic inflos is only implemented "//&
           "for the chemistry case."
        print*,"Boundary treatment skipped."
        return
      endif

      call der_onesided_4_slice(rho_full,sgn,grad_rho(:,:,3),nnn,3)
      call der_onesided_4_slice(pp,sgn,grad_pp(:,:,3),nnn,3)
      call der_onesided_4_slice(f,sgn,iux,dui_dxj(:,:,1,3),nnn,3)
      call der_onesided_4_slice(f,sgn,iuy,dui_dxj(:,:,2,3),nnn,3)
      call der_onesided_4_slice(f,sgn,iuz,dui_dxj(:,:,3,3),nnn,3)

      if (nxgrid /= 1) then

          do i=m1,m2
          call der_pencil(1,f(:,i,nnn,iux),tmp11(:,i))
          call der_pencil(1,f(:,i,nnn,iuy),tmp21(:,i))
          call der_pencil(1,f(:,i,nnn,iuz),tmp31(:,i))
          call der_pencil(1,rho_full(:,i,nnn),tmp1_rho(:,i))
          call der_pencil(1,pp(:,i,nnn),tmp1_pp(:,i))
          call der_pencil(1,mom1(:,i,nnn),dmom1_dx(:,i))
          call der_pencil(1,rhoE_pU(:,i,nnn,1),drhoE_pU(:,i,1))
         enddo


      else
        tmp31=0
        tmp21=0
        tmp11=0
        tmp1_rho=0
        tmp1_pp=0
        dmom1_dx=0
        drhoE_pU(:,:,1)=0
      endif
        dui_dxj(:,:,1,1)=tmp11(l1:l2,m1:m2)
        dui_dxj(:,:,2,1)=tmp21(l1:l2,m1:m2)
        dui_dxj(:,:,3,1)=tmp31(l1:l2,m1:m2)
        grad_rho(:,:,1)=tmp1_rho(l1:l2,m1:m2)
        grad_pp(:,:,1)=tmp1_pp(l1:l2,m1:m2)

      if (nygrid /= 1) then

          do i=l1,l2
          call der_pencil(2,f(i,:,nnn,iux),tmp12(i,:))
          call der_pencil(2,f(i,:,nnn,iuy),tmp22(i,:))
          call der_pencil(2,f(i,:,nnn,iuz),tmp32(i,:))
          call der_pencil(2,rho_full(i,:,nnn),tmp2_rho(i,:))
          call der_pencil(2,pp(i,:,nnn),tmp2_pp(i,:))
          call der_pencil(2,mom2(i,:,nnn),dmom2_dy(i,:))
          call der_pencil(2,rhoE_pU(i,:,nnn,2),drhoE_pU(i,:,2))
         enddo
      else
        tmp32=0
        tmp22=0
        tmp12=0
        tmp2_rho=0
        tmp2_pp=0
        dmom2_dy=0
        drhoE_pU(:,:,2)=0
      endif
        dui_dxj(:,:,1,2)=tmp12(l1:l2,m1:m2)
        dui_dxj(:,:,2,2)=tmp22(l1:l2,m1:m2)
        dui_dxj(:,:,3,2)=tmp32(l1:l2,m1:m2)
        grad_rho(:,:,2)=tmp2_rho(l1:l2,m1:m2)
        grad_pp(:,:,2)=tmp2_pp(l1:l2,m1:m2)

!
      select case (topbot)
      case ('bot')
        N_5=KK*(cs20_ar(l1:l2,m1:m2)/gamma0(l1:l2,m1:m2)*&
            rho0(l1:l2,m1:m2)-p_inf(1:nx,1:ny,1))
        N_1 = (f(l1:l2,m1:m2,nnn,iuz) - cs0_ar(l1:l2,m1:m2))*&
            (grad_pp(:,:,3)- rho0(l1:l2,m1:m2)*cs0_ar(l1:l2,m1:m2)*dui_dxj(:,:,3,3))
      case ('top')
        N_1=KK*(cs20_ar(l1:l2,m1:m2)/gamma0(l1:l2,m1:m2)*&
          rho0(l1:l2,m1:m2)-p_inf(1:nx,1:ny,nz))
        N_5 = (f(l1:l2,m1:m2,nnn,iuz) + cs0_ar(l1:l2,m1:m2))*&
            (grad_pp(:,:,3)+ rho0(l1:l2,m1:m2)*cs0_ar(l1:l2,m1:m2)*dui_dxj(:,:,3,3))
      endselect
!
        N_2 = f(l1:l2,m1:m2,nnn,iuz)  &
           *(cs20_ar(l1:l2,m1:m2)*grad_rho(:,:,3)-grad_pp(:,:,3))
        N_3 = f(l1:l2,m1:m2,nnn,iuz)*dui_dxj(:,:,1,3)
        N_4 = f(l1:l2,m1:m2,nnn,iuz)*dui_dxj(:,:,2,3)
!

        df(l1:l2,m1:m2,nnn,irho_tmp) = &
            drho_prefac*(N_2+0.5*(N_5 + N_1))
        df(l1:l2,m1:m2,nnn,iuz) = -1./&
          (2.*rho0(l1:l2,m1:m2)*cs0_ar(l1:l2,m1:m2))*(N_5 - N_1)
        df(l1:l2,m1:m2,nnn,iux) = -N_3
        df(l1:l2,m1:m2,nnn,iuy) = -N_4
        df(l1:l2,m1:m2,nnn,ilnTT) = -1./&
         (rho0(l1:l2,m1:m2)*cs20_ar(l1:l2,m1:m2))*(-N_2 &
        +0.5*(gamma0(l1:l2,m1:m2)-1.)*(N_5+N_1))
!NNNNNNNNN

       if ((nxgrid /= 1) .or. (nygrid /= 1)) then

         T_1_x(:,:)=-dmom1_dx(l1:l2,m1:m2)/rho0(l1:l2,m1:m2)
         T_1_y(:,:)=-dmom2_dy(l1:l2,m1:m2)/rho0(l1:l2,m1:m2)
         T_2_x(:,:)=-f(l1:l2,m1:m2,nnn,iux)*dui_dxj(:,:,1,1) &
                  -grad_pp(:,:,1)/rho0(l1:l2,m1:m2)
         T_2_y(:,:)=-f(l1:l2,m1:m2,nnn,iuy)*dui_dxj(:,:,1,2)

         T_3_x(:,:)=-f(l1:l2,m1:m2,nnn,iux)*dui_dxj(:,:,2,1)
         T_3_y(:,:)=-f(l1:l2,m1:m2,nnn,iuy)*dui_dxj(:,:,2,2)  &
                   -grad_pp(:,:,2)/rho0(l1:l2,m1:m2)
         T_4_x(:,:)=-f(l1:l2,m1:m2,nnn,iux)*dui_dxj(:,:,3,1)
         T_4_y(:,:)=-f(l1:l2,m1:m2,nnn,iuy)*dui_dxj(:,:,3,2)
         T_5_x(:,:)=+drho_prefac(:,:)*(gamma0(l1:l2,m1:m2)-1.)*gamma0(l1:l2,m1:m2) &
                   *(drhoE_pU(l1:l2,m1:m2,1))
         T_5_y(:,:)=+drho_prefac(:,:)*(gamma0(l1:l2,m1:m2)-1.)*gamma0(l1:l2,m1:m2) &
                   *(drhoE_pU(l1:l2,m1:m2,2))

         df(l1:l2,m1:m2,nnn,irho_tmp) = df(l1:l2,m1:m2,nnn,irho_tmp) + T_1_x+T_1_y
         df(l1:l2,m1:m2,nnn,iux) =      df(l1:l2,m1:m2,nnn,iux)      + T_2_x+T_2_y
         df(l1:l2,m1:m2,nnn,iuy) =      df(l1:l2,m1:m2,nnn,iuy)      + T_3_x+T_3_y
         df(l1:l2,m1:m2,nnn,iuz) =      df(l1:l2,m1:m2,nnn,iuz)      + T_4_x+T_4_y
         df(l1:l2,m1:m2,nnn,ilnTT) =    df(l1:l2,m1:m2,nnn,ilnTT)    + T_5_x+T_5_y
! if ((nxgrid /= 1) .and. (nygrid /= 1)) then
 ! if ((nxgrid /= 1) .and. (nygrid == 1)) then

!
! Corner points
!

       if (nxgrid /= 1) then
         L_1=(f(l1:l2,m1:m2,nnn,iux) - cs0_ar(l1:l2,m1:m2))&
            *(grad_pp(:,:,1)-rho0(l1:l2,m1:m2)*cs0_ar(l1:l2,m1:m2)*dui_dxj(:,:,1,1))
          L_2=f(l1:l2,m1:m2,nnn,iux)*(cs20_ar(l1:l2,m1:m2) &
                  *grad_rho(:,:,1)-grad_pp(:,:,1))
          L_3=f(l1:l2,m1:m2,nnn,iux)*dui_dxj(:,:,2,1)
          L_4=f(l1:l2,m1:m2,nnn,iux)*dui_dxj(:,:,3,1)
          L_5=(f(l1:l2,m1:m2,nnn,iux) + cs0_ar(l1:l2,m1:m2))&
            *(grad_pp(:,:,1)+ rho0(l1:l2,m1:m2)*cs0_ar(l1:l2,m1:m2)*dui_dxj(:,:,1,1))

         if (x(l1)==xyz0(1)) then
          L_1(nx,:)=KK(nx,:)*(cs20_ar(l2,m1:m2)/gamma0(l2,m1:m2)*&
           rho0(l2,m1:m2)-p_inf(nx,:,nnn-3))
         endif
         if (x(l2)==Lxyz(1)) then
          L_5(1,:)=KK(1,:)*(cs20_ar(l1,m1:m2)/gamma0(l1,m1:m2)*&
            rho0(l1,m1:m2)-p_inf(1,:,nnn-3))
         endif

       else
        L_1=0; L_2=0; L_3=0; L_4=0; L_5=0
       endif

       if (nygrid /= 1)  then
           M_1=(f(l1:l2,m1:m2,nnn,iuy) - cs0_ar(l1:l2,m1:m2))&
            *(grad_pp(:,:,2)-rho0(l1:l2,m1:m2)*cs0_ar(l1:l2,m1:m2)*dui_dxj(:,:,2,2))
          M_3=f(l1:l2,m1:m2,nnn,iuy)*dui_dxj(:,:,1,2)
          M_2=f(l1:l2,m1:m2,nnn,iuy) &
           *(cs20_ar(l1:l2,m1:m2)*grad_rho(:,:,2)-grad_pp(:,:,2))
          M_4=f(l1:l2,m1:m2,nnn,iuy)*dui_dxj(:,:,3,2)
          M_5=(f(l1:l2,m1:m2,nnn,iuy) + cs0_ar(l1:l2,m1:m2))&
            *(grad_pp(:,:,2)+ rho0(l1:l2,m1:m2)*cs0_ar(l1:l2,m1:m2)*dui_dxj(:,:,2,2))

         if (y(m2)==Lxyz(2)) then
          M_1(:,ny)=KK(:,ny)*(cs20_ar(l1:l2,m2)/gamma0(l1:l2,m2)*&
            rho0(l1:l2,m2)-p_inf(:,ny,nnn-3))
         endif

         if (y(m1)==xyz0(2)) then
          M_5(:,1)=KK(:,1)*(cs20_ar(l1:l2,m1)/gamma0(l1:l2,m1)*&
            rho0(l1:l2,m1)-p_inf(:,1,nnn-3))
         endif

      else
        M_1=0; M_2=0; M_3=0; M_4=0; M_5=0
       endif

        do j=1,2
         if (j==1) then
           ll=l1; lll=1
           if (x(l1)==xyz0(1))   lcorner_x=.true.
         elseif (j==2) then
           ll=l2; lll=nx
           if (x(l2)==Lxyz(1))   lcorner_x=.true.
         endif
        if (lcorner_x)  then
         df(ll,m1:m2,nnn,irho_tmp) = &
           + drho_prefac(lll,:)*(L_2(lll,:)+0.5*(L_5(lll,:) + L_1(lll,:))) &
           + T_1_y(lll,:) &
           + drho_prefac(lll,:)*(N_2(lll,:)+0.5*(N_5(lll,:) + N_1(lll,:)))
         df(ll,m1:m2,nnn,iux) =  -1./&
           (2.*rho0(ll,m1:m2)*cs0_ar(ll,m1:m2))*(L_5(lll,:) - L_1(lll,:)) &
           +T_2_y(lll,:)-N_3(lll,:)
         df(ll,m1:m2,nnn,iuy) =  - L_3(lll,:) +T_3_y(lll,:)-N_4(lll,:)
         df(ll,m1:m2,nnn,iuz) =  - L_4(lll,:)  + T_4_y(lll,:) - 1./&
            (2.*rho0(ll,m1:m2)*cs0_ar(ll,m1:m2))*(N_5(lll,:) - N_1(lll,:))
         df(ll,m1:m2,nnn,ilnTT) = drho_prefac(lll,:)*(-L_2(lll,:) &
           + 0.5*(gamma0(ll,m1:m2)-1.)*(L_5(lll,:)+L_1(lll,:)))  &
           + T_5_y(lll,:) &
           + drho_prefac(lll,:)*(-N_2(lll,:) &
           + 0.5*(gamma0(ll,m1:m2)-1.)*(N_5(lll,:)+N_1(lll,:)))
         lcorner_x=.false.
       endif
       enddo




        do i=1,2
         if (i==1) then
          mm=m1; mmm=1
          if (y(m1)==xyz0(2))   lcorner_y=.true.
         elseif (i==2) then
          mm=m2; mmm=ny
          if (y(m2)==Lxyz(2))  lcorner_y=.true.
         endif

        if (lcorner_y)  then
         df(l1:l2,mm,nnn,irho_tmp) = T_1_x(:,mmm)&
          + drho_prefac(:,mmm)*(M_2(:,mmm)+0.5*(M_5(:,mmm) + M_1(:,mmm))) &
          + drho_prefac(:,mmm)*(N_2(:,mmm)+0.5*(N_5(:,mmm) &
          + N_1(:,mmm)))
         df(l1:l2,mm,nnn,iux) =  T_2_x(:,mmm) -M_3(:,mmm)-N_3(:,mmm)
         df(l1:l2,mm,nnn,iuy) =  T_3_x(:,mmm) - 1./&
           (2.*rho0(l1:l2,mm)*cs0_ar(l1:l2,mm))*(M_5(:,mmm) - M_1(:,mmm))-N_4(:,mmm)
         df(l1:l2,mm,nnn,iuz) = T_4_x(:,mmm)- M_4(:,mmm)  - 1./&
             (2.*rho0(l1:l2,mm)*cs0_ar(l1:l2,mm))*(N_5(:,mmm) - N_1(:,mmm))
          df(l1:l2,mm,nnn,ilnTT) =T_5_x(:,mmm) &
           +drho_prefac(:,mmm)*(-M_2(:,mmm) &
           +0.5*(gamma0(l1:l2,mm)-1.)*(M_5(:,mmm)+M_1(:,mmm))) &
           +drho_prefac(:,mmm)*(-N_2(:,mmm) &
           +0.5*(gamma0(l1:l2,mm)-1.)*(N_5(:,mmm)+N_1(:,mmm)))
          lcorner_y=.false.
        endif
       enddo


         do i=1,2
         if (i==1) then
          mm=m1; mmm=1
          if (y(m1)==xyz0(2))   lcorner_y=.true.
         elseif (i==2) then
          mm=m2; mmm=ny
          if (y(m2)==Lxyz(2))  lcorner_y=.true.
         endif

         do j=1,2
         if (j==1) then
           ll=l1; lll=1
           if (x(l1)==xyz0(1))   lcorner_x=.true.
         elseif (j==2) then
           ll=l2; lll=nx
           if (x(l2)==Lxyz(1))   lcorner_x=.true.
         endif

   !     if ((x(l1)==xyz0(1))   .or. (x(l2)==Lxyz(1)))   lcorner_x=.true.
   !     if ((y(m1)==xyz0(2))   .or. (y(m2)==Lxyz(2)))   lcorner_y=.true.

        if ((lcorner_x) .and. (lcorner_y)) then
         df(ll,mm,nnn,irho_tmp) = &
           drho_prefac(lll,mmm)*(L_2(lll,mmm)+0.5*(L_5(lll,mmm) + L_1(lll,mmm))) &
         + drho_prefac(lll,mmm)*(M_2(lll,mmm)+0.5*(M_5(lll,mmm) + M_1(lll,mmm))) &
         + drho_prefac(lll,mmm)*(N_2(lll,mmm)+0.5*(N_5(lll,mmm) + N_1(lll,mmm)))

         df(ll,mm,nnn,iux) =  -1./&
           (2.*rho0(ll,mm)*cs0_ar(ll,mm))*(L_5(lll,mmm) - L_1(lll,mmm)) &
           -M_3(lll,mmm)-N_3(lll,mmm)
         df(ll,mm,nnn,iuy) =  - L_3(lll,mmm) - 1./&
           (2.*rho0(ll,mm)*cs0_ar(ll,mm))*(M_5(lll,mmm) - M_1(lll,mmm))-N_4(lll,mmm)
         df(ll,mm,nnn,iuz) =  - L_4(lll,mmm) - M_4(lll,mmm)  - 1./&
            (2.*rho0(ll,mm)*cs0_ar(ll,mm))*(N_5(lll,mmm) - N_1(lll,mmm))
         df(ll,mm,nnn,ilnTT) = drho_prefac(lll,mmm)*(-L_2(lll,mmm) &
           +0.5*(gamma0(ll,mm)-1.)*(L_5(lll,mmm)+L_1(lll,mmm)))  &
           +drho_prefac(lll,mmm)*(-M_2(lll,mmm) &
           +0.5*(gamma0(ll,mm)-1.)*(M_5(lll,mmm)+M_1(lll,mmm))) &
           +drho_prefac(lll,mmm)*(-N_2(lll,mmm) &
           +0.5*(gamma0(ll,mm)-1.)*(N_5(lll,mmm)+N_1(lll,mmm)))
          lcorner_x=.false.
          lcorner_y=.false.
       endif
       enddo
       enddo

     endif

      if (nchemspec>1) then
       do k=1,nchemspec
         call der_onesided_4_slice(f,sgn,ichemspec(k),dYk_dz(l1:l2,m1:m2),nnn,3)

         do i=m1,m2
          call der_pencil(1,f(:,i,nnn,ichemspec(k)),dYk_dx(:,i))
         enddo
         do i=l1,l2
          call der_pencil(2,f(i,:,nnn,ichemspec(k)),dYk_dy(i,:))
         enddo


          df(l1:l2,m1:m2,nnn,ichemspec(k))=&
              -f(l1:l2,m1:m2,nnn,iux)*dYk_dx(l1:l2,m1:m2) &
             -f(l1:l2,m1:m2,nnn,iuy)*dYk_dy(l1:l2,m1:m2) &
              -f(l1:l2,m1:m2,nnn,iuz)*dYk_dz(l1:l2,m1:m2) &
               +RHS_Y(l1:l2,m1:m2,nnn,k)

       ! if (lfilter) then
         do i=l1,l2
         do j=m1,m2
           if ((f(i,j,nnn,ichemspec(k))+df(i,j,nnn,ichemspec(k))*dt)<-1e-25 ) then
             df(i,j,nnn,ichemspec(k))=-1e-25*dt
           endif
           if ((f(i,j,nnn,ichemspec(k))+df(i,j,nnn,ichemspec(k))*dt)>1.) then
             f(i,j,nnn,ichemspec(k))=1.*dt
           endif
         enddo
         enddo
       ! endif
       enddo
      endif
!
!

    endsubroutine bc_nscbc_nref_subout_z
!***********************************************************************
  subroutine no_nscbc(f,topbot,val)
!
!   nscbc case
!   subsonic inflow boundary conditions
!   now it is 2D case (should be corrected for 3D case)
!    ux,uy,T are fixed, drho  is calculated
!
!   16-nov-08/natalia: coded.
!
      use EquationOfState, only: cs0, cs20
      use Deriv, only: der_onesided_4_slice, der_pencil
      use Chemistry
!
      real, dimension (mx,my,mz,mfarray) :: f
      character (len=3) :: topbot
      real, dimension (mcom), optional :: val
      integer :: lll, sgn,k
      real :: u_t, T_t, lnrho_t
      real, dimension(nchemspec) :: YYi
!
      intent(inout) :: f
!
      logical :: non_zero_transveral_velo
      integer :: dir, dir2, dir3, igrid, jgrid, imin, imax, jmin, jmax
      real, allocatable, dimension(:,:,:) :: u_in
      real, allocatable, dimension(:,:) :: scalar_profile
!
! reading data from file
!
      if (inlet_from_file) then

        dir=1; dir2=2; dir3=3
        igrid=ny; jgrid=nz
        imin=m1; imax=m2
        jmin=n1; jmax=n2
        non_zero_transveral_velo=.false.
!
        allocate(u_in(igrid,jgrid,3))
!
        call find_velocity_at_inlet(u_in,non_zero_transveral_velo,&
              Lx_in,nx_in,u_t,dir,m1_in,m2_in,n1_in,n2_in,imin,imax,jmin,jmax,&
              igrid,jgrid,scalar_profile)
      endif
!
!    if (.not.present(val)) call stop_it(&
!          'bc_nscbc_subin_x: you must specify fbcx)')
!
      u_t=val(iux)
      T_t=val(ilnTT)
      lnrho_t=val(ilnrho)
      do k=1,nchemspec
       YYi(k)=val(ichemspec(k))
      enddo
!
      select case (topbot)
      case ('bot')
        lll = l1
        sgn = 1
      case ('top')
        lll = l2
        sgn = -1
      case default
         print*, "bc_no_nscbc: ", topbot, " should be `top' or `bot'"
      endselect
!
        do k=1,nchemspec
         f(lll,m1:m2,n1:n2,ichemspec(k))=YYi(k)
        enddo
         f(lll,m1:m2,n1:n2,ilnTT)  = T_t
!         f(lll,m1:m2,n1:n2,ilnrho) = lnrho_t
!
        if (inlet_from_file) then
         f(lll,m1:m2,n1:n2,iux) = u_t + u_in(:,:,1)
         f(lll,m1:m2,n1:n2,iuy) = u_in(:,:,2)
         f(lll,m1:m2,n1:n2,iuz) = u_in(:,:,3)
        else
         f(lll,m1:m2,n1:n2,iux) = u_t
         f(lll,m1:m2,n1:n2,iuy) = 0.
         f(lll,m1:m2,n1:n2,iuz) = 0.
        endif
!
      if (allocated(u_in)) deallocate(u_in)
!
 endsubroutine no_nscbc
!***********************************************************************
    subroutine final_velocity_profile(f,dir,topbot,imin,imax,jmin,jmax,igrid,jgrid)
!
!  Create the final inlet velocity profile.
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, allocatable, dimension(:,:)   :: u_profile
      integer :: j,kkk,jjj,bound
      integer, dimension(10) :: stat
      real :: rad_2
      integer, intent(in) :: dir,imin,imax,jmin,jmax,igrid,jgrid
      character (len=3) :: topbot
!
      select case (topbot)
        case ('bot')
          if (dir==1) bound = l1
          if (dir==2) bound = m1
          if (dir==3) bound = n1
        case ('top')
          if (dir==1) bound = l2
          if (dir==2) bound = m2
          if (dir==3) bound = n2
        case default
           print*, "final_velocity_profile: ", topbot, " should be `top' or `bot'"
      endselect
!
!  Allocate allocatables
!
      stat=0
      allocate(u_profile(igrid,jgrid  ),STAT=stat(1))
      if (maxval(stat)>0) &
          call stop_it("Couldn't allocate memory for all vars in find_velocity_at_inlet")
!
! Define velocity profile at inlet
!
        do j=1,ninit
          select case (velocity_profile(j))
!
          case ('nothing')
            if (lroot .and. it==1 .and. j == 1 .and. lfirst) &
                print*,'velocity_profile: nothing'
!
          case ('gaussian')
            if (lroot .and. it==1 .and. lfirst) &
                print*,'velocity_profile: gaussian'
            do jjj=imin,imax
              do kkk=jmin,jmax
                if (dir==1) then
                  rad_2=((y(jjj)-jet_center(1))**2+(z(kkk)-jet_center(1))**2)
                elseif (dir==2) then
                  rad_2=((x(jjj)-jet_center(1))**2+(z(kkk)-jet_center(1))**2)
                elseif (dir==3) then
                  rad_2=((x(jjj)-jet_center(1))**2+(y(kkk)-jet_center(1))**2)
                endif
                  u_profile(jjj-imin+1,kkk-jmin+1)=exp(-rad_2/sigma**2)
!print*,u_profile(jjj-imin+1,kkk-jmax+1),rad_2,sigma
              enddo
            enddo
          end select
!
        enddo
!
          if (dir==1) then
            f(bound,m1:m2,n1:n2,iux) = f(bound,m1:m2,n1:n2,iux)*u_profile
            f(bound,m1:m2,n1:n2,iuy) = f(bound,m1:m2,n1:n2,iuy)*u_profile
            f(bound,m1:m2,n1:n2,iuz) = f(bound,m1:m2,n1:n2,iuz)*u_profile
          elseif (dir==2) then
            f(l1:l2,bound,n1:n2,iux) = f(l1:l2,bound,n1:n2,iux)*u_profile
            f(l1:l2,bound,n1:n2,iuy) = f(l1:l2,bound,n1:n2,iuy)*u_profile
            f(l1:l2,bound,n1:n2,iuz) = f(l1:l2,bound,n1:n2,iuz)*u_profile
          elseif (dir==3) then
            f(l1:l2,m1:m2,bound,iux) = f(l1:l2,m1:m2,bound,iux)*u_profile
            f(l1:l2,m1:m2,bound,iuy) = f(l1:l2,m1:m2,bound,iuy)*u_profile
            f(l1:l2,m1:m2,bound,iuz) = f(l1:l2,m1:m2,bound,iuz)*u_profile
          endif
!
!
!  Deallocate
!
        deallocate(u_profile)
!
      endsubroutine final_velocity_profile
!***********************************************************************
endmodule NSCBC
