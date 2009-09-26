! $Id: boundcond.f90 11648 2009-09-09 09:14:01Z nils.e.haugen $
!
!  Module for NSCBC (Navier-Stokes Characteristic Boundary Conditions).  
!  NSCBC is an alternative way of imposing (time-dependent) boundary 
!  conditions through solving differential equations on the boundaries.
!
! 2009.09.09 (Nils Erland) : Moved all NSCBC stuff from boundcond.f90 to 
!                            this module.
!
module NSCBC
!
  use Cdata
  use Cparam
  use Messages
  use Mpicomm
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
  character(len=40) :: turb_inlet_dir=''
  real :: nscbc_sigma_out = 1.,nscbc_sigma_in = 1., p_infty=1.
  logical :: inlet_from_file=.false.
  logical :: first_NSCBC=.true.
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
!
namelist /NSCBC_init_pars/  &
    nscbc_bc, nscbc_sigma_in, nscbc_sigma_out, p_infty, inlet_from_file,&
    turb_inlet_dir

namelist /NSCBC_run_pars/  &
    nscbc_bc, nscbc_sigma_in, nscbc_sigma_out, p_infty, inlet_from_file,&
    turb_inlet_dir
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
      use Chemistry, only: bc_nscbc_subin_x,bc_nscbc_nref_subout_x,&
          bc_nscbc_nref_subout_y

      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      character (len=nscbc_len), dimension(3) :: bc12
      character (len=3) :: topbot
      character (len=60) :: turbfile
      integer j,k,direction,ip_ok,ip_test
      real, dimension(mcom) :: valx,valy,valz
      logical :: proc_at_inlet
      integer :: ipx_in, ipy_in, ipz_in, iproc_in, nprocx_in, nprocy_in, nprocz_in
      character (len=120) :: directory_in
      character (len=5) :: chproc_in

      intent(inout) :: f
      intent(inout) :: df
      intent(in)    :: j
!
      proc_at_inlet=.false.
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
                if ((topbot=='bot'.and.ipx==0).or.&
                    (topbot=='top'.and.ipx==nprocx-1)) then
                  proc_at_inlet=.true.
                  ipx_in=0
                  nprocx_in=1
                endif
              elseif (j==2) then
                if ((topbot=='bot'.and.ipy==0).or.&
                    (topbot=='top'.and.ipy==nprocy-1)) then
                  proc_at_inlet=.true.
                  ipy_in=0
                  nprocy_in=1
                endif
              elseif (j==3) then
                if ((topbot=='bot'.and.ipz==0).or.&
                    (topbot=='top'.and.ipz==nprocz-1)) then
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
!             print*,'datadir=',datadir
              call chn(iproc_in,chproc_in)
              call safe_character_assign(directory_in,&
                  trim(turb_inlet_dir)//'/data/proc'//chproc_in)
!              print*,'directory_in=',directory_in
              call safe_character_assign(turbfile,&
                  trim(directory_in)//'/var.dat')
!              print*,'turbfile=',turbfile
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
!print*,'Lx_in,Ly_in,Lz_in=',Lx_in,Ly_in,Lz_in
!print*,'l1_in,m1_in,n1_in=',l1_in,m1_in,n1_in
!print*,'l2_in,m2_in,n2_in=',l2_in,m2_in,n2_in
!
!  Check size of arrays
!            
!!$              if (j==1) then
!!$                if ((my_in .ne. my) .or. (mz_in .ne. mz)) then
!!$                  print*,'my_in,my,mz_in,mz=',my_in,my,mz_in,mz
!!$                  call stop_it ('The turbulent inlet grid does not conform!')
!!$                endif
!!$                if ((Ly_in .ne. Lxyz(2)/nprocy) .or. (Lz_in .ne. Lxyz(3)/nprocz)) then
!!$                  print*,'Lx_in,Ly_in,Lz_in,Lxyz=',Lx_in,Ly_in,Lz_in,Lxyz
!!$                  call stop_it ('The turbulent inlet grid does not conform!')
!!$                endif                
!!$              elseif (j==2) then
!!$                if ((mx_in .ne. mx) .or. (mz_in .ne. mz)) then
!!$                  print*,'mx_in,mx,mz_in,mz=',mx_in,mx,mz_in,mz
!!$                  call stop_it ('The turbulent inlet grid does not conform !')
!!$                endif
!!$                if ((Lx_in .ne. Lxyz(1)/nprocx) .or. (Lz_in .ne. Lxyz(3)/nprocz)) then
!!$                  print*,'Lx_in,Ly_in,Lz_in,Lxyz=',Lx_in,Ly_in,Lz_in,Lxyz
!!$                  call stop_it ('The turbulent inlet grid does not conform!')
!!$                endif                
!!$              elseif (j==3) then
!!$                if ((my_in .ne. my) .or. (mx_in .ne. mx)) then
!!$                  print*,'my_in,my,mx_in,mx=',my_in,my,mx_in,mx
!!$                  call stop_it ('The turbulent inlet grid does not conform!')
!!$                endif
!!$                if ((Lx_in .ne. Lxyz(1)/nprocx) .or. (Ly_in .ne. Lxyz(2)/nprocz)) then
!!$                  print*,'Lx_in,Ly_in,Lz_in,Lxyz=',Lx_in,Ly_in,Lz_in,Lxyz
!!$                  call stop_it ('The turbulent inlet grid does not conform!')
!!$                endif                
!!$              endif
            endif            
          endif
        endif
!
!  Check if this is a physical boundary
!
        if (ip_test==ip_ok) then
          select case(bc12(j))
          case('part_ref_outlet')
!   Partially reflecting outlet.
            if (j==1) then 
              call bc_nscbc_prf_x(f,df,topbot,.false.)
            elseif (j==2) then 
              call bc_nscbc_prf_y(f,df,topbot,.false.)
            elseif (j==3) then 
              call fatal_error("nscbc_boundtreat_xyz",'bc_nscbc_prf_z is not yet implemented')
            endif
          case('part_ref_inlet')
!   Partially reflecting inlet, ie. impose a velocity u_t.
            if (j==1) then 
              direction = 1
              call bc_nscbc_prf_x(f,df,topbot,.true.,linlet=.true.,&
                  u_t=valx(direction))
            elseif (j==2) then 
              direction = 2
              call bc_nscbc_prf_y(f,df,topbot,.true.,linlet=.true.,&
                  u_t=valy(direction))
            elseif (j==3) then 
              direction = 3
              call fatal_error("nscbc_boundtreat_xyz",'bc_nscbc_prf_z is not yet implemented')
            endif
          case('ref_inlet')
!   Partially reflecting inlet, ie. impose a velocity u_t.
            if (j==1) then 
              direction = 1
              call bc_nscbc_prf_x(f,df,topbot,.false.,linlet=.true.,&
                  u_t=valx(direction))
            elseif (j==2) then 
              direction = 2
              call bc_nscbc_prf_y(f,df,topbot,.false.,linlet=.true.,&
                  u_t=valy(direction))
            elseif (j==3) then 
              direction = 3
              call fatal_error("nscbc_boundtreat_xyz",'bc_nscbc_prf_z is not yet implemented')
            endif
          case('subsonic_inflow')
! Subsonic inflow 
            if (j==1) then
              call bc_nscbc_subin_x(f,df,topbot,valx)
            elseif (j==2) then
            endif
          case('subson_nref_outflow')
            if (j==1) then
              call bc_nscbc_nref_subout_x(f,df,topbot,nscbc_sigma_out)
            elseif (j==2) then
              call bc_nscbc_nref_subout_y(f,df,topbot,nscbc_sigma_out)
            endif
          case('')
!   Do nothing.
          case('none')
            print*,'nscbc_boundtreat_xyz: doing nothing!'
!   Do nothing.
          case default
            call fatal_error("nscbc_boundtreat_xyz",'You must specify nscbc bouncond!')
          endselect
        endif
      end do
!
    endsubroutine
!***********************************************************************
    subroutine bc_nscbc_prf_x(f,df,topbot,non_reflecting_inlet,linlet,u_t)
!
!   Calculate du and dlnrho at a partially reflecting outlet/inlet normal to 
!   x-direction acc. to LODI relations. Uses a one-sided finite diff. stencil.
!
!   7-jul-08/arne: coded.
!  25-nov-08/nils: extended to work in multiple dimensions and with cross terms
!                  i.e. not just the LODI equations.
!
      use EquationOfState, only: cs0, cs20
      use Deriv, only: der_onesided_4_slice, der_pencil, der2_pencil
      use Chemistry
      use Viscosity

      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      character (len=3) :: topbot
      logical, optional :: linlet
      logical :: llinlet, non_reflecting_inlet
      real, optional :: u_t
      real, dimension(ny,nz) :: dlnrho_dx
      real, dimension(ny,nz) :: rho0, L_1, L_3, L_4, L_5,parallell_term_uz
      real, dimension(ny,nz) :: parallell_term_rho,dlnrho_dy,dlnrho_dz
      real, dimension(ny,nz) :: parallell_term_ux,d2u1_dy2,d2u1_dz2
      real, dimension(ny,nz) :: d2u2_dy2,d2u2_dz2,d2u3_dy2,d2u3_dz2
      real, dimension(ny,nz) :: prefac1, prefac2,parallell_term_uy
      real, dimension(ny,nz,3) :: grad_rho, u_in
      real, dimension(ny,nz,3,3) :: dui_dxj
      real, dimension (my,mz) :: cs0_ar,cs20_ar
      real, dimension (my,mz) :: tmp22,tmp12,tmp2_lnrho,tmp33,tmp13,tmp3_lnrho
      real, dimension (my,mz) :: tmp23,tmp32
      real, dimension (ny) :: tmpy
      real, dimension (nz) :: tmpz
      real :: Mach,KK,nu
      integer lll,i
      integer sgn
      real :: shift, grid_shift, weight, round
      integer :: iround,lowergrid,uppergrid

      intent(inout) :: f
      intent(inout) :: df

      llinlet = .false.
      if (present(linlet)) llinlet = linlet
      if (llinlet.and..not.present(u_t)) call stop_it(&
           'bc_nscbc_prf_x: when using linlet=T, you must also specify u_t)')
      select case(topbot)
      case('bot')
        lll = l1
        sgn = 1
      case('top')
        lll = l2
        sgn = -1
      case default
        print*, "bc_nscbc_prf_x: ", topbot, " should be `top' or `bot'"
      endselect
!
!  Find density
!
      if (ldensity_nolog) then
        rho0 = f(lll,m1:m2,n1:n2,irho)
      else
        rho0 = exp(f(lll,m1:m2,n1:n2,ilnrho))
      endif
!
!  Get viscoity
!
      call getnu(nu)
!
!  Set arrays for the speed of sound and for the speed of sound squared (is it
!  really necessarry to have both arrays?) 
!  Set prefactors to be used later.
!
      if (leos_idealgas) then
        cs20_ar(m1:m2,n1:n2)=cs20
        cs0_ar(m1:m2,n1:n2)=cs0
        prefac1 = -1./(2.*cs20)
        prefac2 = -1./(2.*rho0*cs0)
      elseif (leos_chemistry) then
        call fatal_error('bc_nscbc_prf_x',&
            'This sub routine is not yet adapted to work with leos_chemsitry!')
      else
        print*,"bc_nscbc_prf_x: leos_idealgas=",leos_idealgas,"."
        print*,"NSCBC boundary treatment only implemented for an ideal gas." 
        print*,"Boundary treatment skipped."
        return
      endif
!
!  Calculate one-sided derivatives in the boundary normal direction
!
      call der_onesided_4_slice(f,sgn,ilnrho,grad_rho(:,:,1),lll,1)
      call der_onesided_4_slice(f,sgn,iux,dui_dxj(:,:,1,1),lll,1)
      call der_onesided_4_slice(f,sgn,iuy,dui_dxj(:,:,2,1),lll,1)
      call der_onesided_4_slice(f,sgn,iuz,dui_dxj(:,:,3,1),lll,1)
!
!  Do central differencing in the directions parallell to the boundary 
!  first in the y-direction......
!
      if (nygrid /= 1) then
        do i=n1,n2
          call der_pencil(2,f(lll,:,i,iuy),tmp22(:,i))
          call der_pencil(2,f(lll,:,i,ilnrho),tmp2_lnrho(:,i))
          call der_pencil(2,f(lll,:,i,iux),tmp12(:,i))
          call der_pencil(2,f(lll,:,i,iuz),tmp32(:,i))
          call der2_pencil(2,f(lll,:,i,iux),tmpy)
          d2u1_dy2(:,i-n1+1)=tmpy(:)
          call der2_pencil(2,f(lll,:,i,iuy),tmpy)
          d2u2_dy2(:,i-n1+1)=tmpy(:)
          call der2_pencil(2,f(lll,:,i,iuz),tmpy)
          d2u3_dy2(:,i-n1+1)=tmpy(:)
        enddo
      else
        tmp32=0
        tmp22=0
        tmp12=0
        tmp2_lnrho=0
        d2u1_dy2=0
        d2u2_dy2=0
        d2u3_dy2=0
      endif
      dui_dxj(:,:,1,2)=tmp12(m1:m2,n1:n2)
      dui_dxj(:,:,2,2)=tmp22(m1:m2,n1:n2)
      dui_dxj(:,:,3,2)=tmp32(m1:m2,n1:n2)
      grad_rho(:,:,2)=tmp2_lnrho(m1:m2,n1:n2)
!
!  .... then in the z-direction
!
      if (nzgrid /= 1) then
        do i=m1,m2
          call der_pencil(3,f(lll,i,:,iuz),tmp33(i,:))
          call der_pencil(3,f(lll,i,:,ilnrho),tmp3_lnrho(i,:))
          call der_pencil(3,f(lll,i,:,iux),tmp13(i,:))
          call der_pencil(3,f(lll,i,:,iuy),tmp23(i,:))
          call der2_pencil(3,f(lll,i,:,iux),tmpz)
          d2u1_dz2(i-m1+1,:)=tmpz(:)
          call der2_pencil(3,f(lll,i,:,iuy),tmpz)
          d2u2_dz2(i-m1+1,:)=tmpz(:)
          call der2_pencil(3,f(lll,i,:,iuz),tmpz)
          d2u3_dz2(i-m1+1,:)=tmpz(:)
        enddo
      else
        tmp33=0
        tmp23=0
        tmp13=0
        tmp3_lnrho=0
        d2u1_dz2=0
        d2u2_dz2=0
        d2u3_dz2=0
      endif
      dui_dxj(:,:,1,3)=tmp13(m1:m2,n1:n2)
      dui_dxj(:,:,2,3)=tmp23(m1:m2,n1:n2)
      dui_dxj(:,:,3,3)=tmp33(m1:m2,n1:n2)
      grad_rho(:,:,3)=tmp3_lnrho(m1:m2,n1:n2)
!
!  Find divergence of rho if we solve for logarithm of rho
!
      if (.not. ldensity_nolog) then
        do i=1,3
          grad_rho(:,:,i)=grad_rho(:,:,i)*rho0
        enddo
      endif
!
!  Find Mach number 
!  (NILS: I do not think this is a good way to determine the Mach
!  number since this is a local Mach number for this processor. Furthermore
!  by determining the Mach number like this we will see the Mach number varying
!  with the phase of an acoustic wave as the wave pass through the boundary.
!  I think that what we really want is a Mach number averaged over the 
!  timescale of several acoustic waves. How could this be done????)
!
        Mach=sum(f(lll,m1:m2,n1:n2,iux)/cs0_ar(m1:m2,n1:n2))/(ny*nz)
!
!  Find the L_i's (which really is the Lodi equations)
!
      if (llinlet) then
!  This L_1 is correct only for p=cs^2*rho
        L_1 = (f(lll,m1:m2,n1:n2,iux) - sgn*cs0_ar(m1:m2,n1:n2))&
            *(cs20_ar(m1:m2,n1:n2)*grad_rho(:,:,1) &
            - sgn*rho0*cs0_ar(m1:m2,n1:n2)*dui_dxj(:,:,1,1))
!
!  Find velocity at inlet
!
        if (inlet_from_file) then
          round=t*u_t/Lx_in
          iround=int(round)
          shift=round-iround
          grid_shift=shift*nx_in
          lowergrid=l1_in+int(grid_shift)
          uppergrid=lowergrid+1
          weight=grid_shift-int(grid_shift)
          u_in(:,:,:)&
              =f_in(lowergrid,m1_in:m2_in,n1_in:n2_in,iux:iuz)*(1-weight)&
              +f_in(uppergrid,m1_in:m2_in,n1_in:n2_in,iux:iuz)*weight
          u_in(:,:,1)=u_in(:,:,1)+u_t

        else
          u_in(:,:,1)=u_t
          u_in(:,:,2)=0.
          u_in(:,:,3)=0.
        endif
        if (non_reflecting_inlet) then
!
!  The inlet in non-reflecting only when nscbc_sigma_in is set to 0, this 
!  might however lead to problems as the inlet velocity will tend to drift 
!  away from the target velocity u_t. This problem should be overcome by 
!  setting a small but non-zero nscbc_sigma_in.
!
          L_3=nscbc_sigma_in*(f(lll,m1:m2,n1:n2,iuy)-u_in(:,:,2))&
              *cs0_ar(m1:m2,n1:n2)/Lxyz(1)
          L_4=nscbc_sigma_in*(f(lll,m1:m2,n1:n2,iuz)-u_in(:,:,3))&
              *cs0_ar(m1:m2,n1:n2)/Lxyz(1)
          L_5 = nscbc_sigma_in*cs20_ar(m1:m2,n1:n2)*rho0&
              *sgn*(f(lll,m1:m2,n1:n2,iux)-u_in(:,:,1))*(1-Mach**2)/Lxyz(1)
        else
          L_3=0
          L_4=0
          L_5 = L_1
        endif
      else
!
!  Find the parameter determining 
!
        KK=nscbc_sigma_out*(1-Mach**2)*cs0/Lxyz(1)
!
!  Find the L_i's
!
        L_1 = KK*(rho0*cs20-p_infty)
        L_3 = f(lll,m1:m2,n1:n2,iux)*dui_dxj(:,:,2,1)
        L_4 = f(lll,m1:m2,n1:n2,iux)*dui_dxj(:,:,3,1)
        L_5 = (f(lll,m1:m2,n1:n2,iux) - sgn*cs0_ar(m1:m2,n1:n2))*&
             (cs20_ar(m1:m2,n1:n2)*grad_rho(:,:,1)&
             - sgn*rho0*cs0_ar(m1:m2,n1:n2)*dui_dxj(:,:,1,1))
      end if
!
!  Add terms due to derivatives parallell to the boundary
!  NILS: Viscous terms in the x direction are missing!
!
!!$      parallell_term_rho &
!!$           =rho0*dui_dxj(:,:,2,2)+f(lll,m1:m2,n1:n2,iuy)*grad_rho(:,:,2)&
!!$           +rho0*dui_dxj(:,:,3,3)+f(lll,m1:m2,n1:n2,iuz)*grad_rho(:,:,3)
!!$      parallell_term_ux &
!!$           =f(lll,m1:m2,n1:n2,iuy)*dui_dxj(:,:,1,2)&
!!$           +f(lll,m1:m2,n1:n2,iuz)*dui_dxj(:,:,1,3)&
!!$           +nu*(d2u1_dy2+d2u1_dz2)
!!$      parallell_term_uy &
!!$           =f(lll,m1:m2,n1:n2,iuy)*dui_dxj(:,:,2,2)&
!!$           +f(lll,m1:m2,n1:n2,iuz)*dui_dxj(:,:,2,3)&
!!$           +cs20_ar(m1:m2,n1:n2)*grad_rho(:,:,2)/rho0&
!!$           +nu*(d2u2_dy2+d2u2_dz2)
!!$      parallell_term_uz &
!!$           =f(lll,m1:m2,n1:n2,iuy)*dui_dxj(:,:,3,2)&
!!$           +f(lll,m1:m2,n1:n2,iuz)*dui_dxj(:,:,3,3)&
!!$           +cs20_ar(m1:m2,n1:n2)*grad_rho(:,:,3)/rho0&
!!$           +nu*(d2u3_dy2+d2u3_dz2)


!print*,'x----------------------------------------------------'
!!$print*,'uy=',f(lll,m1:m2,n1:n2,iuy)
!!$print*,'uz=',f(lll,m1:m2,n1:n2,iuz)
!!$print*,'dui_dxj(:,:,1,2)=',dui_dxj(:,:,1,2)
!!$print*,'dui_dxj(:,:,1,3)=',dui_dxj(:,:,1,3)
!!$print*,'parallell_term_ux=',parallell_term_ux
!!$print*,maxval(f(lll,m1:m2,n1:n2,iuy)),maxval(dui_dxj(:,:,1,2))
!!$print*,maxval(f(lll,m1:m2,n1:n2,iuz)),maxval(dui_dxj(:,:,1,3))
!!$print*,minval(f(lll,m1:m2,n1:n2,iuy)),minval(dui_dxj(:,:,1,2))
!!$print*,minval(f(lll,m1:m2,n1:n2,iuz)),minval(dui_dxj(:,:,1,3))
!!$print*,minval(parallell_term_rho),maxval(parallell_term_rho)
!!$print*,minval(parallell_term_ux),maxval(parallell_term_ux)
!!$print*,minval(parallell_term_uy),maxval(parallell_term_uy)
!!$print*,minval(parallell_term_uz),maxval(parallell_term_uz)


!
!  NILS: Currently the implementation with the parallell terms does not 
!  NILS: seem to work. Set at parallell terms to zero for now.
!
      parallell_term_rho=0
      parallell_term_ux=0
      parallell_term_uy=0
      parallell_term_uz=0
!
!  Find the evolution equations at the boundary
!
      select case(topbot)
     ! NB: For 'top' L_1 plays the role of L5 and L_5 the role of L1
      case('bot')
        df(lll,m1:m2,n1:n2,ilnrho) = prefac1*(L_1 + L_5)-parallell_term_rho
        if (llinlet) then
          df(lll,m1:m2,n1:n2,iux) = prefac2*( L_5 - L_1)-parallell_term_ux
        else
          df(lll,m1:m2,n1:n2,iux) = prefac2*( L_1 - L_5)-parallell_term_ux
        endif
        df(lll,m1:m2,n1:n2,iuy) = -L_3-parallell_term_uy
        df(lll,m1:m2,n1:n2,iuz) = -L_4-parallell_term_uz
      case('top')
        df(lll,m1:m2,n1:n2,ilnrho) = prefac1*(L_1 + L_5)-parallell_term_rho
        if (llinlet) then
          df(lll,m1:m2,n1:n2,iux) = prefac2*( L_1 - L_5)-parallell_term_ux
        else
          df(lll,m1:m2,n1:n2,iux) = prefac2*(-L_1 + L_5)-parallell_term_ux
        endif
        df(lll,m1:m2,n1:n2,iuy) = -L_3-parallell_term_uy
        df(lll,m1:m2,n1:n2,iuz) = -L_4-parallell_term_uz
      endselect
!
!  Check if we are solving for logrho or rho
!
      if (.not. ldensity_nolog) then
        df(lll,m1:m2,n1:n2,irho)=df(lll,m1:m2,n1:n2,irho)/rho0
      endif
!
! Impose required variables at the boundary
!
      if (llinlet) then
        if (.not. non_reflecting_inlet) then
          f(lll,m1:m2,n1:n2,iux) = u_in(:,:,1)
          f(lll,m1:m2,n1:n2,iuy) = u_in(:,:,2)
          f(lll,m1:m2,n1:n2,iuz) = u_in(:,:,3)
        endif
      endif

    endsubroutine bc_nscbc_prf_x
!***********************************************************************
    subroutine bc_nscbc_prf_y(f,df,topbot,non_reflecting_inlet,linlet,u_t)
!
!   Calculate du and dlnrho at a partially reflecting outlet/inlet normal to 
!   y-direction acc. to LODI relations. Uses a one-sided finite diff. stencil.
!
!   7-jul-08/arne: coded.
!  25-nov-08/nils: extended to work in multiple dimensions and with cross terms
!                  i.e. not just the LODI equations.
!
      use EquationOfState, only: cs0, cs20
      use Deriv, only: der_onesided_4_slice, der_pencil, der2_pencil
      use Chemistry
      use Viscosity

      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      character (len=3) :: topbot
      logical, optional :: linlet
      logical :: llinlet, non_reflecting_inlet
      real, optional :: u_t
      real, dimension(nx,nz) :: dlnrho_dx,dlnrho_dy,dlnrho_dz
      real, dimension(nx,nz) :: rho0, L_1, L_3, L_4, L_5,parallell_term_uz
      real, dimension(nx,nz) :: parallell_term_rho
      real, dimension(nx,nz) :: parallell_term_ux,d2u1_dx2,d2u1_dz2
      real, dimension(nx,nz) :: d2u2_dx2,d2u2_dz2,d2u3_dx2,d2u3_dz2
      real, dimension(nx,nz) :: prefac1, prefac2,parallell_term_uy
      real, dimension(nx,nz,3) :: grad_rho, u_in
      real, dimension(nx,nz,3,3) :: dui_dxj
      real, dimension (mx,mz) :: cs0_ar,cs20_ar
      real, dimension (mx,mz) :: tmp22,tmp12,tmp2_lnrho,tmp33,tmp13,tmp3_lnrho
      real, dimension (mx,mz) :: tmp23,tmp32,tmp21,tmp31,tmp11,tmp1_lnrho
      real, dimension (nx) :: tmpx
      real, dimension (nz) :: tmpz
      real :: Mach,KK,nu
      integer lll,i
      integer sgn
      real :: shift, grid_shift, weight, round
      integer :: iround,lowergrid,uppergrid

      intent(inout) :: f
      intent(inout) :: df


      llinlet = .false.
      if (present(linlet)) llinlet = linlet
      if (llinlet.and..not.present(u_t)) call stop_it(&
           'bc_nscbc_prf_y: when using linlet=T, you must also specify u_t)')
      select case(topbot)
      case('bot')
        lll = m1
        sgn = 1
      case('top')
        lll = m2
        sgn = -1
      case default
        print*, "bc_nscbc_prf_y: ", topbot, " should be `top' or `bot'"
      endselect
!
!  Find density
!
      if (ldensity_nolog) then
        rho0 = f(l1:l2,lll,n1:n2,irho)
      else
        rho0 = exp(f(l1:l2,lll,n1:n2,ilnrho))
      endif
!
!  Get viscoity
!
      call getnu(nu)
!
!  Set arrays for the speed of sound and for the speed of sound squared (is it
!  really necessarry to have both arrays?) 
!  Set prefactors to be used later.
!
      if (leos_idealgas) then
        cs20_ar(l1:l2,n1:n2)=cs20
        cs0_ar(l1:l2,n1:n2)=cs0
        prefac1 = -1./(2.*cs20)
        prefac2 = -1./(2.*rho0*cs0)
      elseif (leos_chemistry) then
        call fatal_error('bc_nscbc_prf_x',&
            'This sub routine is not yet adapted to work with leos_chemsitry!')
      else
        print*,"bc_nscbc_prf_y: leos_idealgas=",leos_idealgas,"."
        print*,"NSCBC boundary treatment only implemented for an ideal gas." 
        print*,"Boundary treatment skipped."
        return
      endif
!
!  Calculate one-sided derivatives in the boundary normal direction
!
      call der_onesided_4_slice(f,sgn,ilnrho,grad_rho(:,:,2),lll,2)
      call der_onesided_4_slice(f,sgn,iux,dui_dxj(:,:,1,2),lll,2)
      call der_onesided_4_slice(f,sgn,iuy,dui_dxj(:,:,2,2),lll,2)
      call der_onesided_4_slice(f,sgn,iuz,dui_dxj(:,:,3,2),lll,2)
!
!  Do central differencing in the directions parallell to the boundary 
!  first in the x-direction......
!
      if (nxgrid /= 1) then
        do i=n1,n2
          call der_pencil(1,f(:,lll,i,iuy),tmp21(:,i))
          call der_pencil(1,f(:,lll,i,ilnrho),tmp1_lnrho(:,i))
          call der_pencil(1,f(:,lll,i,iux),tmp11(:,i))
          call der_pencil(1,f(:,lll,i,iuz),tmp31(:,i))
          call der2_pencil(1,f(:,lll,i,iux),tmpx)
          d2u1_dx2(:,i-n1+1)=tmpx
          call der2_pencil(1,f(:,lll,i,iuy),tmpx)
          d2u2_dx2(:,i-n1+1)=tmpx
          call der2_pencil(1,f(:,lll,i,iuz),tmpx)
          d2u3_dx2(:,i-n1+1)=tmpx
        enddo
      else
        tmp31=0
        tmp21=0
        tmp11=0
        tmp1_lnrho=0
        d2u1_dx2=0
        d2u2_dx2=0
        d2u3_dx2=0
      endif
      dui_dxj(:,:,3,1)=tmp31(l1:l2,n1:n2)
      dui_dxj(:,:,2,1)=tmp21(l1:l2,n1:n2)
      dui_dxj(:,:,1,1)=tmp11(l1:l2,n1:n2)
      grad_rho(:,:,1)=tmp1_lnrho(l1:l2,n1:n2)
!
!  .... then in the z-direction
!
      if (nzgrid /= 1) then
        do i=l1,l2
          call der_pencil(3,f(i,lll,:,iuz),tmp33(i,:))
          call der_pencil(3,f(i,lll,:,ilnrho),tmp3_lnrho(i,:))
          call der_pencil(3,f(i,lll,:,iux),tmp13(i,:))
          call der_pencil(3,f(i,lll,:,iuy),tmp23(i,:))
          call der2_pencil(3,f(i,lll,:,iux),tmpz)
          d2u1_dz2(i,:)=tmpz
          call der2_pencil(3,f(i,lll,:,iuy),tmpz)
          d2u2_dz2(i,:)=tmpz
          call der2_pencil(3,f(i,lll,:,iuz),tmpz)
          d2u3_dz2(i,:)=tmpz
        enddo
      else
        tmp33=0
        tmp23=0
        tmp13=0
        tmp3_lnrho=0
        d2u1_dz2=0
        d2u2_dz2=0
        d2u3_dz2=0
      endif
      dui_dxj(:,:,3,3)=tmp33(l1:l2,n1:n2)
      dui_dxj(:,:,2,3)=tmp23(l1:l2,n1:n2)
      dui_dxj(:,:,1,3)=tmp13(l1:l2,n1:n2)
      grad_rho(:,:,3)=tmp3_lnrho(l1:l2,n1:n2)
!
!  Find divergence of rho if we solve for logarithm of rho
!
      if (.not. ldensity_nolog) then
        do i=1,3
          grad_rho(:,:,i)=grad_rho(:,:,i)*rho0
        enddo
      endif
!
!  Find Mach number 
!  (NILS: I do not think this is a good way to determine the Mach
!  number since this is a local Mach number for this processor. Furthermore
!  by determining the Mach number like this we will see the Mach number varying
!  with the phase of an acoustic wave as the wave pass through the boundary.
!  I think that what we really want is a Mach number averaged over the 
!  timescale of several acoustic waves. How could this be done????)
!
      Mach=sum(f(l1:l2,lll,n1:n2,iuy)/cs0_ar(l1:l2,n1:n2))/(nx*nz)
!
!  Find the L_i's (which really are the Lodi equations)
!
      if (llinlet) then
        L_1 = (f(l1:l2,lll,n1:n2,iuy) - sgn*cs0_ar(l1:l2,n1:n2))*&
            (cs20_ar(l1:l2,n1:n2)*grad_rho(:,:,2) &
            - sgn*rho0*cs0_ar(l1:l2,n1:n2)*dui_dxj(:,:,2,2))
!
!  Find velocity at inlet
!
        if (inlet_from_file) then
!yyyy
          round=t*u_t/Ly_in
          iround=int(round)
          shift=round-iround
          grid_shift=shift*ny_in
          lowergrid=m1_in+int(grid_shift)
          uppergrid=lowergrid+1
!!$          print*,'lowergrid=',lowergrid
!!$          print*,'t=',t
!!$          print*,'u_t=',u_t
!!$          print*,'Ly_in=',Ly_in
!!$          print*,'round=',round
!!$          print*,'shift=',shift
!!$          print*,'grid_shift=',grid_shift
          weight=grid_shift-int(grid_shift)
!!$if (lroot) then
!!$print*,'lowergrid=',lowergrid
!!$print*,'t=',t
!!$print*,'u_t=',u_t
!!$print*,'Ly_in=',Ly_in
!!$print*,'round=',round
!!$print*,'iround=',iround
!!$print*,'shift=',shift
!!$print*,'grid_shift=',grid_shift
!!$endif

          u_in(:,:,:)&
              =f_in(l1_in:l2_in,lowergrid,n1_in:n2_in,iux:iuz)*(1-weight)&
              +f_in(l1_in:l2_in,uppergrid,n1_in:n2_in,iux:iuz)*weight
          u_in(:,:,2)=u_in(:,:,2)+u_t
        else
          u_in(:,:,1)=0.
          u_in(:,:,2)=u_t
          u_in(:,:,3)=0.
        endif
        if (non_reflecting_inlet) then
!
!  The inlet is non-reflecting only when nscbc_sigma_in is set to 0, this 
!  might however lead to problems as the inlet velocity will tend to drift 
!  away from the target velocity u_t. This problem should be overcome by 
!  setting a small but non-zero nscbc_sigma_in.
!
          L_3=nscbc_sigma_in*(f(l1:l2,lll,n1:n2,iux)-u_in(:,:,1))&
              *cs0_ar(l1:l2,n1:n2)/Lxyz(2)
          L_4=nscbc_sigma_in*(f(l1:l2,lll,n1:n2,iuz)-u_in(:,:,3))&
              *cs0_ar(l1:l2,n1:n2)/Lxyz(2)
          L_5 = nscbc_sigma_in*cs20_ar(l1:l2,n1:n2)*rho0&
              *sgn*(f(l1:l2,lll,n1:n2,iuy)-u_in(:,:,2))*(1-Mach**2)/Lxyz(2)
        else
          L_3=0
          L_4=0
          L_5 = L_1
        endif
      else
!
!  Find the parameter determining 
!
        KK=nscbc_sigma_out*(1-Mach**2)*cs0/Lxyz(2)
!
!  Find the L_i's
!
        L_1 = KK*(rho0*cs20-p_infty)
        L_3 = f(l1:l2,lll,n1:n2,iuy)*dui_dxj(:,:,1,2)
        L_4 = f(l1:l2,lll,n1:n2,iuy)*dui_dxj(:,:,3,2)
        L_5 = (f(l1:l2,lll,n1:n2,iuy) - sgn*cs0_ar(l1:l2,n1:n2))*&
             (cs20_ar(l1:l2,n1:n2)*grad_rho(:,:,2)&
             -sgn*rho0*cs0_ar(l1:l2,n1:n2)*dui_dxj(:,:,2,2))
      end if
!
!  Add terms due to derivatives parallell to the boundary
!
!!$      parallell_term_rho &
!!$          =rho0*dui_dxj(:,:,1,1)+f(l1:l2,lll,n1:n2,iux)*grad_rho(:,:,1)&
!!$          +rho0*dui_dxj(:,:,3,3)+f(l1:l2,lll,n1:n2,iuz)*grad_rho(:,:,3)
!!$      parallell_term_ux &
!!$           =f(l1:l2,lll,n1:n2,iux)*dui_dxj(:,:,1,1)&
!!$           +f(l1:l2,lll,n1:n2,iuz)*dui_dxj(:,:,1,3)&
!!$           +cs20_ar(l1:l2,n1:n2)*grad_rho(:,:,1)/rho0&
!!$           +nu*(d2u1_dx2+d2u1_dz2)
!!$      parallell_term_uy &
!!$           =f(l1:l2,lll,n1:n2,iux)*dui_dxj(:,:,2,1)&
!!$           +f(l1:l2,lll,n1:n2,iuz)*dui_dxj(:,:,2,3)&
!!$           +nu*(d2u2_dx2+d2u2_dz2)
!!$      parallell_term_uz &
!!$           =f(l1:l2,lll,n1:n2,iux)*dui_dxj(:,:,3,1)&
!!$           +f(l1:l2,lll,n1:n2,iuz)*dui_dxj(:,:,3,3)&
!!$           +cs20_ar(l1:l2,n1:n2)*grad_rho(:,:,3)/rho0&
!!$           +nu*(d2u3_dx2+d2u3_dz2)

!print*,'y----------------------------------------------------'
!!$print*,'ux=',f(l1:l2,lll,n1:n2,iux)
!!$print*,'uz=',f(l1:l2,lll,n1:n2,iuz)
!!$print*,'dui_dxj(:,:,2,1)=',dui_dxj(:,:,2,1)
!!$print*,'dui_dxj(:,:,2,3)=',dui_dxj(:,:,2,3)
!!$print*,'parallell_term_uy=',parallell_term_uy
!!$print*,maxval(f(l1:l2,lll,n1:n2,iux)),maxval(dui_dxj(:,:,2,1))
!!$print*,maxval(f(l1:l2,lll,n1:n2,iuz)),maxval(dui_dxj(:,:,2,3))
!!$print*,minval(f(l1:l2,lll,n1:n2,iux)),minval(dui_dxj(:,:,2,1))
!!$print*,minval(f(l1:l2,lll,n1:n2,iuz)),minval(dui_dxj(:,:,2,3))
!!$print*,minval(parallell_term_rho),maxval(parallell_term_rho)
!!$print*,minval(parallell_term_ux),maxval(parallell_term_ux)
!!$print*,minval(parallell_term_uy),maxval(parallell_term_uy)
!!$print*,minval(parallell_term_uz),maxval(parallell_term_uz)

!!$
!!$
!!$
      parallell_term_rho=0
      parallell_term_ux=0
      parallell_term_uy=0
      parallell_term_uz=0
!
!  Find the evolution equations at the boundary
!
      select case(topbot)
      ! NB: For 'top' L_1 plays the role of L5 and L_5 the role of L1
      case('bot')
        df(l1:l2,lll,n1:n2,ilnrho) = prefac1*(L_5 + L_1)-parallell_term_rho
        if (llinlet) then
          df(l1:l2,lll,n1:n2,iuy) = prefac2*( L_5 - L_1)-parallell_term_uy
        else
          df(l1:l2,lll,n1:n2,iuy) = prefac2*(-L_5 + L_1)-parallell_term_uy
        endif
        df(l1:l2,lll,n1:n2,iux) = -L_3-parallell_term_ux
        df(l1:l2,lll,n1:n2,iuz) = -L_4-parallell_term_uz
      case('top')
        df(l1:l2,lll,n1:n2,ilnrho) = prefac1*(L_1 + L_5)-parallell_term_rho
        if (llinlet) then
          df(l1:l2,lll,n1:n2,iuy) = prefac2*(-L_5 + L_1)-parallell_term_uy
        else
          df(l1:l2,lll,n1:n2,iuy) = prefac2*( L_5 - L_1)-parallell_term_uy
        endif
        df(l1:l2,lll,n1:n2,iux) = -L_3-parallell_term_ux
        df(l1:l2,lll,n1:n2,iuz) = -L_4-parallell_term_uz
      endselect
!
!  Check if we are solving for logrho or rho
!
      if (.not. ldensity_nolog) then
        df(l1:l2,lll,n1:n2,irho)=df(l1:l2,lll,n1:n2,irho)/rho0
      endif
!
! Impose required variables at the boundary
!
      if (llinlet) then
        if (.not. non_reflecting_inlet) then
          f(l1:l2,lll,n1:n2,iux) = u_in(:,:,1)
          f(l1:l2,lll,n1:n2,iuy) = u_in(:,:,2)
          f(l1:l2,lll,n1:n2,iuz) = u_in(:,:,3)
        endif
      endif
!
    endsubroutine bc_nscbc_prf_y
!***********************************************************************
    subroutine read_NSCBC_init_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
      integer :: i 

      if (present(iostat)) then
        read(unit,NML=NSCBC_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=NSCBC_init_pars,ERR=99)
      endif
!
      do i=1,3
        if (nscbc_bc(i) /= '') lnscbc = .true.
      enddo
!
      if (lnscbc) call parse_nscbc(nscbc_bc,nscbc_bc1,nscbc_bc2)
!
99    return
    endsubroutine read_NSCBC_init_pars
!***********************************************************************
    subroutine read_NSCBC_run_pars(unit,iostat)

      use Sub, only : rdim

      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
      integer :: i,stat 
      logical :: exist
      character (len=130) :: file

      if (present(iostat)) then
        read(unit,NML=NSCBC_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=NSCBC_run_pars,ERR=99)
      endif
!
      do i=1,3
        if (nscbc_bc(i) /= '') lnscbc = .true.
      enddo
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
!print*,'turb_inlet_dir=',turb_inlet_dir
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
    endsubroutine read_NSCBC_run_pars
!***********************************************************************
    subroutine write_NSCBC_init_pars(unit)
      integer, intent(in) :: unit

      write(unit,NML=NSCBC_init_pars)

    endsubroutine write_NSCBC_init_pars
!***********************************************************************
    subroutine write_NSCBC_run_pars(unit)
      integer, intent(in) :: unit

      write(unit,NML=NSCBC_run_pars)

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
      end subroutine NSCBC_clean_up
!***********************************************************************
  end module NSCBC
