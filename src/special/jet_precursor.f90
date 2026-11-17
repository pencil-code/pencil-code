! $Id$
!
!  Special module for a z-periodic ("temporal") jet precursor: the z-averaged
!  velocity is relaxed towards a tabulated axisymmetric profile - the axial
!  component to U(r) read from jet_table_file, the lateral components to zero -
!  and the velocity is damped beyond a cylindrical radius (a sponge).  With
!  periodic boundaries in all directions the turbulence of the jet's shear
!  layer develops freely while the mean profile is held, which gives a
!  statistically steady jet cross-section; the time series of its middle
!  plane (slices) is the inlet of a spatially developing jet run with
!  NSCBC turb_inlet_mode='slices'.
!
!  Forcing terms added to du/dt (per pencil):
!    - (<u_z>_z(x,y) - U(r))/tau_relax         r = sqrt((x-xc)^2 + (y-yc)^2)
!    - <u_x>_z/tau_relax,  - <u_y>_z/tau_relax  (if lrelax_lateral)
!    - dampfac * step(r, rdamp, wdamp) * u       (if dampfac > 0)
!  where <.>_z is the average over z (the periodic direction), computed here
!  from the f array every sub-step (special_after_boundary), so the module does
!  not depend on hydro's lcalc_uumeanxy.
!
!  jet_table_file: first line n, then n lines  r  U(r)  (code units); linear
!  interpolation in r, constant outside the table.
!
!  Not to be confused with special/jet.f90, which sets up a spatially
!  developing coaxial jet (initial condition and inlet boundary conditions
!  with turbulence from a file, used by samples/2d-tests/cellular-flame).
!  This module has no inlet: the jet is periodic in z and only forced.
!
!  11-nov-26/nils: moved here from hydro.f90 (uuprof='jet_table', ldampuext_cyl)
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of special_dummies.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lspecial = .true.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
!***************************************************************
!
module Special
!
  use Cdata
  use Quiet
  use Sub, only: step, finalize_aver
  use Messages, only: fatal_error, not_implemented
!
  implicit none
!
  include '../special.h'
!
  character (len=fnlen) :: jet_table_file='jet_profile.dat'
  real, dimension(2) :: jet_table_center=(/0.,0./)
  real :: tau_relax=0.0            ! relaxation time of the z-mean; <= 0: no relaxation
  logical :: lrelax_lateral=.true. ! relax the z-mean of u_x, u_y to zero as well
  real :: dampfac=0.0              ! damping rate of the sponge; <= 0: no sponge
  real :: rdamp=impossible, wdamp=0.0   ! sponge beyond rdamp, step width wdamp
!
  real, dimension(:,:), allocatable :: uz_target_xy      ! U(r) on the (x,y) plane (nx,my)
  real, dimension(:,:,:), allocatable :: uumxy_loc       ! z-mean of u (mx,my,3)
!
  namelist /special_run_pars/ &
      jet_table_file, jet_table_center, tau_relax, lrelax_lateral, dampfac, rdamp, wdamp
!
  contains
!***********************************************************************
    subroutine initialize_special(f)
!
      real, contiguous, dimension(:,:,:,:) :: f
!
      if (lyang) call not_implemented('initialize_special','jet_precursor on the Yin-Yang grid')
      if (tau_relax>0.0) then
        if (.not.lequidist(3)) call fatal_error('initialize_special', &
            'the z-average of jet_precursor requires an equidistant z grid')
        if (lrun) call read_jet_table
        if (.not.allocated(uumxy_loc)) allocate(uumxy_loc(mx,my,3))
        uumxy_loc=0.0
      endif
      if (dampfac>0.0 .and. rdamp==impossible) call fatal_error('initialize_special', &
          'dampfac > 0 needs rdamp')
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_special
!***********************************************************************
    subroutine pencil_criteria_special
!
      if (dampfac>0.0) then
        lpenc_requested(i_uu)=.true.
        lpenc_requested(i_rcyl_mn)=.true.
      endif
!
    endsubroutine pencil_criteria_special
!***********************************************************************
    subroutine read_special_run_pars(iomsg)
!
      use File_io, only: parallel_unit
!
      character(len=*), intent(out) :: iomsg
      integer :: iostat
!
      read(parallel_unit, NML=special_run_pars, IOSTAT=iostat, IOMSG=iomsg)
      if (iostat==0) iomsg=""
!
    endsubroutine read_special_run_pars
!***********************************************************************
    subroutine write_special_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=special_run_pars)
!
    endsubroutine write_special_run_pars
!***********************************************************************
    subroutine read_jet_table
!
!  Read the axisymmetric target profile U(r) and put it on the local (x,y) plane.
!
      use File_io, only: file_exists
!
      integer, parameter :: lun=94
      integer :: ntab, i, l, m
      real, dimension(:), allocatable :: rtab, utab
      real :: rad, w
!
      if (.not.file_exists(jet_table_file)) &
          call fatal_error('read_jet_table','no such file: '//trim(jet_table_file))
      open(lun,file=jet_table_file,status='old',action='read')
      read(lun,*) ntab
      allocate(rtab(ntab),utab(ntab))
      do i=1,ntab
        read(lun,*) rtab(i), utab(i)
      enddo
      close(lun)
!
      if (.not.allocated(uz_target_xy)) allocate(uz_target_xy(nx,my))
      do m=1,my
        do l=l1,l2
          rad=sqrt((x(l)-jet_table_center(1))**2+(y(m)-jet_table_center(2))**2)
          if (rad<=rtab(1)) then
            uz_target_xy(l-l1+1,m)=utab(1)
          elseif (rad>=rtab(ntab)) then
            uz_target_xy(l-l1+1,m)=utab(ntab)
          else
            i=1
            do while (rtab(i+1)<rad)
              i=i+1
            enddo
            w=(rad-rtab(i))/(rtab(i+1)-rtab(i))
            uz_target_xy(l-l1+1,m)=(1.0-w)*utab(i)+w*utab(i+1)
          endif
        enddo
      enddo
      if (lroot) print*,'read_jet_table: ',ntab,' rows from '//trim(jet_table_file)// &
          ', U(r=0) = ',utab(1),', U(r_max) = ',utab(ntab),' (r_max = ',rtab(ntab),')'
      deallocate(rtab,utab)
!
    endsubroutine read_jet_table
!***********************************************************************
    subroutine special_after_boundary(f)
!
!  The z-average of the velocity, every sub-step (as hydro's lcalc_uumeanxy).
!
      real, contiguous, dimension(:,:,:,:) :: f
!
      integer :: j
!
      if (tau_relax<=0.0) return
      do j=1,3
        uumxy_loc(:,:,j)=sum(f(:,:,n1:n2,iuu+j-1),3)/nzgrid_eff
      enddo
      call finalize_aver(nprocz,3,uumxy_loc)
!
    endsubroutine special_after_boundary
!***********************************************************************
    subroutine special_calc_hydro(f,df,p)
!
      real, contiguous, dimension(:,:,:,:) :: f
      real, contiguous, dimension(:,:,:,:) :: df
      type(pencil_case), intent(in) :: p
!
      real, dimension(nx) :: pdamp
      real :: tau1
      integer :: i
!
      if (tau_relax>0.0) then
        tau1=1.0/tau_relax
        df(l1:l2,m,n,iuz)=df(l1:l2,m,n,iuz)-tau1*(uumxy_loc(l1:l2,m,3)-uz_target_xy(:,m))
        if (lrelax_lateral) then
          df(l1:l2,m,n,iux)=df(l1:l2,m,n,iux)-tau1*uumxy_loc(l1:l2,m,1)
          df(l1:l2,m,n,iuy)=df(l1:l2,m,n,iuy)-tau1*uumxy_loc(l1:l2,m,2)
        endif
      endif
!
      if (dampfac>0.0) then
        pdamp=step(p%rcyl_mn,rdamp,wdamp)
        do i=1,3
          df(l1:l2,m,n,iux+i-1)=df(l1:l2,m,n,iux+i-1)-dampfac*pdamp*p%uu(:,i)
        enddo
      endif
!
      call keep_compiler_quiet(f)
!
    endsubroutine special_calc_hydro
!***********************************************************************
!************        DO NOT DELETE THE FOLLOWING       **************
!********************************************************************
!**  This is an automatically generated include file that creates  **
!**  copies dummy routines from nospecial.f90 for any Special      **
!**  routines not implemented in this file                         **
!**                                                                **
    include '../special_dummies.inc'
!********************************************************************
endmodule Special
