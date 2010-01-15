! $Id: internal_flow.f90 12795 2010-01-03 14:03:57Z ajohan@strw.leidenuniv.nl $

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lspecial = .true.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
!***************************************************************

module Special

  use Cdata
  use Cparam
  use Messages
  use Sub, only: keep_compiler_quiet

  implicit none

  include '../special.h'

  !
  ! Slice precalculation buffers
  !

  integer :: dummy
  character(len=24) :: initspecial='nothing'

!!  character, len(50) :: initcustom

! input parameters
  namelist /jet_init_pars/ &
      initspecial
  ! run parameters
  namelist /jet_run_pars/  &
       dummy

  contains

!***********************************************************************
    subroutine register_special()
!
!  Configure pre-initialised (i.e. before parameter read) variables
!  which should be know to be able to evaluate
!
!  6-oct-03/tony: coded
!
      if (lroot) call svn_id( &
           "$Id: internal_flow.f90 12795 2010-01-03 14:03:57Z ajohan@strw.leidenuniv.nl $")
!
    endsubroutine register_special
!***********************************************************************
    subroutine initialize_special(f)
!
!  called by run.f90 after reading parameters, but before the time loop
!
!  06-oct-03/tony: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
!!
!!  Initialize any module variables which are parameter dependent
!!
!
! DO NOTHING
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_special
!***********************************************************************
    subroutine init_special(f)
!
!  initialise special condition; called from start.f90
!  06-oct-2003/tony: coded
!
      use Mpicomm
      use Sub
      use Initcond
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: i,j,jjj,kkk
      real, dimension(3) :: velo,tmp
      real, dimension(2) :: radius,theta,jet_center
      real :: radius_mean, velocity_ratio,An,u_t,rad
      logical :: non_zero_transveral_velo
!
      intent(inout) :: f
      !
      ! Select case
      !
      select case (initspecial)
      case ('nothing'); if (lroot) print*,'init_special: nothing'
      case ('coaxial_jet')      
        u_t=5.
        velocity_ratio=3.3
        velo(1)=u_t
        velo(2)=velo(1)*velocity_ratio
        velo(3)=0.04*velo(2)
        radius(1)=0.0182
        radius(2)=radius(1)*2.
        radius_mean=(radius(1)+radius(2))/2.
        theta(1)=radius(1)/13.
        theta(2)=radius(2)/20.
        jet_center(1)=0
        jet_center(2)=0
!
! Set velocity profiles
!          
        do jjj=1,ny
          do kkk=1,nz
            rad=sqrt(&
                (y(jjj+m1-1)-jet_center(1))**2+&
                (z(kkk+n1-1)-jet_center(2))**2)
            ! Add mean velocity profile
            if (rad < radius_mean) then
              f(:,jjj+m1-1,kkk+n1-1,1)&
                  =(velo(1)+velo(2))/2&
                  +(velo(2)-velo(1))/2*tanh((rad-radius(1))/(2*theta(1)))
            else
              f(:,jjj+m1-1,kkk+n1-1,1)&
                  =(velo(2)+velo(3))/2&
                  +(velo(3)-velo(2))/2*tanh((rad-radius(2))/(2*theta(2)))
            endif
          enddo
        enddo
        f(:,:,:,iuy:iuz)=0        
      case default
        !
        !  Catch unknown values
        !
        if (lroot) print*,'init_special: No such value for initspecial: ', &
             trim(initspecial)
        call stop_it("")
      endselect
      !
    endsubroutine init_special
!***********************************************************************
    subroutine pencil_criteria_special()
!
!  All pencils that this special module depends on are specified here.
!
!  18-07-06/tony: coded
!
    endsubroutine pencil_criteria_special
!***********************************************************************
    subroutine pencil_interdep_special(lpencil_in)
!
!  Interdependency among pencils provided by this module are specified here.
!
!  18-07-06/tony: coded
!
      logical, dimension(npencils) :: lpencil_in
!
      call keep_compiler_quiet(lpencil_in)
!
    endsubroutine pencil_interdep_special
!***********************************************************************
    subroutine calc_pencils_special(f,p)
!
!  Calculate Hydro pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!   24-nov-04/tony: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      intent(in) :: f
      intent(inout) :: p
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)
!
    endsubroutine calc_pencils_special
!***********************************************************************
    subroutine dspecial_dt(f,df,p)
!
!  calculate right hand side of ONE OR MORE extra coupled PDEs
!  along the 'current' Pencil, i.e. f(l1:l2,m,n) where
!  m,n are global variables looped over in equ.f90
!
!  Due to the multi-step Runge Kutta timestepping used one MUST always
!  add to the present contents of the df array.  NEVER reset it to zero.
!
!  several precalculated Pencils of information are passed if for
!  efficiency.
!
!   06-oct-03/tony: coded
!
      use Diagnostics
      use Mpicomm
      use Sub
      use Deriv, only: der_pencil
      use Viscosity, only: getnu
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (3) :: meanx_oo
      real, dimension (3) :: meanx_uu
      real, dimension (nx,3) :: ufluct
      real, dimension (nx) :: ufluct2
      type (pencil_case) :: p
      integer :: i,j
      real, dimension (my) :: tmp,du_mean_dy
      real :: tau_tmp,nu
!
      intent(in) :: f,p
      intent(inout) :: df
!
!
    endsubroutine dspecial_dt
!***********************************************************************
    subroutine read_special_init_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat

 
      if (present(iostat)) then
        read(unit,NML=jet_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=jet_init_pars,ERR=99)
      endif

99    return
    endsubroutine read_special_init_pars
!***********************************************************************
    subroutine write_special_init_pars(unit)
!
      integer, intent(in) :: unit

      write(unit,NML=jet_init_pars)

    endsubroutine write_special_init_pars
!***********************************************************************
    subroutine read_special_run_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat

      if (present(iostat)) then
        read(unit,NML=jet_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=jet_run_pars,ERR=99)
      endif

99    return
    endsubroutine read_special_run_pars
!***********************************************************************
    subroutine write_special_run_pars(unit)
!
      integer, intent(in) :: unit

      write(unit,NML=jet_run_pars)

    endsubroutine write_special_run_pars
!***********************************************************************
    subroutine rprint_special(lreset,lwrite)
!
!  reads and registers print parameters relevant to special
!
!   06-oct-03/tony: coded
!
      use Diagnostics
      use Sub
!
      integer :: iname
      logical :: lreset,lwr
      logical, optional :: lwrite
!
!!$      lwr = .false.
!!$      if (present(lwrite)) lwr=lwrite
!!$!
!!$!  reset everything in case of reset
!!$!  (this needs to be consistent with what is defined above!)
!!$!
!!$      if (lreset) then
!!$        idiag_turbint=0
!!$        idiag_tau_w=0
!!$        idiag_uxm_central=0
!!$      endif
!!$!
!!$      do iname=1,nname
!!$        call parse_name(iname,cname(iname),cform(iname),'turbint',idiag_turbint)
!!$        call parse_name(iname,cname(iname),cform(iname),'tau_w',idiag_tau_w)
!!$        call parse_name(iname,cname(iname),cform(iname),'uxm_central',idiag_uxm_central)
!!$      enddo
!!$!
!!$!  write column where which magnetic variable is stored
!!$      if (lwr) then
!!$        write(3,*) 'i_turbint=',idiag_turbint
!!$        write(3,*) 'i_tau_w=',idiag_tau_w
!!$        write(3,*) 'i_uxm_central=',idiag_uxm_central
!!$      endif
!
    endsubroutine rprint_special
!***********************************************************************
    subroutine get_slices_special(f,slices)
!
!  Write slices for animation of special variables.
!
!  26-jun-06/tony: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
      !
      !  Loop over slices
      !
!!$      select case (trim(slices%name))
!!$        !
!!$        !  Vorticity (derived variable)
!!$        !
!!$      case ('oo_meanx')
!!$        if (slices%index == 3) then
!!$          slices%ready = .false.
!!$        else
!!$          slices%index = slices%index+1
!!$          slices%xy=>oo_xy_meanx(:,:,slices%index)
!!$          if (slices%index < 3) slices%ready = .true.
!!$        endif
!!$      case ('uu_meanx')
!!$        if (slices%index >= 3) then
!!$          slices%ready = .false.
!!$        else
!!$          slices%index = slices%index+1
!!$          slices%xy=uu_xy_meanx(:,:,slices%index)
!!$          if (slices%index < 3) slices%ready = .true.
!!$        endif
!!$      endselect
!
    endsubroutine get_slices_special
!***********************************************************************
    subroutine calc_lspecial_pars(f)
!
!  Mean flow velocitites
!
!  14-mar-08/nils: coded
!
      use Sub
      use Mpicomm, only: mpireduce_sum, mpibcast_real
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension(nygrid,3) :: mean_u_tmp
      real :: faq
      integer :: j,k
!!$!
!!$!  calculate mean of velocity in xz planes
!!$!
!!$      if (lvideo.and.lfirst .or. ldiagnos) then
!!$        mean_u_tmp=0
!!$        faq=nxgrid*nzgrid
!!$        do j=m1,m2
!!$          do k=1,3
!!$            mean_u_tmp(j+ny*ipy-nghost,k)=sum(f(l1:l2,j,n1:n2,k+iux-1))/faq
!!$          enddo
!!$        enddo        
!!$        do k=1,3
!!$          call mpireduce_sum(mean_u_tmp(:,k),mean_u(:,k),nygrid)
!!$          call mpibcast_real(mean_u(:,k),nygrid)
!!$        enddo
!!$      endif
!
    endsubroutine calc_lspecial_pars
!***********************************************************************
    subroutine special_calc_density(f,df,p)
!
!   calculate a additional 'special' term on the right hand side of the
!   entropy equation.
!
!   Some precalculated pencils of data are passed in for efficiency
!   others may be calculated directly from the f array
!
!   06-oct-03/tony: coded
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p

!!
!!  SAMPLE IMPLEMENTATION
!!     (remember one must ALWAYS add to df)
!!
!!
!!  df(l1:l2,m,n,ilnrho) = df(l1:l2,m,n,ilnrho) + SOME NEW TERM
!!
!!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(p)
!
    endsubroutine special_calc_density
!***********************************************************************
    subroutine special_calc_hydro(f,df,p)
!
!   calculate a additional 'special' term on the right hand side of the
!   entropy equation.
!
!   Some precalculated pencils of data are passed in for efficiency
!   others may be calculated directly from the f array
!
!   06-oct-03/tony: coded
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p

!!
!!  SAMPLE IMPLEMENTATION
!!     (remember one must ALWAYS add to df)
!!
!!
!NILS      if (m>=18.and.m<=20) then
!NILS        df(l1:l2,m,n,iux) = df(l1:l2,m,n,iux) - f(l1:l2,m,n,iux)*100
!NILS        df(l1:l2,m,n,iuy) = df(l1:l2,m,n,iuy) - f(l1:l2,m,n,iuy)*100
!NILS        df(l1:l2,m,n,iuz) = df(l1:l2,m,n,iuz) - f(l1:l2,m,n,iuz)*100
!NILS      endif
!!
!!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(p)
!
    endsubroutine special_calc_hydro
!***********************************************************************
    subroutine special_calc_magnetic(f,df,p)
!
!   calculate a additional 'special' term on the right hand side of the
!   entropy equation.
!
!   Some precalculated pencils of data are passed in for efficiency
!   others may be calculated directly from the f array
!
!   06-oct-03/tony: coded
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p

!!
!!  SAMPLE IMPLEMENTATION
!!     (remember one must ALWAYS add to df)
!!
!!
!!  df(l1:l2,m,n,iux) = df(l1:l2,m,n,iux) + SOME NEW TERM
!!  df(l1:l2,m,n,iuy) = df(l1:l2,m,n,iuy) + SOME NEW TERM
!!  df(l1:l2,m,n,iuz) = df(l1:l2,m,n,iuz) + SOME NEW TERM
!!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(p)
!
    endsubroutine special_calc_magnetic
!!***********************************************************************
    subroutine special_calc_entropy(f,df,p)
!
!   calculate a additional 'special' term on the right hand side of the
!   entropy equation.
!
!   Some precalculated pencils of data are passed in for efficiency
!   others may be calculated directly from the f array
!
!   06-oct-03/tony: coded
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p

!!
!!  SAMPLE IMPLEMENTATION
!!     (remember one must ALWAYS add to df)
!!
!!
!!  df(l1:l2,m,n,ient) = df(l1:l2,m,n,ient) + SOME NEW TERM
!!
!!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(p)
!
    endsubroutine special_calc_entropy
!***********************************************************************
    subroutine special_boundconds(f,bc)
!
!   Some precalculated pencils of data are passed in for efficiency
!   others may be calculated directly from the f array
!
!   2008-06-19/nils: coded
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      type (boundary_condition) :: bc
!!$!
!!$      select case (bc%bcname)
!!$      case ('poi')
!!$        select case (bc%location)
!!$        case (iBC_X_TOP)
!!$          call bc_poi_x(f,-1,'top',iux,REL=.true.,val=bc%value1)
!!$        case (iBC_X_BOT)
!!$          call bc_poi_x(f,-1,'bot',iux,REL=.true.,val=bc%value1)
!!$        end select
!!$      end select
!
    endsubroutine special_boundconds
!***********************************************************************
    subroutine special_before_boundary(f)
!
!   Possibility to modify the f array before the boundaries are
!   communicated.
!
!   Some precalculated pencils of data are passed in for efficiency
!   others may be calculated directly from the f array
!
!   06-jul-06/tony: coded
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine special_before_boundary
!***********************************************************************

!********************************************************************
!************        DO NOT DELETE THE FOLLOWING       **************
!********************************************************************
!**  This is an automatically generated include file that creates  **
!**  copies dummy routines from nospecial.f90 for any Special      **
!**  routines not implemented in this file                         **
!**                                                                **
    include '../special_dummies.inc'
!********************************************************************

endmodule Special

