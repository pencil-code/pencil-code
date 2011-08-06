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
  logical :: first_time=.true. 
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
  integer :: nghost_in
  integer :: m1_in  
  integer :: n1_in
  integer :: l2_in
  integer :: m2_in
  integer :: n2_in
  real :: Lx_in
  real :: Ly_in
  real :: Lz_in
  character (len=60) :: turbfile
  character(len=40) :: turb_inlet_dir='' 
  logical :: proc_at_inlet
  integer :: ipx_in, ipy_in, ipz_in, iproc_in, nprocx_in, nprocy_in, nprocz_in
  character (len=120) :: directory_in
  character (len=5) :: chproc_in
  real, dimension(2) :: radius=(/0.0182,0.0364/)
  real, dimension(2) :: momentum_thickness=(/0.014,0.0182/)
  real, dimension(2) :: jet_center=(/0.,0./)
  real :: u_t=5.,velocity_ratio=3.3




!!  character, len(50) :: initcustom

! input parameters
  namelist /jet_init_pars/ &
      initspecial,turb_inlet_dir,u_t,velocity_ratio,radius,momentum_thickness,&
      jet_center
  ! run parameters
  namelist /jet_run_pars/  &
       turb_inlet_dir

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
    subroutine initialize_special(f,lstarting)
!
!  called by run.f90 after reading parameters, but before the time loop
!
!  19-jan-10/nils: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      logical :: lstarting
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(lstarting)
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
      real :: radius_mean,An,rad
      logical :: non_zero_transveral_velo
!
      intent(inout) :: f
      !
      ! Select case
      !
      select case (initspecial)
      case ('nothing'); if (lroot) print*,'init_special: nothing'
      case ('coaxial_jet')      
        if (lroot) print*,'init_special: coaxial_jet'
        velo(1)=u_t
        velo(2)=velo(1)*velocity_ratio
        velo(3)=0.04*velo(2)
        radius_mean=(radius(1)+radius(2))/2.
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
                  +(velo(2)-velo(1))/2*tanh((rad-radius(1))/&
                  (2*momentum_thickness(1)))
            else
              f(:,jjj+m1-1,kkk+n1-1,1)&
                  =(velo(2)+velo(3))/2&
                  +(velo(3)-velo(2))/2*tanh((rad-radius(2))/&
                  (2*momentum_thickness(2)))
            endif
          enddo
        enddo
        f(:,:,:,iuy:iuz)=0   

        case ('single_jet')
          if (lroot) print*,'init_special: single_jet'
          velo(1)=u_t
          velo(2)=velo(1)/velocity_ratio
!
! Set velocity profiles
!
          do jjj=1,ny 
            do kkk=1,nz
              rad=sqrt(&
                  (y(jjj+m1-1)-jet_center(1))**2+&
                  (z(kkk+n1-1)-jet_center(1))**2)
              !Add velocity profile
              f(:,jjj+m1-1,kkk+n1-1,1)&
                  =velo(1)*(1-tanh((rad-radius(1))/&
                  momentum_thickness(1)))*0.5+velo(2)
            enddo
          enddo
          f(:,:,:,iuy:iuz)=0

        case ('single_laminar_wall_jet')
          if (lroot) print*,'init_special: single_laminar_wall_jet'
!
! Set velocity profiles
!
!!$          do jjj=1,ny 
!!$            do kkk=1,nz
!!$              rad=sqrt(&
!!$                  (y(jjj+m1-1)-jet_center(1))**2+&
!!$                  (z(kkk+n1-1)-jet_center(1))**2)
!!$              !Add velocity profile
!!$              if (rad < radius(1)) then 
!!$                f(:,jjj+m1-1,kkk+n1-1,1)=u_t*(1-(rad/radius(1))**2)
!!$              else
!!$                f(:,jjj+m1-1,kkk+n1-1,1)=0.
!!$              endif
!!$            enddo
!!$          enddo



          do jjj=1,my 
            do kkk=1,mz
              rad=sqrt(&
                  (y(jjj)-jet_center(1))**2+&
                  (z(kkk)-jet_center(1))**2)
              !Add velocity profile
              if (rad < radius(1)) then 
                f(:,jjj,kkk,1)=u_t*(1-(rad/radius(1))**2)
              else
                f(:,jjj,kkk,1)=0.
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
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      type (boundary_condition) :: bc
      character (len=3) :: topbot
!


      topbot='top'
      if (bc%location==-1) topbot='bot'

      select case (bc%bcname)
      case ('tur')
        select case (bc%location)
        case (iBC_X_TOP)
          call bc_turb(f,bc%value1,'top',1,bc%ivar)
        case (iBC_X_BOT)
          call bc_turb(f,bc%value1,'bot',1,bc%ivar)
          bc%done=.true.
        end select
      case ('wi')
        call bc_wi_x(f,+1,topbot,bc%ivar,val=bc%value1)
        bc%done=.true.
      case ('wip')
        call bc_wip_x(f,+1,topbot,bc%ivar,val=bc%value1)
        bc%done=.true.
!!$      case ('wo')
!!$        call bc_wo_x(f,+1,topbot,bc%ivar,val=bc%value1)
!!$        bc%done=.true.
      end select
!
    endsubroutine special_boundconds
!***********************************************************************
    subroutine bc_turb(f,u_t,topbot,j,ivar)
!
! Use a prerun simulation as inlet condition
!
! 2010-10-15/nils: coded
!      
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f 
      real :: shift, grid_shift, weight, round
      integer :: iround,lower,upper,ii,j,imin,imax,ivar
      character (len=3) :: topbot
      real :: T_t,u_t
!
! Read data from file only initially
! At later times this is stored in processor memory.
!
      if (first_time) then
        call read_turbulent_data(topbot,j)
      end if
!
! Set all ghost points at the same time
!
      if (ivar==iux) then 
!
! Set which ghost points to update
!
        if (topbot=='bot') then
          imin=1
          imax=3
        else
          imin=l1+1
          imax=mx
        endif
!
! Set data at ghost points
!
        if (Lx_in == 0) call fatal_error('bc_nscbc_prf_x',&
            'Lx_in=0. Check that the precisions are the same.')
        round=t*u_t/Lx_in
        iround=int(round)
        shift=round-iround
        grid_shift=shift*nx_in
        lower=l1_in+int(grid_shift)
        upper=lower+1
        weight=grid_shift-int(grid_shift)
        f(imin:imax,m1:m2,n1:n2,iux:iuz)&
            =f_in(lower-3:lower-1,m1_in:m2_in,n1_in:n2_in,iux:iuy)*(1-weight)&
            +f_in(upper-3:upper-1,m1_in:m2_in,n1_in:n2_in,iux:iuy)*weight
        f(imin:imax,m1:m2,n1:n2,ilnrho)&
            =f_in(lower-3:lower-1,m1_in:m2_in,n1_in:n2_in,ilnrho)*(1-weight)&
            +f_in(upper-3:upper-1,m1_in:m2_in,n1_in:n2_in,ilnrho)*weight
!
! Add mean flow velocity on top of the turubulence
!      
        f(imin:imax,m1:m2,n1:n2,iux)=f(imin:imax,m1:m2,n1:n2,iux)+u_t
      endif
!
    end subroutine bc_turb
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
    subroutine read_turbulent_data(topbot,j)
!
!   This subroutine will read data from pre-run isotropic box
!
!   19-jan-10/nils: coded
!
      use General, only: safe_character_assign, chn
      use Sub, only : rdim
!
      character (len=3), intent(in) :: topbot
      integer, intent(in) :: j
      integer :: i,stat 
      logical :: exist
      character (len=fnlen) :: file
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
        call fatal_error('bc_turb','Could not find file!')
      endif
!
! Allocate array for data to be used at the inlet.
! For now every processor reads all the data - this is clearly an overkill,
! but we leav it like this during the development of this feature.
!
      allocate( f_in(mx_in,my_in,mz_in,nv_in),STAT=stat)
      if (stat>0) call fatal_error('bc_turb',&
          "Couldn't allocate memory for f_in ")
      allocate( x_in(mx_in),STAT=stat)
      if (stat>0) call fatal_error('bc_turb',&
          "Couldn't allocate memory for x_in ")
      allocate( y_in(my_in),STAT=stat)
      if (stat>0) call fatal_error('bc_turb',&
          "Couldn't allocate memory for y_in ")
      allocate( z_in(mz_in),STAT=stat)
      if (stat>0) call fatal_error('bc_turb',&
          "Couldn't allocate memory for z_in ")
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
          if ((topbot=='bot'.and.lfirst_proc_x).or.&
              (topbot=='top'.and.llast_proc_x)) then
            proc_at_inlet=.true.
            ipx_in=0
            nprocx_in=1
          endif
        elseif (j==2) then
          if ((topbot=='bot'.and.lfirst_proc_y).or.&
              (topbot=='top'.and.llast_proc_y)) then
            proc_at_inlet=.true.
            ipy_in=0
            nprocy_in=1
          endif
        elseif (j==3) then
          if ((topbot=='bot'.and.lfirst_proc_z).or.&
              (topbot=='top'.and.llast_proc_z)) then
            proc_at_inlet=.true.
            ipz_in=0
            nprocz_in=1
          endif
        else
          call fatal_error("bc_turb_x",'No such direction!')
        endif
        iproc_in=ipz_in*nprocy_in*nprocx_in+ipy_in*nprocx_in+ipx_in
      else
        call fatal_error("bc_turb_x",&
            'lprocz_slowest=F not implemeted for inlet from file!')
      endif
!
!  Read data only if required, i.e. if we are at a processor handling inlets
!
      if (proc_at_inlet) then
        print*,'datadir=',datadir
        call chn(iproc_in,chproc_in)
        call safe_character_assign(directory_in,&
            trim(turb_inlet_dir)//'/data/proc'//chproc_in)
        print*,'directory_in=',directory_in
        call safe_character_assign(turbfile,&
            trim(directory_in)//'/var.dat')
        print*,'turbfile=',turbfile
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
        first_time=.false.
        close(1)
      endif
!
    end subroutine read_turbulent_data
!***********************************************************************
    subroutine bc_wip_x(f,sgn,topbot,j,rel,val)
!
!
!  23-may-10/nils+marianne: coded
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      real :: val
      integer :: sgn,i,j
      logical, optional :: rel
      logical :: relative

      real :: radius, rad,z0,y0,yp0,zp0,dyp,minrad,distance,yp1
      integer :: jj,kk,nrad,irad,nb_lines,nb_colums,iz,iy
      real, allocatable, dimension(:,:,:) :: center
!
! Allocate array
!
      dyp=0.01
      radius=0.025
      distance=(dyp+radius)*2.
      yp0=xyz0(2)+radius*1.1
      yp1=xyz1(2)-radius*1.1
      nb_lines =floor((yp1-yp0)/(2*radius+dyp))
      nb_colums=max(floor(Lxyz(3)/(2*radius+dyp)),1)
      allocate(center(nb_lines,nb_colums,2))
!
      
      if (Lxyz(3) > 2*radius) then
        zp0=xyz0(3)+radius*1.3
      else
        zp0=xyz0(3)
      endif
      nrad=Lxyz(2)/(2*radius+dyp)
      z0=0.
!
! Find the center position of all the holes and put in array
!
      do iy=1,Nb_lines
        do iz=1,Nb_colums
          if (mod(iz,2)==0) then
            center(iy,iz,1)=yp0-distance/sqrt(2.)+(iy-1)*distance
          else
            center(iy,iz,1)=yp0+distance*(iy-1)                  
          endif
          center(iy,iz,2)=zp0-distance/sqrt(2.)*(iz-1)
        enddo
      enddo
!
!  Loop over all grid points
!
        do jj=1,my
          do kk=1,mz
            minrad=impossible
            do iy=1,Nb_lines
              do iz=1,Nb_colums
                rad=sqrt((y(jj)-center(iy,iz,1))**2+(z(kk)-center(iy,iz,2))**2)
                if (rad < minrad) minrad=rad
              enddo
            enddo
!
!  Zero derivative for density
!
            if (j == ilnrho) then
              do i=1,nghost
                f(l1-i,jj,kk,j)= f(l1+i,jj,kk,j)
              enddo
!
!  Constant temperature
!
            elseif (j == ilnTT) then
              f(l1,jj,kk,j)=val
              do i=1,nghost
                f(l1-i,jj,kk,j)=2*val-f(l1+i,jj,kk,j)
              enddo
            else
!
! Check if we are inside the radius if the inlet
!
              if (minrad > radius) then
            ! Solid wall
                do i=1,nghost
                  if (j <= iuz) then
                    f(l1,jj,kk,j)=0.
                    f(l1-i,jj,kk,j)=-f(l1+i,jj,kk,j)
                  else
                    f(l1-i,jj,kk,j)= f(l1+i,jj,kk,j)
                  endif
                enddo
              else
              ! Inlet           
                if (j == iux) then
                  f(l1,jj,kk,j)=val*(1-(minrad/radius)**2)
                  do i=1,nghost
                    f(l1-i,jj,kk,j)=2*val*(1-(minrad/radius)**2)-f(l1+i,jj,kk,j)
                  enddo
                elseif (j == iuy .or. j == iuz) then
                  f(l1,jj,kk,j)=0.
                  do i=1,nghost
                    f(l1-i,jj,kk,j)=-f(l1+i,jj,kk,j)
                  enddo
                else
                  do i=1,nghost
                    f(l1-i,jj,kk,j)=2*val-f(l1+i,jj,kk,j)
                  enddo
                endif
              endif
            endif
         enddo
        enddo
!
  endsubroutine bc_wip_x
!***********************************************************************
    subroutine bc_wi_x(f,sgn,topbot,j,rel,val)
!
!
!  23-may-10/nils+marianne: coded
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      real  :: val
      integer :: sgn,i,j
      logical, optional :: rel
      logical :: relative
      real :: radius, rad,z0,y0
      integer :: jj,kk
!
!
!
      radius=0.015
      z0=0.
      y0=0.
!
!  Loop over all grid points
!
        do jj=1,my
          do kk=1,mz
            rad=sqrt((y(jj)-y0)**2+(z(kk)-z0)**2)
!
!  Zero derivative for density
!
            if (j == ilnrho) then
              do i=1,nghost
                f(l1-i,jj,kk,j)= f(l1+i,jj,kk,j)
              enddo
!
!  Constant temperature
!
            elseif (j == ilnTT) then
              f(l1,jj,kk,j)=val
              do i=1,nghost
                f(l1-i,jj,kk,j)=2*val-f(l1+i,jj,kk,j)
              enddo
            else
!
! Check if we are inside the radius if the inlet
!
              if (rad > radius) then
            ! Solid wall
                do i=1,nghost
                  if (j <= iuz) then
                    f(l1,jj,kk,j)=0.
                    f(l1-i,jj,kk,j)=-f(l1+i,jj,kk,j)
                  else
                    f(l1-i,jj,kk,j)= f(l1+i,jj,kk,j)
                  endif
                enddo
              else
              ! Inlet           
                if (j == iux) then
                  f(l1,jj,kk,j)=val*(1-(rad/radius)**2)
                  do i=1,nghost
                    f(l1-i,jj,kk,j)=2*val*(1-(rad/radius)**2)-f(l1+i,jj,kk,j)
                  enddo
                elseif (j == iuy .or. j == iuz) then
                  f(l1,jj,kk,j)=0.
                  do i=1,nghost
                    f(l1-i,jj,kk,j)=-f(l1+i,jj,kk,j)
                  enddo
                else
                  do i=1,nghost
                    f(l1-i,jj,kk,j)=2*val-f(l1+i,jj,kk,j)
                  enddo
                endif
              endif
            endif
         enddo
        enddo
!
  endsubroutine bc_wi_x
!!$!***********************************************************************
!!$    subroutine bc_wo_x(f,sgn,topbot,j,rel,val)
!!$!
!!$!  23-may-10/nils: coded
!!$!
!!$      use EquationOfState
!!$      use chemistry
!!$!
!!$      character (len=3) :: topbot
!!$      real, dimension (mx,my,mz,mfarray) :: f
!!$      real :: val
!!$      integer :: sgn,i,j
!!$      logical, optional :: rel
!!$      logical :: relative
!!$
!!$      real :: radius, rad,z0,y0,T0,P0,r,mu
!!$      integer :: jj,kk
!!$!
!!$!
!!$!
!!$      radius=0.0
!!$      z0=0.
!!$      y0=0.
!!$!
!!$! First of all we set one-sided derivatives for all variables over the full
!!$! boundary. This will then mostly be overwritten later.
!!$
!!$!NILS: COmmented out the below statements for now in order to make the 
!!$!NILS: code compile after I moved this from the boundary module.
!!$!NILS: will have to duplicate whatever bc i need here since this 
!!$!NILS: module is used by the boundary module.....
!!$!
!!$!      call bc_onesided_x(f,topbot,j)
!!$!      call bc_extrap_2_1(f,topbot,j)
!!$!
!!$!  Loop over all grid points
!!$!
!!$        do jj=1,my
!!$          do kk=1,mz
!!$            rad=sqrt((y(jj)-y0)**2+(z(kk)-z0)**2)
!!$!
!!$!  Zero derivative for density
!!$!
!!$            if (j == ilnrho) then
!!$              do i=1,nghost
!!$                f(l2+i,jj,kk,j)= f(l2-i,jj,kk,j)
!!$              enddo
!!$            else
!!$!
!!$! Check if we are inside the radius if the inlet
!!$!
!!$              if (rad > radius) then
!!$                ! Solid wall
!!$                do i=1,nghost
!!$                  if (j <= iuz) then
!!$                    f(l2,jj,kk,j)=0.
!!$                    f(l2+i,jj,kk,j)=-f(l2-i,jj,kk,j)
!!$                  elseif (j == ilnTT) then
!!$                    f(l2,jj,kk,j)=val
!!$                    f(l2+i,jj,kk,j)=2*val-f(l2-i,jj,kk,j)
!!$                  else
!!$                    f(l2+i,jj,kk,j)= f(l2-i,jj,kk,j)
!!$                  endif
!!$                enddo
!!$              else
!!$                ! Outlet           
!!$                if (j == ilnTT) then
!!$                  call getmu(f,mu,l2,jj,kk)
!!$                  r=Rgas*mu
!!$                  P0=1.013e6
!!$                  T0=min(P0/(exp(f(l2,jj,kk,ilnrho))*r),2500.)
!!$                 f(l2,jj,kk,j)=log(T0)
!!$                  do i=1,nghost
!!$                    f(l2+i,jj,kk,j)=2*log(T0)-f(l2-i,jj,kk,j)
!!$                  enddo
!!$                else
!!$                  ! Do nothing because one sided conditions has already
!!$                  ! been set in the top of the rutine
!!$                endif
!!$              endif
!!$            endif
!!$         enddo
!!$        enddo
!!$!
!!$  endsubroutine bc_wo_x
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

