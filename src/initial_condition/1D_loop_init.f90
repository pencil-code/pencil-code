! $Id$
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: linitial_condition = .true.
!
!***************************************************************
module InitialCondition
!
  use Cdata
  use Messages
  use Sub, only: keep_compiler_quiet
!
  implicit none
!
  include '../initial_condition.h'
!
  character (len=labellen) :: stratitype='nothing'
  character (len=labellen) :: lnrho_init='nothing'
  character (len=labellen) :: lnTT_init='nothing'
  character (len=labellen) :: aa_init='nothing'
  logical :: set_lnTT_first =.true.
  real :: rho_init=0.
  real :: T0=6000.,T1=1e6,z0_tanh=4e6,width_tanh=1e6
  character (len=labellen) :: direction='z'
  real, dimension(4) :: mpoly_special = (/1.3,1000.,-1.04,500./)
  real, dimension(3) :: zpoly = (/0.,3.,5./)
!
  namelist /initial_condition_pars/ &
      lnrho_init,lnTT_init,stratitype,rho_init,direction, &
      set_lnTT_first,T0,T1,z0_tanh,width_tanh,mpoly_special,zpoly
!
contains
!***********************************************************************
    subroutine register_initial_condition()
!
!  Register variables associated with this module; likely none.
!
!  04-sep-10/bing: coded
!
      if (lroot) call svn_id( &
          "$Id$")
!
    endsubroutine register_initial_condition
!***********************************************************************
    subroutine initialize_initial_condition(f)
!
!  Initialize any module variables which are parameter dependent.
!
!  14-dec-10/bing: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      if (iproc==0) then
        write(*,*) "-------------------------------------------------------------"
        write(*,*) "Parameters to be set in run.in:"
        write(*,*)
        write(*,*) "Kpara=",2e-11 /unit_density/unit_velocity**3./ &
            unit_length*unit_temperature**3.5
!
        write(*,*) "hcond_grad=",1e9*dxmax**3*unit_temperature/unit_velocity**3
        write(*,*) "-------------------------------------------------------------"
      endif
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_initial_condition
!***********************************************************************
    subroutine read_initial_condition_pars(unit,iostat)
!
!  04-sep-10/bing: coded
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=initial_condition_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=initial_condition_pars,ERR=99)
      endif
!
 99  return
!
    endsubroutine read_initial_condition_pars
!***********************************************************************
    subroutine write_initial_condition_pars(unit)
!
!  04-sep-10/bing: coded
!
      integer, intent(in) :: unit
!
      write(unit,NML=initial_condition_pars)
!
    endsubroutine write_initial_condition_pars
!***********************************************************************
    subroutine initial_condition_lnrho(f)
!
!  Initialize logarithmic density.
!
!  04-sep-10/bing: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      select case (lnTT_init)
      case ('nothing')
        ! do nothing
      case ('prof_lnTT_z','prof_lnTT_loop')
        call setup_profiles(f)
      case ('tanh_z','tanh_loop')
        call setup_tanh(f)
      case default
        call fatal_error('initial_condition', &
            'no such value for lnTT_init')
      endselect
!
      select case (lnrho_init)
      case ('nothing')
        ! do nothing
      case ('prof_lnrho_z','prof_lnrho_loop')
        call setup_profiles(f)
      case ('hydrostatic')
        call hydrostatic_lnTT(f)
      case default
        call fatal_error('initial_condition_lnrho', &
            'no such value for lnTT_init')
      endselect
!
      call write_stratification_dat(f)
!
    endsubroutine initial_condition_lnrho
!***********************************************************************
    subroutine setup_profiles(f)
!
!  Read and set vertical profiles for initial temperature and density.
!  Initial temperature profile is given in ln(T) [K] over z [Mm]
!  Initial density profile is given in ln(rho) [kg/m^3] over z [Mm]
!
!  04-sep-10/bing: coded
!
      use EquationOfState, only: get_cp1,gamma,gamma_m1,cs20
      use Mpicomm, only: mpibcast_int, mpibcast_real, stop_it_if_any
      use Messages, only: warning
      use Syscalls, only: file_exists, file_size
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      real :: cp1=1.,ztmp
      integer :: lend,lend_b8,ierr
      integer :: i,j
      integer, parameter :: unit=12
!
! file location settings
      character (len=*), parameter :: lnrho_dat = 'prof_lnrho.dat'
      character (len=*), parameter :: lnT_dat = 'prof_lnT.dat'
!
      integer :: prof_nz
      real, dimension (:), allocatable :: prof_z, prof_lnrho, prof_lnTT
      real, dimension (mx) :: profile_x
!
      logical :: lread_lnrho=.false., lread_lnTT=.false.
!
      inquire(IOLENGTH=lend) 1.0
      inquire(IOLENGTH=lend_b8) 1.0d0
!
      lread_lnTT=((lnTT_init=='prof_lnTT_z').or.(lnTT_init=='prof_lnTT_loop'))
      lread_lnrho=((lnrho_init=='prof_lnrho_z').or.(lnrho_init=='prof_lnrho_loop'))
!
! read density profile for interpolation
      if (lread_lnrho) then
!
! file access is only done on the MPI root rank
        if (lroot) then
          if (.not. file_exists (lnrho_dat)) call stop_it_if_any ( &
              .true., 'setup_special: file not found: '//trim(lnrho_dat))
! find out, how many data points our profile file has
          prof_nz = (file_size (lnrho_dat) - 2*2*4) / (lend*4 * 2)
        endif
        call stop_it_if_any(.false.,'')
        call mpibcast_int (prof_nz,1)
!
        allocate (prof_z(prof_nz), prof_lnrho(prof_nz), stat=ierr)
!
        if (lroot) then
          open (unit,file=lnrho_dat,form='unformatted',status='unknown', &
              recl=lend*prof_nz)
          read (unit,iostat=ierr) prof_lnrho
          read (unit,iostat=ierr) prof_z
          if (ierr /= 0) call stop_it_if_any(.true.,'setup_special: '// &
              'Error reading stratification file: "'//trim(lnrho_dat)//'"')
          close (unit)
        endif
        call stop_it_if_any(.false.,'')
!
        call mpibcast_real (prof_lnrho,prof_nz)
        call mpibcast_real (prof_z,prof_nz)
!
! convert from logarithmic SI to Pencil units
        prof_lnrho = prof_lnrho - alog(real(unit_density))
!
! convert z coordinates from [Mm] to Pencil units
        prof_z = prof_z / unit_length
!
! interpolate density profile to Pencil grid
        if (lnrho_init=='prof_lnrho_z') then
          do j = l1-nghost, l2+nghost
            ztmp = Lxyz(1) / pi *sin( x(j)/Lxyz(1)*pi )
            if (ztmp < prof_z(1) ) then
              profile_x(j) = prof_lnrho(1)
            elseif (z(j) >= prof_z(prof_nz)) then
              profile_x(j) = prof_lnrho(prof_nz)
            else
              do i = 1, prof_nz-1
                if ((ztmp >= prof_z(i)) .and. (ztmp < prof_z(i+1))) then
                  ! linear interpolation: y = m*(x-x1) + y1
                  profile_x(j) = (prof_lnrho(i+1)-prof_lnrho(i)) / &
                      (prof_z(i+1)-prof_z(i)) * (ztmp-prof_z(i)) + prof_lnrho(i)
                  exit
                endif
              enddo
            endif
          enddo
        elseif (lnrho_init=='prof_lnrho_loop') then
          do j = l1-nghost, l2+nghost
            if (x(j) < prof_z(1) ) then
              profile_x(j) = prof_lnrho(1)
            elseif (x(j) >= prof_z(prof_nz)) then
              profile_x(j) = prof_lnrho(prof_nz)
            else
              do i = 1, prof_nz-1
                if ((x(j) >= prof_z(i)) .and. (x(j) < prof_z(i+1))) then
                  ! linear interpolation: y = m*(x-x1) + y1
                  profile_x(j) = (prof_lnrho(i+1)-prof_lnrho(i)) / &
                      (prof_z(i+1)-prof_z(i)) * (x(j)-prof_z(i)) + prof_lnrho(i)
                  exit
                endif
              enddo
            endif
          enddo
        endif
!
        if (allocated (prof_z)) deallocate (prof_z)
        if (allocated (prof_lnrho)) deallocate (prof_lnrho)
!
        if (ldensity_nolog) then
          f(:,:,:,irho)=spread(spread(exp(profile_x),2,my),3,mz)
        else
         f(:,:,:,ilnrho)=spread(spread(profile_x,2,my),3,mz)
        endif
      endif
!
! read temperature profile for interpolation
      if (lread_lnTT) then
!
! file access is only done on the MPI root rank
        if (lroot) then
          if (.not. file_exists (lnT_dat)) call stop_it_if_any ( &
              .true., 'setup_special: file not found: '//trim(lnT_dat))
! find out, how many data points our profile file has
          prof_nz = (file_size (lnT_dat) - 2*2*4) / (lend*8/lend_b8 * 2)
        endif
        call stop_it_if_any(.false.,'')
        call mpibcast_int(prof_nz,1)
!
        allocate (prof_z(prof_nz), prof_lnTT(prof_nz), stat=ierr)
!
        if (lroot) then
          open (unit,file=lnT_dat,form='unformatted',status='unknown', &
              recl=lend*prof_nz)
          read (unit,iostat=ierr) prof_lnTT
          read (unit,iostat=ierr) prof_z
          if (ierr /= 0) call stop_it_if_any(.true.,'setup_special: '// &
              'Error reading stratification file: "'//trim(lnT_dat)//'"')
          close (unit)
        endif
        call stop_it_if_any(.false.,'')
!
        call mpibcast_real(prof_lnTT,prof_nz)
        call mpibcast_real(prof_z,prof_nz)
!
! convert from logarithmic SI to Pencil units
        prof_lnTT = prof_lnTT - alog(real(unit_temperature))
!
! convert z coordinates from [Mm] to Pencil units
        prof_z = prof_z / unit_length
!
! interpolate temperature profile to Pencil grid
!
        if (lnTT_init=='prof_lnTT_z') then
          do j = l1-nghost, l2+nghost
            ztmp = Lxyz(1) / pi *sin(x(j)/Lxyz(1)*pi )
            if (ztmp < prof_z(1) ) then
              profile_x(j) = prof_lnTT(1)
            elseif (ztmp >= prof_z(prof_nz)) then
              profile_x(j) = prof_lnTT(prof_nz)
            else
              do i = 1, prof_nz-1
                if ((ztmp >= prof_z(i)) .and. (ztmp < prof_z(i+1))) then
                  ! linear interpolation: y = m*(x-x1) + y1
                  profile_x(j) = (prof_lnTT(i+1)-prof_lnTT(i)) / &
                      (prof_z(i+1)-prof_z(i)) * (ztmp-prof_z(i)) + prof_lnTT(i)
                  exit
                endif
              enddo
            endif
          enddo
        else if (lnTT_init=='prof_lnTT_loop') then
          do j = l1-nghost, l2+nghost
            if (x(j) < prof_z(1) ) then
              profile_x(j) = prof_lnTT(1)
            elseif (x(j) >= prof_z(prof_nz)) then
              profile_x(j) = prof_lnTT(prof_nz)
            else
              do i = 1, prof_nz-1
                if ((x(j) >= prof_z(i)) .and. (x(j) < prof_z(i+1))) then
                  ! linear interpolation: y = m*(x-x1) + y1
                  profile_x(j) = (prof_lnTT(i+1)-prof_lnTT(i)) / &
                      (prof_z(i+1)-prof_z(i)) * (x(j)-prof_z(i)) + prof_lnTT(i)
                  exit
                endif
              enddo
            endif
          enddo
        endif
!
        if (allocated (prof_z)) deallocate (prof_z)
        if (allocated (prof_lnTT)) deallocate (prof_lnTT)
!
        if (ltemperature) then
          if (ltemperature_nolog) then
            f(:,:,:,iTT)=spread(spread(exp(profile_x),2,my),3,mz)
          else
            f(:,:,:,ilnTT)=spread(spread(profile_x,2,my),3,mz)
          endif
        else if (lthermal_energy) then
          if (leos) call get_cp1(cp1)
          f(:,:,:,ieth)=spread(spread(exp(profile_x),2,my),3,mz)
          if (ldensity_nolog) then
            f(:,:,:,ieth)=f(:,:,:,ieth)*f(:,:,:,irho)/(gamma*cp1)
          else
            f(:,:,:,ieth)=f(:,:,:,ieth)*exp(f(:,:,:,ilnrho))/(gamma*cp1)
          endif
        else if (lentropy .and. (.not. pretend_lnTT)) then
          if (leos) call get_cp1(cp1)
          f(:,:,:,iss) = (log(gamma_m1/cs20/cp1)+spread(spread(profile_x,2,my),3,mz)- &
              gamma_m1*(f(:,:,:,ilnrho)-log(rho_init))) / cp1 /gamma
        else
          call fatal_error('setup_profiles', &
              'Not implemented for current set of thermodynamic variables.')
        endif
      endif
!
    endsubroutine setup_profiles
!***********************************************************************
    subroutine hydrostatic_lnTT(f)
!
! Solves the hydrostatic equilibrium using the
! temperature from the f-array.
! Integration is done using the trapezoidal rule.
!
      use EquationOfState, only: gamma,gamma_m1,get_cp1,lnrho0
      use Gravity, only: get_xgravity
      use Mpicomm, only: mpisend_real,mpirecv_real
!
      real, dimension (mx,my,mz,mfarray) :: f
      real :: konst,cp1=1.,lnrho_0,int
      real, dimension(mx) :: TT,lnTT,xgrav
      integer :: i,ii
!
      if (nygrid/=1 .or. nzgrid/=1) call fatal_error('hydrostatic_lnTT', &
          'only for nygrid=nzgrid=1')
!
      if (leos) call get_cp1(cp1)
      call get_xgravity(xgrav)
!
      konst = gamma*cp1/gamma_m1
!
      lnrho_0=alog(rho_init)
      f(l1,:,:,ilnrho)= lnrho_0
!
      if (ltemperature) then
        if (ltemperature_nolog) then
          TT = f(:,m1,n1,iTT)
          lnTT = alog(f(:,m1,n1,iTT))
        else
          TT = exp(f(:,m1,n1,ilnTT))
          lnTT = f(:,m1,n1,ilnTT)
        endif
      else
        call fatal_error('hydrostatic_lnTT', &
            'only for ltemperature=T up to now')
      endif
!
      do ii=0,nprocx-1
        if (ipx==ii) then
          do i=l1+1,l2+nghost
            int = 0.5 * (x(i)-x(i-1)) * (xgrav(i)*(1/TT(i-1)+1/TT(i)))
            f(i,:,:,ilnrho)=f(i-1,:,:,ilnrho)-lnTT(i)+ &
                lnTT(i-1)+konst*int
          enddo
          if (ipx < nprocx-1) call mpisend_real(f(l2+1,m1,n1,ilnrho),1,iproc+1,iproc)
        elseif (ipx==ii+1 .and. ipx<=nprocx-1) then
          call mpirecv_real(lnrho_0,1,iproc-1,iproc-1)
          f(l1,:,:,ilnrho) = lnrho_0
        endif
      enddo
!
!  Fill the lower most ghost celĺs. Can be overriden by standard
!  boundary conditions.
!
      if (ipx==0) then
        do i=l1-1,1,-1
          int = 0.5 * (x(i)-x(i+1)) * (xgrav(i)*(1/TT(i+1)+1/TT(i)))
          f(i,:,:,ilnrho)=f(i+1,:,:,ilnrho)-lnTT(i)+ &
              lnTT(i+1)+konst*int
        enddo
      endif
!
    endsubroutine hydrostatic_lnTT
!***********************************************************************
    subroutine write_stratification_dat(f)
!
!  Writes the initial density temperature stratification into each
!  proc subfolder.
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      integer :: unit=12,lend
      real, dimension (mx) :: xwrite_density=0.,xwrite_energy=0.
      real :: dummy=1.
!
      character (len=*), parameter :: filename='/strat.dat'
!
      inquire(IOLENGTH=lend) dummy
!
      if (ldensity) xwrite_density=f(:,m1,n1,ilnrho)
      if (lentropy) xwrite_energy=f(:,m1,n1,iss)
      if (ltemperature) then
        if (ltemperature_nolog) then
          xwrite_energy=f(:,m1,n1,iTT)
        else
          xwrite_energy=f(:,m1,n1,ilnTT)
        endif
      endif
        !
      open(unit,file=trim(directory_snap)//filename, &
          form='unformatted',status='unknown',recl=lend*mx)
      write(unit) x
      write(unit) xwrite_density
      write(unit) xwrite_energy
      close(unit)
!
    endsubroutine write_stratification_dat
!***********************************************************************
    function lin_inpol(f2,f1,x2,x1,yn)
!
      real :: f2,f1,x2,x1,yn
      real :: lin_inpol
!
      lin_inpol = (f2-f1)/(x2-x1)*(yn-x1) + f1
!
    endfunction lin_inpol
!***********************************************************************\
    subroutine setup_tanh(f)
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (mx) :: TT,z_SI,TT_var
      integer :: i,j
!
      z_SI = Lxyz(1) / pi *sin( x/Lxyz(1)*pi ) *unit_length
!
      TT = (T1-T0)*(0.5*tanh((z_SI-z0_tanh)/width_tanh)+0.5)+T0
!
      if (ltemperature) then
        if (ltemperature_nolog) then
          TT_var = TT / unit_temperature
        else
          TT_var = alog(TT /real(unit_temperature))
        endif
      else
        TT_var = impossible
        call fatal_error('setup_tanh','only works for ltemperature=T')
      endif
!
      do i=1,my
        do j=1,mz
          f(:,i,j,ilnTT) = TT_var
        enddo
      enddo
!
    endsubroutine setup_tanh
!***********************************************************************
    subroutine piecewice_poly(f)
!
      use EquationOfState, only: gamma, gamma_m1, get_cp1
      use Gravity, only: gravz
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
      real :: Ttop,T2, T1, T0, cp1=1, temp
      real :: lnrhotop, lnrho2, lnrho1, lnrho0, ztop
      real :: lnrhobot,zbot,Tbot
      real, dimension(4) :: beta
      integer :: i
!
      if (leos) call get_cp1(cp1)
!
!  Top boundary values.
!
      ztop=xyz0(3)+Lxyz(3)
      zbot=xyz0(3)
!
!  Temperature gradients.
!
      beta = cp1*gravz/(mpoly_special+1.)*gamma/gamma_m1
!
!
      T0 = 6000./unit_temperature
      lnrho0 = alog(3e-4/real(unit_density))
!
      Tbot = T0 - beta(1)*(zpoly(1)-zbot)
      T1   = T0 + beta(2)*(zpoly(2)-zpoly(1))
      T2   = T1 + beta(3)*(zpoly(3)-zpoly(2))
      Ttop = T2 + beta(4)*(ztop-zpoly(3))
!
!
      lnrhobot =  lnrho0+mpoly_special(1)*log(Tbot/T0)
      lnrho1   =  lnrho0+mpoly_special(2)*log(T1/T0)
      lnrho2   =  lnrho1+mpoly_special(3)*log(T2/T1)
      lnrhotop =  lnrho2+mpoly_special(4)*log(Ttop/T2)
!
      if (iproc==0) then
      print*,'########################################'
      print*,'Beta',beta
      print*,"TTbot",Tbot*unit_temperature
      print*,"TT1",T0*unit_temperature
      print*,"TT2",T1*unit_temperature
      print*,"TT3",T2*unit_temperature
      print*,"TTtop",Ttop*unit_temperature
      print*,"rhobot",exp(lnrhobot)*unit_density
      print*,"rho0",exp(lnrho0)*unit_density
      print*,"rho1",exp(lnrho1)*unit_density
      print*,"rho2",exp(lnrho2)*unit_density
      print*,"rhotop",exp(lnrhotop)*unit_density
      print*,'########################################'
      endif
      do  i=1,mz
!
        if (z(i) < zpoly(1)) then
          Temp = Tbot + beta(1)*(z(i)-zbot)
          f(:,:,i,ilnTT)=log(Temp)
          f(:,:,i,ilnrho)=lnrhobot+mpoly_special(1)*log(Temp/Tbot)
!
        elseif (z(i) >=zpoly(1) .and. z(i) < zpoly(2)) then
          Temp = T0 + beta(2)*(z(i)-zpoly(1))
          f(:,:,i,ilnTT)=log(Temp)
          f(:,:,i,ilnrho)=lnrho0+mpoly_special(2)*log(Temp/T0)
!
        elseif (z(i) >= zpoly(2) .and. z(i) <zpoly(3)) then
          Temp = T1 + beta(3)*(z(i)-zpoly(2))
          f(:,:,i,ilnTT)=log(Temp)
          f(:,:,i,ilnrho)=lnrho1+mpoly_special(3)*log(Temp/T1)
!
        elseif (z(i) >= zpoly(3)) then
          Temp = T2 + beta(4)*(z(i)-zpoly(3))
          f(:,:,i,ilnTT)=log(Temp)
          f(:,:,i,ilnrho)=lnrho2+mpoly_special(4)*log(Temp/T2)
!
        endif
!
      enddo
!
  endsubroutine piecewice_poly
!***********************************************************************
!
!********************************************************************
!************        DO NOT DELETE THE FOLLOWING       **************
!********************************************************************
!**  This is an automatically generated include file that creates  **
!**  copies dummy routines from noinitial_condition.f90 for any    **
!**  InitialCondition routines not implemented in this file        **
!**                                                                **
    include '../initial_condition_dummies.inc'
!********************************************************************
  endmodule InitialCondition
