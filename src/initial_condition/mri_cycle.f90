!  
!
!  Initial condition (magnetic field, velocity) 
!   for the MRI cycle described in 
!
!  "Periodic magnetorotational dynamo action as a prototype 
!   of nonlinear magnetic field generation in shear flows"
!   by Herault, Rincon, Cossu, Lesur, Ogilvie & Longaretti
!   (arXiv:1109.1811)
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
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Mpicomm, only: stop_it
  use Messages
!
  implicit none
!
  include '../initial_condition.h'
!
  real :: kx,ky,kz0,ampl_norm=0.3,power_kz=0.
  real,dimension(:),allocatable,save :: ampl_kz
  integer :: nikz=1
  character (len=labellen) :: kzspec='none'
  real :: ampl_kz0
!
  namelist /initial_condition_pars/ &
    nikz,ampl_norm,power_kz,ampl_kz0,kzspec
!
  contains
!***********************************************************************
    subroutine register_initial_condition()
!
!  Configure pre-initialised (i.e. before parameter read) variables
!  which should be know to be able to evaluate
!
!  07-oct-09/wlad: coded
!
!  Identify CVS/SVN version information.
!
      if (lroot) call svn_id( " ")
!
    endsubroutine register_initial_condition
!***********************************************************************
    subroutine initial_condition_uu(f)
!
!  Initialize the velocity field.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      integer :: i,m,n,ikz
      real :: energy_kz
      real :: kx,ky,kz,kvec
!
      allocate(ampl_kz(0:nikz-1))
!
      select case(kzspec)
      case ('flat')
        if (lroot) print*,'setting flat spectrum for kz amplitudes'
        do i=0,nikz-1
          ampl_kz(i)= 1.
        enddo
      case('power-law')
        if (lroot) then
          print*,'setting power-law spectrum for kz amplitudes'
          print*, 'power_kz=',power_kz
        endif
        do i=1,nikz-1
          ampl_kz(i)= ((2*pi/Lz)*float(i))**(-power_kz)
        enddo
        ampl_kz(0)=ampl_kz0
      case default
        if (lroot) print*, 'No such value for kzspec'
        call stop_it('mri_cycle: initial_condition_uu')
      endselect
!
! normalize
!
        energy_kz=0.5*sum(ampl_kz*ampl_kz)
        ampl_kz= ampl_norm*qshear*(ampl_kz/sqrt(energy_kz))
!
      kx=-2*pi/Lx
      ky=2*pi/Ly
      do m=1,my;do n=1,mz
        do ikz=0,nikz-1
          kz=(2.*pi/Lz)*ikz
          kz0=kz
          if (ikz.eq.0) kz0=1.
          kvec=sqrt(kx*kx+ky*ky+kz*kz)
          f(:,m,n,iux) = f(:,m,n,iux) + ampl_kz(ikz)*(kx/kvec) &
              *cos(kx*x+ky*y(m)+kz*z(n))
          f(:,m,n,iuy) = f(:,m,n,iuy) + ampl_kz(ikz)*(ky/kvec) &
              *cos(kx*x+ky*y(m)+kz*z(n))
          f(:,m,n,iuz) = f(:,m,n,iuz) + ampl_kz(ikz)*(kz/kvec) &
              *cos(kx*x+ky*y(m)+kz*z(n))
        enddo
      enddo;enddo
!
    endsubroutine initial_condition_uu
!***********************************************************************
    subroutine initial_condition_lnrho(f)
!
!
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f

!
    endsubroutine initial_condition_lnrho
!***********************************************************************
    subroutine initial_condition_aa(f)
!
!  Initialize the magnetic vector potential. Constant plasma 
!  beta magnetic field. 
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      integer :: i,m,n,ikz
      real :: kx,ky,kz,kvec
!
!  
!
      kx=-2*pi/Lx
      ky=2*pi/Ly
      do m=1,my;do n=1,mz
!
! full spectrum (in kz) of leading (in x,y) shearing-wave packets
!
        do ikz=0,nikz-1
          kz=(2*pi/Lz)*ikz
          kz0=kz
          if (ikz.eq.0) kz0=1.
          kvec=sqrt(kx*kx+ky*ky+kz*kz)
          f(:,m,n,iax) = f(:,m,n,iax) + ampl_kz(ikz) &
              *(sqrt(ky*ky+kz*kz)/(kz0*ky))*sin(kx*x+ky*y(m)+kz*z(n))
          f(:,m,n,iay) = f(:,m,n,iay) + ampl_kz(ikz) &
              *(sqrt(kx*kx+kz*kz)/(kz0*kx))*sin(kx*x+ky*y(m)+kz*z(n))
          f(:,m,n,iaz) = f(:,m,n,iaz) + ampl_kz(ikz) &
              *(sqrt(ky*ky+kx*kx)/(kx *ky))*sin(kx*x+ky*y(m)+kz*z(n))
        enddo
!
! fundamental mode B_{0x}
!
        kz = 2*pi/Lz
        f(:,m,n,iay) = f(:,m,n,iay) - 0.046*qshear*sin(kz*z(n))/kz
!
! modulation B_{mod y}
!
        kx = 2*pi/Lx
        f(:,m,n,iax) = f(:,m,n,iax) + 0.110*qshear*sin(kz*z(n))/kz *cos(kx*x)
      enddo;enddo
!
    endsubroutine initial_condition_aa
!***********************************************************************
    subroutine read_initial_condition_pars(unit,iostat)
!
!  07-may-09/wlad: coded
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
99    return
!
    endsubroutine read_initial_condition_pars
!***********************************************************************
    subroutine write_initial_condition_pars(unit)
!     
      integer, intent(in) :: unit
!
      write(unit,NML=initial_condition_pars)
!
    endsubroutine write_initial_condition_pars
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
