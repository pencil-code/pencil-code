! Isolated, localized, entropy perturbations ('blobs'; also see 'blob',
! 'blob_hs', and 'blob_radeq' in entropy.f90)
!
! 04-Nov-2026/Kishore: added
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
  use Quiet
  use Cdata
!
  implicit none
!
  include '../initial_condition.h'
!
  real, dimension(ninit) :: ampl_ss=0, radius_ss=1, center1_x=0, center1_y=0, center1_z=0
  character (len=labellen), dimension(ninit) :: blob_type=''
!
  namelist /initial_condition_pars/ ampl_ss, radius_ss, center1_x, center1_y, &
  center1_z, blob_type
!
  real, dimension(:,:), pointer :: reference_state
!
  contains
!***********************************************************************
    subroutine initialize_initial_condition(f)
!
!  Initialize any module variables which are parameter dependent.
!
      use SharedVariables, only: get_shared_variable
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
!  Get reference_state. Requires that density is initialized before energy.
!
      if (lreference_state) call get_shared_variable('reference_state',reference_state)
!
    endsubroutine initialize_initial_condition
!***********************************************************************
    subroutine initial_condition_all(f,profiles)
!
!  Initializes all the f arrays in one call.  This subroutine is called last.
!
      use Messages, only: fatal_error
      use EquationOfState, only: get_gamma_etc

      real, dimension (mx,my,mz,mfarray), optional, intent(inout):: f
      real, dimension (:,:),              optional, intent(out)  :: profiles
!
      real :: cp,cp1
      integer :: j
!
      do j=1,ninit
        select case (trim(blob_type(j)))
!
          case ('cos^2')
            if (lroot) print*,'initial_condition_all: adding '//trim(blob_type(j))//' blob'
            call get_gamma_etc(cp=cp); cp1=1./cp
            call blob_compact(ampl_ss(j),f,iss,radius_ss(j),center1_x(j),center1_y(j),center1_z(j))
            call blob_compact(-cp1*ampl_ss(j),f,ilnrho,radius_ss(j),center1_x(j),center1_y(j),center1_z(j))
            if (ldensity_nolog) then
                f(:,:,:,irho) = exp(f(:,:,:,ilnrho))
              if (lreference_state) then
                do m=1,my; do n=1,mz;
                      f(l1:l2,m,n,irho) = f(l1:l2,m,n,irho) - reference_state(:,iref_rho)
                enddo; enddo
              endif
            endif
!
          case ('')
            !do nothing
!
          case default
            call fatal_error('initial_condition_all', 'no such blob profile: '//trim(blob_type(j)))
!
        end select
      enddo
!
    endsubroutine initial_condition_all
!***********************************************************************
    subroutine read_initial_condition_pars(iomsg)
!
      use File_io, only: parallel_unit
!
      character(LEN=*), intent(out) :: iomsg
!
      integer :: iostat
!
      read(parallel_unit, NML=initial_condition_pars, IOSTAT=iostat, IOMSG=iomsg)
      if (iostat==0) iomsg=""
!
    endsubroutine read_initial_condition_pars
!***********************************************************************
    subroutine write_initial_condition_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=initial_condition_pars)
!
    endsubroutine write_initial_condition_pars
!********************************************************************
    subroutine blob_compact(ampl,f,i,radius,xblob,yblob,zblob)
!
!  Single blob, compactly supported
!
!  03-Mov-2026/Kishore: coded modifying sub.f90/blob
!
      integer, intent(in) :: i
      real, contiguous, dimension(:,:,:,:), intent(inout) :: f
      real, intent(in) :: ampl,radius
      real, optional, intent(in) :: xblob,yblob,zblob
!
      real, dimension (mx) :: delx
      real, dimension (my) :: dely
      real, dimension (mz) :: delz
      real :: x01=0.,y01=0.,z01=0., fact, delr
      integer :: l,m,n
!
      if (present(xblob)) x01=xblob
      if (present(yblob)) y01=yblob
      if (present(zblob)) z01=zblob
!
      if (ampl==0) then
        return
      endif
!
      fact=pi/(2*radius)
!
      delx=x-x01
      dely=y-y01
      delz=z-z01
!
      do n=1,mz; do m=1,my; do l=1,mx
        delr = min(radius, sqrt(delx(l)**2 + dely(m)**2 + delz(n)**2))
        f(l,m,n,i) = f(l,m,n,i) + ampl*cos(delr*fact)**2
      enddo; enddo; enddo
!
    endsubroutine blob_compact
!***********************************************************************
!********************************************************************
!************        DO NOT DELETE THE FOLLOWING       **************
!********************************************************************
!**  This is an automatically generated include file that creates  **
!**  copies dummy routines from noinitial_condition.f90 for any    **
!**  InitialCondition routines not implemented in this file        **
!**                                                                **
    include '../initial_condition_dummies.inc'
!********************************************************************!
endmodule InitialCondition
