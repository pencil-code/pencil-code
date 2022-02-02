!  Initial condition (density, magnetic field, velocity) 
!  for a field created externally in pencil code format.
!  The file must be in an appropriate format with the right
!  order of indices and match with start.in.
!  USE WITH CAUTION!
!
!  14-jun-17/simon
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
  use Messages
!
  implicit none
!
  include '../initial_condition.h'
!
  character (len=50) :: inFile='var.dat'
!
  namelist /initial_condition_pars/ &
      inFile
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
      if (lroot) call svn_id( &
           "$Id: iucaa_logo.f90 19193 2012-06-30 12:55:46Z wdobler $")
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
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_uu
!***********************************************************************
    subroutine initial_condition_lnrho(f)
!
!  Initialize logarithmic density. init_lnrho 
!  will take care of converting it to linear 
!  density if you use ldensity_nolog
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_lnrho
!***********************************************************************
    subroutine initial_condition_aa(f)
!
!  Initialize the magnetic vector potential.
!
!  15-jun-17/simon: coded
!
!  Load the magnetic field from a vtk file created by GLEMuR
!  and compute the magnetic vector potential from it.
!
!  Created 2017-06-15 by Simon Candelaresi (Iomsn)
!
      use Poisson
      use Sub

      real, dimension (mx,my,mz,mfarray) :: f
      real :: bb_mean_x, bb_mean_y, bb_mean_z
      integer :: l, j, ju, proc
      real, dimension (mx,my,mz,mfarray) :: cc  ! This is phi for poisson.f90
      real, dimension (nx,ny,nz,3) :: curl_cc ! This will be equivalent to the vector potential A.
      ! The next 2 variables are used for the uncurling.
      real, dimension (nx,ny,nz,3) :: jj, tmpJ 
!
!  Read the externally created varfile.
!
! open the binary file for this proc
      write(inFile, '(A3, I0)') 'VAR', (ipz + ipy*nprocz + ipx*nprocy*nprocz)
      write(*,*) inFile
      open(0, file = inFile, form = 'unformatted', action = 'read', access = 'stream')
      read(0) f(:,:,:,iax:iaz)
      close(0)
!
! compute the mean magnetic fluxes through the three directions
      bb_mean_x = sum(f(l1:l2,m1:m2,n1:n2,iax))/size(f(l1:l2,m1:m2,n1:n2,iax))
      bb_mean_y = sum(f(l1:l2,m1:m2,n1:n2,iay))/size(f(l1:l2,m1:m2,n1:n2,iay))
      bb_mean_z = sum(f(l1:l2,m1:m2,n1:n2,iaz))/size(f(l1:l2,m1:m2,n1:n2,iaz))
!
!  Transform the magnetic field into a vector potential
!
!     Compute curl(B) = J for the Poisson solver
      do m=m1,m2
        do n=n1,n2
          call curl(f,iaa,jj(:,m-nghost,n-nghost,:))
        enddo
      enddo
      tmpJ = -jj
!
!     Use the Poisson solver to solve \nabla^2 A = -J for A
      do j=1,3
        call inverse_laplacian(tmpJ(:,:,:,j))
      enddo
!
!     Overwrite the f-array with the correct vector potential A
      do j=1,3
        ju=iaa-1+j
        f(l1:l2,m1:m2,n1:n2,ju) = tmpJ(:,:,:,j)
      enddo

! !  Use the Poisson solver to solve \nabla^2 cc = -B for cc
!       cc(:,:,:,iax:iaz) = -f(:,:,:,iax:iaz)
!       do j=1,3
!         call inverse_laplacian(cc(:,:,:,iax+j-1))
!       enddo
! ! 
! !  Compute curl(cc) = A
!       do m=m1,m2
!          do n=n1,n2
!             call curl(cc,iaa,curl_cc(:,m-nghost,n-nghost,:))
!          enddo
!       enddo
!       
!       f(l1:l2,m1:m2,n1:n2,iax:iaz) = curl_cc(:,:,:,iax:iaz)
!
! !  Overwrite the f-array with the correct vector potential A
!       do j=1,3
!           ju=iaa-1+j
!           f(l1:l2,m1:m2,n1:n2,ju) = tmp(:,:,:,j)
!       enddo
!
! !     Add a background field
!       do l=1,mx
!         do m=1,my
!           f(l,m,:,iax) = f(l,m,:,iax) - y(m)*bb_mean_z/2.
!           f(l,m,:,iay) = f(l,m,:,iay) + x(l)*bb_mean_z/2.
!         enddo
!       enddo
!       do l=1,mx
!         do m=1,my
!           f(l,m,:,iay) = f(l,m,:,iay) - z(:)*bb_mean_x/2.
!           f(l,m,:,iaz) = f(l,m,:,iaz) + y(l)*bb_mean_x/2.
!         enddo
!       enddo
!       do l=1,mx
!         do m=1,my
!           f(l,m,:,iaz) = f(l,m,:,iaz) - x(m)*bb_mean_y/2.
!           f(l,m,:,iax) = f(l,m,:,iax) + z(:)*bb_mean_y/2.
!         enddo
!       enddo
!
    endsubroutine initial_condition_aa
!***********************************************************************
    subroutine read_initial_condition_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=initial_condition_pars, IOSTAT=iostat)
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
!***********************************************************************
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
