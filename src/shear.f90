! $Id: shear.f90,v 1.19 2004-07-08 09:07:49 ajohan Exp $

!  This modules deals with all aspects of shear; if no
!  shear is invoked, a corresponding replacement dummy
!  routine is used instead which absorbs all the calls to the
!  shear relevant subroutines listed in here.
!  Shear can either be given relative to Omega (using qshear),
!  or in absolute fashion via the parameters Sshear.

module Shear

  use Sub
  use Cdata

  implicit none

  namelist /shear_init_pars/ &
       qshear,Sshear

  namelist /shear_run_pars/ &
       qshear,Sshear

  integer :: i_dtshear=0

  contains

!***********************************************************************
    subroutine register_shear()
!
!  Initialise variables
!
!  2-july-02/nils: coded
!
      use Mpicomm
!
      logical, save :: first=.true.
!
      if (.not. first) call stop_it('register_shear called twice')
      first = .false.
!
      lshear = .true.
!
!  identify version number
!
      if (lroot) call cvs_id( &
           "$Id: shear.f90,v 1.19 2004-07-08 09:07:49 ajohan Exp $")
!
    endsubroutine register_shear
!***********************************************************************
    subroutine initialize_shear()
!
!  21-nov-02/tony: coded
!  08-jul-04/anders: Sshear calculated whenever qshear /= 0

!  calculate shear flow velocity; if qshear is given then Sshear=-qshear*Omega
!  is calculated. Otherwise Sshear keeps its value from the input list.
!
      if (qshear /= 0.) Sshear=-qshear*Omega
      if (ip <= 12) print*,'initialize_shear: Sshear,qshear=',Sshear,qshear

    endsubroutine initialize_shear
!***********************************************************************
    subroutine shearing(f,df)
!
!  Calculates the shear terms, -uy0*df/dy (shearing sheat approximation)
!
!  2-jul-02/nils: coded
!  6-jul-02/axel: runs through all nvar variables; added timestep check
! 16-aug-02/axel: use now Sshear which is calculated in param_io.f90
! 20-aug-02/axel: added magnetic stretching term
!
      use Cparam
      use Deriv
      use Hydro, only:theta
!
      integer :: j,k
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx) :: uy0,dfdy
!
      intent(in)  :: f
!
!  print identifier
!
      if (headtt.or.ldebug) print*,'shearing: Sshear,qshear=',Sshear,qshear
!
!  add shear term, -uy0*df/dy, for all variables
!
      uy0=Sshear*x(l1:l2)
      do j=1,nvar
        call der(f,j,dfdy,2)
        df(l1:l2,m,n,j)=df(l1:l2,m,n,j)-uy0*dfdy
      enddo
!
! Taking care of the fact that the Coriolis force changes when 
! we have got shear. The rest of the Coriolis force is calculated 
! in hydro.
!
      if (lhydro) then
        if (theta==0) then
          df(l1:l2,m,n,iuy)=df(l1:l2,m,n,iuy)-Sshear*f(l1:l2,m,n,iux)
        else
          if (headtt) print*,'Sure you want Sshear with finite theta??'
          df(l1:l2,m,n,iuy)=df(l1:l2,m,n,iuy)&
               -Sshear*cos(theta*pi/180.)*f(l1:l2,m,n,iux)
        endif
      endif
!
!  Loop over dust layers
!
      if (ldustvelocity) then
        do k=1,ndustspec
!
!  Correct Coriolis force term for all dust layers 
!
          if (theta==0) then
            df(l1:l2,m,n,iudy(k)) = df(l1:l2,m,n,iudy(k)) &
                - Sshear*f(l1:l2,m,n,iudx(k))
          else
            if (headtt) print*,'Sure you want Sshear with finite theta??'
            df(l1:l2,m,n,iudy(k)) = df(l1:l2,m,n,iudy(k)) &
                - Sshear*cos(theta*pi/180.)*f(l1:l2,m,n,iudx(k))
          endif
!
!  End loop over dust layers
!
        enddo
      endif
!
!  Magnetic stretching term
!
      if (lmagnetic) then
        df(l1:l2,m,n,iax)=df(l1:l2,m,n,iax)-Sshear*f(l1:l2,m,n,iay)
      endif
!
!  take shear into account for calculating time step
!
      if (lfirst.and.ldt) advec_shear=abs(uy0*dy_1(m))
!
!  Calculate shearing related diagnostics
!
      if (ldiagnos) then
        if (i_dtshear/=0) call max_mn_name(advec_shear/cdt,i_dtshear,l_dt=.true.)
      endif
!
    end subroutine shearing
!***********************************************************************
    subroutine advance_shear
!
!  advance shear distance, deltay, using dt. Using t instead introduces
!  significant errors when nt = t/dt exceeds ~100,000 steps.
!  This formulation works also when Sshear is changed during the run.
!
! 18-aug-02/axel: incorporated from nompicomm.f90
!
      use Cdata
      use Mpicomm, only: stop_it
!
!  Works currently only when Sshear is not positive
!
      if (Sshear>0.) then
        if(lroot) print*,'Note: must use non-positive values of Sshear'
        call stop_it("")
      endif
!
!  Make sure deltay is in the range 0 <= deltay < Ly (assuming Sshear<0).
!
      deltay=deltay-Sshear*Lx*dt
      deltay=deltay-int(deltay/Ly)*Ly
!
!  print identifier
!
      if (headtt.or.ldebug) print*,'advance_shear: deltay=',deltay
!
    end subroutine advance_shear
!***********************************************************************
    subroutine rprint_shear(lreset,lwrite)
!
!  reads and registers print parameters relevant to shearing
!
!   2-jul-04/tobi: adapted from entropy
!
      use Cdata
      use Sub
!
      integer :: iname
      logical :: lreset,lwr
      logical, optional :: lwrite
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
!  reset everything in case of reset
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        i_dtshear=0
      endif
!
!  iname runs through all possible names that may be listed in print.in
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'dtshear',i_dtshear)
      enddo
!
!  write column where which magnetic variable is stored
!
      if (lwr) then
        write(3,*) 'i_dtshear=',i_dtshear
      endif
!
    endsubroutine rprint_shear
!***********************************************************************
  end module Shear
