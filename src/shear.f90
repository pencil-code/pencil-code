! $Id: shear.f90,v 1.34 2006-06-25 08:15:51 ajohan Exp $

!  This modules deals with all aspects of shear; if no
!  shear is invoked, a corresponding replacement dummy
!  routine is used instead which absorbs all the calls to the
!  shear relevant subroutines listed in here.
!  Shear can either be given relative to Omega (using qshear),
!  or in absolute fashion via the parameters Sshear.

module Shear

  use Sub
  use Cdata
  use Messages

  implicit none

  real, dimension (nz) :: uy0_extra, duy0dz_extra
  real :: eps_vshear=0.0
  logical :: luy0_extra=.false.

  include 'shear.h'

  namelist /shear_init_pars/ &
       qshear,Sshear,deltay,eps_vshear,Omega

  namelist /shear_run_pars/ &
       qshear,Sshear,deltay,eps_vshear,Omega

  integer :: idiag_dtshear=0

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
           "$Id: shear.f90,v 1.34 2006-06-25 08:15:51 ajohan Exp $")
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
      if (qshear/=0.0) Sshear=-qshear*Omega
      if (lroot .and. ip<=12) &
          print*,'initialize_shear: Sshear,qshear=',Sshear,qshear
!
!  Possible to add extra rotation profile.
!
      if (eps_vshear/=0.0) then
        if (lroot) print*, 'initialize_shear: eps_vshear=', eps_vshear
        uy0_extra=Sshear*eps_vshear*Lx*cos(2*pi/Lz*z(n1:n2))
        duy0dz_extra=-Sshear*eps_vshear*Lx*sin(2*pi/Lz*z(n1:n2))*2*pi/Lz
        luy0_extra=.true.
      endif
!
    endsubroutine initialize_shear
!***********************************************************************
    subroutine read_shear_init_pars(unit,iostat)
!
!
!    
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!      
      if (present(iostat)) then
        read(unit,NML=shear_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=shear_init_pars,ERR=99)
      endif
!      
99    return
!
    endsubroutine read_shear_init_pars
!***********************************************************************
    subroutine write_shear_init_pars(unit)
!    
!
!
      integer, intent(in) :: unit
!
      write(unit,NML=shear_init_pars)
!      
    endsubroutine write_shear_init_pars
!***********************************************************************
    subroutine read_shear_run_pars(unit,iostat)
!
!
!    
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=shear_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=shear_run_pars,ERR=99)
      endif
!      
99    return
!
    endsubroutine read_shear_run_pars
!***********************************************************************
    subroutine write_shear_run_pars(unit)
!    
      integer, intent(in) :: unit
!
      write(unit,NML=shear_run_pars)
!
    endsubroutine write_shear_run_pars
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
      use Cdata
      use Deriv
!--   use Testfield, only:ntestfield
      integer, parameter :: ntestfield=36
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
      if (headtt.or.ldebug) print*, 'shearing: Sshear,qshear=', Sshear, qshear
!
!  add shear term, -uy0*df/dy, for all variables
!
      uy0=Sshear*x(l1:l2)
!
!  Add extra rotation profile.
!
      if (luy0_extra) uy0=uy0+uy0_extra(n-nghost)
!        
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
        if (luy0_extra) df(l1:l2,m,n,iuy) = df(l1:l2,m,n,iuy) &
            -f(l1:l2,m,n,iuz)*duy0dz_extra(n-nghost)
      endif
!
!  Loop over dust species
!
      if (ldustvelocity) then
        do k=1,ndustspec
!
!  Correct Coriolis force term for all dust species
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
!  End loop over dust species
!
        enddo
      endif
!
!  Magnetic stretching term
!
      if (lmagnetic) then
        df(l1:l2,m,n,iax)=df(l1:l2,m,n,iax)-Sshear*f(l1:l2,m,n,iay)
        if (luy0_extra) df(l1:l2,m,n,iaz)=df(l1:l2,m,n,iaz)-duy0dz_extra(n-nghost)*f(l1:l2,m,n,iay)
      endif
!
!  Testfield stretching term
!
      if (ltestfield) then
        do j=iaatest,iaatest+ntestfield-1,3
          df(l1:l2,m,n,j)=df(l1:l2,m,n,j)-Sshear*f(l1:l2,m,n,j+1)
        enddo
      endif
!
!  take shear into account for calculating time step
!
      if (lfirst.and.ldt) advec_shear=abs(uy0*dy_1(m))
!
!  Calculate shearing related diagnostics
!
      if (ldiagnos) then
        if (idiag_dtshear/=0) &
            call max_mn_name(advec_shear/cdt,idiag_dtshear,l_dt=.true.)
      endif
!
    end subroutine shearing
!***********************************************************************
    subroutine advance_shear(dt_shear)
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
      real :: dt_shear
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
      deltay=deltay-Sshear*Lx*dt_shear
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
        idiag_dtshear=0
      endif
!
!  iname runs through all possible names that may be listed in print.in
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'dtshear',idiag_dtshear)
      enddo
!
!  write column where which magnetic variable is stored
!
      if (lwr) then
        write(3,*) 'i_dtshear=',idiag_dtshear
      endif
!
    endsubroutine rprint_shear
!***********************************************************************
  endmodule Shear
