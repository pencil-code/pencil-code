! $Id: nohydro_file.f90,v 1.13 2003-10-24 13:17:31 dobler Exp $

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
!***************************************************************

module Hydro

  use Cparam
  use Density

  implicit none

  integer :: dummyuu           ! We cannot define empty namelists
  namelist /hydro_init_pars/ dummyuu
  namelist /hydro_run_pars/  dummyuu

  ! other variables (needs to be consistent with reset list below)
  integer :: i_u2m=0,i_um2=0,i_oum=0,i_o2m=0
  integer :: i_urms=0,i_umax=0,i_orms=0,i_omax=0
  integer :: i_Marms=0,i_Mamax=0

  contains

!***********************************************************************
    subroutine register_hydro()
!
!  Initialise variables which should know that we solve the hydro
!  equations: iuu, etc; increase nvar accordingly.
!
!  6-nov-01/wolf: coded
!
      use Cdata
      use Mpicomm, only: lroot,stop_it
      use Sub
!
      logical, save :: first=.true.
!
      if (.not. first) call stop_it('register_hydro called twice')
      first = .false.
!
      lhydro = .false.
!
!  identify version number (generated automatically by CVS)
!
      if (lroot) call cvs_id( &
           "$RCSfile: nohydro_file.f90,v $", &
           "$Revision: 1.13 $", &
           "$Date: 2003-10-24 13:17:31 $")
!
    endsubroutine register_hydro
!***********************************************************************
    subroutine initialize_hydro()
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  24-nov-02/tony: coded 
!
!  do nothing
!
    endsubroutine initialize_hydro
!***********************************************************************
    subroutine init_uu(f,xx,yy,zz)
!
!  initialise uu and lnrho; called from start.f90
!  Should be located in the Hydro module, if there was one.
!
!   7-jun-02/axel: adapted from hydro
!
      use Cdata
      use Sub
!
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (mx,my,mz) :: xx,yy,zz
!
      if(ip==0) print*,f,xx,yy,zz  !(keep compiler quiet)
    endsubroutine init_uu
!***********************************************************************
    subroutine duu_dt(f,df,uu,glnrho,divu,rho1,u2,uij)
!
!  velocity evolution, dummy routine
!  This routine is used in kinematic dynamo calculations;
!  allow for different prescribed velocity fields.
!
!   7-jun-02/axel: adapted from hydro
!
      use Cdata
      use Magnetic
      use General
      use Sub
      use IO
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx,3,3) :: uij
      real, dimension (nx,3) :: uu,glnrho
      real, dimension (nx) :: divu,u2,rho1
!
      real, save, dimension (nx,ny,nz,3) :: uuu
      logical, save :: first=.true.
      integer :: j
!
      if (kinflow=='file') then
        if (first) then
          print*,'duu_dt: read array uu.dat'
          call input_array(trim(directory)//'/uu.dat',uuu,nx,ny,nz,3)
          first=.false.
        else
          do j=1,3
            uu(:,j)=uuu(:,m-nghost,n-nghost,j)
          enddo
        endif
      elseif (kinflow=='ABC') then
        if (headtt) print*,'ABC flow'
        uu(:,1)=ABC_A*sin(kz_aa*z(n))    +ABC_C*cos(ky_aa*y(m))
        uu(:,2)=ABC_B*sin(kx_aa*x(l1:l2))+ABC_A*cos(kz_aa*z(n))
        uu(:,3)=ABC_C*sin(ky_aa*y(m))    +ABC_B*cos(kx_aa*x(l1:l2))
      else
        if (headtt) print*,'uu=0'
        uu=0.
      endif
!
!  maximum squared avection speed (for timestep)
!
      if (lfirst.and.ldt) then
        call dot2_mn(uu,u2)
        maxadvec2=amax1(maxadvec2,u2)
      endif
!
!  Calculate maxima and rms values for diagnostic purposes
!  (The corresponding things for magnetic fields etc happen inside magnetic etc)
!  The length of the timestep is not known here (--> moved to prints.f90)
!
      if (ldiagnos) then
        call dot2_mn(uu,u2)
        if (i_u2m/=0) call sum_mn_name(u2,i_u2m)
        if (i_um2/=0) call max_mn_name(u2,i_um2)
      endif
!
      if(ip==0) print*,f,df,glnrho,divu,rho1,u2  !(keep compiler quiet)
    endsubroutine duu_dt
!***********************************************************************
    subroutine rprint_hydro(lreset,lwrite)
!
!  reads and registers print parameters relevant for hydro part
!
!   8-jun-02/axel: adapted from hydro
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
        i_u2m=0;i_um2=0;i_oum=0;i_o2m=0
      endif
!
!  iname runs through all possible names that may be listed in print.in
!
      if(lroot.and.ip<14) print*,'run through parse list'
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'u2m',i_u2m)
        call parse_name(iname,cname(iname),cform(iname),'um2',i_um2)
!       call parse_name(iname,cname(iname),cform(iname),'o2m',i_o2m)
!       call parse_name(iname,cname(iname),cform(iname),'oum',i_oum)
      enddo
!
!  write column where which magnetic variable is stored
!
      if (lwr) then
        write(3,*) 'i_u2m=',i_u2m
        write(3,*) 'i_um2=',i_um2
        write(3,*) 'i_o2m=',i_o2m
        write(3,*) 'i_oum=',i_oum
        write(3,*) 'i_urms=',i_urms
        write(3,*) 'i_umax=',i_umax
        write(3,*) 'i_orms=',i_orms
        write(3,*) 'i_omax=',i_omax
        write(3,*) 'nname=',nname
        write(3,*) 'iuu=',iuu
        write(3,*) 'iux=',iux
        write(3,*) 'iuy=',iuy
        write(3,*) 'iuz=',iuz
      endif
!
    endsubroutine rprint_hydro
!***********************************************************************
    subroutine calc_mflow
!
!  dummy routine 
!          
!  19-jul-03/axel: adapted from hydro
!   
    endsubroutine calc_mflow
!***********************************************************************

endmodule Hydro
