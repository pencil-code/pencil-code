! $Id: nohydro.f90,v 1.31 2004-07-05 22:19:50 theine Exp $

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

  real :: othresh=0.,othresh_per_orms=0.,orms=0.,othresh_scl=1.
  real :: nu_turb=0.,nu_turb0=0.,tau_nuturb=0.,nu_turb1=0.
  logical :: lcalc_turbulence_pars=.false.
  real :: kep_cutoff_pos_ext= huge(1.0),kep_cutoff_width_ext=0.0
  real :: kep_cutoff_pos_int=-huge(1.0),kep_cutoff_width_int=0.0
  real :: u_out_kep=0.0

  integer :: dummyuu           ! We cannot define empty namelists
  namelist /hydro_init_pars/ dummyuu
  namelist /hydro_run_pars/  dummyuu

  real :: theta=0.
  real :: Hp,cs_ave,alphaSS,ul0,tl0,eps_diss,teta,ueta,tl01,teta1

  ! other variables (needs to be consistent with reset list below)
  integer :: i_u2m=0,i_um2=0,i_oum=0,i_o2m=0
  integer :: i_uxpt=0,i_uypt=0,i_uzpt=0
  integer :: i_dtu=0,i_dtv=0,i_urms=0,i_umax=0,i_uzrms=0,i_uzmax=0
  integer :: i_orms=0,i_omax=0
  integer :: i_ux2m=0, i_uy2m=0, i_uz2m=0
  integer :: i_uxuym=0, i_uxuzm=0, i_uyuzm=0
  integer :: i_ruxm=0,i_ruym=0,i_ruzm=0,i_rumax=0
  integer :: i_uxmz=0,i_uymz=0,i_uzmz=0,i_umx=0,i_umy=0,i_umz=0
  integer :: i_uxmxy=0,i_uymxy=0,i_uzmxy=0
  integer :: i_Marms=0,i_Mamax=0
  integer :: i_divu2m=0,i_epsK=0
  integer :: i_urmphi=0,i_upmphi=0,i_uzmphi=0,i_u2mphi=0,i_oumphi=0

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
           "$Id: nohydro.f90,v 1.31 2004-07-05 22:19:50 theine Exp $")
!
    endsubroutine register_hydro
!***********************************************************************
    subroutine initialize_hydro(f,lstarting)
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  24-nov-02/tony: coded 
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      logical :: lstarting
!      
      if (ip == 0) print*,f,lstarting  !(keep compiler quiet)
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
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz) :: xx,yy,zz
!
      if(ip==0) print*,f,xx,yy,zz  !(keep compiler quiet)
    endsubroutine init_uu
!***********************************************************************
    subroutine duu_dt(f,df,uu,glnrho,divu,rho1,u2,uij,shock,gshock)
!
!  velocity evolution, dummy routine
!  This routine is used in kinematic dynamo calculations;
!  allow for different prescribed velocity fields.
!
!   7-jun-02/axel: adapted from hydro
!
      use Cdata
      use Magnetic
      use Sub
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx,3,3) :: uij
      real, dimension (nx,3) :: uu,oo,glnrho,gshock
      real, dimension (nx) :: u2,uz,divu,o2,ou,rho1,shock
!
      if (kinflow=='ABC') then
        if (headtt) print*,'ABC flow'
        uu(:,1)=ABC_A*sin(kz_aa*z(n))    +ABC_C*cos(ky_aa*y(m))
        uu(:,2)=ABC_B*sin(kx_aa*x(l1:l2))+ABC_A*cos(kz_aa*z(n))
        uu(:,3)=ABC_C*sin(ky_aa*y(m))    +ABC_B*cos(kx_aa*x(l1:l2))
      elseif (kinflow=='roberts') then
        if (headtt) print*,'Glen Roberts flow; kx_aa,ky_aa=',kx_aa,ky_aa
        uu(:,1)=+sin(kx_aa*x(l1:l2))*cos(ky_aa*y(m))
        uu(:,2)=-cos(kx_aa*x(l1:l2))*sin(ky_aa*y(m))
        uu(:,3)=+sin(kx_aa*x(l1:l2))*sin(ky_aa*y(m))*sqrt(2.)
      else
        if (headtt) print*,'uu=0'
        uu=0.
      endif
!
!  uu/dx for timestep
!
      if (lfirst.and.ldt) advec_uu=abs(uu(:,1))*dx_1(l1:l2)+ &
                                   abs(uu(:,2))*dy_1(  m  )+ &
                                   abs(uu(:,3))*dz_1(  n  )
      if (headtt.or.ldebug) print*,'duu_dt: max(advec_uu) =',maxval(advec_uu)
!
!  Calculate maxima and rms values for diagnostic purposes
!  (The corresponding things for magnetic fields etc happen inside magnetic etc)
!  The length of the timestep is not known here (--> moved to prints.f90)
!
      if (ldiagnos) then
        call dot2_mn(uu,u2)
        if (i_urms/=0)  call sum_mn_name(u2,i_urms,lsqrt=.true.)
        if (i_umax/=0)  call max_mn_name(u2,i_umax,lsqrt=.true.)
        if (i_uzrms/=0) call sum_mn_name(uu(:,3)**2,i_uzrms,lsqrt=.true.)
        if (i_uzmax/=0) call max_mn_name(uu(:,3)**2,i_uzmax,lsqrt=.true.)
        if (i_u2m/=0)   call sum_mn_name(u2,i_u2m)
        if (i_um2/=0)   call max_mn_name(u2,i_um2)
      endif
!
      if(ip==0) print*,f,df,glnrho,divu,rho1,u2,uij,shock,gshock
                                               !(keep compiler quiet)
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
        i_u2m=0; i_um2=0; i_oum=0; i_o2m=0
        i_uxpt=0; i_uypt=0; i_uzpt=0
        i_dtu=0; i_dtv=0; i_urms=0; i_umax=0; i_uzrms=0; i_uzmax=0;
        i_orms=0; i_omax=0
        i_ruxm=0; i_ruym=0; i_ruzm=0; i_rumax=0
        i_ux2m=0; i_uy2m=0; i_uz2m=0
        i_uxuym=0; i_uxuzm=0; i_uyuzm=0
        i_umx=0; i_umy=0; i_umz=0
        i_Marms=0; i_Mamax=0
        i_divu2m=0; i_epsK=0
        i_urmphi=0; i_upmphi=0; i_uzmphi=0; i_u2mphi=0; i_oumphi=0
      endif
!
!  iname runs through all possible names that may be listed in print.in
!
      if(lroot.and.ip<14) print*,'run through parse list'
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'u2m',i_u2m)
        call parse_name(iname,cname(iname),cform(iname),'um2',i_um2)
        call parse_name(iname,cname(iname),cform(iname),'o2m',i_o2m)
        call parse_name(iname,cname(iname),cform(iname),'oum',i_oum)
        call parse_name(iname,cname(iname),cform(iname),'dtu',i_dtu)
        call parse_name(iname,cname(iname),cform(iname),'dtv',i_dtv)
        call parse_name(iname,cname(iname),cform(iname),'urms',i_urms)
        call parse_name(iname,cname(iname),cform(iname),'umax',i_umax)
        call parse_name(iname,cname(iname),cform(iname),'uzrms',i_uzrms)
        call parse_name(iname,cname(iname),cform(iname),'uzmax',i_uzmax)
        call parse_name(iname,cname(iname),cform(iname),'ux2m',i_ux2m)
        call parse_name(iname,cname(iname),cform(iname),'uy2m',i_uy2m)
        call parse_name(iname,cname(iname),cform(iname),'uz2m',i_uz2m)
        call parse_name(iname,cname(iname),cform(iname),'uxuym',i_uxuym)
        call parse_name(iname,cname(iname),cform(iname),'uxuzm',i_uxuzm)
        call parse_name(iname,cname(iname),cform(iname),'uyuzm',i_uyuzm)
        call parse_name(iname,cname(iname),cform(iname),'orms',i_orms)
        call parse_name(iname,cname(iname),cform(iname),'omax',i_omax)
        call parse_name(iname,cname(iname),cform(iname),'ruxm',i_ruxm)
        call parse_name(iname,cname(iname),cform(iname),'ruym',i_ruym)
        call parse_name(iname,cname(iname),cform(iname),'ruzm',i_ruzm)
        call parse_name(iname,cname(iname),cform(iname),'rumax',i_rumax)
        call parse_name(iname,cname(iname),cform(iname),'umx',i_umx)
        call parse_name(iname,cname(iname),cform(iname),'umy',i_umy)
        call parse_name(iname,cname(iname),cform(iname),'umz',i_umz)
        call parse_name(iname,cname(iname),cform(iname),'Marms',i_Marms)
        call parse_name(iname,cname(iname),cform(iname),'Mamax',i_Mamax)
        call parse_name(iname,cname(iname),cform(iname),'divu2m',i_divu2m)
        call parse_name(iname,cname(iname),cform(iname),'epsK',i_epsK)
        call parse_name(iname,cname(iname),cform(iname),'uxpt',i_uxpt)
        call parse_name(iname,cname(iname),cform(iname),'uypt',i_uypt)
        call parse_name(iname,cname(iname),cform(iname),'uzpt',i_uzpt)
      enddo
!
!  write column where which hydro variable is stored
!
      if (lwr) then
        write(3,*) 'i_u2m=',i_u2m
        write(3,*) 'i_um2=',i_um2
        write(3,*) 'i_o2m=',i_o2m
        write(3,*) 'i_oum=',i_oum
        write(3,*) 'i_dtu=',i_dtu
        write(3,*) 'i_dtv=',i_dtv
        write(3,*) 'i_urms=',i_urms
        write(3,*) 'i_umax=',i_umax
        write(3,*) 'i_uzrms=',i_uzrms
        write(3,*) 'i_uzmax=',i_uzmax
        write(3,*) 'i_ux2m=',i_ux2m
        write(3,*) 'i_uy2m=',i_uy2m
        write(3,*) 'i_uz2m=',i_uz2m
        write(3,*) 'i_uxuym=',i_uxuym
        write(3,*) 'i_uxuzm=',i_uxuzm
        write(3,*) 'i_uyuzm=',i_uyuzm
        write(3,*) 'i_orms=',i_orms
        write(3,*) 'i_omax=',i_omax
        write(3,*) 'i_ruxm=',i_ruxm
        write(3,*) 'i_ruym=',i_ruym
        write(3,*) 'i_ruzm=',i_ruzm
        write(3,*) 'i_rumax=',i_rumax
        write(3,*) 'i_umx=',i_umx
        write(3,*) 'i_umy=',i_umy
        write(3,*) 'i_umz=',i_umz
        write(3,*) 'i_Marms=',i_Marms
        write(3,*) 'i_Mamax=',i_Mamax
        write(3,*) 'i_divu2m=',i_divu2m
        write(3,*) 'i_epsK=',i_epsK
        write(3,*) 'nname=',nname
        write(3,*) 'iuu=',iuu
        write(3,*) 'iux=',iux
        write(3,*) 'iuy=',iuy
        write(3,*) 'iuz=',iuz
        write(3,*) 'i_uxpt=',i_uxpt
        write(3,*) 'i_uypt=',i_uypt
        write(3,*) 'i_uzpt=',i_uzpt
        write(3,*) 'i_uxmz=',i_uxmz
        write(3,*) 'i_uymz=',i_uymz
        write(3,*) 'i_uzmz=',i_uzmz
        write(3,*) 'i_uxmxy=',i_uxmxy
        write(3,*) 'i_uymxy=',i_uymxy
        write(3,*) 'i_uzmxy=',i_uzmxy
        write(3,*) 'i_urmphi=',i_urmphi
        write(3,*) 'i_upmphi=',i_upmphi
        write(3,*) 'i_uzmphi=',i_uzmphi
        write(3,*) 'i_u2mphi=',i_u2mphi
        write(3,*) 'i_oumphi=',i_oumphi
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
    subroutine calc_turbulence_pars(f)
!
!  dummy routine
!
!  18-may-04/anders: adapted from hydro
!
      real, dimension (mx,my,mz,mvar+maux) :: f
!
      if(ip==0) print*,f  !(keep compiler quiet)
!      
    endsubroutine calc_turbulence_pars
!***********************************************************************

endmodule Hydro
