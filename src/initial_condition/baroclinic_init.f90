! $Id$
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: linitial_condition = .true.
!
!***************************************************************
!
!
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
  real :: rho_rms=0.05,xmid=0.0
  integer :: xmodes=10,ymodes=10,zmodes=0
  logical :: lunstratified=.false.,lgaussian_distributed_noise=.true.
!
  namelist /initial_condition_pars/ &
      xmodes,ymodes,zmodes,rho_rms,xmid,lunstratified,&
      lgaussian_distributed_noise
!
  contains
!***********************************************************************
    subroutine register_initial_condition()
!
!  Register variables associated with this module; likely none.
!
!  07-may-09/wlad: coded
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
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_initial_condition
!***********************************************************************
    subroutine initial_condition_uu(f)
!
!  Initialize the velocity field.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
!  SAMPLE IMPLEMENTATION
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_uu
!***********************************************************************
    subroutine initial_condition_lnrho(f)
!
!  Initialize logarithmic density. init_lnrho will take care of
!  converting it to linear density if you use ldensity_nolog.
!
!  07-may-09/wlad: coded
!
      use General, only: random_number_wrapper
      use Mpicomm
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (nx,ny,nz) :: lump_of_sines
      real :: Lx,Ly,Lz,d0,phase,xi,yi,zi,nw1,fac
      real :: fmeantmp_rho2,fmeantmp_rho
      real :: fmean_rho2,fmean_rho
      real :: unnormalized_rho_rms
      real :: normalization_factor
      integer :: i,mm,nn,ll,mm1,ll1,nn1,itmp,irho
!
      Lx=Lxyz(1) ; Ly=Lxyz(2) ; Lz=Lxyz(3)
!
      d0=0.2*Lx
!
      if (lroot) print*,'domain size=',Lx,Ly,Lz
!
!  save the original linear density
!  on the free spot of shock, which isn't set yet
!
      if (lshock) then
        itmp=ishock
        f(l1:l2,m1:m2,n1:n2,itmp)=f(l1:l2,m1:m2,n1:n2,ilnrho)
      else
        if (lroot) then
          print*,'This is a baroclinic run. You expect '
          print*,'many vortices to form in this run, do you not?'
          print*,'These beasts usually get supersonic, which is '
          print*,'actually what prevents them from growing too much,'
          print*,'as then they shock and dissipate. Yet, you are'
          print*,'NOT using shock viscosity. I recommend you to stop'
          print*,'and switch SHOCK=shock_highorder in src/Makefile.local'
        endif
        call fatal_error("","")
      endif
!
!  Work with LINEAR density. As this is start time, there is no confusion
!  linear-log. All initial conditions are to be coded in log. We will
!  construct the density here as linear, then pass to log
!
      irho=ilnrho
!
!  Have to set it to something first, otherwise the sines
!  can lead to negative densities.
!
      f(l1:l2,m1:m2,n1:n2,irho)=1.
      lump_of_sines=0.
!
      if (lroot) then
        print*,'initial_condition_lnrho: Calculating the finite-amplitude '
        print*,'                         pertubations. It will loop the whole '
        print*,'                         grid for every mode. May take a '
        print*,'                         little while.'
      endif
!
      do ll=-xmodes,xmodes ; do mm=0,ymodes ; do nn=-zmodes,zmodes
!
        if (lroot) call random_number_wrapper(phase)
        call mpibcast_real(phase,1)
!
        do i=1,nx ; do m=1,ny ; do n=1,nz
          ll1=i+l1-1 ; xi=x(ll1)
          mm1=m+m1-1 ; yi=y(mm1)
          nn1=n+n1-1 ; zi=z(nn1)
!
          lump_of_sines(i,m,n)=lump_of_sines(i,m,n) + &
             sin(2*pi*(ll*xi/Lx + mm*yi/Ly + nn*zi/Lz+ phase))
!
        enddo;enddo;enddo !end grid loop
      enddo;enddo;enddo !end modes loop
!
!  Now construct the density and fill in the normalization
!  constants needed to get a rms equal to the input rho_rms.
!
      fmeantmp_rho2=0.
      fmeantmp_rho=0.
      nw1=1./(nxgrid*nygrid*nzgrid*1.0)
!
      do n=1,nz
        nn1=n+n1-1
        do m=1,ny
          mm1=m+m1-1
          do i=1,nx
            ll1=i+l1-1 ; xi=x(ll1)
!
            if (lgaussian_distributed_noise) then
              fac=exp(-(.5*(xi-xmid)/d0)**2)
            else
              fac=1
            endif
!
            f(ll1,mm1,nn1,irho)=f(ll1,mm1,nn1,irho) + &
                 lump_of_sines(i,m,n)*fac
          enddo
          fmeantmp_rho2=fmeantmp_rho2+nw1*sum(f(l1:l2,mm1,nn1,irho)**2)
          fmeantmp_rho =fmeantmp_rho +nw1*sum(f(l1:l2,mm1,nn1,irho))
        enddo
      enddo !end grid loop
!
!  Sum the normalization constants over processors, and perform
!  the normalization.
!
      call mpireduce_sum(fmeantmp_rho2,fmean_rho2)
      call mpireduce_sum(fmeantmp_rho,fmean_rho)
      call mpibcast_real(fmean_rho2,1)
      call mpibcast_real(fmean_rho,1)
!
      unnormalized_rho_rms=sqrt(fmean_rho2-fmean_rho**2)
      normalization_factor=rho_rms/unnormalized_rho_rms
!
!  Assumes rho0=1.
!
      f(l1:l2,m1:m2,n1:n2,irho)=1.+&
          normalization_factor*(f(l1:l2,m1:m2,n1:n2,irho)-1.)
!
      if (lroot) then
        print*,'max density (linear): ',maxval(f(l1:l2,m1:m2,n1:n2,irho))
        print*,'min density (linear): ',minval(f(l1:l2,m1:m2,n1:n2,irho))
        print*,'rms density (linear): ',&
            unnormalized_rho_rms*normalization_factor
      endif
!
!  convert to log before finishing.
!
      f(l1:l2,m1:m2,n1:n2,ilnrho)=&
          f(l1:l2,m1:m2,n1:n2,itmp)+log(f(l1:l2,m1:m2,n1:n2,irho))
!
!  keep the original stratification in the ishock slot since it
!  will be needed when setting the entropy
!
    endsubroutine initial_condition_lnrho
!***********************************************************************
    subroutine initial_condition_ss(f)
!
!  Initialize entropy.
!
!  07-may-09/wlad: coded
!
      use EquationOfState, only: gamma,gamma_m1,gamma1,cs20,rho0,lnrho0
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (nx) :: lnrho,lnTT,TT,rho
      real :: cp,cv,cp1,lnTT0,pp0,TT0
      integer :: irho
!
!  Get the density and use a constant pressure entropy condition
!
      cp=1.
      cp1=1/cp
      cv=gamma1*cp
!
      TT0 = cs20*cp1/gamma_m1 ; lnTT0=log(TT0)
      pp0=(cp-cv)*TT0*rho0
!
      irho=ilnrho
!
      do m=m1,m2;do n=n1,n2
!
        if (nzgrid==1.or.&
            (nzgrid/=1.and.lunstratified)) then
          !unstratified case, everything is fresh
          rho = f(l1:l2,m,n,ilnrho)
          TT=(pp0/((cp-cv)*rho))/TT0
        else
          !stratified case, there's a previous
          !density stratification, stored in ishock
          rho = f(l1:l2,m,n,ilnrho)/exp(f(l1:l2,m,n,ishock))
          TT=1/rho !now this is just the temperature perturbation
        endif
!
        lnTT = log(TT)
        lnrho=log(rho)
!
!  assumes cp=1
!
        f(l1:l2,m,n,iss)=f(l1:l2,m,n,iss) + &
            cv*(lnTT-gamma_m1*lnrho)
!
      enddo;enddo
!
!  Revert the original variable (ishock) to zero.
!
      f(:,:,:,ishock)=0.
!
    endsubroutine initial_condition_ss
!********************************************************************
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
