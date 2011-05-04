! $Id: noinitial_condition.f90 15726 2010-12-22 20:34:26Z ccyang@ucolick.org $
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: linitial_condition = .true.
!
!***************************************************************
!
! To use your additional initial condition code, edit the
! Makefile.local in the src directory under the run directory in which
! you wish to use your initial condition. Add a line that says e.g.
!
!    INITIAL_CONDITION =   initial_condition/mhs_equilibrium
!
! Here mhs_equilibrium is replaced by the filename of your new file,
! not including the .f90 extension.
!
! This module is based on Tony's special module.
!
module InitialCondition
!
  use Cdata
  use Cparam
  use Messages
  use Sub, only: keep_compiler_quiet
!
  implicit none
!
  include '../initial_condition.h'
  integer :: mk
  real :: kav
  double precision, allocatable,dimension (:) :: kkx,kky,kkz
  integer :: no_of_modes
  real :: relhel=0.,ampluu,amplaa
  real :: scale_kvectorx=1.,scale_kvectory=1.,scale_kvectorz=1.
  logical :: lheluu=.false.,lhelaa=.false.,lscale_kvector_fac,lscale_kvector_tobox
  logical :: lkvec_allocated=.false.
  logical, dimension (3) :: extent
!
!!  integer :: dummy
!
  namelist /initial_condition_pars/  &
  no_of_modes,relhel,lheluu,lhelaa,ampluu,amplaa,scale_kvectorx,scale_kvectory,&
  scale_kvectorz
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
         "$Id: noinitial_condition.f90 15726 2010-12-22 20:34:26Z ccyang@ucolick.org $")
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
      call read_inputk()
! make lhelaa etc consistent
      if((amplaa.ne.0).and. (.not.lhelaa)) & 
          call inevitably_fatal_error ('initialize_initial_condition:', & 
              'amplaa nonzero but lhelaa zero')
      if((ampluu.ne.0).and. (.not.lheluu)) & 
          call inevitably_fatal_error ('initialize_initial_condition:', & 
              'ampluu nonzero but lheluu zero')
!
    endsubroutine initialize_initial_condition
!***********************************************************************
    subroutine initial_condition_all(f)
!
!  Initializes all the f arrays in one call.  This subroutine is called last.
!
!  21-dec-10/ccyang: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
!  SAMPLE IMPLEMENTATION
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_all
!***********************************************************************
    subroutine initial_condition_uu(f)
!
!  Initialize the velocity field.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      if (lheluu) call generate_random_hel (f,iuu,ampluu)
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
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
!  SAMPLE IMPLEMENTATION
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_lnrho
!***********************************************************************
    subroutine initial_condition_ss(f)
!
!  Initialize entropy.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
!  SAMPLE IMPLEMENTATION
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_ss
!***********************************************************************
    subroutine initial_condition_aa(f)
!
!  Initialize the magnetic vector potential.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      if (lhelaa) call generate_random_hel(f,iaa,amplaa)
!
    endsubroutine initial_condition_aa
!***********************************************************************
    subroutine initial_condition_aatest(f)
!
!  Initialize testfield.
!
!  04-may-11/dhruba: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_aatest
!***********************************************************************
    subroutine initial_condition_uutest(f)
!
!  Initialize testflow.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_uutest
!***********************************************************************
    subroutine initial_condition_lncc(f)
!
!  Initialize passive scalar.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_lncc
!***********************************************************************
    subroutine initial_condition_chiral(f)
!
!  Initialize chiral.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_chiral
!***********************************************************************
    subroutine initial_condition_chemistry(f)
!
!  Initialize chemistry.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_chemistry
!***********************************************************************
    subroutine initial_condition_uud(f)
!
!  Initialize dust fluid velocity.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_uud
!***********************************************************************
    subroutine initial_condition_nd(f)
!
!  Initialize dust fluid density.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_nd
!***********************************************************************
    subroutine initial_condition_uun(f)
!
!  Initialize neutral fluid velocity.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_uun
!***********************************************************************
    subroutine initial_condition_lnrhon(f)
!
!  Initialize neutral fluid density.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_lnrhon
!***********************************************************************
    subroutine initial_condition_ecr(f)
!
!  Initialize cosmic rays.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_ecr
!***********************************************************************
    subroutine initial_condition_fcr(f)
!
!  Initialize cosmic ray flux.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_fcr
!***********************************************************************
    subroutine initial_condition_solid_cells(f)
!
!  Initialize solid cells.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_solid_cells
!***********************************************************************
    subroutine initial_condition_cctest(f)
!
!  Initialize testscalar.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_cctest
!***********************************************************************
    subroutine initial_condition_xxp(f,fp)
!
!  Initialize particles' positions.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (:,:), intent(inout) :: fp
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(fp)
!
    endsubroutine initial_condition_xxp
!***********************************************************************
    subroutine initial_condition_vvp(f,fp)
!
!  Initialize particles' velocities.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (:,:), intent(inout) :: fp
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(fp)
!
    endsubroutine initial_condition_vvp
!***********************************************************************
    subroutine read_initial_condition_pars(unit,iostat)
!
!  07-may-09/wlad: coded
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
!  Uncomment these lines back in when you turn 
!  this file into an initial condition, so it 
!  is able to read the namelist.  
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
!  07-may-09/wlad: coded
!
      integer, intent(in) :: unit
!
!  Uncomment this line back in when you turn 
!  this file into an initial condition, so it 
!  is able to write the namelist.  
!
      write(unit,NML=initial_condition_pars)
!
    endsubroutine write_initial_condition_pars
!***********************************************************************
    subroutine read_inputk
!
!  read input file containg k vectors
!
!  04-may-11/dhruba: coded
!
      logical :: lkini_dot_dat_exists
      integer :: ik
!--------------------------------------------
      inquire(FILE="kini.dat", EXIST=lkini_dot_dat_exists)
      if (lkini_dot_dat_exists) then
          open(9,file='kini.dat',status='old')
          read(9,*) mk,kav
          allocate(kkx(mk),kky(mk),kkz(mk))
          lkvec_allocated=.true.
          read(9,*) (kkx(ik),ik=1,mk)
          read(9,*) (kky(ik),ik=1,mk)
          read(9,*) (kkz(ik),ik=1,mk)
          close(9)
        else
          call inevitably_fatal_error ('read_inputk:', & 
              'you must give an input kini.dat file')
        endif
        extent(1)=nx/=1
        extent(2)=ny/=1
        extent(3)=nz/=1
!--------------------------------------------
    endsubroutine read_inputk
!***********************************************************************
    subroutine generate_random_hel (f,ivar,amp)
!
!  generates initial condition which is a sum of random beltrami
! waves. Copied from forcing_hel in forcing module with minor
! modifications. 
!
!  04-may-11/dhruba: coded
!
      use Sub
      use General, only : random_number_wrapper
      real, dimension (mx,my,mz,mfarray) :: f
      integer, intent(in) :: ivar
      real, intent (in) :: amp
      complex, dimension (mx) :: fx
      complex, dimension (my) :: fy
      complex, dimension (mz) :: fz
      real, dimension (3) :: coef1,coef2
      real, dimension(3) :: e1,e2,ee,kk      
      real, dimension (2) :: fran
      real :: phase,phi,kx,ky,kz,k2,k,force_ampl,pi_over_Lx,norm,ffnorm
      real :: ex,ey,ez,kde,sig,fact,kex,key,kez,kkex,kkey,kkez
      integer :: ikk,ik,j
! Randomly selects nk number of modes from the mk read from the file
! kini.dat. Then adds their contribution with random phases.  
      do ikk=1,no_of_modes
        call random_number_wrapper(fran)
        phase=pi*(2*fran(1)-1.)
        ik=mk*(.9999*fran(2))+1
! scale the kvectors if demanded
        if (lscale_kvector_fac) then
          kx=kkx(ik)*scale_kvectorx
          ky=kky(ik)*scale_kvectory
          kz=kkz(ik)*scale_kvectorz
          pi_over_Lx=0.5
        elseif (lscale_kvector_tobox) then
          kx=kkx(ik)*(2.*pi/Lxyz(1))
          ky=kky(ik)*(2.*pi/Lxyz(2))
          kz=kkz(ik)*(2.*pi/Lxyz(3))
          pi_over_Lx=pi/Lxyz(1)
        else
          kx=kkx(ik)
          ky=kky(ik)
          kz=kkz(ik)
          pi_over_Lx=0.5
        endif
!
!  compute k^2 and output wavenumbers
!
        k2=kx**2+ky**2+kz**2
        k=sqrt(k2)
!
!  Find e-vector:
!  Start with old method (not isotropic) for now.
!  Pick e1 if kk not parallel to ee1. ee2 else.
!
        if ((ky==0).and.(kz==0)) then
          ex=0; ey=1; ez=0
        else
          ex=1; ey=0; ez=0
        endif
 !
 !  Isotropize ee in the plane perp. to kk by
 !  (1) constructing two basis vectors for the plane perpendicular
 !      to kk, and
 !  (2) choosing a random direction in that plane (angle phi)
 !  Need to do this in order for the forcing to be isotropic.
 !
        kk = (/kx, ky, kz/)
        ee = (/ex, ey, ez/)
        call cross(kk,ee,e1)
        call dot2(e1,norm); e1=e1/sqrt(norm) ! e1: unit vector perp. to kk
        call cross(kk,e1,e2)
        call dot2(e2,norm); e2=e2/sqrt(norm) ! e2: unit vector perp. to kk, e1
        call random_number_wrapper(phi); phi = phi*2*pi
        ee = cos(phi)*e1 + sin(phi)*e2
        ex=ee(1); ey=ee(2); ez=ee(3)
!
!  k.e
!
        call dot(kk,ee,kde)
!
!  k x e
!
        kex=ky*ez-kz*ey
        key=kz*ex-kx*ez
        kez=kx*ey-ky*ex
!
!  k x (k x e)
!
        kkex=ky*kez-kz*key
        kkey=kz*kex-kx*kez
        kkez=kx*key-ky*kex
!
! set the amplitude
!
        ffnorm=sqrt(1.+relhel**2)*k*sqrt(k2-kde**2)/sqrt(kav)
        fact=amp/ffnorm
!
        fx=exp(cmplx(0.,kx*x+phase))*fact
        fy=exp(cmplx(0.,ky*y))
        fz=exp(cmplx(0.,kz*z))
!
!  prefactor; treat real and imaginary parts separately (coef1 and coef2),
!  so they can be multiplied by different profiles below.
!
        coef1(1)=k*kex; coef2(1)=relhel*kkex
        coef1(2)=k*key; coef2(2)=relhel*kkey
        coef1(3)=k*kez; coef2(3)=relhel*kkez
!
        do n=n1,n2;do m=m1,m2
          do j=1,3
            if (extent(j)) then
              f(l1:l2,m,n,ivar+j-1)= f(l1:l2,m,n,ivar+j-1) + &
                  fact*real(cmplx(coef1(j),coef2(j))*fx(l1:l2)*fy(m)*fz(n))
            endif
          enddo
        enddo;enddo
!
      enddo ! loop over ikk
!
    endsubroutine generate_random_hel
!***********************************************************************
    subroutine initial_condition_clean_up
      if (lkvec_allocated) deallocate(kkx,kky,kkz)
    endsubroutine initial_condition_clean_up
!***********************************************************************
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
