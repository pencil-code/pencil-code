! $Id$
!
!  This module fetchs the position of the massive n-body particles 
!  and applies extra shock dissipation around them. This allows for
!  less shock dissipation to be used elsewhere in the quiescent 
!  regions of the simulation box. 
!  
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lspecial = .true.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
!***************************************************************

!-------------------------------------------------------------------
!
! HOW TO USE THIS FILE
! --------------------
!
! The rest of this file may be used as a template for your own
! special module.  Lines which are double commented are intended
! as examples of code.  Simply fill out the prototypes for the
! features you want to use.
!
! Save the file with a meaningful name, eg. geo_kws.f90 and place
! it in the $PENCIL_HOME/src/special directory.  This path has
! been created to allow users to optionally check their contributions
! in to the Pencil-Code CVS repository.  This may be useful if you
! are working on/using the additional physics with somebody else or
! may require some assistance from one of the main Pencil-Code team.
!
! To use your additional physics code edit the Makefile.local in
! the src directory under the run directory in which you wish to
! use your additional physics.  Add a line with all the module
! selections to say something like:
!
!    SPECIAL=special/localshock
!
! Where nstar it replaced by the filename of your new module
! upto and not including the .f90
!
!--------------------------------------------------------------------

module Special

  use Cdata
  use EquationOfState
  use Cparam
  use Messages
  use Sub, only: keep_compiler_quiet

  implicit none

  include '../special.h'

! input parameters
  
  real :: rmask=0.,rmask2,rmask1,rmask12,eta_shock_local=0.
  real :: diffrho_shock_local=0.,nu_shock_local=0.
  integer :: dummy
! global arrays
  real, dimension(nx,ny,nz)   :: shock_mask
  real, dimension(nx,ny,nz,3) :: gshock_mask
! "internal" pencils that don't need to be cast into p%
  real, dimension(nx)   :: shock_masked
  real, dimension(nx,3) :: gshock_masked
!
!  start parameters
!
  namelist /special_init_pars/ dummy!maxspar
!
!  run parameters
!
  namelist /special_run_pars/ &
      eta_shock_local,diffrho_shock_local,nu_shock_local,rmask
!
!
! Keep some over used pencils
!
!!
!! Declare any index variables necessary for main or
!!
!!   integer :: iSPECIAL_VARIABLE_INDEX=0
!
!! other variables (needs to be consistent with reset list below)
!

!
! hydro diagnostics
!

!
! magnetic diagnostics
!

!
! 1D average diagnostics
! 

!
  contains

!***********************************************************************
    subroutine register_special()
!
!  Configure pre-initialised (i.e. before parameter read) variables
!  which should be know to be able to evaluate
!
!
!  6-oct-03/tony: coded
!
      if (lroot) call svn_id( &
          "$Id$")
!
    endsubroutine register_special
!***********************************************************************
    subroutine initialize_special(f,lstarting)
!
!  called by run.f90 after reading parameters, but before the time loop
!
!  06-oct-03/tony: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      logical :: lstarting
!
!  Initialize any module variables which are parameter dependant
!
      rmask2=rmask**2
      rmask1=1./rmask
      rmask12=rmask1**2
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(lstarting)
!
    endsubroutine initialize_special
!***********************************************************************
    subroutine init_special(f)
!
!  initialise special condition; called from start.f90
!  06-oct-2003/tony: coded
!
      real, dimension (mx,my,mz,mvar+maux) :: f
!
      intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine init_special
!***********************************************************************
    subroutine pencil_interdep_special(lpencil_in)
!
!  Interdependency among pencils provided by this module are specified here.
!
!  18-07-06/tony: coded                                                         
!
      logical, dimension(npencils) :: lpencil_in
!
      call keep_compiler_quiet(lpencil_in)
!
    endsubroutine pencil_interdep_special
!***********************************************************************
    subroutine pencil_criteria_special()
!
!  All pencils that this special module depends on are specified here.
!
!  18-07-06/tony: coded
!
      if (ldensity.or.lhydro.or.lmagnetic) then 
        lpenc_requested(i_gshock)=.true.
        lpenc_requested(i_shock)=.true.
      endif
!      
      if (ldensity) then 
        if (ldensity_nolog) then 
          lpenc_requested(i_grho)=.true.
          lpenc_requested(i_del2rho)=.true.
          if (lhydro) lpenc_requested(i_glnrho)=.true.
        else
          lpenc_requested(i_glnrho)=.true.
          lpenc_requested(i_del2lnrho)=.true.
          lpenc_requested(i_glnrho2)=.true.
        endif
!
        if (lhydro) then
          lpenc_requested(i_divu)=.true.
          lpenc_requested(i_graddivu)=.true.
          !lpenc_diagnos(i_diffus_total)=.true.
        endif
      endif
!
      if (lmagnetic) then 
        lpenc_requested(i_del2a)=.true.
        lpenc_requested(i_diva)=.true.
      endif
!
    endsubroutine pencil_criteria_special
!***********************************************************************
    subroutine dspecial_dt(f,df,p)
!
!  calculate right hand side of ONE OR MORE extra coupled PDEs
!  along the 'current' Pencil, i.e. f(l1:l2,m,n) where
!  m,n are global variables looped over in equ.f90
!
!  Due to the multi-step Runge Kutta timestepping used one MUST always
!  add to the present contents of the df array.  NEVER reset it to zero.
!
!  several precalculated Pencils of information are passed if for
!  efficiency.
!
!   06-oct-03/tony: coded
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      intent(in) :: f,p
      intent(inout) :: df
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(p)
!
    endsubroutine dspecial_dt
!***********************************************************************
    subroutine get_slices_special(f,slices)
!
!  Write slices for animation of special variables.
!
!  26-jun-06/tony: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(slices%ready)
!
    endsubroutine get_slices_special
!***********************************************************************
    subroutine read_special_init_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat

      if (present(iostat)) then
        read(unit,NML=special_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=special_init_pars,ERR=99)
      endif

99    return
    endsubroutine read_special_init_pars
!***********************************************************************
    subroutine write_special_init_pars(unit)
      integer, intent(in) :: unit

      write(unit,NML=special_init_pars)

    endsubroutine write_special_init_pars
!***********************************************************************
    subroutine read_special_run_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat

      if (present(iostat)) then
        read(unit,NML=special_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=special_run_pars,ERR=99)
      endif

99    return
    endsubroutine read_special_run_pars
!***********************************************************************
    subroutine write_special_run_pars(unit)
      integer, intent(in) :: unit

      write(unit,NML=special_run_pars)

    endsubroutine write_special_run_pars
!***********************************************************************
    subroutine rprint_special(lreset,lwrite)
!
!  reads and registers print parameters relevant to special
!
!   06-oct-03/tony: coded
!
      use Sub
!
!  define diagnostics variable
!
      integer :: iname,inamer
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
      endif
!
      do iname=1,nname
!        call parse_name(iname,cname(iname),cform(iname),'ur2m',idiag_ur2m)
      enddo
!
      do inamer=1,nnamer
!        call parse_name(inamer,cnamer(inamer),cform(inamer),'urupmr',idiag_urupmr)
      enddo
!
!  write column where which special variable is stored
!
      if (lwr) then
        !hydro
!        write(3,*) 'i_urm=',idiag_urm
      endif
!
    endsubroutine rprint_special
!***********************************************************************
    subroutine calc_pencils_special(f,p)
!   06-oct-03/tony: coded
!
      use Mpicomm
      use General, only: spline
!
      real, dimension(mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
      integer :: i,mg,ng
!
      mg=m-m1+1 ; ng=n-n1+1
      shock_masked=p%shock*shock_mask(:,mg,ng)
      do i=1,3
        gshock_masked(:,i)=&
             shock_mask(:,mg,ng)*p%gshock(:,i) + p%shock*gshock_mask(:,mg,ng,i)
      enddo
!
    endsubroutine calc_pencils_special
!***********************************************************************
    subroutine calc_lspecial_pars(f)
!
!  dummy routine
!
!  15-jan-08/axel: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine calc_lspecial_pars
!***********************************************************************
    subroutine special_calc_density(f,df,p)
!
      use Sub
!
      real, dimension (mx,my,mz,mvar+maux), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p 
      real, dimension (nx) :: fdiff,gshockgrho,gshockglnrho
!
!  Only calculate the shock if a value of lshock_local
!  in this pencil is true
!
      if (headtt) print*,'special_calc_density: add shock diffusion'
!
      if (ldensity_nolog) then
        call dot_mn(gshock_masked,p%grho,gshockgrho)
        fdiff=diffrho_shock_local*&
             (shock_masked*p%del2rho + gshockgrho)
      else
        call dot_mn(gshock_masked,p%glnrho,gshockglnrho)
        fdiff=diffrho_shock_local*&
             (shock_masked*(p%del2lnrho+p%glnrho2) + gshockglnrho)
      endif
!
      df(l1:l2,m,n,ilnrho) = df(l1:l2,m,n,ilnrho) + fdiff
!
    endsubroutine special_calc_density
!***********************************************************************
    subroutine special_calc_hydro(f,df,p)
!
!   16-jul-06/wlyra: coded
!
      use Sub
!
      real, dimension (mx,my,mz,mvar+maux), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
      real, dimension(nx,3) :: tmp,tmp2,fvisc
!
      if (headtt) print*,'special_calc_hydro: add shock viscosity'
!
      if (ldensity) then
        call multsv(p%divu,p%glnrho,tmp2)
        tmp=tmp2 + p%graddivu
        call multsv(nu_shock_local*shock_masked,tmp,tmp2)
        call multsv_add(tmp2,nu_shock_local*p%divu,gshock_masked,tmp)
        fvisc=tmp
!
        df(l1:l2,m,n,iux:iuz) = df(l1:l2,m,n,iux:iuz) + fvisc
      endif
!
    endsubroutine special_calc_hydro
!***********************************************************************
    subroutine special_calc_magnetic(f,df,p)
!
!   06-oct-03/tony: coded
!
     real, dimension (mx,my,mz,mvar+maux), intent(in) :: f
     real, dimension (mx,my,mz,mvar), intent(inout) :: df
     type (pencil_case), intent(in) :: p
     real, dimension(nx,3) :: fres
     integer :: i
!
     if (headtt) print*,'special_calc_magnetic: add shock resistivity'
!
     do i=1,3
       fres(:,i) = eta_shock_local*&
            (shock_masked*p%del2a(:,i)+p%diva*gshock_masked(:,i))
     enddo
!
     df(l1:l2,m,n,iax:iaz) = df(l1:l2,m,n,iax:iaz) + fres
!
    endsubroutine special_calc_magnetic
!***********************************************************************
    subroutine special_calc_entropy(f,df,p)
!
      real, dimension (mx,my,mz,mvar+maux), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
!
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(p)
!
    endsubroutine special_calc_entropy
!***********************************************************************
    subroutine special_calc_particles(fp)
!      
      real, dimension(mpar_loc,mpvar) :: fp
!
    endsubroutine special_calc_particles
!***********************************************************************
    subroutine special_calc_particles_nbody(fsp)
!
!  Calculate the shock mask and its gradient. It's better to 
!  do it analytically than to use grad on the shock_mask array.
!  
!  10-dec-08/wlad: coded
!
      real, dimension(nspar,mpvar) :: fsp
      real, dimension(nx,3) :: tmp
      real, dimension(nx) :: xc,rp2
      real :: e1,e2,e3
      integer :: ks,i,mg,ng
!
      if (headtt) print*,'special_calc_particles_nbody: calculate mask'
!
      shock_mask(:,:,:)=0. ; gshock_mask(:,:,:,:)=0.
      do ks=1,nspar
        e1=fsp(ks,1) ; e2=fsp(ks,2) ; e3=fsp(ks,3) 
!
!  if the particle is inside the box... 
!
        if ((e1.ge.xyz0(1)).and.(e1.le.xyz1(1)).and.&
            (e2.ge.xyz0(2)).and.(e2.le.xyz1(2)).and.&
            (e3.ge.xyz0(3)).and.(e3.le.xyz1(3))) then
!
!  ... calculate rp2.
!
          xc=x(l1:l2)
          do m=m1,m2
          do n=n1,n2
            if (lcartesian_coords) then 
              rp2=(xc-e1)**2+(y(m)-e2)**2+(z(n)-e3)**2
            elseif (lcylindrical_coords) then 
              rp2=xc**2+e1**2 - 2*xc*e1*cos(y(m)-e2) + (z(n)-e3)**2
            elseif (lspherical_coords) then 
              rp2=xc**2 + e1**2 - 2*xc*e1*&
                  (cos(y(m))*cos(e2)+sin(y(m))*sin(e2)*cos(z(n)-e3))
            else
              call fatal_error("special_calc_particles_nbody",&
                  "invalid coordinate system")
            endif
!
            mg=m-m1+1 ; ng=n-n1+1
!
!  Shock-mask: gaussian 
!  rmask12=1./rmask**2
!
            shock_mask(:,mg,ng)=shock_mask(:,mg,ng)+exp(-.5*rp2*rmask12)
            call get_gradshock(shock_mask(:,mg,ng),e1,e2,e3,tmp)
            gshock_mask(:,mg,ng,:)=gshock_mask(:,mg,ng,:)+tmp
!
          enddo
          enddo
        endif
      enddo
!
    endsubroutine special_calc_particles_nbody
!***********************************************************************
    subroutine get_gradshock(fshock,e1,e2,e3,gradshock)
!
!  Analitical grad-shock: works only for the gaussian mask
!
!  10-dec-08/wlad: coded
!
      real, dimension(nx,3) :: gradshock
      real, dimension(nx) :: fshock,base
      real :: e1,e2,e3
!
      if (headtt.and.(m==m1).and.(n==n1)) &
           print*,'calculate shock mask gradient'
!
      base=-rmask12*fshock
!
      if (lcartesian_coords) then
        gradshock(:,1)=base*(x(l1:l2)-e1)
        gradshock(:,2)=base*(y(  m  )-e2)
        gradshock(:,3)=base*(z(  n  )-e3)
      elseif (lcylindrical_coords) then 
        gradshock(:,1)=base*(x(l1:l2)-e1*cos(y(m)-e2))
        gradshock(:,2)=base*(e1*sin(y(m)-e2))
        gradshock(:,3)=base*(z(n)-e3)
      elseif (lspherical_coords) then
        gradshock(:,1)=base*(x(l1:l2)-e1*&
             (cos(y(m))*cos(e2) + sin(y(m))*sin(e2)*cos(z(n)-e3)))
        gradshock(:,2)=base*e1*&
             (sin(y(m))*cos(e2) - cos(y(m))*sin(e2)*cos(z(n)-e3))
        gradshock(:,3)=base*e1*sin(e2)*sin(z(n)-e3)
      else
        call fatal_error("get_gradshock","invalid coordinate system")
      endif
!
    endsubroutine get_gradshock
!***********************************************************************
    subroutine special_before_boundary(f)
!
!   Possibility to modify the f array before the boundaries are
!   communicated.
!
!   Some precalculated pencils of data are passed in for efficiency
!   others may be calculated directly from the f array
!
!   Called from equ, but before the evolving loop that calls the 
!   dynamical equations.
!
!   06-jul-06/tony: coded
!
      use Mpicomm
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine special_before_boundary
!**********************************************************************
!**********************************************************************
!************        DO NOT DELETE THE FOLLOWING       **************
!********************************************************************
!**  This is an automatically generated include file that creates  **
!**  copies dummy routines from nospecial.f90 for any Special      **
!**  routines not implemented in this file                         **
!**                                                                **
    include '../special_dummies.inc'
!********************************************************************

endmodule Special

