! $Id$
!
!  This modules deals with all aspects of testfield fields; if no
!  testfield fields are invoked, a corresponding replacement dummy
!  routine is used instead which absorbs all the calls to the
!  testfield relevant subroutines listed in here.
!
!  NOTE: since the fall of 2007 we have been using the routine
!  testfield_z.f90, not this one! For more information, please
!  contact Axel Brandenburg <brandenb@nordita.org>
!
!  Alex Hubbard and Matthias Rheinhardt have then developed the
!  present module from inactive/testfield.f90 rather than the
!  inactive/testfield_xz.f90 that also exists.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: ltestfield = .true.
! CPARAM logical, parameter :: ltestfield_xz = .true.
! CPARAM logical, parameter :: ltestfield_z = .false.
! CPARAM logical, parameter :: ltestfield_xy = .false.
!
!***************************************************************
!
module Testfield
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
  use Testfield_general
!
  implicit none
!
  include '../testfield.h'
!
! run parameters
!
  logical :: lflucts_with_xyaver=.false.
  real    :: ktestfield_x=1., ktestfield_z=1., xx0=0., zz0=0.
!
  namelist /testfield_run_pars/       &
           B_ext,                     &
           reinitialize_aatest,       &
           lsoca,                     &
           etatest,                   &
           itestfield,                &
           ktestfield_x,              &
           ktestfield_z,              &
           xx0,                       &
           zz0,                       &
           daainit,                   &
           linit_aatest,              &
           rescale_aatest,            &
           lflucts_with_xyaver
!
! diagnostic variables
!
  integer, parameter :: n_cdiags = 62, len_cdiags = 8
  character(LEN=len_cdiags), dimension(n_cdiags) :: cdiags = &
  (/ 'alp11   ','alp21   ','alp31   ','alp12   ','alp22   ','alp32   ','alp13   ','alp23   ','alp33   ',&   ! DIAG_DOC: $\alpha_{ij}$       
     'eta111  ','eta211  ','eta311  ','eta121  ','eta221  ','eta321  ','eta131  ','eta231  ','eta331  ',&   ! DIAG_DOC: $\eta_{ijk}$
     'eta113  ','eta213  ','eta313  ','eta123  ','eta223  ','eta323  ','eta133  ','eta233  ','eta333  ',&
     'alp11cc ','alp11cs ','alp11sc ','alp11ss ','eta123cc','eta123cs','eta123sc','eta123ss',           &   ! DIAG_DOC: $\alpha_{11,\rm hh},$ 
                                                                                                            ! Diag_DOC: $\eta_{123,\rm hh}, {\rm h}={\rm c,s}$
     'E11     ','E21     ','E31     ','E12     ','E22     ','E32     ','E13     ','E23     ','E33     ',&   ! DIAG_DOC: ${\cal E}^j_i$
     'E14     ','E24     ','E34     ','E15     ','E25     ','E35     ','E16     ','E26     ','E36     ',&   
     'E17     ','E27     ','E37     ','E18     ','E28     ','E38     ','E19     ','E29     ','E39     '  /) 
!
  integer, dimension(n_cdiags):: idiags=0, idiags_z=0, idiags_xz=0
  integer, parameter :: idiag_base_end=27, idiag_Eij_start=36, idiag_Eij_end=idiag_Eij_start+27-1
!
  integer, dimension(4) :: idiag_alp11h, idiag_eta123h            
  equivalence(idiags(idiag_base_end+1),idiag_alp11h), (idiags(idiag_base_end+5),idiag_eta123h)      ! alias names for selected diagnostics
!
!  work variables
!
  real, dimension(nx)              :: cx,sx                     !GPU => DEVICE
  real, dimension(nz)              :: cz,sz
  real, dimension(nx,nz,3,3)       :: Minv                      !GPU => DEVICE
  real, dimension(nx,nz,3,njtest)  :: uxbtestm
  logical, dimension(idiag_base_end):: twod_need_1d, twod_need_2d
  logical, dimension(2)            :: needed2d
!
  contains
!
!***********************************************************************
    subroutine initialize_testfield(f)
!
!  Perform any post-parameter-read initialization
!
!   2-jun-05/axel: adapted from magnetic
!  27-jun-13/MR  : set itestfield='1' as it is only implemented case
!                  set lcalc_uumeanxz=.true., itestfield now string     
!   9-sep-13/MR  : intro'd use of initialize_testfield_general
!  20-oct-13/MR  : call to calc_means_hydro intro'd (needed when used with hydro_kinematic)
!
      use Hydro, only: lcalc_uumeanxz, calc_means_hydro
!
      real, dimension(mx,my,mz,mfarray), intent(inout):: f
!
      call initialize_testfield_general(f)
!
      !!!if (reinitialize_aatest) f(:,:,:,iaatest:iaatest+ntestfield-1)=0.  !!!TBC
!
! calculate inverse matrix for determination of the turbulent coefficients
!
      call calc_inverse_matrix(x(l1:l2),z(n1:n2),ktestfield_x,ktestfield_z,xx0,zz0,Minv,cx,sx,cz,sz)
!
      if (lhydro_kinematic) call calc_means_hydro(f)
!
!  write testfield information to a file (for convenient post-processing)
!
      if (lroot) then
        open(1,file=trim(datadir)//'/testfield_info.dat',STATUS='unknown')
        write(1,'(a,i1)') 'lsoca=',merge(1,0,lsoca)
        write(1,'(a,a)') 'itestfield=',itestfield
        write(1,'(a,2(f3.0))') 'ktestfield_x,z=', ktestfield_x, ktestfield_z
        close(1)
      endif
!
    endsubroutine initialize_testfield
!***********************************************************************
    subroutine read_testfield_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=testfield_run_pars, IOSTAT=iostat)
!
    endsubroutine read_testfield_run_pars
!***********************************************************************
    subroutine write_testfield_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=testfield_run_pars)
!
    endsubroutine write_testfield_run_pars
!***********************************************************************
    subroutine daatest_dt(f,df,p)
!
!  testfield evolution
!  calculate da^(q)/dt=uxB^(q)+eta*del2A^(q), where q=1,...,9
!
!   3-jun-05/axel: coded
!  27-jun-13/MR  : correct calculation of uufluct intro'd
!                  moved calculation of xy-averaged quantities to 
!                  calc_coefficients, completed
!  20-oct-13/MR  : cases for itestfield='linear','1-alt' added
!
      use Hydro, only: uumxz
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      intent(in)     :: f,p
      intent(inout)  :: df
!
      integer :: nl

      nl=n-n1+1
      select case (itestfield)
        case ('1')         ; call rhs_daatest(f,df,p,uumxz(l1:l2,n,:),uxbtestm(:,nl,:,:),set_bbtest)
        case ('1-alt')     ; call rhs_daatest(f,df,p,uumxz(l1:l2,n,:),uxbtestm(:,nl,:,:),set_bbtest_alt)
        case ('alt-II')    ; call rhs_daatest(f,df,p,uumxz(l1:l2,n,:),uxbtestm(:,nl,:,:),set_bbtest_alt_II)
        case ('2')         ; call rhs_daatest(f,df,p,uumxz(l1:l2,n,:),uxbtestm(:,nl,:,:),set_bbtest2)
        case ('3')         ; call rhs_daatest(f,df,p,uumxz(l1:l2,n,:),uxbtestm(:,nl,:,:),set_bbtest3)
        case ('4','linear'); call rhs_daatest(f,df,p,uumxz(l1:l2,n,:),uxbtestm(:,nl,:,:),set_bbtest4)
        case default       ; call fatal_error('daatest_dt','undefined itestfield')
      endselect
!
    endsubroutine daatest_dt
!***********************************************************************
    subroutine get_slices_testfield(f,slices)
!
!  Dummy routine
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(slices%ready)
!
    endsubroutine get_slices_testfield
!***********************************************************************
    subroutine testfield_before_boundary(f)
!
!  Actions to take before boundary conditions are set.
!
!    4-oct-18/axel+nishant: adapted from testflow
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine testfield_before_boundary
!***********************************************************************
    subroutine testfield_after_boundary(f)
!
!  calculate <uxb>, which is needed when lsoca=.false.
!
!  21-jan-06/axel: coded
!  27-sep-13/MR  : corrected fac, simplified communication,
!                  added call to calc_pencils_hydro to enable work
!                  with hydro_kinematic (restricted to p%uu)
!  20-oct-13/MR  : call of calc_coefficients corrected
!  30-oct-13/MR  : added parameter nygrid to call of calc_coefficients
!
      use Sub, only: curl, cross_mn, finalize_aver
      use Diagnostics, only: ysum_mn_name_xz_npar,xysum_mn_name_z_npar
      use Hydro, only: calc_pencils_hydro
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
!
      real, dimension (nx,3) :: btest,uxbtest
      integer :: jtest,j, nl
      logical :: headtt_save
      real :: fac
      type (pencil_case) :: p
      logical, dimension(npencils) :: lpenc_loc
!
      logical :: need_output
!
!  In this routine we will reset headtt after the first pencil,
!  so we need to reset it afterwards.
!
      headtt_save=headtt
      fac=1./nygrid
      need_output = (ldiagnos .and. needed2d(1)) .or. &
                    (l2davgfirst .and. needed2d(2))
!
!  do each of the 9 test fields at a time
!  but exclude redundancies, e.g. if the averaged field lacks x extent.
!  Note: the same block of lines occurs again further up in the file.
!
      uxbtestm=0.
      if (.not.lsoca .or. need_output) then
!
        lpenc_loc = .false.; lpenc_loc(i_uu)=.true.
!
        do jtest=1,njtest
!
          iaxtest=iaatest+3*(jtest-1)
          iaztest=iaxtest+2

          do n=n1,n2
!
            nl=n-n1+1
!
            do m=m1,m2
!
              call curl(f,iaxtest,btest)
              call calc_pencils_hydro(f,p,lpenc_loc)
              call cross_mn(p%uu,btest,uxbtest)
!
!  without SOCA, the alpha tensor is anisotropic even for the standard
!  Roberts flow. To check that this is because of the averaging that
!  enters, we allow ourselves to turn it off with lflucts_with_xyaver=.true.
!  It is off by default.
!
              do j=1,3
                if (lflucts_with_xyaver) then
                  uxbtestm(:,nl,j,jtest)=spread(sum( &
                    uxbtestm(:,nl,j,jtest)+fac*uxbtest(:,j),1),1,nx)/nx
                else
                  uxbtestm(:,nl,j,jtest)= &
                    uxbtestm(:,nl,j,jtest)+fac*uxbtest(:,j)
                endif
              enddo
              headtt=.false.
            enddo
          enddo
        enddo
        call finalize_aver(nprocy,2,uxbtestm)
      endif
!
!  reset headtt
!
      headtt=headtt_save
!
      if (need_output) call calc_coefficients(idiags(1:idiag_base_end),idiags_z(1:idiag_base_end),idiags_xz(1:idiag_base_end), &
                                              idiags(idiag_Eij_start:idiag_Eij_end),idiags_z(idiag_Eij_start:idiag_Eij_end),   &
                                              idiags_xz(idiag_Eij_start:idiag_Eij_end), &
                                              idiags(idiag_base_end+1:idiag_base_end+4),idiags(idiag_base_end+5:idiag_base_end+8),&
                                              uxbtestm,Minv,ysum_mn_name_xz_npar,xysum_mn_name_z_npar,  &
                                              twod_need_1d,twod_need_2d,needed2d,ny)
!
    endsubroutine testfield_after_boundary
!***********************************************************************
    subroutine set_bbtest(bbtest,jtest)
!
!  set testfield
!
!   3-jun-05/axel: coded
!
      real, dimension (nx,3) :: bbtest
      integer :: jtest
!
      intent(in)  :: jtest
      intent(out) :: bbtest
!
      integer :: nl
!
!  set bbtest for each of the 9 cases
!
      nl = n-n1+1
!
      select case (jtest)
!
      case (1); bbtest(:,1)=cx*cz(nl); bbtest(:,2)=0.; bbtest(:,3)=0.
      case (2); bbtest(:,1)=sx*cz(nl); bbtest(:,2)=0.; bbtest(:,3)=0.
      case (3); bbtest(:,1)=cx*sz(nl); bbtest(:,2)=0.; bbtest(:,3)=0.
!
      case (4); bbtest(:,1)=0.; bbtest(:,2)=cx*cz(nl); bbtest(:,3)=0.
      case (5); bbtest(:,1)=0.; bbtest(:,2)=sx*cz(nl); bbtest(:,3)=0.
      case (6); bbtest(:,1)=0.; bbtest(:,2)=cx*sz(nl); bbtest(:,3)=0.
!
      case (7); bbtest(:,1)=0.; bbtest(:,2)=0.; bbtest(:,3)=cx*cz(nl)
      case (8); bbtest(:,1)=0.; bbtest(:,2)=0.; bbtest(:,3)=sx*cz(nl)
      case (9); bbtest(:,1)=0.; bbtest(:,2)=0.; bbtest(:,3)=cx*sz(nl)
!
      case default; bbtest=0.
!
      endselect
!
    endsubroutine set_bbtest
!***********************************************************************
    subroutine set_bbtest_alt(bbtest,jtest)
!
!  tries of alternative testfields
!
!  23-oct-13/MR: coded
!  30-oct-13/MR: added another alternative
!  
      real, dimension (nx,3) :: bbtest
      integer :: jtest
!
      intent(in)  :: jtest
      intent(out) :: bbtest
!
      integer :: nl
!
!  set bbtest for each of the 9 cases
!
      nl = n-n1+1
!
      select case (jtest)
!
      case (1); bbtest(:,1)=cx*cz(nl)+sx*sz(nl); bbtest(:,2)=0.; bbtest(:,3)=0.
      case (2); bbtest(:,1)=sx*cz(nl)          ; bbtest(:,2)=0.; bbtest(:,3)=0.
!      case (2); bbtest(:,1)=sx*cz(nl)+cx*sz(nl); bbtest(:,2)=0.; bbtest(:,3)=0.
      case (3); bbtest(:,1)=cx*sz(nl);           bbtest(:,2)=0.; bbtest(:,3)=0.
!      case (3); bbtest(:,1)=cx*sz(nl)+sx*sz(nl); bbtest(:,2)=0.; bbtest(:,3)=0.
!
      case (4); bbtest(:,1)=0.; bbtest(:,2)=cx*cz(nl)+sx*sz(nl); bbtest(:,3)=0.
      case (5); bbtest(:,1)=0.; bbtest(:,2)=sx*cz(nl)          ; bbtest(:,3)=0.
!      case (5); bbtest(:,1)=0.; bbtest(:,2)=sx*cz(nl)+cx*sz(nl); bbtest(:,3)=0.
      case (6); bbtest(:,1)=0.; bbtest(:,2)=cx*sz(nl);           bbtest(:,3)=0.
!      case (6); bbtest(:,1)=0.; bbtest(:,2)=cx*sz(nl)+sx*sz(nl); bbtest(:,3)=0.
!
      case (7); bbtest(:,1)=0.; bbtest(:,2)=0.; bbtest(:,3)=cx*cz(nl)+sx*sz(nl)
      case (8); bbtest(:,1)=0.; bbtest(:,2)=0.; bbtest(:,3)=sx*cz(nl)
!      case (8); bbtest(:,1)=0.; bbtest(:,2)=0.; bbtest(:,3)=sx*cz(nl)+cx*sz(nl)
      case (9); bbtest(:,1)=0.; bbtest(:,2)=0.; bbtest(:,3)=cx*sz(nl)
!      case (9); bbtest(:,1)=0.; bbtest(:,2)=0.; bbtest(:,3)=cx*sz(nl)+sx*sz(nl)
!
      case default; bbtest=0.
!
      endselect
!
    endsubroutine set_bbtest_alt
!***********************************************************************
    subroutine set_bbtest_alt_II(bbtest,jtest)
!
!  tries of alternative testfields
!
!  23-oct-13/MR: coded
!  30-oct-13/MR: added another alternative
!  
      real, dimension (nx,3) :: bbtest
      integer :: jtest
!
      intent(in)  :: jtest
      intent(out) :: bbtest
!
      integer :: nl
!
!  set bbtest for each of the 9 cases
!
      nl = n-n1+1
!
      select case (jtest)
!
      case (1); bbtest(:,1)=cx*cz(nl)+sx*sz(nl); bbtest(:,2)=0.; bbtest(:,3)=0.
      case (2); bbtest(:,1)=sx*cz(nl)+cx*sz(nl); bbtest(:,2)=0.; bbtest(:,3)=0.
      case (3); bbtest(:,1)=cx*sz(nl)+sx*sz(nl); bbtest(:,2)=0.; bbtest(:,3)=0.
!
      case (4); bbtest(:,1)=0.; bbtest(:,2)=cx*cz(nl)+sx*sz(nl); bbtest(:,3)=0.
      case (5); bbtest(:,1)=0.; bbtest(:,2)=sx*cz(nl)+cx*sz(nl); bbtest(:,3)=0.
      case (6); bbtest(:,1)=0.; bbtest(:,2)=cx*sz(nl)+sx*sz(nl); bbtest(:,3)=0.
!
      case (7); bbtest(:,1)=0.; bbtest(:,2)=0.; bbtest(:,3)=cx*cz(nl)+sx*sz(nl)
      case (8); bbtest(:,1)=0.; bbtest(:,2)=0.; bbtest(:,3)=sx*cz(nl)+cx*sz(nl)
      case (9); bbtest(:,1)=0.; bbtest(:,2)=0.; bbtest(:,3)=cx*sz(nl)+sx*sz(nl)
!
      case default; bbtest=0.
!
      endselect
!
    endsubroutine set_bbtest_alt_II
!***********************************************************************
    subroutine set_bbtest2(bbtest,jtest)
!
!  set alternative testfield
!
!  10-jun-05/axel: adapted from set_bbtest
!  27-aug-13/MR: removed unneeded calculations
!
      real, dimension (nx,3) :: bbtest
      integer :: jtest
!
      intent(in)  :: jtest
      intent(out) :: bbtest
!
      real, dimension (nx) :: xz
      real :: szl
      integer :: nl
!
      nl = n-n1+1
      szl=sz(nl)
      xz=cx*cz(nl)
!
!  set bbtest for each of the 9 cases
!
      select case (jtest)
      case (1); bbtest(:,1)=xz ; bbtest(:,2)=0. ; bbtest(:,3)=0.
      case (2); bbtest(:,1)=0. ; bbtest(:,2)=xz ; bbtest(:,3)=0.
      case (3); bbtest(:,1)=0. ; bbtest(:,2)=0. ; bbtest(:,3)=xz
      case (4); bbtest(:,1)=szl; bbtest(:,2)=0. ; bbtest(:,3)=0.
      case (5); bbtest(:,1)=0. ; bbtest(:,2)=szl; bbtest(:,3)=0.
      case (6); bbtest(:,1)=0. ; bbtest(:,2)=0. ; bbtest(:,3)=szl
      case (7); bbtest(:,1)=sx ; bbtest(:,2)=0. ; bbtest(:,3)=0.
      case (8); bbtest(:,1)=0. ; bbtest(:,2)=sx ; bbtest(:,3)=0.
      case (9); bbtest(:,1)=0. ; bbtest(:,2)=0. ; bbtest(:,3)=sx
      case default; bbtest=0.
      endselect
!
    endsubroutine set_bbtest2
!***********************************************************************
    subroutine set_bbtest4(bbtest,jtest)
!
!  set testfield using constant and linear functions
!
!  15-jun-05/axel: adapted from set_bbtest3
!  27-aug-13/MR: made zz scalar
!  20-oct-13/MR: changed order of testfields
!
      real, dimension (nx,3) :: bbtest
      integer :: jtest
!
      intent(in)  :: jtest
      intent(out) :: bbtest
!
      real, dimension (nx) :: xx
      real :: zz
!
      xx=x(l1:l2)
      zz=z(n)
!
!  set bbtest for each of the 9 cases
!
      select case (jtest)
      case (1); bbtest(:,1)=1.; bbtest(:,2)=0.; bbtest(:,3)=0.
      case (2); bbtest(:,1)=xx; bbtest(:,2)=0.; bbtest(:,3)=0.
      case (3); bbtest(:,1)=zz; bbtest(:,2)=0.; bbtest(:,3)=0.
!
      case (4); bbtest(:,1)=0.; bbtest(:,2)=1.; bbtest(:,3)=0.
      case (5); bbtest(:,1)=0.; bbtest(:,2)=xx; bbtest(:,3)=0.
      case (6); bbtest(:,1)=0.; bbtest(:,2)=zz; bbtest(:,3)=0.
!
      case (7); bbtest(:,1)=0.; bbtest(:,2)=0.; bbtest(:,3)=1.
      case (8); bbtest(:,1)=0.; bbtest(:,2)=0.; bbtest(:,3)=xx
      case (9); bbtest(:,1)=0.; bbtest(:,2)=0.; bbtest(:,3)=zz
!
      case default; bbtest=0.
      endselect
!
    endsubroutine set_bbtest4
!***********************************************************************
    subroutine set_bbtest3(bbtest,jtest)
!
!  set alternative testfield
!
!  10-jun-05/axel: adapted from set_bbtest
!  27-aug-13/MR: removed unneeded calculations
!
      real, dimension (nx,3) :: bbtest
      integer :: jtest
!
      intent(in)  :: jtest
      intent(out) :: bbtest
!
      integer :: nl
!
!  set bbtest for each of the 9 cases
!
      nl=n-n1+1

      select case (jtest)
      case (1); bbtest(:,1)=1.    ; bbtest(:,2)=0.    ; bbtest(:,3)=0.
      case (2); bbtest(:,1)=0.    ; bbtest(:,2)=1.    ; bbtest(:,3)=0.
      case (3); bbtest(:,1)=0.    ; bbtest(:,2)=0.    ; bbtest(:,3)=1.
      case (4); bbtest(:,1)=sz(nl); bbtest(:,2)=0.    ; bbtest(:,3)=0.
      case (5); bbtest(:,1)=0.    ; bbtest(:,2)=sz(nl); bbtest(:,3)=0.
      case (6); bbtest(:,1)=0.    ; bbtest(:,2)=0.    ; bbtest(:,3)=cz(nl)
      case (7); bbtest(:,1)=cx    ; bbtest(:,2)=0.    ; bbtest(:,3)=0.
      case (8); bbtest(:,1)=0.    ; bbtest(:,2)=sx    ; bbtest(:,3)=0.
      case (9); bbtest(:,1)=0.    ; bbtest(:,2)=0.    ; bbtest(:,3)=sx
      case default; bbtest=0.
      endselect
!
    endsubroutine set_bbtest3
!***********************************************************************
    subroutine rprint_testfield(lreset,lwrite)
!
!  reads and registers print parameters relevant for testfield fields
!
!   3-jun-05/axel: adapted from rprint_magnetic
!  26-feb-13/MR  : output of ntestfield in index.pro added
!   6-mar-13/MR  : alternative parse_name used
!
      use Diagnostics, only: parse_name
!
      logical :: lreset
      logical, optional :: lwrite
!
      integer :: iname,inamez,inamexz,i
!
!  reset everything in case of RELOAD
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        idiags=0
        idiags_z=0
        idiags_xz=0
      endif
!
!  check for those quantities that we want to evaluate online
!
!  volume averages
!
      do iname=1,nname
        do i=1,size(idiags)
          call parse_name(iname,cname,cform,cdiags(i),idiags(i))
        enddo
      enddo
!
!  xy-averages
!
      do inamez=1,nnamez
        do i=1,size(idiags_z)
          call parse_name(inamez,cnamez,cformz,trim(cdiags(i))//'z',idiags_z(i))
        enddo
      enddo
!
!  y-averages
!
      do inamexz=1,nnamexz
        do i=1,size(idiags_xz)
          call parse_name(inamexz,cnamexz,cformxz,trim(cdiags(i))//'xz',idiags_xz(i))
        enddo
      enddo
!
      needed2d = diagnos_interdep(idiags(1:idiag_base_end),idiags_z(1:idiag_base_end), &
                                  idiags_xz(1:idiag_base_end),twod_need_1d,twod_need_2d)
!
    endsubroutine rprint_testfield
!***********************************************************************
endmodule Testfield
