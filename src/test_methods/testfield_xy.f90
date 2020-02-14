! $Id: testfield_xy.f90 20934 2013-08-22 19:27:54Z mreinhardt@nordita.org $
!
!  This modules deals with all aspects of testfield fields; if no
!  testfield fields are invoked, a corresponding replacement dummy
!  routine is used instead which absorbs all the calls to the
!  testfield relevant subroutines listed in here.
!
!  Derived from testfield_xz by P. Käpylä, A. Brandenburg and 
!  M. Rheinhardt
!   27-aug-13/pete: adapted from testfield_xz.f90
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: ltestfield_xy = .true.
! CPARAM logical, parameter :: ltestfield_z = .false.
! CPARAM logical, parameter :: ltestfield_xz  = .false.
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
  real    :: ktestfield_x=1., ktestfield_y=1., xx0=0., yy0=0.
!
  namelist /testfield_run_pars/       &
           B_ext,                     &
           reinitialize_aatest,       &
           lsoca,                     &
           etatest,                   &
           itestfield,                &
           ktestfield_x,              &
           ktestfield_y,              &
           xx0,                       &
           yy0,                       &
           daainit,                   &
           linit_aatest,              &
           rescale_aatest,            &
           lflucts_with_xyaver
!
! diagnostic variables
!
  integer, parameter :: n_cdiags = 62, len_cdiags = 8
  character(LEN=len_cdiags), dimension(n_cdiags) :: cdiags = &
  (/ 'alp11   ','alp21   ','alp31   ','alp12   ','alp22   ','alp32   ','alp13   ','alp23   ','alp33   ',&   
! DIAG_DOC: $\alpha_{ij}$       
     'eta111  ','eta211  ','eta311  ','eta121  ','eta221  ','eta321  ','eta131  ','eta231  ','eta331  ',&   
! DIAG_DOC: $\eta_{ijk}$
     'eta112  ','eta212  ','eta312  ','eta122  ','eta222  ','eta322  ','eta132  ','eta232  ','eta332  ',&
     'alp11cc ','alp11cs ','alp11sc ','alp11ss ','eta122cc','eta122cs','eta122sc','eta122ss'           ,&   
! DIAG_DOC: $\alpha_{11,\rm hh},$ 
                                                                                                            
! Diag_DOC: $\eta_{122,\rm hh}, {\rm h}={\rm c,s}$
     'E11     ','E21     ','E31     ','E12     ','E22     ','E32     ','E13     ','E23     ','E33     ',& 
! DIAG_DOC: ${\cal E}^i_j$
     'E14     ','E24     ','E34     ','E15     ','E25     ','E35     ','E16     ','E26     ','E36     ',&   
     'E17     ','E27     ','E37     ','E18     ','E28     ','E38     ','E19     ','E29     ','E39     ' /) 
!
  integer, dimension(n_cdiags):: idiags=0, idiags_x=0, idiags_xy=0
  integer, parameter :: idiag_base_end=27, idiag_Eij_start=36, idiag_Eij_end=idiag_Eij_start+27-1
!
  integer, dimension(4) :: idiag_alp11h, idiag_eta122h            
  equivalence(idiags(idiag_base_end+1),idiag_alp11h), (idiags(idiag_base_end+5),idiag_eta122h)      
! alias names for selected diagnostics
!
!  work variables
!
  real, dimension (nx)              :: cx,sx
  real, dimension (ny)              :: cy,sy
  real, dimension (nx,ny,3,3)       :: Minv
  real, dimension (nx,ny,3,njtest)  :: uxbtestm
  logical, dimension(idiag_base_end):: twod_need_1d, twod_need_2d
  logical, dimension(2)             :: needed2d
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
!                  set lcalc_uumeanxy=.true., itestfield now string     
!  20-oct-13/MR  : calculation of means in hydro triggered
! 
      use Hydro, only: lcalc_uumeanxy, calc_means_hydro
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call initialize_testfield_general(f)
           
      !!!if (reinitialize_aatest) f(:,:,:,iaatest:iaatest+ntestfield-1)=0.  !!! TBC
!
! calculate inverse matrix for determination of the turbulent coefficients
!
      call calc_inverse_matrix(x(l1:l2),y(m1:m2),ktestfield_x,ktestfield_y,xx0,yy0,Minv,cx,sx,cy,sy)
!
      if (lhydro_kinematic) call calc_means_hydro(f)
!
!  write testfield information to a file (for convenient post-processing)
!
      if (lroot) then
        open(1,file=trim(datadir)//'/testfield_info.dat',STATUS='unknown')
        write(1,'(a,i1)') 'lsoca='  ,merge(1,0,lsoca)
        write(1,'(a,a)') 'itestfield=',itestfield
        write(1,'(a,2(f3.0))') 'ktestfield_x,y=', ktestfield_x, ktestfield_y
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
!
      use Hydro, only: uumxy
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      intent(in)    :: f,p
      intent(inout) :: df
!
      integer :: ml
!
        ml=m-m1+1
        select case (itestfield)
          case ('1','1-alt') ; call rhs_daatest(f,df,p,uumxy(l1:l2,m,:),uxbtestm(:,ml,:,:),set_bbtest )
          case ('2')         ; call rhs_daatest(f,df,p,uumxy(l1:l2,m,:),uxbtestm(:,ml,:,:),set_bbtest2)
          case ('3')         ; call rhs_daatest(f,df,p,uumxy(l1:l2,m,:),uxbtestm(:,ml,:,:),set_bbtest3)
          case ('4','linear'); call rhs_daatest(f,df,p,uumxy(l1:l2,m,:),uxbtestm(:,ml,:,:),set_bbtest4)
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
!  04-oct-13/MR  : removed p from parameter list, introd calculation of
!                  hydro pencils (restricted); simplified communication
!
      use Sub, only: curl, cross_mn, finalize_aver
      use Hydro, only:  calc_pencils_hydro
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
!
      real, dimension (nx,3) :: btest,uxbtest
      integer :: jtest,j,ml
      logical :: headtt_save
      real :: fac
      type (pencil_case) :: p
      logical, dimension(npencils) :: lpenc_loc
!
      logical :: need_output
!
!  In this routine we will reset headtt after the first pencil,
!  so we need to restore it afterwards.
!
      headtt_save=headtt
      fac=1./nzgrid
      need_output = (ldiagnos .and. needed2d(1)) .or. &
                    (l2davgfirst .and. needed2d(2))
!
!  do each of the 9 test fields at a time
!  but exclude redundancies, e.g. if the averaged field lacks x extent.
!  Note: the same block of lines occurs again further up in the file.
!
      uxbtestm=0.
!
      if (.not.lsoca .or. need_output) then

        lpenc_loc = .false.; lpenc_loc(i_uu)=.true.
!
        do jtest=1,njtest
!
          iaxtest=iaatest+3*(jtest-1)
          iaztest=iaxtest+2

          do m=m1,m2
!
            ml=m-m1+1
!
            do n=n1,n2
!
              call curl(f,iaxtest,btest)
              call calc_pencils_hydro(f,p,lpenc_loc)
!
! U x btest = (\mean{U} + u') x btest
!
              call cross_mn(p%uu,btest,uxbtest)
!
!  without SOCA, the alpha tensor is anisotropic even for the standard
!  Roberts flow. To check that this is because of the averaging that
!  enters, we allow ourselves to turn it off with lflucts_with_xyaver=.true.
!  It is off by default.
!
              do j=1,3
                if (lflucts_with_xyaver) then    !!! TBC
                  uxbtestm(:,ml,j,jtest)=spread(sum( &
                    uxbtestm(:,ml,j,jtest)+fac*uxbtest(:,j),1),1,nx)/nx
                else
                  uxbtestm(:,ml,j,jtest)= &
                    uxbtestm(:,ml,j,jtest)+fac*uxbtest(:,j)
                endif
              enddo
              headtt=.false.
            enddo
          enddo
        enddo
!      
        call finalize_aver(nprocz,3,uxbtestm)
!
      endif
!
!  reset headtt
!
      headtt=headtt_save
!
      if (need_output) call calc_coeffs
!
    endsubroutine testfield_after_boundary
!***********************************************************************
    subroutine calc_coeffs
!
!  interface to enable use of of calc_coefficients originally developed for testfield_xz
!
!  28-aug-13/MR: coded
!  02-dec-13/MR: fixed bugs: missing mapping in twod_need_1d, twod_need_2d added, wrong 'backmapping' removed
!
      use Diagnostics, only: zsum_mn_name_xy_mpar_scal,yzsum_mn_name_x_mpar
!
      integer, dimension(idiag_base_end) :: idiags_map
      integer :: i,j
      real, dimension (nx,ny,3,njtest)  :: tmp

!  consider rotation of coordinate system about x axis by pi/2, that is: x --> x, y --> z, z --> -y;
!
!      turbulent coefficients calculated here:
!
!       1            2          3          4          5          6          7          8          9
!     'alp11   ','alp21   ','alp31   ','alp12   ','alp22   ','alp32   ','alp13   ','alp23   ','alp33   ' 
!       10           11         12         13         14         15         16         17         18     
!     'eta111  ','eta211  ','eta311  ','eta121  ','eta221  ','eta321  ','eta131  ','eta231  ','eta331  '
!       19           20         21         22         23         24         25         26         27
!     'eta112  ','eta212  ','eta312  ','eta122  ','eta222  ','eta322  ','eta132  ','eta232  ','eta332 

!      correspondingly in testfield_xz:
!
!       1            2          3          4          5          6          7          8          9
!      'alp11   ','alp21   ','alp31   ','alp12   ','alp22   ','alp32   ','alp13   ','alp23   ','alp33   '
!       10           11         12         13         14         15         16         17         18     
!      'eta111  ','eta211  ','eta311  ','eta121  ','eta221  ','eta321  ','eta131  ','eta231  ','eta331  '
!       19           20         21         22         23         24         25         26         27
!      'eta113  ','eta213  ','eta313  ','eta123  ','eta223  ','eta323  ','eta133  ','eta233  ','eta333  '

!     hence the following mapping by idiags_map, where a negative sign indicates that the sign of the coefficient must
!     be inverted:
!
!     running index   1   2   3     4    5    6     7    8   9    10   11   12    13  14   15   16   17  18
!
      idiags_map = (/ 1,  3, -2,    7,   9,  -8,   -4,  -6,  5,   10,  12, -11,   16, 18, -17, -13, -15, 14, &
!
!     running index  19   20   21   22  23   24   25   26  27
!
!!!                     19, 21, -20,-25, 27,-26, 22,-24, 23 /)
                     19,  21, -20,  25, 27, -26, -22, -24, 23 /)      !tb checked
!
      tmp=uxbtestm(:,:,(/1,3,2/),:); tmp(:,:,2,:)=-tmp(:,:,2,:)
!      tmp=uxbtestm(:,:,(/1,3,2/),:); tmp(:,:,3,:)=-tmp(:,:,3,:)

      call calc_coefficients( idiags(abs(idiags_map)),idiags_x(abs(idiags_map)),idiags_xy(abs(idiags_map)), &
                              idiags(idiag_Eij_start:idiag_Eij_end),idiags_x(idiag_Eij_start:idiag_Eij_end),   &
                              idiags_xy(idiag_Eij_start:idiag_Eij_end), &
                              idiag_alp11h, idiag_eta122h, &
                              tmp,Minv,zsum_mn_name_xy_mpar_scal,yzsum_mn_name_x_mpar, &
                              twod_need_1d(abs(idiags_map)),twod_need_2d(abs(idiags_map)),needed2d,nz )
!
!  sign inversion if necessary
!
      if (ldiagnos .and. needed2d(1)) then
        where(idiags(1:idiag_base_end)/=0) fname(idiags(1:idiag_base_end)) =  sign(1.,float(idiags_map)) &
                                                                             *fname(idiags(1:idiag_base_end))
        do i=1,nprocz; do j=1,nz
          where(idiags_x(1:idiag_base_end)/=0) &
            fnamex(j,i,idiags_x(1:idiag_base_end)) = sign(1.,float(idiags_map))*fnamex(j,i,idiags_x(1:idiag_base_end))
        enddo; enddo
      endif
      
      if (l2davgfirst .and. needed2d(2)) then
        do i=1,ny; do j=1,nx
          where(idiags_xy(1:idiag_base_end)/=0) &
            fnamexy(idiags_xy(1:idiag_base_end),j,i) = sign(1.,float(idiags_map))*fnamexy(idiags_xy(1:idiag_base_end),j,i)
        enddo; enddo
      endif

    endsubroutine calc_coeffs
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
      integer :: ml
!
!  set bbtest for each of the 9 cases
!
      ml = m-m1+1
!
      select case (jtest)
!
      case (1); bbtest(:,1)=cx*cy(ml); bbtest(:,2)=0.; bbtest(:,3)=0.
      case (2); bbtest(:,1)=sx*cy(ml); bbtest(:,2)=0.; bbtest(:,3)=0.
      case (3); bbtest(:,1)=cx*sy(ml); bbtest(:,2)=0.; bbtest(:,3)=0.
!
      case (4); bbtest(:,1)=0.; bbtest(:,2)=0.; bbtest(:,3)=-cx*cy(ml)
      case (5); bbtest(:,1)=0.; bbtest(:,2)=0.; bbtest(:,3)=-sx*cy(ml)
      case (6); bbtest(:,1)=0.; bbtest(:,2)=0.; bbtest(:,3)=-cx*sy(ml)
!
      case (7); bbtest(:,1)=0.; bbtest(:,2)=cx*cy(ml); bbtest(:,3)=0.
      case (8); bbtest(:,1)=0.; bbtest(:,2)=sx*cy(ml); bbtest(:,3)=0.
      case (9); bbtest(:,1)=0.; bbtest(:,2)=cx*sy(ml); bbtest(:,3)=0.
!
      case default; bbtest=0.
!
      endselect
!
    endsubroutine set_bbtest
!***********************************************************************
    subroutine set_bbtest2(bbtest,jtest)
!
!  set alternative testfield
!
!  10-jun-05/axel: adapted from set_bbtest
!
      real, dimension (nx,3) :: bbtest
      integer :: jtest
!
      intent(in)  :: jtest
      intent(out) :: bbtest
!
      real, dimension (nx) :: xy
      real :: syl

      integer :: ml
!
      ml = m-m1+1
      syl=sy(ml)
      xy=cx*cy(ml)
!
!  set bbtest for each of the 9 cases
!
      select case (jtest)
      case (1); bbtest(:,1)=xy ; bbtest(:,2)=0. ; bbtest(:,3)=0.
      case (2); bbtest(:,1)=0. ; bbtest(:,2)=xy ; bbtest(:,3)=0.
      case (3); bbtest(:,1)=0. ; bbtest(:,2)=0. ; bbtest(:,3)=xy
      case (4); bbtest(:,1)=syl; bbtest(:,2)=0. ; bbtest(:,3)=0.
      case (5); bbtest(:,1)=0. ; bbtest(:,2)=syl; bbtest(:,3)=0.
      case (6); bbtest(:,1)=0. ; bbtest(:,2)=0. ; bbtest(:,3)=syl
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
!  20-oct-13/MR  : changed order of testfields
!
      real, dimension (nx,3) :: bbtest
      integer :: jtest
!
      intent(in) :: jtest
      intent(out):: bbtest
!
      real, dimension (nx) :: xx
      real :: yy

      xx=x(l1:l2)
      yy=y(m)
!
!  set bbtest for each of the 9 cases
!
      select case (jtest)
      case (1); bbtest(:,1)=1.; bbtest(:,2)=0.; bbtest(:,3)=0.
      case (2); bbtest(:,1)=xx; bbtest(:,2)=0.; bbtest(:,3)=0.
      case (3); bbtest(:,1)=yy; bbtest(:,2)=0.; bbtest(:,3)=0.
!
      case (4); bbtest(:,1)=0.; bbtest(:,2)=0.; bbtest(:,3)=-1.
      case (5); bbtest(:,1)=0.; bbtest(:,2)=0.; bbtest(:,3)=-xx
      case (6); bbtest(:,1)=0.; bbtest(:,2)=0.; bbtest(:,3)=-yy
!
      case (7); bbtest(:,1)=0.; bbtest(:,2)=1.; bbtest(:,3)=0.
      case (8); bbtest(:,1)=0.; bbtest(:,2)=xx; bbtest(:,3)=0.
      case (9); bbtest(:,1)=0.; bbtest(:,2)=yy; bbtest(:,3)=0.
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
!
      real, dimension (nx,3) :: bbtest
      integer :: jtest
!
      intent(in)  :: jtest
      intent(out) :: bbtest
!
      integer :: ml
!
      ml=m-m1+1
!
!  set bbtest for each of the 9 cases
!
      select case (jtest)
      case (1); bbtest(:,1)=1.    ; bbtest(:,2)=0.    ; bbtest(:,3)=0.
      case (2); bbtest(:,1)=0.    ; bbtest(:,2)=1.    ; bbtest(:,3)=0.
      case (3); bbtest(:,1)=0.    ; bbtest(:,2)=0.    ; bbtest(:,3)=1.
      case (4); bbtest(:,1)=sy(ml); bbtest(:,2)=0.    ; bbtest(:,3)=0.
      case (5); bbtest(:,1)=0.    ; bbtest(:,2)=sy(ml); bbtest(:,3)=0.
      case (6); bbtest(:,1)=0.    ; bbtest(:,2)=0.    ; bbtest(:,3)=cy(ml)
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
      integer :: iname,inamex,inamexy,i
!
!  reset everything in case of RELOAD
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        idiags=0
        idiags_x=0
        idiags_xy=0
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
!  yz-averages
!
      do inamex=1,nnamex
        do i=1,size(idiags_x)
          call parse_name(inamex,cnamex,cformx,trim(cdiags(i))//'x',idiags_x(i))
        enddo
      enddo
!
!  z-averages
!
      do inamexy=1,nnamexy
        do i=1,size(idiags_xy)
          call parse_name(inamexy,cnamexy,cformxy,trim(cdiags(i))//'xy',idiags_xy(i))
        enddo
      enddo
!
      needed2d = diagnos_interdep(idiags(1:idiag_base_end),idiags_x(1:idiag_base_end), &
                                  idiags_xy(1:idiag_base_end),twod_need_1d,twod_need_2d)
!
    endsubroutine rprint_testfield
!***********************************************************************
endmodule Testfield
