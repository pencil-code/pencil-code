! $Id$
!
!  This module contains useful subroutines.
!
!  Rules:
!
!    - Please do not put very specific subroutines here. If a subroutine
!      is only needed by a single module, put it directly in that module.
!
!    - Please DO NOT use large arrays or global arrays
!      [e.g. of size (mx,my,mz) or (nxgrid,nygrid,nzgrid)]
!
module Sub
!
  use Cdata
  use Cparam
  use Messages
!
  implicit none
!
  private
!
  public :: step,stepdown,der_stepdown
  public :: ylm, ylm_other
  public :: kronecker_delta
!
  public :: identify_bcs, parse_bc, parse_bc_rad, parse_bc_radg
  public :: inverse_parse_bc
!
  public :: poly
  public :: blob, vecout
  public :: cubic_step, cubic_der_step, quintic_step, quintic_der_step, erfunc
  public :: sine_step, interp1
  public :: hypergeometric2F1
  public :: gamma_function
  public :: one_minus_exp
!
  public :: get_nseed
  public :: get_where
!
  public :: grad, grad5, div, div_mn, curl, curli, curl_mn, curl_other
  public :: curl_horizontal
  public :: div_other
  public :: gij, g2ij, gij_etc
  public :: gijl_symmetric
  public :: bij_tilde
  public :: symmetrise3x3_ut2lt
  public :: der_step
  public :: der6_step
  public :: u_dot_grad, h_dot_grad, h_dot_grad_vec
  public :: u_dot_grad_alt
  public :: nou_dot_grad_scl
  public :: u_dot_grad_mat
  public :: del2, del2v, del2v_etc, del2fj,d2fi_dxj,del2fi_dxjk
  public :: del2m3x3_sym
  public :: del4v, del4, del2vi_etc, del4graddiv
  public :: del6v, del6, del6_other, del6fj, del6fjv, del6_strict
  public :: gradf_upw1st, doupwind
  public :: matrix2linarray, linarray2matrix
  public :: ScalarTripleProduct
  public :: det3X3mat,Inv2_3X3mat
!
  public :: dot, dot2, dot_mn, dot_mn_sv, dot_mn_sm, dot2_mn, dot_add, dot_sub, dot2fj
  public :: dot_mn_vm, dot_mn_vm_trans, div_mn_2tensor, trace_mn
  public :: dyadic2, dyadic2_other
  public :: cross, cross_mn, cross_mixed
  public :: sum_mn, max_mn
  public :: multsv, multsv_add, multsv_mn, multsv_mn_add
  public :: multvs, multvv_mat
  public :: multmm_sc
  public :: multm2, multm2_mn, multm2_sym, multm2_sym_mn
  public :: multmv, multmv_mn, multmv_transp
  public :: mult_matrix
!
  public :: read_line_from_file, control_file_exists
  public :: noform
!
  public :: update_snaptime, read_snaptime
  public :: shift_dt, set_dt
  public :: parse_shell
  public :: get_radial_distance, power_law
!
  public :: max_for_dt,unit_vector
!
  public :: write_dx_general, rdim
  public :: write_xprof, write_zprof, remove_prof
!
  public :: tensor_diffusion_coef
!
  public :: smooth_kernel, despike
  public :: smooth, get_smooth_kernel
!
  public :: ludcmp, lubksb
  public :: bspline_basis, bspline_interpolation, bspline_precondition
!
  public :: gij_psi, gij_psi_etc
  public :: xlocation, ylocation, zlocation, location_in_proc, position
  public :: register_report_aux
  public :: fourier_single_mode
  public :: remove_mean,global_mean
  public :: insert
  public :: find_max_fvec, find_rms_fvec, find_xyrms_fvec
  public :: finalize_aver
  public :: meanyz
  public :: calc_all_diff_fluxes
  public :: periodic_fold_back
  public :: lower_triangular_index
!
  public :: vec_dot_3tensor
  public :: traceless_strain,calc_sij2
  public :: calc_del6_for_upwind
!
  public :: remove_mean_value
  public :: stagger_to_base_interp_1st, stagger_to_base_interp_3rd
!
  interface poly                ! Overload the `poly' function
    module procedure poly_0
    module procedure poly_1
    module procedure poly_3
  endinterface
!
  interface grad                 ! Overload the `grad' function
    module procedure grad_main   ! grad of an 'mvar' variable
    module procedure grad_other  ! grad of another field (mx,my,mz)
  endinterface
!
  interface del2                 ! Overload the `del2' function
    module procedure del2_main
    module procedure del2_other
  endinterface
!
  interface cross
    module procedure cross_mn
    module procedure cross_0
    module procedure cross_mixed
  endinterface
!
  interface u_dot_grad
    module procedure u_dot_grad_scl
    module procedure u_dot_grad_vec
  endinterface
!
  interface u_dot_grad_alt
    module procedure u_dot_grad_scl_alt
    module procedure u_dot_grad_vec_alt
  endinterface
!
  interface h_dot_grad
    module procedure h_dot_grad_scl
    module procedure h_dot_grad_vec
  endinterface
!
  interface dot
    module procedure dot_mn_sv
    module procedure dot_mn
    module procedure dot_0
  endinterface
!
  interface dot2
    module procedure dot2_mn
    module procedure dot2_0
  endinterface
!
  interface dot_add
    module procedure dot_mn_add
  endinterface
!
  interface dot_sub
    module procedure dot_mn_sub
  endinterface
!
  interface multsv
    module procedure multsv_mn
  endinterface
!
  interface multsv_add
    module procedure multsv_add_mn
  endinterface
!
  interface multvs
    module procedure multvs_mn
  endinterface
!
  interface multvv_mat
    module procedure multvv_mat_mn
  endinterface
!
  interface multmm_sc
    module procedure multmm_sc_mn
  endinterface
!
  interface multm2
    module procedure multm2_mn
  endinterface
!
  interface multm2_sym
    module procedure multm2_sym_mn
  endinterface
!
  interface multmv_transp
    module procedure multmv_mn_transp
  endinterface
!
  interface multmv
    module procedure multmv_mn
  endinterface
!
  interface max_for_dt
    module procedure max_for_dt_nx_nx
    module procedure max_for_dt_1_nx
    module procedure max_for_dt_1_1_1_nx
  endinterface
!
  interface step
    module procedure step_scalar
    module procedure step_vector
  endinterface
!
  interface cubic_step
    module procedure cubic_step_pt
    module procedure cubic_step_mn
  endinterface
!
  interface cubic_der_step
    module procedure cubic_der_step_pt
    module procedure cubic_der_step_mn
  endinterface
!
  interface quintic_step
    module procedure quintic_step_pt
    module procedure quintic_step_mn
  endinterface
!
  interface quintic_der_step
    module procedure quintic_der_step_pt
    module procedure quintic_der_step_mn
  endinterface
!
  interface sine_step
    module procedure sine_step_pt
    module procedure sine_step_mn
  endinterface
!
  interface power_law
     module procedure power_law_pt
     module procedure power_law_mn
  endinterface
!
  interface insert                 ! Overload the 'insert' function
    module procedure insert_carray
    module procedure insert_carray_mult
    module procedure insert_rarray
  endinterface
!
  interface finalize_aver
    module procedure finalize_aver_1D
    module procedure finalize_aver_2D
    module procedure finalize_aver_3D
    module procedure finalize_aver_4D
  endinterface finalize_aver
!
  interface meanyz
     module procedure meanyz_s
     module procedure meanyz_v
  endinterface
!
!  extended intrinsic operators to do some scalar/vector pencil arithmetic
!  Tobi: Array valued functions do seem to be slower than subroutines,
!        hence commented out for the moment.
!
!  public :: operator(*),operator(+),operator(/),operator(-)
!
!  interface operator (*)
!    module procedure pencil_multiply1
!    module procedure pencil_multiply2
!  endinterface
!
!  interface operator (+)
!    module procedure pencil_add1
!    module procedure pencil_add2
!  endinterface
!
!  interface operator (/)
!    module procedure pencil_divide1
!    module procedure pencil_divide2
!  endinterface
!
!  interface operator (-)
!    module procedure pencil_subtract1
!    module procedure pencil_subtract2
!  endinterface
!
!ajwm Commented pending a C replacement
!  INTERFACE getenv
!    SUBROUTINE GETENV (VAR, VALUE)
!      CHARACTER(LEN=*) VAR, VALUE
!    endsubroutine
!  END INTERFACE
!
  real, dimension(7,7,7), parameter :: smth_kernel = reshape((/ &
 6.03438e-15,9.07894e-11,1.24384e-08,5.46411e-08,1.24384e-08,9.07894e-11,5.03438e-15,9.07894e-11,2.21580e-07,9.14337e-06,&
 2.69243e-05,9.14337e-06,2.21580e-07,9.07894e-11, 1.24384e-08, 9.14337e-06, 0.000183649, 0.000425400, 0.000183649, 9.14337e-06,&
 1.24384e-08,5.46411e-08,2.69243e-05,0.000425400, 0.000909623, 0.000425400, 2.69243e-05, 5.46411e-08, 1.24384e-08, 9.14337e-06,&
 0.000183649,0.000425400,0.000183649,9.14337e-06, 1.24384e-08, 9.07894e-11, 2.21580e-07, 9.14337e-06, 2.69243e-05, 9.14337e-06,&
 2.21580e-07,9.07894e-11,5.03438e-15,9.07894e-11, 1.24384e-08, 5.46411e-08, 1.24384e-08, 9.07894e-11, 5.03438e-15, 9.07894e-11,&
 2.21580e-07,9.14337e-06,2.69243e-05,9.14337e-06, 2.21580e-07, 9.07894e-11, 2.21580e-07, 7.31878e-05, 0.000909623, 0.00179548, &
 0.000909623,7.31878e-05,2.21580e-07,9.14337e-06, 0.000909623, 0.00550289, 0.00854438, 0.00550289, 0.000909623, 9.14337e-06,   &
 2.69243e-05, 0.00179548, 0.00854438,0.0122469, 0.00854438,      0.00179548, 2.69243e-05, 9.14337e-06, 0.000909623, 0.00550289,&
  0.00854438, 0.00550289,0.000909623,9.14337e-06, 2.21580e-07,  7.31878e-05, 0.000909623, 0.00179548, 0.000909623, 7.31878e-05,&
 2.21580e-07,9.07894e-11,2.21580e-07,9.14337e-06, 2.69243e-05, 9.14337e-06, 2.21580e-07, 9.07894e-11, 1.24384e-08, 9.14337e-06,&
 0.000183649,0.000425400,0.000183649,9.14337e-06, 1.24384e-08, 9.14337e-06, 0.000909623, 0.00550289, 0.00854438, 0.00550289,   &
 0.000909623,9.14337e-06,0.000183649,0.00550289, 0.0162043,    0.0197919, 0.0162043, 0.00550289, 0.000183649, 0.000425400,     &
  0.00854438,  0.0197919,  0.0223153,0.0197919, 0.00854438,        0.000425400, 0.000183649, 0.00550289, 0.0162043, 0.0197919, &
   0.0162043, 0.00550289,0.000183649,9.14337e-06, 0.000909623,  0.00550289, 0.00854438, 0.00550289, 0.000909623, 9.14337e-06,  &
 1.24384e-08,9.14337e-06,0.000183649,0.000425400, 0.000183649, 9.14337e-06, 1.24384e-08, 5.46411e-08, 2.69243e-05, 0.000425400,&
 0.000909623,0.000425400,2.69243e-05,5.46411e-08, 2.69243e-05, 0.00179548, 0.00854438, 0.0122469, 0.00854438, 0.00179548,      &
 2.69243e-05,0.000425400, 0.00854438,0.0197919, 0.0223153,      0.0197919, 0.00854438, 0.000425400, 0.000909623, 0.0122469,    &
   0.0223153,  0.0232260,  0.0223153,0.0122469, 0.000909623,       0.000425400, 0.00854438, 0.0197919, 0.0223153, 0.0197919,   &
  0.00854438,0.000425400,2.69243e-05,0.00179548, 0.00854438,   0.0122469, 0.00854438, 0.00179548, 2.69243e-05, 5.46411e-08,    &
 2.69243e-05,0.000425400,0.000909623,0.000425400, 2.69243e-05, 5.46411e-08, 1.24384e-08, 9.14337e-06, 0.000183649, 0.000425400,&
 0.000183649,9.14337e-06,1.24384e-08,9.14337e-06, 0.000909623, 0.00550289, 0.00854438, 0.00550289, 0.000909623, 9.14337e-06,   &
 0.000183649, 0.00550289,  0.0162043,0.0197919, 0.0162043,        0.00550289, 0.000183649, 0.000425400, 0.00854438, 0.0197919, &
   0.0223153,  0.0197919, 0.00854438,0.000425400, 0.000183649,    0.00550289, 0.0162043, 0.0197919, 0.0162043, 0.00550289,     &
 0.000183649,9.14337e-06,0.000909623,0.00550289, 0.00854438,   0.00550289, 0.000909623, 9.14337e-06, 1.24384e-08, 9.14337e-06, &
 0.000183649,0.000425400,0.000183649,9.14337e-06, 1.24384e-08, 9.07894e-11, 2.21580e-07, 9.14337e-06, 2.69243e-05, 9.14337e-06,&
 2.21580e-07,9.07894e-11,2.21580e-07,7.31878e-05, 0.000909623, 0.00179548, 0.000909623, 7.31878e-05, 2.21580e-07, 9.14337e-06, &
 0.000909623, 0.00550289, 0.00854438,0.00550289, 0.000909623,    9.14337e-06, 2.69243e-05, 0.00179548, 0.00854438, 0.0122469,  &
  0.00854438, 0.00179548,2.69243e-05,9.14337e-06, 0.000909623,  0.00550289, 0.00854438, 0.00550289, 0.000909623, 9.14337e-06,  &
 2.21580e-07,7.31878e-05,0.000909623,0.00179548, 0.000909623,  7.31878e-05, 2.21580e-07, 9.07894e-11, 2.21580e-07, 9.14337e-06,&
 2.69243e-05,9.14337e-06,2.21580e-07,9.07894e-11, 5.03438e-15, 9.07894e-11, 1.24384e-08, 5.46411e-08, 1.24384e-08, 9.07894e-11,&
 5.03438e-15,9.07894e-11,2.21580e-07,9.14337e-06, 2.69243e-05, 9.14337e-06, 2.21580e-07, 9.07894e-11, 1.24384e-08, 9.14337e-06,&
 0.000183649,0.000425400,0.000183649,9.14337e-06, 1.24384e-08, 5.46411e-08, 2.69243e-05, 0.000425400, 0.000909623, 0.000425400,&
 2.69243e-05,5.46411e-08,1.24384e-08,9.14337e-06, 0.000183649, 0.000425400, 0.000183649, 9.14337e-06, 1.24384e-08, 9.07894e-11,&
 2.21580e-07,9.14337e-06,2.69243e-05,9.14337e-06, 2.21580e-07, 9.07894e-11, 5.03438e-15, 9.07894e-11, 1.24384e-08, 5.46411e-08,&
 1.24384e-08,9.07894e-11,5.03438e-15 /), (/ 7, 7, 7 /))
!
  contains
!
!***********************************************************************
    subroutine max_mn(a,res)
!
!  successively calculate maximum of a, which is supplied at each call.
!  Start from scratch if lfirstpoint=.true.
!
!   1-apr-01/axel+wolf: coded
!
      real, dimension (nx), intent(in) :: a
      real, intent(inout) :: res
!
      if (lfirstpoint) then
        res=maxval(a)
      else
        res=max(res,maxval(a))
      endif
!
    endsubroutine max_mn
!***********************************************************************
    subroutine mean_mn(a,res)
!
!  Successively calculate mean of a, which is supplied at each call.
!  Start from zero if lfirstpoint=.true.
!
!  17-dec-01/wolf: coded
!  20-jun-07/dhruba:adapted for spherical polar coordinate system
!
      real, dimension (nx) :: a
      real :: res
      integer :: isum
!
      if (lfirstpoint) then
        if (lspherical_coords) then
          res = 0.
          do isum=l1,l2
            res = res+x(isum)*x(isum)*sinth(m)*a(isum)*a(isum)
          enddo
        else
          res=sum(dble(a))     ! sum at double precision to improve accuracy
        endif
      else
        if (lspherical_coords) then
          do isum=l1,l2
            res = res+x(isum)*x(isum)*sinth(m)*a(isum)*a(isum)
          enddo
        else
          res=res+sum(dble(a))
        endif
      endif
!
      if (lcylindrical_coords) &
          call fatal_error('mean_mn','not implemented for cylindrical')
!
    endsubroutine mean_mn
!***********************************************************************
    subroutine meanyz_s(f,iif,mean,lexp)
!
!  Calculates mean of variable at subscript iif of f over all x coordinate surfaces 
!  by proper integration.
!
!  7-jun-15/MR: coded
!
      use General, only: loptest

      real, dimension (mx,my,mz,mfarray) :: f
      integer :: iif
      real, dimension(nx) :: mean
      logical, optional :: lexp

      integer :: nn, ll, lll

      do ll=l1,l2
        lll=ll-l1+1
        mean(lll)=0.
        do nn=n1,n2
          if (loptest(lexp)) then
            mean(lll) = mean(lll) + dVol_z(nn)*sum(dVol_y(m1:m2)*exp(f(ll,m1:m2,nn,iif)))
          else
            mean(lll) = mean(lll) + dVol_z(nn)*sum(dVol_y(m1:m2)*f(ll,m1:m2,nn,iif))
          endif
        enddo
        mean(lll) = mean(lll)/Area_yz
      enddo

      call finalize_aver(nprocyz,23,mean)

    endsubroutine meanyz_s
!***********************************************************************
    subroutine meanyz_v(f,iif1,iif2,mean,lexp)
!
!  Calculates mean of variables at subscripts iif1 through iif2 of f over all x coordinate surfaces 
!  by proper integration.
!
!  7-jun-15/MR: coded
!
      use General, only: loptest

      real, dimension (mx,my,mz,mfarray) :: f
      integer :: iif1,iif2
      real, dimension(nx,iif2-iif1+1) :: mean
      logical, optional :: lexp

      integer :: nn, ll, lll, iif, iil

      do iif=iif1,iif2
        iil=iif-iif1+1
        do ll=l1,l2
          lll=ll-l1+1
          mean(lll,iil)=0.
          do nn=n1,n2
            if (loptest(lexp)) then
              mean(lll,iil) = mean(lll,iil) + dVol_z(nn)*sum(dVol_y(m1:m2)*exp(f(ll,m1:m2,nn,iif)))
            else
              mean(lll,iil) = mean(lll,iil) + dVol_z(nn)*sum(dVol_y(m1:m2)*f(ll,m1:m2,nn,iif))
            endif
          enddo
          mean(lll,iil) = mean(lll,iil)/Area_yz
        enddo
      enddo

      call finalize_aver(nprocyz,23,mean)

    endsubroutine meanyz_v
!***********************************************************************
    subroutine rms_mn(a,res)
!
!  Successively calculate rms of a, which is supplied at each call.
!  Start from zero if lfirstpoint=.true.
!
!   1-apr-01/axel+wolf: coded
!
      real, dimension (nx) :: a
      real :: res
      integer :: isum
!
      if (lfirstpoint) then
        if (lspherical_coords) then
          res = 0.
          do isum=l1,l2
            res = res+x(isum)*x(isum)*sinth(m)*a(isum)*a(isum)
          enddo
        else
          res=sum(a**2)
        endif
      else
        if (lspherical_coords) then
          do isum=l1,l2
            res = res+x(isum)*x(isum)*sinth(m)*a(isum)*a(isum)
          enddo
        else
          res=res+sum(a**2)
        endif
      endif
!
      if (lcylindrical_coords) &
          call fatal_error('rms_mn','not implemented for cylindrical')
!
    endsubroutine rms_mn
!***********************************************************************
    subroutine rms2_mn(a2,res)
!
!  Successively calculate rms of a, with a2=a^2 being supplied at each
!  call.
!  Start from zero if lfirstpoint=.true.
!
!   1-apr-01/axel+wolf: coded
!
      real, dimension (nx) :: a2
      real :: res
      integer :: isum
!
      if (lfirstpoint) then
        if (lspherical_coords) then
          res = 0.
          do isum=l1,l2
            res = res+x(isum)*x(isum)*sinth(m)*a2(isum)
          enddo
        else
          res=sum(a2)
        endif
      else
        if (lspherical_coords) then
          do isum=l1,l2
            res = res+x(isum)*x(isum)*sinth(m)*a2(isum)
          enddo
        else
          res=res+sum(a2)
        endif
      endif
!
      if (lcylindrical_coords) &
          call fatal_error('rms2_mn','not implemented for cylindrical')
!
    endsubroutine rms2_mn
!***********************************************************************
    subroutine sum_mn(a,res)
!
!  Successively calculate the sum over all points of a, which is supplied
!  at each call.
!  Start from zero if lfirstpoint=.true.
!
!   1-apr-01/axel+wolf: coded
!
      real, dimension (nx) :: a
      real :: res
      integer :: isum
!
      if (lfirstpoint) then
        if (lspherical_coords) then
          res = 0.
          do isum=l1,l2
            res = res+x(isum)*x(isum)*sinth(m)*a(isum)
          enddo
        else
          res=sum(a)
        endif
      else
        if (lspherical_coords) then
          do isum=l1,l2
            res = res+x(isum)*x(isum)*sinth(m)*a(isum)
          enddo
        else
          res=res+sum(a)
        endif
      endif
!
      if (lcylindrical_coords) &
          call fatal_error('sum_mn','not implemented for cylindrical')
!
    endsubroutine sum_mn
!***********************************************************************
    subroutine dot_mn(a,b,c,ladd)
!
!  Dot product, c=a.b, on pencil arrays
!
!   3-apr-01/axel+gitta: coded
!  24-jun-08/MR: ladd added for incremental work
!
      use General, only: loptest
!
      real, dimension (:,:) :: a,b
      real, dimension (:) :: c
!
      logical, optional :: ladd
!
      intent(in) :: a,b,ladd
      intent(inout) :: c
!
      integer :: i
      logical :: l0
!
      l0 = .not.loptest(ladd)

      do i=1,size(a,2)
        if (l0) then
          c=a(:,i)*b(:,i)
          l0=.false.
        else
          c=c+a(:,i)*b(:,i)
        endif
      enddo
!
    endsubroutine dot_mn
!***********************************************************************
    subroutine vec_dot_3tensor(a,b,c)
!
!  Dot product of a vector with 3 tensor,
!   c_ij = a_k b_ijk
!
!  28-aug-08/dhruba : coded
!
      real, dimension (nx,3) :: a
      real, dimension (nx,3,3) :: c
      real, dimension (nx,3,3,3) :: b
      integer :: i,j,k
!
      intent(in) :: a,b
      intent(out) :: c
!
      do i=1,3
        do j=1,3
          c(:,i,j) = 0.
          do k=1,3
            c(:,i,j) = c(:,i,j) + a(:,k)*b(:,i,j,k)
          enddo
        enddo
      enddo
!
    endsubroutine vec_dot_3tensor
!***********************************************************************
    subroutine contract_jk3(a,c)
!
!  Contracts the jk of a_ijk.
!
!  20-aug-08/dhruba: coded
!
      real, dimension (nx,3,3,3) :: a
      real, dimension (nx,3) :: c
      integer :: i,j,k
!
      intent(in) :: a
      intent(out) :: c
!
      c=0
      do i=1,3; do j=1,3; do k=1,3
        c(:,i)=c(:,i)+a(:,i,j,k)
      enddo; enddo; enddo
!
    endsubroutine contract_jk3
!***********************************************************************
    subroutine matrix2linarray(mm,aa)
!
! converts a 3X3 matrix to an array of length 9
!
!18-sep-15/dhruba: coded
!
      real,dimension(3,3),intent(in) :: mm
      real,dimension(9),intent(out) :: aa
!
      integer :: i,j,ij
      ij=0
      do i=1,3
        do j=1,3
          ij=ij+1
          aa(ij) = mm(i,j)
        enddo
      enddo
!
    endsubroutine matrix2linarray
!***********************************************************************
    subroutine linarray2matrix(aa,mm)
!
! converts a 3X3 matrix to an array of length 9
!
!18-sep-15/dhruba: coded
!
      real,dimension(9),intent(in) :: aa
      real,dimension(3,3),intent(out) :: mm
!
      integer :: i,j,ij
      ij=0
      do i=1,3
        do j=1,3
          ij=ij+1
          mm(i,j) = aa(ij)
       enddo
      enddo
!
    endsubroutine linarray2matrix
!***********************************************************************
    subroutine Inv2_3X3mat(A,QQ)
!
! calcualates the second invariant of a 3X3 matrix directly
!
      real,dimension(3,3), intent(in) :: A
      real, intent(out) :: QQ
!
      QQ =   A(1,1)*A(2,2)  &
            + A(2,2)*A(3,3)  &
            + A(1,1)*A(3,2)  &
            - A(1,2)*A(2,1)  &
            - A(2,3)*A(3,2)  &
            - A(1,3)*A(3,1)
      
    endsubroutine Inv2_3X3mat
!***********************************************************************
    subroutine ScalarTripleProduct(A,B,C,product)
      real, dimension(3), intent(in) :: A,B,C
      real :: product
      real,dimension(3,3) :: Mat
      Mat(:,1) = A
      Mat(:,2) = B
      Mat(:,3) = C
      call det3X3mat(Mat,product)
!        
    endsubroutine ScalarTripleProduct
!***********************************************************************
    subroutine det3X3mat(A,det)
!
! calculates determinant of a 3X3 matrix directly
!
      real,dimension(3,3), intent(in) :: A
      real, intent(out) :: det

      det =   A(1,1)*A(2,2)*A(3,3)  &
            - A(1,1)*A(2,3)*A(3,2)  &
            - A(1,2)*A(2,1)*A(3,3)  &
            + A(1,2)*A(2,3)*A(3,1)  &
            + A(1,3)*A(2,1)*A(3,2)  &
            - A(1,3)*A(2,2)*A(3,1)

    endsubroutine det3X3mat
!***********************************************************************
    subroutine dot_mn_sv(a,b,c)
!
!  Dot product, c=a.b, between non-pencilized vector and pencil array.
!
!  10-oct-06/axel: coded
!
      real, dimension (3)    :: a
      real, dimension (nx,3) :: b
      real, dimension (nx)   :: c
!
      intent(in) :: a,b
      intent(out) :: c
!
      c=a(1)*b(:,1)+a(2)*b(:,2)+a(3)*b(:,3)
!
    endsubroutine dot_mn_sv
!***********************************************************************
    subroutine dot_mn_sm(a,b,c)
!
!  Dot product, c=a.b, between non-pencilized vector and pencil matrix.
!
!  10-oct-06/axel: coded
!
      real, dimension (3)      :: a
      real, dimension (nx,3,3) :: b
      real, dimension (nx,3)   :: c
      integer :: i
!
      intent(in) :: a,b
      intent(out) :: c
!
      do i=1,3
        c(:,i)=a(1)*b(:,i,1)+a(2)*b(:,i,2)+a(3)*b(:,i,3)
      enddo
!
    endsubroutine dot_mn_sm
!***********************************************************************
    subroutine dot_mn_vm(a,b,c)
!
!  Dot product, c=a.b, between pencil vector and pencil matrix.
!
!  10-oct-06/axel: coded
!
      real, dimension (nx,3)   :: a
      real, dimension (nx,3,3) :: b
      real, dimension (nx,3)   :: c
      integer :: i
!
      intent(in) :: a,b
      intent(out) :: c
!
      do i=1,3
        c(:,i)=a(:,1)*b(:,i,1)+a(:,2)*b(:,i,2)+a(:,3)*b(:,i,3)
      enddo
!
    endsubroutine dot_mn_vm
!***********************************************************************
    subroutine dot_mn_vm_trans(a,b,c)
!
!  Dot product, c=a.b, between pencil vector and pencil matrix.
!  I think the name of dot_mn_vm is not right and should have been transposed.
!
!  10-oct-06/axel: coded
!
      real, dimension (nx,3)      :: a
      real, dimension (nx,3,3) :: b
      real, dimension (nx,3)   :: c
      integer :: i
!
      intent(in) :: a,b
      intent(out) :: c
!
      do i=1,3
        c(:,i)=a(:,1)*b(:,1,i)+a(:,2)*b(:,2,i)+a(:,3)*b(:,3,i)
      enddo
!
    endsubroutine dot_mn_vm_trans
!***********************************************************************
    subroutine dot_0(a,b,c)
!
!  Dot product, c=a.b, of two simple 3-d arrays.
!
!  11-mar-04/wolf: coded
!
      real, dimension (:) :: a,b
      real :: c
!
      intent(in) :: a,b
      intent(out) :: c
!
      c = dot_product(a,b)
!
    endsubroutine dot_0
!***********************************************************************
    subroutine dot2_mn(a,b,fast_sqrt,precise_sqrt)
!
!  Dot product with itself, to calculate max and rms values of a vector.
!  FAST_SQRT is only correct for ~1e-18 < |a| < 1e18 (for single precision);
!  PRECISE_SQRT works for full range.
!
!  29-sep-97/axel: coded
!   1-apr-01/axel: adapted for cache-efficient sub-array formulation
!  25-jun-05/bing: added optional args for calculating |a|
!
      real, dimension (nx,3) :: a
      real, dimension (nx) :: b,a_max
      logical, optional :: fast_sqrt,precise_sqrt
      logical :: fast_sqrt1,precise_sqrt1
!
      intent(in) :: a,fast_sqrt,precise_sqrt
      intent(out) :: b
!
!     ifc treats these variables as SAVE so we need to reset
      if (present(fast_sqrt)) then
        fast_sqrt1=fast_sqrt
      else
        fast_sqrt1=.false.
      endif
!
      if (present(precise_sqrt)) then
        precise_sqrt1=precise_sqrt
      else
        precise_sqrt1=.false.
      endif
!
!  Rescale by factor a_max before taking sqrt.
!  In single precision this increases the dynamic range from 1e18 to 1e36.
!  To avoid division by zero when calculating a_max, we add tini.
!
      if (precise_sqrt1) then
        a_max=tini+maxval(abs(a),dim=2)
        b=(a(:,1)/a_max)**2+(a(:,2)/a_max)**2+(a(:,3)/a_max)**2
        b=a_max*sqrt(b)
      else
        b=a(:,1)**2+a(:,2)**2+a(:,3)**2
        if (fast_sqrt1) b=sqrt(b)
      endif
!
    endsubroutine dot2_mn
!***********************************************************************
    subroutine dot2_0(a,b)
!
!  Dot product, c=a.b, of two simple 3-d arrays.
!
!  11-mar-04/wolf: coded
!
      real, dimension (:) :: a
      real :: b
!
      intent(in) :: a
      intent(out) :: b
!
      b = dot_product(a,a)
!
    endsubroutine dot2_0
!***********************************************************************
    subroutine dot_mn_add(a,b,c)
!
!  Dot product, add to previous value.
!
!  11-nov-02/axel: adapted from dot_mn
!
      real, dimension (nx,3) :: a,b
      real, dimension (nx) :: c
!
      intent(in) :: a,b
      intent(inout) :: c
!
      c=c+a(:,1)*b(:,1)+a(:,2)*b(:,2)+a(:,3)*b(:,3)
!
    endsubroutine dot_mn_add
!***********************************************************************
    subroutine dot_mn_sub(a,b,c)
!
!  Dot product, subtract from previous value.
!
!  21-jul-03/axel: adapted from dot_mn_sub
!
      real, dimension (nx,3) :: a,b
      real, dimension (nx) :: c
!
      intent(in) :: a,b
      intent(inout) :: c
!
      c=c-(a(:,1)*b(:,1)+a(:,2)*b(:,2)+a(:,3)*b(:,3))
!
    endsubroutine dot_mn_sub
!***********************************************************************
    subroutine dot2fj(a,vec,b)
!
!  Dot product with itself, multiplied by anisotropic factor.
!
!  20-dec-12/wlad: adapted from dot2_mn
!
      real, dimension (nx,3) :: a
      real, dimension (3)    :: vec
      real, dimension (nx)   :: b
      integer :: j
!
      intent(in) :: a,vec
      intent(out) :: b
!
      b=0.
      do j=1,3
         b=b+vec(j)*a(:,j)**2
      enddo
!
    endsubroutine dot2fj
!***********************************************************************
    subroutine dyadic2(a,b)
!
!  Dyadic product of a penciled vector with itself.
!
!  24-jan-09/axel: coded
!
      real, dimension (nx,3) :: a
      real, dimension (nx,3,3) :: b
!
      intent(in) :: a
      intent(out) :: b
!
!  diagonal components
!
      b(:,1,1)=a(:,1)**2
      b(:,2,2)=a(:,2)**2
      b(:,3,3)=a(:,3)**2
!
!  upper off-diagonal components
!
      b(:,1,2)=a(:,1)*a(:,2)
      b(:,1,3)=a(:,1)*a(:,3)
      b(:,2,3)=a(:,2)*a(:,3)
!
!  lower off-diagonal components
!
      b(:,2,1)=b(:,1,2)
      b(:,3,1)=b(:,1,3)
      b(:,3,2)=b(:,2,3)
!
    endsubroutine dyadic2
!***********************************************************************
    function dyadic2_other(a) result(b)
!
!  Dyadic product of a vector with itself for non-pencil data.
!
!  26-may-19/MR: aped from dyadic2
!
      real, dimension (3) :: a
      real, dimension (3,3) :: b
!
      intent(in) :: a
!
!  diagonal components
!
      b(1,1)=a(1)**2
      b(2,2)=a(2)**2
      b(3,3)=a(3)**2
!
!  upper off-diagonal components
!
      b(1,2)=a(1)*a(2)
      b(1,3)=a(1)*a(3)
      b(2,3)=a(2)*a(3)
!
!  lower off-diagonal components
!
      b(2,1)=b(1,2)
      b(3,1)=b(1,3)
      b(3,2)=b(2,3)
!
    endfunction dyadic2_other
!***********************************************************************
    subroutine trace_mn(a,b)
!
!  Trace of a matrix.
!
!   3-apr-01/axel+gitta: coded
!
      real, dimension (nx,3,3) :: a
      real, dimension (nx) :: b
!
      intent(in) :: a
      intent(out) :: b
!
      b=a(:,1,1)+a(:,2,2)+a(:,3,3)
!
    endsubroutine trace_mn
!***********************************************************************
    subroutine multvv_mat_mn(a,b,c)
!
!  Vector multiplied with vector, gives matrix.
!
!   21-dec-01/nils: coded
!   16-jul-02/nils: adapted from pencil_mpi
!
      real, dimension (nx,3) :: a,b
      real, dimension (nx,3,3) :: c
      integer :: i,j
!
      do i=1,3; do j=1,3
        c(:,i,j)=a(:,j)*b(:,i)
      enddo; enddo
!
    endsubroutine multvv_mat_mn
!***********************************************************************
    subroutine multmm_sc_mn(a,b,c)
!
!  Matrix multiplied with matrix, gives scalar.
!
!   21-dec-01/nils: coded
!   16-jul-02/nils: adapted from pencil_mpi
!
      real, dimension (nx,3,3) :: a,b
      real, dimension (nx) :: c
      integer :: i,j
!
      c=0.0
      do i=1,3; do j=1,3
        c=c+a(:,i,j)*b(:,i,j)
      enddo; enddo
!
    endsubroutine multmm_sc_mn
!***********************************************************************
    subroutine mult_matrix(a,b,c)
!
!  Matrix multiplication of two pencil variables.
!
      real, dimension (nx,3,3) :: a,b
      real, dimension (nx,3,3) :: c
      integer :: ix
!
!  24-feb-11/dhruba: using the f90 command matmul
!
      do ix=1,nx
        c(ix,:,:)=matmul(a(ix,:,:),b(ix,:,:))
      enddo
!
    endsubroutine mult_matrix
!***********************************************************************
    subroutine multm2_mn(a,b)
!
!  Matrix squared, gives scalar.
!
!  11-nov-02/axel: adapted from multmm_sc_mn
!
      real, dimension (nx,3,3) :: a
      real, dimension (nx) :: b
      integer :: i,j
!
      b=0.0
      do i=1,3; do j=1,3
        b=b+a(:,i,j)**2
      enddo; enddo
!
    endsubroutine multm2_mn
!***********************************************************************
    subroutine multm2_sym_mn(a,b)
!
!  Symmetric matrix squared, gives scalar.
!
!  24-aug-2011/Bourdin.KIS: adapted from multm2_mn
!
      real, dimension (nx,3,3), intent(in) :: a
      real, dimension (nx), intent(out) :: b
!
      integer :: i, j
!
      b = a(:,1,1)**2
      do i = 2, 3
        b = b + a(:,i,i)**2
        do j = 1, i-1
          b = b + 2 * a(:,i,j)**2
        enddo
      enddo
!
    endsubroutine multm2_sym_mn
!***********************************************************************
    subroutine multmv_mn(a,b,c,ladd)
!
!  Matrix multiplied with vector, gives vector.
!
!  C_i = A_{i,j} B_j
!
!   3-apr-01/axel+gitta: coded
!  24-jun-08/MR: ladd added for incremental work
!
      use General, only: loptest
!
      real, dimension (nx,3,3) :: a
      real, dimension (nx,3) :: b,c
      real, dimension (nx) :: tmp
      integer :: i,j
      logical, optional :: ladd
!
      intent(in) :: a,b,ladd
      intent(out) :: c
!
      do i=1,3
!
        j=1
        tmp=a(:,i,j)*b(:,j)
        do j=2,3
          tmp=tmp+a(:,i,j)*b(:,j)
        enddo
!
        if (loptest(ladd)) then
          c(:,i)=c(:,i)+tmp
        else
          c(:,i)=tmp
        endif
!
      enddo
!
    endsubroutine multmv_mn
!***********************************************************************
    subroutine multmv_mn_transp(a,b,c,ladd)
!
!  Transposed matrix multiplied with vector, gives vector.
!  Could have called multvm_mn, but this may not be clear enough.
!
!  C_i = A_{j,i} B_j
!
!  21-jul-03/axel: adapted from multmv_mn
!  24-jun-08/MR: ladd added for incremental work
!
      use General, only: loptest
!
      real, dimension (nx,3,3) :: a
      real, dimension (nx,3) :: b,c
      real, dimension (nx) :: tmp
      integer :: i,j
      logical, optional :: ladd
!
      intent(in) :: a,b,ladd
      intent(inout) :: c
!
      do i=1,3
        j=1
        tmp=a(:,j,i)*b(:,j)
        do j=2,3
          tmp=tmp+a(:,j,i)*b(:,j)
        enddo
!
        if (loptest(ladd)) then
          c(:,i)=c(:,i)+tmp
        else
          c(:,i)=tmp
        endif
!
      enddo
!
    endsubroutine multmv_mn_transp
!***********************************************************************
    subroutine multsv_mn(a,b,c)
!
!  Vector multiplied with scalar, gives vector.
!
!  22-nov-01/nils erland: coded
!  10-oct-03/axel: a is now the scalar (now consistent with old routines)
!  24-jun-08/MR: ladd added for incremental work
!  28-feb-10/bing: removed ladd keyword, use multsv_mn_add instead
!
      intent(in) :: a,b
      intent(out) :: c
!
      real, dimension (nx,3) :: b,c
      real, dimension (nx) :: a
      integer :: i
!
      do i=1,3
        c(:,i)=a*b(:,i)
      enddo
!
    endsubroutine multsv_mn
!***********************************************************************
    subroutine multsv_mn_add(a,b,c)
!
!  Vector multiplied with scalar, gives vector.
!
!  22-nov-01/nils erland: coded
!  10-oct-03/axel: a is now the scalar (now consistent with old routines)
!  24-jun-08/MR: ladd added for incremental work
!
      intent(in) :: a,b
      intent(inout) :: c
!
      real, dimension (nx,3) :: b,c
      real, dimension (nx) :: a
      integer :: i
!
      do i=1,3
        c(:,i)=c(:,i)+a*b(:,i)
      enddo
!
    endsubroutine multsv_mn_add
!***********************************************************************
    subroutine multsv_add_mn(a,b,c,d)
!
!  Multiply scalar with a vector and add to another vector.
!
!  29-oct-97/axel: coded
!
      real, dimension (nx,3) :: a,c,d
      real, dimension (nx) :: b
      integer :: j
!
      intent(in) :: a,b,c
      intent(out) :: d
!
      do j=1,3
        d(:,j)=a(:,j)+b*c(:,j)
      enddo
!
    endsubroutine multsv_add_mn
!***********************************************************************
    subroutine multvs_mn(a,b,c)
!
!  Vector pencil multiplied with scalar pencil, gives vector pencil.
!
!  22-nov-01/nils erland: coded
!
      real, dimension (nx,3) :: a, c
      real, dimension (nx) :: b
      integer :: i
!
      do i=1,3
        c(:,i)=a(:,i)*b(:)
      enddo
!
    endsubroutine multvs_mn
!***********************************************************************
    subroutine cross_mn(a,b,c)
!
!  Cross product, c = a x b, for pencil variables.
!  Previously called crossp.
!
      real, dimension (nx,3) :: a,b,c
!
      intent(in) :: a,b
      intent(out) :: c
!
      c(:,1)=a(:,2)*b(:,3)-a(:,3)*b(:,2)
      c(:,2)=a(:,3)*b(:,1)-a(:,1)*b(:,3)
      c(:,3)=a(:,1)*b(:,2)-a(:,2)*b(:,1)
!
    endsubroutine cross_mn
!***********************************************************************
    subroutine cross_mixed(a,b,c)
!
!  Cross product, c = a x b, of a pencil and a non-pencil variable.
!
!  17-apr-2015/MR: coded
!
      real, dimension (nx,3) :: a,c
      real, dimension (3) :: b
!
      intent(in) :: a,b
      intent(out) :: c
!
      c(:,1)=a(:,2)*b(3)-a(:,3)*b(2)
      c(:,2)=a(:,3)*b(1)-a(:,1)*b(3)
      c(:,3)=a(:,1)*b(2)-a(:,2)*b(1)
!
    endsubroutine cross_mixed
!***********************************************************************
    subroutine cross_0(a,b,c)
!
!  Cross product, c = a x b, for simple 3-d vectors (independent of position).
!
      real, dimension (3) :: a,b,c
!
      intent(in) :: a,b
      intent(out) :: c
!
      c(1)=a(2)*b(3)-a(3)*b(2)
      c(2)=a(3)*b(1)-a(1)*b(3)
      c(3)=a(1)*b(2)-a(2)*b(1)
!
    endsubroutine cross_0
!***********************************************************************
    subroutine gij(f,k,g,nder)
!
!  Calculate gradient of a vector, return matrix.
!
!   3-apr-01/axel+gitta: coded
!
      use Deriv, only: der,der2,der3,der4,der5,der6
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx,3,3) :: g
      real, dimension (nx) :: tmp
      integer :: i,j,k,k1,nder
!
      intent(in) :: f,k
      intent(out) :: g
!
      k1=k-1
      do i=1,3; do j=1,3
        if (nder == 1) then
          call der(f,k1+i,tmp,j)
        elseif (nder == 2) then
          call der2(f,k1+i,tmp,j)
        elseif (nder == 3) then
          call der3(f,k1+i,tmp,j)
        elseif (nder == 4) then
          call der4(f,k1+i,tmp,j)
        elseif (nder == 5) then
          call der5(f,k1+i,tmp,j)
        elseif (nder == 6) then
          call der6(f,k1+i,tmp,j)
        endif
        g(:,i,j)=tmp
      enddo; enddo
!
    endsubroutine gij
!***********************************************************************
    subroutine gijl_symmetric(f,k,gijl)
!
!  Calculate gradient of a (symmetric) second rank matrix, return 3rd rank
!  matrix
!
!  18-aug-08/dhruba: coded
!
      use Deriv, only: der
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx,3,3,3) :: gijl
      real, dimension (nx,3,3,3) :: tmpg
      real, dimension (nx) :: tmp
      integer :: i,j,l,k,k1
!
      intent(in) :: f,k
      intent(out) :: gijl
!
      k1=k-1
      do i=1,3
        do j=i,3
          k1=k1+1
          do l=1,3
            call der(f,k1,tmp,l)
            tmpg(:,i,j,l) = tmp
          enddo
        enddo
      enddo
      do l=1,3
        call symmetrise3x3_ut2lt(tmpg(:,:,:,l))
      enddo
      gijl=tmpg
!
    endsubroutine gijl_symmetric
!***********************************************************************
    subroutine grad_main(f,k,g)
!
!  Calculate gradient of a scalar, get vector.
!
!  29-sep-97/axel: coded
!
      use Deriv, only: der
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx,3) :: g
      integer :: k
!
      intent(in) :: f,k
      intent(out) :: g
!
      call der(f,k,g(:,1),1)
      call der(f,k,g(:,2),2)
      call der(f,k,g(:,3),3)
!
    endsubroutine grad_main
!***********************************************************************
    subroutine grad_other(f,g)
!
!  For non 'mvar' variable calculate gradient of a scalar, get vector
!
!  26-nov-02/tony: coded
!
      use Deriv, only: der
!
      real, dimension (mx,my,mz) :: f
      real, dimension (nx,3) :: g
      real, dimension (nx) :: tmp
!
      intent(in) :: f
      intent(out) :: g
!
!  Uses overloaded der routine.
!
      call der(f,tmp,1); g(:,1)=tmp
      call der(f,tmp,2); g(:,2)=tmp
      call der(f,tmp,3); g(:,3)=tmp
!
    endsubroutine grad_other
!***********************************************************************
    subroutine grad5(f,k,g)
!
!  Calculate 5th order gradient of a scalar, get vector.
!
!  03-jun-07/anders: adapted
!
      use Deriv, only: der5
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx,3) :: g
      real, dimension (nx) :: tmp
      integer :: k
!
      intent(in) :: f,k
      intent(out) :: g
!
      call der5(f,k,tmp,1); g(:,1)=tmp
      call der5(f,k,tmp,2); g(:,2)=tmp
      call der5(f,k,tmp,3); g(:,3)=tmp
!
    endsubroutine grad5
!***********************************************************************
    subroutine div(f,k,g,ldiff_fluxes)
!
!  Calculate divergence of vector, get scalar.
!
!  13-dec-01/nils: coded
!  16-jul-02/nils: adapted from pencil_mpi
!  31-aug-07/wlad: adapted for cylindrical and spherical coords
!  28-sep-15/Joern+MR: adapted to use for slope-limited diffusive 
!                      flux given on a staggered grid.
!                      Here the vector field v is stored in
!                      f(*,*,*,k:k+dimensionality-1), that is, 
!                      the components in the degenerate directions
!                      *are missing*, and the existing components are
!                      always stored in ascending order,
!                      e.g., [v_x,v_z] or [v_y,v_z].
!  16-nov-16/MR: modifications for non-Cartesian coordinates.
!
      use Deriv, only: der
      use General, only: loptest
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: k
      real, dimension (nx) :: g
      logical, optional :: ldiff_fluxes
!
      intent(in)  :: f,k,ldiff_fluxes
      intent(out) :: g
!
      integer :: k1,i
      real, dimension(nx) :: tmp
      integer, save :: indr=0, indth=0
      logical, save :: s0=.true.
!
      k1=k-1
!
      if (loptest(ldiff_fluxes)) then
        g=0.
        do i=1,dimensionality
          call der_4th_stag(f,k1+i,tmp,dim_mask(i))
!          call der_2nd_stag(f,k1+i,tmp,dim_mask(i))
          g=g+tmp
        enddo

        if (s0) then
          s0=.false. 
          if (dim_mask(1)==1) indr=1
          if (dim_mask(2)==2) indth=2
          if (dim_mask(1)==2) indth=1
        endif
!
        if (lspherical_coords) then
          if (indr>0) g=g+r1_mn*2.*f(l1:l2,m,n,k1+indr)
          if (indth>0) g=g+r1_mn*cotth(m)*f(l1:l2,m,n,k1+indth)
        endif
!
        if (lcylindrical_coords.and.indr>0) g=g+rcyl_mn1*f(l1:l2,m,n,k1+indr)
!
      else
        call der(f,k1+1,tmp,1)
        g=tmp
        call der(f,k1+2,tmp,2)
        g=g+tmp
        call der(f,k1+3,tmp,3)
        g=g+tmp

        if (lspherical_coords) &
          g=g+r1_mn*(2.*f(l1:l2,m,n,k1+1)+cotth(m)*f(l1:l2,m,n,k1+2))
        if (lcylindrical_coords) g=g+rcyl_mn1*f(l1:l2,m,n,k1+1)
!
      endif
!
    endsubroutine div
!***********************************************************************
    subroutine der_2nd(f,k,df,j)

      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(nx), intent(out) :: df
      integer, intent(in) :: j, k
!
      if (j==1) then
        if (nxgrid/=1) then
          df=(f(l1+1:l2+1,m,n,k)-f(l1-1:l2-1,m,n,k))/(2.*dx) 
        else
          df=0.
          if (ip<=5) print*, 'der_2nd: Degenerate case in x-direction'
        endif
      elseif (j==2) then
        if (nygrid/=1) then
          df=(f(l1:l2,m+1,n,k)-f(l1:l2,m-1,n,k))/(2.*dy)
        else
          df=0.
          if (ip<=5) print*, 'der_2nd: Degenerate case in y-direction'
        endif
      elseif (j==3) then
        if (nzgrid/=1) then
          df=(f(l1:l2,m,n+1,k)-f(l1:l2,m,n-1,k))/(2.*dz)
        else
          df=0.
          if (ip<=5) print*, 'der_2nd: Degenerate case in z-direction'
        endif
      endif

    endsubroutine der_2nd
!***********************************************************************
    subroutine der_4th_stag(f,k,df,j)
!
!  Calculates 1st order derivative by a 4th order difference scheme from
!  data given on a grid shifted by half a grid step w.r.t. the point looked at.
!  Only valid for equidistant grids!
!
!  30-sep-15/MR: coded
!   4-feb-16/MR: checked again
!  16-nov-16/MR: modifications for non-Cartesian coordinates.
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(nx), intent(out) :: df
      integer, intent(in) :: j, k
!
      if (j==1) then
        if (nxgrid/=1) then
          df=( -    (f(l1+1:l2+1,m,n,k)-f(l1-2:l2-2,m,n,k))     &
               +27.*(f(l1  :l2  ,m,n,k)-f(l1-1:l2-1,m,n,k)) )/(24.*dx)
        else
          df=0.
          if (ip<=5) print*, 'der_4th_stag: Degenerate case in x-direction'
        endif
      elseif (j==2) then
        if (nygrid/=1) then
          df=( -    (f(l1:l2,m+1,n,k)-f(l1:l2,m-2,n,k))    &
               +27.*(f(l1:l2,m  ,n,k)-f(l1:l2,m-1,n,k)) )/(24.*dy) 
          if (lspherical_coords  ) df = df * r1_mn
          if (lcylindrical_coords) df = df * rcyl_mn1
        else
          df=0.
          if (ip<=5) print*, 'der_4th_stag: Degenerate case in y-direction'
        endif
      elseif (j==3) then
        if (nzgrid/=1) then
          df=( -    (f(l1:l2,m,n+1,k)-f(l1:l2,m,n-2,k))    &
               +27.*(f(l1:l2,m,n  ,k)-f(l1:l2,m,n-1,k)) )/(24.*dz)
          if (lspherical_coords) df = df * r1_mn * sin1th(m)
        else
          df=0.
          if (ip<=5) print*, 'der_4th_stag: Degenerate case in z-direction'
        endif
      endif

    endsubroutine der_4th_stag
!***********************************************************************
    subroutine der_2nd_stag(f,k,df,j)
!
!  Calculates 1st order derivative by a 2nd order difference scheme from
!  data given on a grid shifted by half a grid step w.r.t. the point looked at.
!  Only valid for equidistant grids!
!
!  23-jun-18/JW: Adapted from der_4ht_stag
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(nx), intent(out) :: df
      integer, intent(in) :: j, k
!
      if (j==1) then
        if (nxgrid/=1) then
          df=( f(l1:l2,m,n,k)-f(l1-1:l2-1,m,n,k))/(2.*dx)
        else
          df=0.
          if (ip<=5) print*, 'der_2nd_stag: Degenerate case in x-direction'
        endif
      elseif (j==2) then
        if (nygrid/=1) then
          df=(f(l1:l2,m  ,n,k)-f(l1:l2,m-1,n,k))/(2*dy)
          if (lspherical_coords  ) df = df * r1_mn
          if (lcylindrical_coords) df = df * rcyl_mn1
        else
          df=0.
          if (ip<=5) print*, 'der_2nd_stag: Degenerate case in y-direction'
        endif
      elseif (j==3) then
        if (nzgrid/=1) then
          df=(f(l1:l2,m,n  ,k)-f(l1:l2,m,n-1,k))/(2.*dz)
          if (lspherical_coords) df = df * r1_mn * sin1th(m)
        else
          df=0.
          if (ip<=5) print*, 'der_2nd_stag: Degenerate case in z-direction'
        endif
      endif

    endsubroutine der_2nd_stag
!***********************************************************************
    subroutine div_other(f,g)
!
      use Deriv, only: der
!
      real, dimension (mx,my,mz,3) :: f
      real, dimension (nx) :: g, tmp
!
      call der(f(:,:,:,1),tmp,1)
      g=tmp
      call der(f(:,:,:,2),tmp,2)
      g=g+tmp
      call der(f(:,:,:,3),tmp,3)
      g=g+tmp
!
      if (lspherical_coords) then
        g=g+2.*r1_mn*f(l1:l2,m,n,1)+r1_mn*cotth(m)*f(l1:l2,m,n,2)
      endif
!
      if (lcylindrical_coords) then
        g=g+rcyl_mn1*f(l1:l2,m,n,1)
      endif
!
    endsubroutine div_other
!***********************************************************************
    subroutine div_mn(aij,b,a)
!
!  Calculate divergence from derivative matrix.
!
!  18-sep-04/axel: coded
!  21-feb-07/axel: corrected spherical coordinates
!  14-mar-07/wlad: added cylindrical coordinates
!
      real, dimension (nx,3,3) :: aij
      real, dimension (nx,3) :: a
      real, dimension (nx) :: b
!
      intent(in) :: aij,a
      intent(out) :: b
!
      b=aij(:,1,1)+aij(:,2,2)+aij(:,3,3)
!
!  Adjustments for other coordinate systems.
!
      if (lspherical_coords) then
        b=b+2.*r1_mn*a(:,1)+r1_mn*cotth(m)*a(:,2)
      endif
!
      if (lcylindrical_coords) then
        b=b+rcyl_mn1*a(:,1)
      endif
!
      if (lpipe_coords) then
        b=b+glnCrossSec*a(:,1)
      endif
!
    endsubroutine div_mn
!***********************************************************************
    subroutine div_mn_2tensor(aijk,bi)
!
!  Calculate divergence from derivative matrix.
!
!  07-aug-10/dhruba: coded
!
      real, dimension (nx,3,3,3) :: aijk
      real, dimension (nx,3) :: bi
!
      intent(in) :: aijk
      intent(out) :: bi
      integer :: i
!
      do i=1,3
        bi(:,i)=aijk(:,i,1,1)+aijk(:,i,2,2)+aijk(:,i,3,3)
      enddo
!
!  Adjustments for spherical coordinate system.
!
      if (lspherical_coords) then
        call fatal_error('div_mn_2tensor','not impelmented in sph-coordinate')
      endif
!
      if (lcylindrical_coords) then
        call fatal_error('div_mn_2tensor','not impelmented in cyl-coordinate')
      endif
!
    endsubroutine div_mn_2tensor
!***********************************************************************
    subroutine curl_mn(aij,b,a,lcovariant_derivative)
!
!  Calculate curl from derivative matrix.
!
!  21-jul-03/axel: coded
!  21-feb-07/axel: corrected spherical coordinates
!  14-mar-07/wlad: added cylindrical coordinates
!  16-jun-16/fred: added option to use covariant derivative bij
!
      use General, only: loptest
!
      real, dimension (nx,3,3), intent (in) :: aij
      real, dimension (:,:), intent (in), optional :: a
      logical, intent (in), optional :: lcovariant_derivative
      real, dimension (nx,3), intent (out) :: b
      integer :: i1=1,i2=2,i3=3,i4=4,i5=5,i6=6,i7=7,a1,a2
!
      b(:,1)=aij(:,3,2)-aij(:,2,3)
      b(:,2)=aij(:,1,3)-aij(:,3,1)
      b(:,3)=aij(:,2,1)-aij(:,1,2)
!
      if (present(a)) then
        a1 = 1
        a2 = size(a,1)
        if (((a2 /= nx) .and. (a2 /= mx)) .or. (size(a,2) /= 3)) &
            call fatal_error('curl_mn','array "a" has wrong size, must be (nx,3) or (mx,3)')
        if (a2 == mx) then
          a1 = l1
          a2 = l2
        endif
      endif
!
!  Adjustments for spherical coordinate system.
!
      if (lspherical_coords) then
        if (.not.loptest(lcovariant_derivative)) then
          if (.not. present(a)) then
            call fatal_error('curl_mn','Need a for spherical curl')
          endif
          b(:,1)=b(:,1)+a(a1:a2,3)*r1_mn*cotth(m)
          b(:,2)=b(:,2)-a(a1:a2,3)*r1_mn
          b(:,3)=b(:,3)+a(a1:a2,2)*r1_mn
        endif
      endif
!
!  Adjustments for cylindrical coordinate system.
!  If we go all the way to the center, we need to put a regularity condition.
!  We do this here currently up to second order, and only for curl_mn.
!
      if (lcylindrical_coords.and.present(a)) then
        if (.not.loptest(lcovariant_derivative)) b(:,3)=b(:,3)+a(a1:a2,2)*rcyl_mn1
        if (rcyl_mn(1)==0.) b(i1,3)=(360.*b(i2,3)-450.*b(i3,3)+400.*b(i4,3) &
                                    -225.*b(i5,3)+72.*b(i6,3)-10.*b(i7,3))/147.
      endif
!
    endsubroutine curl_mn
!***********************************************************************
    subroutine curl_horizontal(f,k,g)
!
!  Calculate curl of a vector, whose z component is given.
!
!   8-oct-09/axel: adapted from
!
      use Deriv, only: der
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx,3) :: g
      real, dimension (nx) :: tmp1,tmp2
      integer :: k
!
      intent(in) :: f,k
      intent(out) :: g
!
      call der(f,k,tmp1,2)
      g(:,1)=tmp1
!
      call der(f,k,tmp2,1)
      g(:,2)=-tmp2
!
!g(:,1)=0.
!g(:,2)=0.
      g(:,3)=0.
!
!  Adjustments for spherical corrdinate system.
!
      if (lspherical_coords) then
        g(:,1)=g(:,1)+f(l1:l2,m,n,k)*r1_mn*cotth(m)
        g(:,2)=g(:,2)-f(l1:l2,m,n,k)*r1_mn
      endif
!
      if (lcylindrical_coords) then
        call fatal_error("curl_horizontal","not implemented for cylindrical coords")
!--     g(:,3)=g(:,3)+f(l1:l2,m,n,k1+2)*rcyl_mn1
      endif
!
    endsubroutine curl_horizontal
!***********************************************************************
    subroutine curl(f, k, g, ignoredx)
!
!  Calculate curl of a vector, get vector.
!
!  12-sep-97/axel: coded
!  10-sep-01/axel: adapted for cache efficiency
!  11-sep-04/axel: began adding spherical coordinates
!  21-feb-07/axel: corrected spherical coordinates
!  14-mar-07/wlad: added cylindrical coordinates
!  20-sep-13/ccyang: added optional argument ignoredx
!
      use Deriv, only: der
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(nx,3), intent(out) :: g
      integer, intent(in) :: k
      logical, intent(in), optional :: ignoredx
!
      real, dimension(nx) :: tmp1, tmp2
      logical :: igdx
      integer :: k1
!
      if (present(ignoredx)) then
        igdx = ignoredx
      else
        igdx = .false.
      endif
!
      k1=k-1
!
      call der(f, k1+3, tmp1, 2, ignoredx=igdx)
      call der(f, k1+2, tmp2, 3, ignoredx=igdx)
      g(:,1)=tmp1-tmp2
!
      call der(f, k1+1, tmp1, 3, ignoredx=igdx)
      call der(f, k1+3, tmp2, 1, ignoredx=igdx)
      g(:,2)=tmp1-tmp2
!
      call der(f, k1+2, tmp1, 1, ignoredx=igdx)
      call der(f, k1+1, tmp2, 2, ignoredx=igdx)
      g(:,3)=tmp1-tmp2
!
!  Adjustments for spherical coordinate system.
!
      if (.not. igdx .and. lspherical_coords) then
        g(:,1)=g(:,1)+f(l1:l2,m,n,k1+3)*r1_mn*cotth(m)
        g(:,2)=g(:,2)-f(l1:l2,m,n,k1+3)*r1_mn
        g(:,3)=g(:,3)+f(l1:l2,m,n,k1+2)*r1_mn
      endif
!
      if (.not. igdx .and. lcylindrical_coords) then
        g(:,3)=g(:,3)+f(l1:l2,m,n,k1+2)*rcyl_mn1
      endif
!
    endsubroutine curl
!***********************************************************************
    subroutine curl_other(f,g)
!
!  Calculate curl of a non-mvar vector, get vector.
!
!  23-june-09/wlad: adapted from curl
!
      use Deriv, only: der
!
      real, dimension (mx,my,mz,3) :: f
      real, dimension (nx,3) :: g
      real, dimension (nx) :: tmp1,tmp2
!
      intent(in) :: f
      intent(out) :: g
!
      call der(f(:,:,:,3),tmp1,2)
      call der(f(:,:,:,2),tmp2,3)
      g(:,1)=tmp1-tmp2
!
      call der(f(:,:,:,1),tmp1,3)
      call der(f(:,:,:,3),tmp2,1)
      g(:,2)=tmp1-tmp2
!
      call der(f(:,:,:,2),tmp1,1)
      call der(f(:,:,:,1),tmp2,2)
      g(:,3)=tmp1-tmp2
!
!  Adjustments for spherical corrdinate system.
!
      if (lspherical_coords) then
        g(:,1)=g(:,1)+f(l1:l2,m,n,3)*r1_mn*cotth(m)
        g(:,2)=g(:,2)-f(l1:l2,m,n,3)*r1_mn
        g(:,3)=g(:,3)+f(l1:l2,m,n,2)*r1_mn
      endif
!
      if (lcylindrical_coords) then
        g(:,3)=g(:,3)+f(l1:l2,m,n,2)*rcyl_mn1
      endif
!
    endsubroutine curl_other
!***********************************************************************
    subroutine curli(f,k,g,i)
!
!  Calculate curl of a vector, get vector.
!
!  22-oct-02/axel+tarek: adapted from curl
!
      use Deriv, only: der
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: g
      real, dimension (nx) :: tmp1,tmp2
      integer :: k,k1,i
!
      intent(in) :: f,k,i
      intent(out) :: g
!
      k1=k-1
!
      select case (i)
!
!  1-component
!
      case (1)
        call der(f,k1+3,tmp1,2)
        call der(f,k1+2,tmp2,3)
!
!  corrections for spherical coordinates
!
        if (lspherical_coords) then
          g=tmp1-tmp2+f(l1:l2,m,n,k1+3)*r1_mn*cotth(m)
        else
          g=tmp1-tmp2
        endif
!
!  2-component
!
      case (2)
        call der(f,k1+1,tmp1,3)
        call der(f,k1+3,tmp2,1)
!
!  corrections for spherical coordinates
!
        if (lspherical_coords) then
          g=tmp1-tmp2-f(l1:l2,m,n,k1+3)*r1_mn
        else
          g=tmp1-tmp2
        endif
!
!  3-component
!
      case (3)
        call der(f,k1+2,tmp1,1)
        call der(f,k1+1,tmp2,2)
!
!  corrections for spherical and cylindrical coordinates
!
        if (lspherical_coords) then
          g=tmp1-tmp2+f(l1:l2,m,n,k1+2)*r1_mn
        elseif (lcylindrical_coords) then
          g=tmp1-tmp2+f(l1:l2,m,n,k1+2)*rcyl_mn1
        else
          g=tmp1-tmp2
        endif
!
      endselect
!
    endsubroutine curli
!***********************************************************************
    subroutine del2_main(f,k,del2f)
!
!  Calculate del2 of a scalar, get scalar.
!
!  12-sep-97/axel: coded
!   7-mar-07/wlad: added cylindrical coordinates
!
      use Deriv, only: der,der2
!
      intent(in) :: f,k
      intent(out) :: del2f
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: del2f,d2fdx,d2fdy,d2fdz,tmp
      integer :: k
!
      call der2(f,k,d2fdx,1)
      call der2(f,k,d2fdy,2)
      call der2(f,k,d2fdz,3)
      del2f=d2fdx+d2fdy+d2fdz
!
      if (lcylindrical_coords) then
        call der(f,k,tmp,1)
        del2f=del2f+tmp*rcyl_mn1
      endif
!
      if (lspherical_coords) then
        call der(f,k,tmp,1)
        del2f=del2f+2.*r1_mn*tmp
        call der(f,k,tmp,2)
        del2f=del2f+cotth(m)*r1_mn*tmp
      endif
!
    endsubroutine del2_main
!***********************************************************************
    subroutine del2_other(f,del2f)
!
!  Calculate del2 of a scalar, get scalar.
!   8-may-09/nils: adapted from del2
!
      use Deriv, only: der,der2
!
      intent(in) :: f
      intent(out) :: del2f
!
      real, dimension (mx,my,mz) :: f
      real, dimension (nx) :: del2f,d2fdx,d2fdy,d2fdz,tmp
!
      call der2(f,d2fdx,1)
      call der2(f,d2fdy,2)
      call der2(f,d2fdz,3)
      del2f=d2fdx+d2fdy+d2fdz
!
      if (lcylindrical_coords) then
        call der(f,tmp,1)
        del2f=del2f+tmp*rcyl_mn1
      endif
!
      if (lspherical_coords) then
        call der(f,tmp,1)
        del2f=del2f+2.*r1_mn*tmp
        call der(f,tmp,2)
        del2f=del2f+cotth(m)*r1_mn*tmp
      endif
!
    endsubroutine del2_other
!***********************************************************************
    subroutine del2v(f,k,del2f,fij,pff)
!
!  Calculate del2 of a vector, get vector.
!
!  28-oct-97/axel: coded
!  15-mar-07/wlad: added cylindrical coordinates
!
      use Deriv, only: der
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, optional, dimension(nx,3,3) :: fij
      real, optional, dimension(nx,3) :: pff
      real, dimension (nx,3) :: del2f
      real, dimension (nx) :: tmp
      integer :: i,k,k1
!
      intent(in) :: f,k
      intent(out) :: del2f
!
!  do the del2 diffusion operator
!
      k1=k-1
      do i=1,3
        call del2(f,k1+i,tmp)
        del2f(:,i)=tmp
      enddo
!
      if (lcylindrical_coords) then
        !del2 already contains the extra term 1/r*d(uk)/dt
        call der(f,k1+2,tmp,2)
        del2f(:,1)=del2f(:,1) -(2*tmp+f(l1:l2,m,n,k1+1))*rcyl_mn2
        call der(f,k1+1,tmp,2)
        del2f(:,2)=del2f(:,2) +(2*tmp-f(l1:l2,m,n,k1+2))*rcyl_mn2
      endif
!
      if (lspherical_coords) then
        if (.not. (present(fij) .and. present(pff))) then
           call fatal_error('del2v', &
               'Cannot do a spherical del2v without aij and aa')
        endif
!
! for r component (factors of line elements are taken care of inside p%uij
        del2f(:,1)= del2f(:,1)+&
            r1_mn*(2.*(fij(:,1,1)-fij(:,2,2)-fij(:,3,3) &
            -r1_mn*pff(:,1)-cotth(m)*r1_mn*pff(:,2) ) &
            +cotth(m)*fij(:,1,2) )
! for theta component
        del2f(:,2)=del2f(:,2)+&
            r1_mn*(2.*(fij(:,2,1)-cotth(m)*fij(:,3,3)&
            +fij(:,1,2) )&
            +cotth(m)*fij(:,2,2)-r1_mn*sin1th(m)*sin1th(m)*pff(:,2) )
! for phi component
        del2f(:,3)=del2f(:,3)+&
            r1_mn*(2.*(fij(:,3,1)+fij(:,1,3)&
            +cotth(m)*fij(:,2,3) ) &
            +cotth(m)*fij(:,3,2)-sin1th(m)*pff(:,3) )
      endif
!
    endsubroutine del2v
!***********************************************************************
    subroutine del2m3x3_sym(f,k,del2f)
!
!  Calculate del2 of a 3x3 symmetric matrix, get matrix
!  23-feb-11/dhruba: coded in a new manner
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx,3,3) :: del2f
      real, dimension (nx) :: tmp
      integer :: i,j,k,k1
!
      intent(in) :: f,k
      intent(out) :: del2f
!
!  simple del2 diffusion operator for each component
!
      k1=k-1
      do i=1,3
        do j=i,3
          k1=k1+1
          call del2(f,k1,tmp)
          del2f(:,i,j) = tmp
        enddo
      enddo
      call symmetrise3x3_ut2lt(del2f)
!
      if (lcylindrical_coords) then
        call fatal_error('del2m', &
            'del2m is not implemented in cylindrical coordiates')
      endif
!
      if (lspherical_coords) then
        call fatal_error('del2m', &
            'del2m is not implemented in spherical coordiates')
      endif
!
    endsubroutine del2m3x3_sym
!***********************************************************************
    subroutine del2fj(f,vec,k,del2f)
!
!  Calculate del2 of a scalar, get scalar, adding anisotropic factor.
!
!  20-dec-12/wlad: adapted from del2_main
!
      use Deriv, only: der,der2
!
      intent(in) :: f,k
      intent(out) :: del2f
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: del2f,d2fdx,d2fdy,d2fdz,tmp
      real, dimension (3) :: vec
      integer :: k
!
      call der2(f,k,d2fdx,1)
      call der2(f,k,d2fdy,2)
      call der2(f,k,d2fdz,3)
      del2f=vec(1)*d2fdx+vec(2)*d2fdy+vec(3)*d2fdz
!
      if (lcylindrical_coords.and.vec(1)/=0.) then
        call der(f,k,tmp,1)
        del2f=del2f+vec(1)*tmp*rcyl_mn1
      endif
!
      if (lspherical_coords) then
        if (vec(1)/=0.) then
          call der(f,k,tmp,1)
          del2f=del2f+vec(1)*2.*r1_mn*tmp
        endif
        if (vec(2)/=0.) then
          call der(f,k,tmp,2)
          del2f=del2f+vec(2)*cotth(m)*r1_mn*tmp
        endif
      endif
!
    endsubroutine del2fj
!***********************************************************************
    subroutine del2fi_dxjk(f,k,del2fkdxij)
!
!  Calculate \partial^2f/\partial x_j\partial x_k of a vector, get a 9 dimensional object
!
      use Deriv, only: der2,derij
!
      intent(in) :: f,k
      intent(out) :: del2fkdxij
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx,3,3,3) :: del2fkdxij
      real, dimension(nx) :: tmp
      integer :: k
      integer :: i,j,kincrement
!
!
      if (.not. lcartesian_coords) &
        call fatal_error('sub','del2fi_dxjk implemented only in cartesian')
!
      do kincrement=0,2
       do i=1,3; do j=i,3
         tmp=0.
         if (i.eq.j) then
           call der2(f,k+kincrement,tmp,i)
           del2fkdxij(:,k+kincrement,i,i)=tmp
         else
           call derij(f,k+kincrement,tmp,i,j)
           del2fkdxij(:,k+kincrement,i,j)=tmp
           del2fkdxij(:,k+kincrement,j,i)=tmp
         endif
       enddo;enddo
     enddo
!
   endsubroutine del2fi_dxjk
!***********************************************************************
    subroutine d2fi_dxj(f,k,d2fidxj)
!
!  Calculate d^2f_i/dx^2_j of a vector, get a six dimensional object
!
      use Deriv, only: derij
!
      intent(in) :: f,k
      intent(out) :: d2fidxj
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx,3,3) :: d2fidxj
      real, dimension(nx,3) :: tmp
      integer :: k
!
      if (.not. lcartesian_coords) &
        call fatal_error('sub','d2fidxj implemented only in cartesian')
!
      call d2f_dxj(f,k,tmp)
      d2fidxj(:,1,:) = tmp
      call d2f_dxj(f,k+1,tmp)
      d2fidxj(:,2,:) = tmp
      call d2f_dxj(f,k+2,tmp)
      d2fidxj(:,3,:) = tmp
!
    endsubroutine d2fi_dxj
!***********************************************************************
    subroutine d2f_dxj(f,k,d2fdxj)
!
!  Calculate d^2f/dx^2_j of a scalar, get a three dimensional object
!
      use Deriv, only: der,der2
!
      intent(in) :: f,k
      intent(out) :: d2fdxj
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx,3) :: d2fdxj
      integer :: k
!
      if (.not. lcartesian_coords) &
        call fatal_error('sub','d2f_dxj implemented only in cartesian')
!
      call der2(f,k,d2fdxj(:,1),1)
      call der2(f,k,d2fdxj(:,2),2)
      call der2(f,k,d2fdxj(:,3),3)
!
    endsubroutine d2f_dxj
!***********************************************************************
    subroutine symmetrise3x3_ut2lt (matrix_ut3x3)
!
! sets the lower triangular values of a matrix to its upper
! triangular values. Does not touch the diagonal. Applies
! to 3x3 matrices (pencil) only.
!
!  23-dhruba-11/dhruba: coded
!
      real, dimension(nx,3,3) :: matrix_ut3x3
      integer :: i,j
!
      do i=1,3
        do j=1,i
          if (i/=j) &
            matrix_ut3x3(:,i,j) = matrix_ut3x3(:,j,i)
        enddo
      enddo
!
    endsubroutine symmetrise3x3_ut2lt
!***********************************************************************
    subroutine del2v_etc(f,k,del2,graddiv,curlcurl,gradcurl)
!
!  Calculates a number of second derivative expressions of a vector
!  outputs a number of different vector fields.
!  gradcurl is not the vector gradient.
!  Surprisingly, calling derij only if graddiv or curlcurl are present
!  does not speed up the code on Mephisto @ 32x32x64.
!
!  12-sep-01/axel: coded
!  15-mar-07/wlad: added cylindrical coordinates
!
      use Deriv, only: der,der2,derij
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx,3,3) :: fjji,fijj
      real, dimension (nx,3,3), optional :: gradcurl
      real, dimension (nx,3), optional :: del2,graddiv,curlcurl
      real, dimension (nx,3) ::  fjik
      real, dimension (nx) :: tmp
      integer :: i,j,k,k1
!
      intent(in) :: f,k
      intent(out) :: del2,graddiv,curlcurl,gradcurl
!
!  calculate f_{i,jj} and f_{j,ji}
!  AJ: graddiv needs diagonal elements from the first tmp (derij only sets
!      off-diagonal elements)
!
      k1=k-1
      do i=1,3
      do j=1,3
        if (present(del2) .or. present(curlcurl) .or. present(gradcurl) .or. &
            present(graddiv)) then
          call der2 (f,k1+i,tmp,  j); fijj(:,i,j)=tmp  ! f_{i,jj}
        endif
        if (present(graddiv) .or. present(curlcurl).or. present(gradcurl)) then
          call derij(f,k1+j,tmp,j,i); fjji(:,i,j)=tmp  ! f_{j,ji}
        endif
      enddo
      enddo
!
!  the diagonal terms have not been set in derij; do this now
!  ** They are automatically set above, because derij   **
!  ** doesn't overwrite the value of tmp for i=j!       **
!
!     do j=1,3
!       fjji(:,j,j)=fijj(:,j,j)
!     enddo
!
!  calculate f_{i,jk} for i /= j /= k
!
      if (present(gradcurl)) then
        call derij(f,k1+1,tmp,2,3)
        fjik(:,1)=tmp
        call derij(f,k1+2,tmp,1,3)
        fjik(:,2)=tmp
        call derij(f,k1+3,tmp,1,2)
        fjik(:,3)=tmp
        if (lcylindrical_coords.or.lspherical_coords) &
          call fatal_error('del2v_etc','gradcurl at del2v_etc not '// &
              'implemented for non-cartesian coordinates')
      endif
!
!  del2
!
      if (present(del2)) then
        do i=1,3
          del2(:,i)=fijj(:,i,1)+fijj(:,i,2)+fijj(:,i,3)
        enddo
        if (lcylindrical_coords) then
          !r-component
          call der(f,k1+2,tmp,2)
          del2(:,1)=del2(:,1) -(2*tmp+f(l1:l2,m,n,k1+1))*rcyl_mn2
          call der(f,k1+1,tmp,1)
          del2(:,1)=del2(:,1) + tmp*rcyl_mn1
          !phi-component
          call der(f,k1+1,tmp,2)
          del2(:,2)=del2(:,2) +(2*tmp-f(l1:l2,m,n,k1+2))*rcyl_mn2
          call der(f,k1+2,tmp,1)
          del2(:,2)=del2(:,2) + tmp*rcyl_mn1
          !z-component
          call der(f,k1+3,tmp,1)
          del2(:,3)=del2(:,3) + tmp*rcyl_mn1
        endif
        if (lspherical_coords) call fatal_error('del2v_etc', &
            'del2 at del2v_etc not implemented for spherical '// &
            'coordinates - use gij_etc')
      endif
!
      if (present(graddiv)) then
        do i=1,3
          graddiv(:,i)=fjji(:,i,1)+fjji(:,i,2)+fjji(:,i,3)
        enddo
        if (lcylindrical_coords) then
           call der(f,k1+1,tmp,1)
           graddiv(:,1)=graddiv(:,1)+tmp*rcyl_mn1 - f(l1:l2,m,n,k1+1)*rcyl_mn2
           call der(f,k1+1,tmp,2)
           graddiv(:,2)=graddiv(:,2)+tmp*rcyl_mn1
           call der(f,k1+1,tmp,3)
           graddiv(:,3)=graddiv(:,3)+tmp*rcyl_mn1
        endif
        if (lspherical_coords) call fatal_error('del2v_etc', &
            'graddiv is implemented in gij_etc for spherical '// &
            'coords - use gij_etc')
      endif
!
!  curlcurl
!
      if (present(curlcurl)) then
        curlcurl(:,1)=fjji(:,1,2)-fijj(:,1,2)+fjji(:,1,3)-fijj(:,1,3)
        curlcurl(:,2)=fjji(:,2,3)-fijj(:,2,3)+fjji(:,2,1)-fijj(:,2,1)
        curlcurl(:,3)=fjji(:,3,1)-fijj(:,3,1)+fjji(:,3,2)-fijj(:,3,2)
        if (lcylindrical_coords) call fatal_error('del2v_etc', &
            'curlcurl not implemented for non-cartesian coordinates')
        if (lspherical_coords) call fatal_error('del2v_etc', &
            'curlcurl not implemented for non-cartesian coordinates - '// &
            'use gij_etc')
      endif
!
!  gradcurl (as tensor)
!
      if (present(gradcurl)) then
         gradcurl(:,1,1) = fjik(:,3)   - fjik(:,2)
         gradcurl(:,1,2) = fjji(:,1,3) - fijj(:,3,1)
         gradcurl(:,1,3) = fijj(:,2,1) - fjji(:,1,2)
!
         gradcurl(:,2,1) = fijj(:,3,2) - fjji(:,2,3)
         gradcurl(:,2,2) = fjik(:,1)   - fjik(:,3)
         gradcurl(:,2,3) = fjji(:,2,1) - fijj(:,1,2)
!
         gradcurl(:,3,1) = fjji(:,3,2) - fijj(:,2,3)
         gradcurl(:,3,2) = fijj(:,1,3) - fjji(:,3,1)
         gradcurl(:,3,3) = fjik(:,2)   - fjik(:,1)
!
         if (lcylindrical_coords.or.lspherical_coords) &
             call fatal_error('del2v_etc','gradcurl not implemented '// &
             'for non-cartesian coordinates')
      endif
!
    endsubroutine del2v_etc
!***********************************************************************
    subroutine del2vi_etc(f,k,ii,del2,graddiv,curlcurl)
!
!  Calculates a number of second derivative expressions of a vector.
!  Outputs a number of different vector fields.
!  Surprisingly, calling derij only if graddiv or curlcurl are present
!  does not speed up the code on Mephisto @ 32x32x64.
!  Just do the ith component
!
!   7-feb-04/axel: adapted from del2v_etc
!
      use Deriv, only: der2,derij
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx,3,3) :: fjji,fijj
      real, dimension (nx), optional :: del2,graddiv,curlcurl
      real, dimension (nx) :: tmp
      integer :: i,j,k,k1,ii
!
      intent(in) :: f,k,ii
      intent(out) :: del2,graddiv,curlcurl
!
!  Do the del2 diffusion operator.
!
      k1=k-1
      do i=1,3
      do j=1,3
        call der2 (f,k1+i,tmp,  j); fijj(:,i,j)=tmp  ! f_{i,jj}
        call derij(f,k1+j,tmp,j,i); fjji(:,i,j)=tmp  ! f_{j,ji}
      enddo
      enddo
!
      if (present(del2)) then
        del2=fijj(:,ii,1)+fijj(:,ii,2)+fijj(:,ii,3)
      endif
!
      if (present(graddiv)) then
        graddiv=fjji(:,ii,1)+fjji(:,ii,2)+fjji(:,ii,3)
      endif
!
      if (present(curlcurl)) then
        select case (ii)
        case (1); curlcurl=fjji(:,1,2)-fijj(:,1,2)+fjji(:,1,3)-fijj(:,1,3)
        case (2); curlcurl=fjji(:,2,3)-fijj(:,2,3)+fjji(:,2,1)-fijj(:,2,1)
        case (3); curlcurl=fjji(:,3,1)-fijj(:,3,1)+fjji(:,3,2)-fijj(:,3,2)
        endselect
      endif
!
      if (lcylindrical_coords.or.lspherical_coords) &
          call fatal_error('del2vi_etc', &
          'not implemented for non-cartesian coordinates')
!
    endsubroutine del2vi_etc
!***********************************************************************
    subroutine del4v(f,k,del4f)
!
!  Calculate del4 of a vector, get vector.
!
!  09-dec-03/nils: adapted from del6v
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx,3) :: del4f
      real, dimension (nx) :: tmp
      integer :: i,k,k1
!
      intent(in) :: f,k
      intent(out) :: del4f
!
!  Exit if this is requested for non-cartesian runs.
!
      if (lcylindrical_coords.or.lspherical_coords) &
          call fatal_error('del4v', &
          'not implemented for non-cartesian coordinates')
!
!  Do the del4 diffusion operator.
!
      k1=k-1
      do i=1,3
        call del4(f,k1+i,tmp)
        del4f(:,i)=tmp
      enddo
!
    endsubroutine del4v
!***********************************************************************
    subroutine del6v(f,k,del6f,lstrict)
!
!  Calculate del6 of a vector, get vector.
!
!  28-oct-97/axel: coded
!  24-apr-03/nils: adapted from del2v
!
      use General, only : loptest
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx,3) :: del6f
      real, dimension (nx) :: tmp
      integer :: i,k,k1
      logical, optional :: lstrict
!
      intent(in) :: f,k
      intent(out) :: del6f
!
!  Do the del6 diffusion operator.
!
      k1=k-1
      do i=1,3
        if (loptest(lstrict)) then
          call del6_strict(f,k1+i,tmp)
        else
          call del6(f,k1+i,tmp)
        endif
        del6f(:,i)=tmp
      enddo
!
!  Exit if this is requested for non-cartesian runs.
!
      !!!!if (lcylindrical_coords.or.lspherical_coords) &
      !!!!    call fatal_error('del6v', &
      !!!!    'not implemented for non-cartesian coordinates')
!
    endsubroutine del6v
!***********************************************************************
    subroutine bij_tilde(f,bb,bijtilde,bij_cov_corr)
!
! Calculates \partial B_[r,\theta,\phi]/ \partial r, B_[r,\theta,\phi]/(r \theta),
! and optionally the correction, needed for covariant derivatives, in bij_cov_corr.
!
! 20-nov-16/MR: coded
!
      use Deriv, only: der2,derij,der
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (nx,3),             intent(in) :: bb
      real, dimension (nx,3,*),           intent(out):: bijtilde
      real, dimension (nx,3,3), optional, intent(out):: bij_cov_corr
 
      real, dimension (nx,3) :: bbr1
      real, dimension (nx,2:3) :: aar2
      real, dimension (nx) :: d2adrdt,tmp,tmp1
      integer :: i

      do i=1,3 
        bbr1(:,i)=bb(:,i)*r1_mn
      enddo
!
!  d B_r/dr
!      
      do i=1,2 
        aar2(:,i+1)=f(l1:l2,m,n,iaa+i)*r2_mn
      enddo
!
      call derij(f,iaz,d2adrdt,1,2)     ! (1/r) d^2 a_phi/dr dtheta

      bijtilde(:,1,1) = d2adrdt - (bbr1(:,1) + cotth(m)*(bbr1(:,2)+aar2(:,3)))
!
!  (1/r) d B_r/d theta
!      
      call der2(f,iaz,tmp,2)            ! (1/r^2) d^2 a_phi/dtheta^2
      bijtilde(:,1,2) = tmp - aar2(:,3)*(cotth(m)*cotth(m) + sin2th(m)) + cotth(m)*bbr1(:,1)
!
!  d B_theta/dr
!      
      call der2(f,iaz,tmp,1)            ! d^2 a_phi/dr^2
      bijtilde(:,2,1) = -tmp + (bbr1(:,2)+2.*aar2(:,3))
!
!  (1/r) d B_theta/d theta
!      
      bijtilde(:,2,2) = -d2adrdt - (bbr1(:,1) - cotth(m)*aar2(:,3))
!
!  d B_phi/dr
!      
      call der2(f,iay,tmp,1)            ! d^2 a_theta/dr^2
      call derij(f,iax,d2adrdt,1,2)     ! (1/r) d^2 a_r/dr dtheta
      call der(f,iax,tmp1,2)            ! (1/r) d a_r/dtheta
      bijtilde(:,3,1) = tmp - d2adrdt + bbr1(:,3) - 2.*(aar2(:,2)-tmp1*r1_mn)
!
!  (1/r) d B_phi/d theta
!      
      call der2(f,iax,tmp,2)            ! (1/r^2) d^2 a_r/dtheta^2
      call derij(f,iay,d2adrdt,1,2)     ! (1/r) d^2 a_theta/dr dtheta
      call der(f,iay,tmp1,2)            ! (1/r) d a_theta/dtheta
      bijtilde(:,3,2) = tmp1*r1_mn - (tmp - d2adrdt)
 
      bijtilde(:,:,3)=0. 
      if (present(bij_cov_corr)) then
        bij_cov_corr(:,:,1) = 0.; bij_cov_corr(:,3,2)=0.
        bij_cov_corr(:,1,2) = -bbr1(:,2)
        bij_cov_corr(:,1,3) = -bbr1(:,3)
        bij_cov_corr(:,2,2) =  bbr1(:,1)
        bij_cov_corr(:,2,3) = -cotth(m)*bbr1(:,3)
        bij_cov_corr(:,3,3) =  cotth(m)*bbr1(:,2)+bbr1(:,1)
      endif
 
    endsubroutine bij_tilde 
!***********************************************************************
    subroutine gij_etc(f,iref,aa,aij,bij,del2,graddiv,lcovariant_derivative)
!
!  Calculate B_i,j = eps_ikl A_l,jk and A_l,kk.
!
!  21-jul-03/axel: coded
!  26-jul-05/tobi: do not calculate both d^2 A/(dx dy) and d^2 A/(dy dx)
!  23-feb-07/axel: added spherical coordinates
!   7-mar-07/wlad: added cylindrical coordinates
!  29-aug-13/MR: made bij optional; added error messages for missing optional parameters
!  15-jun-16/fred: bij optional covariant derivative terms; spherical/cylindrical
!
      use Deriv, only: der2,derij
      use General, only: loptest
!
      real, dimension (mx,my,mz,mfarray), intent (in) :: f
      integer, intent (in) :: iref
      logical, intent (in), optional :: lcovariant_derivative
      real, dimension (nx,3), intent (in), optional :: aa
      real, dimension (nx,3,3), intent (in), optional :: aij
      real, dimension (nx,3,3), intent (out), optional :: bij
      real, dimension (nx,3), intent (out), optional :: del2,graddiv
!
!  Locally used variables.
!
      real, dimension (nx,3,3,3) :: d2A
      integer :: iref1,i,j
!
!  Reference point of argument.
!
      iref1=iref-1
!
!  Calculate all (mixed and non-mixed) second derivatives
!  of the vector potential (A_k,ij).
!
!  Do not calculate both d^2 A_k/(dx dy) and d^2 A_k/(dy dx).
!  (This wasn't spotted by me but by a guy from SGI...)
!  Note: for non-cartesian coordinates there are different correction terms,
!  see below.
!
      do i=1,3
        do j=1,3
          call der2(f,iref1+i,d2A(:,j,j,i),j)
        enddo
        call derij(f,iref1+i,d2A(:,2,3,i),2,3); d2A(:,3,2,i)=d2A(:,2,3,i)
        call derij(f,iref1+i,d2A(:,3,1,i),3,1); d2A(:,1,3,i)=d2A(:,3,1,i)
        call derij(f,iref1+i,d2A(:,1,2,i),1,2); d2A(:,2,1,i)=d2A(:,1,2,i)
      enddo
!
!  Corrections for spherical polars from swapping mixed derivatives:
!  Psi_{,theta^ r^} = Psi_{,r^ theta^} - Psi_{,\theta^}/r
!  Psi_{,phi^ r^} = Psi_{,r^ phi^} - Psi_{,\phi^}/r
!  Psi_{,phi^ theta^} = Psi_{,theta^ phi^} - Psi_{,\phi^}*r^{-1}*cot(theta)
!
      if (lspherical_coords) then
        if (.not.present(aij)) &
          call fatal_error('gij_etc', 'aij needed for spherical coordinates')
        do i=1,3
          d2A(:,2,1,i)=d2A(:,2,1,i)-aij(:,i,2)*r1_mn
          d2A(:,3,1,i)=d2A(:,3,1,i)-aij(:,i,3)*r1_mn
          d2A(:,3,2,i)=d2A(:,3,2,i)-aij(:,i,3)*r1_mn*cotth(m)
        enddo
      endif
!
!  For cylindrical, only
!  Psi_{,phi^ pom^} = Psi_{,pom^ phi^} - Psi_{,\phi^}/pom .
!
      if (lcylindrical_coords) then
        if (.not.present(aij)) &
          call fatal_error('gij_etc', 'aij needed for cylindrical coordinates')
        do i=1,3
          d2A(:,2,1,i)=d2A(:,2,1,i)-aij(:,i,2)*rcyl_mn1
        enddo
      endif
!
!  Calculate optionally b_i,j = eps_ikl A_l,kj,
!  del2_i = A_i,jj and graddiv_i = A_j,ji .
!
      if (present(bij)) then
!
        bij(:,1,:)=d2A(:,2,:,3)-d2A(:,3,:,2)
        bij(:,2,:)=d2A(:,3,:,1)-d2A(:,1,:,3)
        bij(:,3,:)=d2A(:,1,:,2)-d2A(:,2,:,1)
!
!  Corrections for spherical coordinates.
!
        if (lspherical_coords) then
          if (.not.present(aa)) &
            call fatal_error('gij_etc', 'aa needed for spherical coordinates')
          bij(:,3,2)=bij(:,3,2)+aij(:,2,2)*r1_mn
          bij(:,2,3)=bij(:,2,3)-aij(:,3,3)*r1_mn
          bij(:,1,3)=bij(:,1,3)+aij(:,3,3)*r1_mn*cotth(m)
          bij(:,3,1)=bij(:,3,1)+aij(:,2,1)*r1_mn         -aa(:,2)*r2_mn
          bij(:,2,1)=bij(:,2,1)-aij(:,3,1)*r1_mn         +aa(:,3)*r2_mn
          bij(:,1,2)=bij(:,1,2)+aij(:,3,2)*r1_mn*cotth(m)-aa(:,3)*r2_mn*sin2th(m)
          if (loptest(lcovariant_derivative)) then
            bij(:,1,1)=bij(:,1,1)+(aij(:,3,2)*r1_mn-aa(:,3)*r2_mn)*cotth(m) ! MR: not correct, but strangely bij(:,1,1)
                                                                            !     does need covariant correction. Why? 
            bij(:,1,2)=bij(:,1,2)+(aij(:,3,1)-aij(:,1,3))*r1_mn+aa(:,3)*r2_mn
            bij(:,1,3)=bij(:,1,3)+(aij(:,1,2)-aij(:,2,1))*r1_mn-aa(:,2)*r2_mn
            bij(:,2,2)=bij(:,2,2)-aij(:,2,3)*r1_mn         +aa(:,3)*r2_mn*cotth(m)
            bij(:,2,3)=bij(:,2,3)+(aij(:,1,2)*r1_mn-aij(:,2,1)*r1_mn-aa(:,2)*r2_mn)&
                                                                         *cotth(m)
            bij(:,3,3)=bij(:,3,3)+(aij(:,3,2)+(aij(:,1,3)-aij(:,3,1))*cotth(m))*r1_mn
          endif
        endif
!
!  Corrections for cylindrical coordinates.
!
        if (lcylindrical_coords) then
          if (.not.present(aa)) &
            call fatal_error('gij_etc', 'aa needed for cylindrical coordinates')
          bij(:,3,2)=bij(:,3,2)+ aij(:,2,2)*rcyl_mn1
!          !bij(:,3,1)=bij(:,3,1)+(aij(:,2,1)+aij(:,1,2))*rcyl_mn1-aa(:,2)*rcyl_mn2
!  FAG:Partial correction to -d2A(:,2,1,1) already applied above +aij(:,i,2)*rcyl_mn1
!
          bij(:,3,1)=bij(:,3,1)+aij(:,2,1)*rcyl_mn1-aa(:,2)*rcyl_mn2
          if (loptest(lcovariant_derivative)) then
            bij(:,1,2)=bij(:,1,2)+(aij(:,3,1)-aij(:,1,3))*rcyl_mn1
            bij(:,2,2)=bij(:,2,2)+(aij(:,3,2)-aij(:,2,3))*rcyl_mn1
            bij(:,3,3)=bij(:,3,3)+aij(:,2,3)*rcyl_mn1   ! MR: correct?
          endif
        endif
      endif
!
!  Calculate del2 and graddiv, if requested.
!
      if (present(graddiv)) then
        graddiv(:,:)=d2A(:,1,:,1)+d2A(:,2,:,2)+d2A(:,3,:,3)
        if (lspherical_coords) then
          if (.not.present(aa)) &
            call fatal_error('gij_etc', 'aa needed for spherical coordinates')
          graddiv(:,1)=graddiv(:,1)+aij(:,1,1)*r1_mn*2+ &
             aij(:,2,1)*r1_mn*cotth(m)-aa(:,2)*r2_mn*cotth(m)-aa(:,1)*r2_mn*2
          graddiv(:,2)=graddiv(:,2)+aij(:,1,2)*r1_mn*2+ &
             aij(:,2,2)*r1_mn*cotth(m)-aa(:,2)*r2_mn*sin2th(m)
          graddiv(:,3)=graddiv(:,3)+aij(:,1,3)*r1_mn*2+ &
             aij(:,2,3)*r1_mn*cotth(m)
        endif
      endif
!
      if (present(del2)) then
        del2(:,:)=d2A(:,1,1,:)+d2A(:,2,2,:)+d2A(:,3,3,:)
        if (lspherical_coords) then
          if (.not.present(aa)) &
            call fatal_error('gij_etc', 'aa needed for spherical coordinates')
          del2(:,1)= del2(:,1)+&
            r1_mn*(2.*(aij(:,1,1)-aij(:,2,2)-aij(:,3,3)&
            -r1_mn*aa(:,1)-cotth(m)*r1_mn*aa(:,2) ) &
            +cotth(m)*aij(:,1,2) )
          del2(:,2)=del2(:,2)+&
            r1_mn*(2.*(aij(:,2,1)-cotth(m)*aij(:,3,3)&
            +aij(:,1,2) )&
            +cotth(m)*aij(:,2,2)-r1_mn*sin2th(m)*aa(:,2) )
          del2(:,3)=del2(:,3)+&
            r1_mn*(2.*(aij(:,3,1)+aij(:,1,3)&
            +cotth(m)*aij(:,2,3) ) &
            +cotth(m)*aij(:,3,2)-r1_mn*sin2th(m)*aa(:,3) )
        endif
        if (lcylindrical_coords) then
          if (.not.present(aa)) &
            call fatal_error('gij_etc', 'aa needed for cylindrical coordinates')
          del2(:,1)= del2(:,1)+&
            rcyl_mn1*(aij(:,1,1)-2*aij(:,2,2))-rcyl_mn2*aa(:,1)
          del2(:,2)=del2(:,2)+&               
            rcyl_mn1*(aij(:,2,1)-2*aij(:,1,2))-rcyl_mn2*aa(:,2)
        endif
      endif
!
    endsubroutine gij_etc
!***********************************************************************
    subroutine g2ij(f,k,g)
!
!  Calculates the Hessian, i.e. all second derivatives of a scalar.
!
!  11-jul-02/axel: coded
!
      use Deriv, only: der2,derij
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx,3,3) :: g
      real, dimension (nx) :: tmp
      integer :: i,j,k
!
      intent(in) :: f,k
      intent(out) :: g
!
!  Run though all 9 possibilities, treat diagonals separately.
!
      do j=1,3
        call der2 (f,k,tmp,j); g(:,j,j)=tmp
        do i=j+1,3
          call derij(f,k,tmp,i,j); g(:,i,j)=tmp; g(:,j,i)=tmp
        enddo
      enddo
!
    endsubroutine g2ij
!***********************************************************************
    subroutine del4(f,k,del4f,ignoredx)
!
!  Calculate del4 (defined here as d^4/dx^4 + d^4/dy^4 + d^4/dz^4, rather
!  than del2^3) of a scalar for hyperdiffusion.
!
!   8-jul-02/wolf: coded
!   9-dec-03/nils: adapted from del6
!
      use Deriv, only: der4
!
      intent(in) :: f,k,ignoredx
      intent(out) :: del4f
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: del4f,d4fdx,d4fdy,d4fdz
      integer :: k
      logical, optional :: ignoredx
      logical :: ignore_dx
!
      if (present(ignoredx)) then
        ignore_dx = ignoredx
      else
        ignore_dx = .false.
!
!  Exit if this is requested for non-cartesian runs.
!
        if (lcylindrical_coords.or.lspherical_coords) &
            call fatal_error('del4', &
            'not implemented for non-cartesian coordinates')
      endif
!
      call der4(f,k,d4fdx,1,ignore_dx)
      call der4(f,k,d4fdy,2,ignore_dx)
      call der4(f,k,d4fdz,3,ignore_dx)
      del4f = d4fdx + d4fdy + d4fdz
!
    endsubroutine del4
!***********************************************************************
    subroutine del6(f,k,del6f,ignoredx)
!
!  Calculate del6 (defined here as d^6/dx^6 + d^6/dy^6 + d^6/dz^6, rather
!  than del2^3) of a scalar for hyperdiffusion. Using INGOREDX
!  calculates something similar to del6, but ignoring the steps dx, dy, dz.
!  Useful for Nyquist filtering, where you just want to remove the
!  Nyquist frequency fully, while retaining the amplitude in small wave
!  numbers.
!
!  8-jul-02/wolf: coded
!  22-jul-11/bing: added ignoredx
!
      use Deriv, only: der6
!
      intent(in) :: f,k,ignoredx
      intent(out) :: del6f
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: del6f,d6fdx,d6fdy,d6fdz
      integer :: k
      logical, optional :: ignoredx
      logical :: ignore_dx
!
      if (present(ignoredx)) then
        ignore_dx = ignoredx
      else
        ignore_dx = .false.
!
!  Exit if this is requested for lspherical_coords run.
!
!!!!        if (lspherical_coords.or.lcylindrical_coords) &
!!!!            call fatal_error('del6', &
!!!!            'not implemented for non-cartesian coordinates')
      endif
!
      call der6(f,k,d6fdx,1,ignore_dx)
      call der6(f,k,d6fdy,2,ignore_dx)
      call der6(f,k,d6fdz,3,ignore_dx)
!
      del6f = d6fdx + d6fdy + d6fdz
!
    endsubroutine del6
!***********************************************************************
    subroutine del6_strict(f,k,del6)
!
!   Calculates del6rho=del2(del2(del2(rho))), with del2=div(grad). The routine is
!   strictly accurate for Cartesian coordinates, and retains all the leading dx1**6
!   terms in cylindrical. The subroutine is small enough and memory-cheap enough
!   that it could be in dlnrhodt. Yet, writing it as a subroutine allows
!   not only for encapsulation but also better documentation.
!
!          d6a    d6a   d6a     / d4d2a    d4d2a    d4d2a    d4d2a    d4d2a    d4d2a  \      d2d2d2a
!   del6 = --- +  --- + --- + 3 | ------ + ------ + ------ + ------ + ------ + ------ | + 6 ---------
!          dx6    dy6   dz6     \ dx4dy2   dx4dz2   dx2dy4   dx2dz4   d4ydz2   dz4dy2 /     dx2dy2dz2
!
!   02-apr-17/wlad: coded
!
      use Deriv, only: der6,der4i2j,der2i2j2k
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(nx) :: del6,tmp
      integer :: k,i,j
!      
      intent(in) :: f,k
      intent(out) :: del6
!
      del6=0.
      do i=1,3
        call der6(f,k,tmp,i)
        del6 = del6 + tmp
        do j=1,3
          if (j/=i) then
            call der4i2j(f,k,tmp,i,j)
            del6 = del6 + 3*tmp
          endif
        enddo
      enddo
      call der2i2j2k(f,k,tmp)
      del6 = del6 + 6*tmp
!
    endsubroutine del6_strict
!***********************************************************************
    subroutine del4graddiv(f,ikk,del4graddivu)
!
!  Calculate del4(grad(div())), which enters in the formulation of strict
!  hyperviscosity. This is strictly accurate for Cartesian, and retains all
!  the leading (dx1**6) terms for polar coordinates. For the x-component of
!  del4(grad(div())), the result of sympy is
!
!               [del4(grad(div(u)))]_x = f(ux) + g(uy) + g(uz)
!
!   where f(ux) = der6x(ux)
!               + 2*(der4x2y(ux)+der4x2z(ux)+der4y2x(ux)+der4z2x(ux))
!               + 4*der2x2y2z(ux)
!
!         g(uy) =   der5x1y(uy) + der1x5y(uy) + 3*der3x3y(uy)
!               + 2*der3x1y2z(uy) + 3*der1x3y2z(uy) + der1x1y4z(uy)
!
!   and similary
!
!         h(uz) =  der5x1z(uz) + der1x5z(uz) + 3*der3x3z(uz)
!               + 2*der3x2y1z(uz) + 3*der1x2y3z(uz) + der1x4y1z(uz)
!
!   Per symmetry, the formulation for the y and z components are
!   identical under the permutation [xyz].
!
!   09-apr-17/wlad: coded
!
      use Deriv, only: der6,der4i2j,der2i2j2k,der5i1j,der3i3j,der3i2j1k,der4i1j1k
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(nx,3) :: del4graddivu
      real, dimension(nx) :: tmp
      integer :: ikk,ki,kj
      integer :: i,j,k
!
      intent(in) :: f,ikk
      intent(out) :: del4graddivu
!
      if (ikk .ne. iuu) call fatal_error("del4graddiv",&
           "del4graddiv only coded for velocity")
!
!  The f(ui) part
!
      do i=1,3
        ki = ikk + (i-1)
        del4graddivu(:,i) = 0.
        call der6(f,ki,tmp,i);          del4graddivu(:,i) =  del4graddivu(:,i) + tmp
        do j=1,3
          if (j/=i) then
            call der4i2j(f,ki,tmp,i,j); del4graddivu(:,i) =  del4graddivu(:,i) + 2*tmp
            call der4i2j(f,ki,tmp,j,i); del4graddivu(:,i) =  del4graddivu(:,i) + 2*tmp
          endif
        enddo
        call der2i2j2k(f,ki,tmp) ;      del4graddivu(:,i) =  del4graddivu(:,i) + 4*tmp
!
!  The g(uj) and h(uk) parts
!
        do j=1,3
          if (j/=i) then
            if ((i==1).and.(j==2)) k=3
            if ((i==1).and.(j==3)) k=2
            if ((i==2).and.(j==1)) k=3
            if ((i==2).and.(j==3)) k=1
            if ((i==3).and.(j==1)) k=2
            if ((i==3).and.(j==2)) k=1
!
            kj = ikk+(j-1)
            call der5i1j(f,kj,tmp,i,j)
            del4graddivu(:,i) =  del4graddivu(:,i) + tmp
            call der5i1j(f,kj,tmp,j,i)
            del4graddivu(:,i) =  del4graddivu(:,i) + tmp
!
            call der3i3j(f,kj,tmp,i,j)
            del4graddivu(:,i) =  del4graddivu(:,i) + 3*tmp
!
            call der3i2j1k(f,kj,tmp,i,k,j)
            del4graddivu(:,i) =  del4graddivu(:,i) + 2*tmp
!
            call der3i2j1k(f,kj,tmp,j,k,i)
            del4graddivu(:,i) =  del4graddivu(:,i) + 3*tmp
!
            call der4i1j1k(f,kj,tmp,k,i,j)
            del4graddivu(:,i) =  del4graddivu(:,i) + tmp
          endif
        enddo
      enddo
!
    endsubroutine del4graddiv
!***********************************************************************
    subroutine del6_other(f,del6f)
!
!  Calculate del6 (defined here as d^6/dx^6 + d^6/dy^6 + d^6/dz^6, rather
!  than del2^3) of a scalar for hyperdiffusion.
!
!  13-jun-05/anders: adapted from del6
!
      use Deriv, only: der6_other
!
      intent(in) :: f
      intent(out) :: del6f
!
      real, dimension (mx,my,mz) :: f
      real, dimension (nx) :: del6f,d6fdx,d6fdy,d6fdz
!
      call der6_other(f,d6fdx,1)
      call der6_other(f,d6fdy,2)
      call der6_other(f,d6fdz,3)
      del6f = d6fdx + d6fdy + d6fdz
!
!  Exit if this is requested for non-cartesian runs.
!
      if (lcylindrical_coords.or.lspherical_coords) &
          call fatal_error('del6_other', &
          'not implemented for non-cartesian coordinates')
!
    endsubroutine del6_other
!***********************************************************************
    subroutine del6fj(f,vec,k,del6f)
!
!  Calculates fj*del6 (defined here as (vecx*d^6/dx^6 + vecy*d^6/dy^6 + vecz*d^6/dz^6)f )
!  needed for hyperdissipation of a scalar (diffrho) with non-cubic cells where
!  the coefficient depends on resolution. Returns scalar.
!
!  30-oct-06/wlad: adapted from del6
!
      use Deriv, only: der6
!
      intent(in) :: f,k
      intent(out) :: del6f
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: del6f,d6fdx,d6fdy,d6fdz
      real, dimension (3) :: vec
      integer :: k
!
      call der6(f,k,d6fdx,1)
      call der6(f,k,d6fdy,2)
      call der6(f,k,d6fdz,3)
      del6f = vec(1)*d6fdx + vec(2)*d6fdy + vec(3)*d6fdz
!
!  Exit if this is requested for non-cartesian runs.
!
      if (lcylindrical_coords.or.lspherical_coords) &
          call fatal_error('del6fj', &
          'not implemented for non-cartesian coordinates')
!
    endsubroutine del6fj
!***********************************************************************
    subroutine del6fjv(f,vec,k,del6f)
!
!  Calculates fj*del6 (defined here as fx*d^6/dx^6 + fy*d^6/dy^6 + fz*d^6/dz^6)
!  needed for hyperdissipation of vectors (visc, res) with non-cubic cells
!  where the coefficient depends on resolution. Returns vector.
!
!  30-oct-06/wlad: adapted from del6v
!
      intent(in) :: f,k
      intent(out) :: del6f
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx,3) :: del6f
      real, dimension (nx) :: tmp
      integer :: i,k,k1
      real, dimension (3) :: vec
!
      k1=k-1
      do i=1,3
        call del6fj(f,vec,k1+i,tmp)
        del6f(:,i)=tmp
      enddo
!
!  Exit if this is requested for non-cartesian runs.
!
      if (lcylindrical_coords.or.lspherical_coords) &
          call fatal_error('del2fjv', &
          'not implemented for non-cartesian coordinates')
!
    endsubroutine del6fjv
!***********************************************************************
    subroutine u_dot_grad_vec(f,k,gradf,uu,ugradf,upwind,ladd)
!
!  u.gradu
!
!  21-feb-07/axel+dhruba: added spherical coordinates
!   7-mar-07/wlad: added cylindrical coordinates
!  24-jun-08/MR: ladd added for incremental work
!
      use General, only: loptest
!
      intent(in) :: f,k,gradf,uu,upwind
      intent(out) :: ugradf
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx,3,3) :: gradf
      real, dimension (nx,3) :: uu,ff,ugradf,grad_f_tmp
      real, dimension (nx) :: tmp
      integer :: j,k
      logical, optional :: upwind,ladd
!
      if (k<1 .or. k>mfarray) then
        call fatal_error('u_dot_grad_vec','variable index is out of bounds')
        return
      endif
!
      do j=1,3
!
        grad_f_tmp = gradf(:,j,:)
        call u_dot_grad_scl(f,k+j-1,grad_f_tmp,uu,tmp,UPWIND=upwind)
        if (loptest(ladd)) then
          ugradf(:,j)=ugradf(:,j)+tmp
        else
          ugradf(:,j)=tmp
        endif
!
      enddo
!
!  Adjustments for spherical coordinate system.
!  The following now works for general u.gradA.
!
      if (lspherical_coords) then
        ff=f(l1:l2,m,n,k:k+2)
        ugradf(:,1)=ugradf(:,1)-r1_mn*(uu(:,2)*ff(:,2)+uu(:,3)*ff(:,3))
        ugradf(:,2)=ugradf(:,2)+r1_mn*(uu(:,2)*ff(:,1)-uu(:,3)*ff(:,3)*cotth(m))
        ugradf(:,3)=ugradf(:,3)+r1_mn*(uu(:,3)*ff(:,1)+uu(:,3)*ff(:,2)*cotth(m))
      endif
!
!  The following now works for general u.gradA.
!
      if (lcylindrical_coords) then
        ff=f(l1:l2,m,n,k:k+2)
        ugradf(:,1)=ugradf(:,1)-rcyl_mn1*(uu(:,2)*ff(:,2))
        ugradf(:,2)=ugradf(:,2)+rcyl_mn1*(uu(:,2)*ff(:,1))
      endif
!
    endsubroutine u_dot_grad_vec
!***********************************************************************
    subroutine u_dot_grad_vec_alt(f,k,gradf,uu,ugradf,iadvec,ladd)
!
!  u.gradu
!
!  21-feb-07/axel+dhruba: added spherical coordinates
!   7-mar-07/wlad: added cylindrical coordinates
!  24-jun-08/MR: ladd added for incremental work
!
      use General, only: loptest
!
      intent(in) :: f,k,gradf,uu,iadvec
      intent(out) :: ugradf
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx,3,3) :: gradf
      real, dimension (nx,3) :: uu,ff,ugradf,grad_f_tmp
      real, dimension (nx) :: tmp
      integer :: j,k,iadvec
      logical, optional :: ladd
      logical :: ladd1
!
      if (k<1 .or. k>mfarray) then
        call fatal_error('u_dot_grad_vec','variable index is out of bounds')
        return
      endif
!
      ladd1=loptest(ladd)
!
      do j=1,3
        grad_f_tmp = gradf(:,j,:)
        call u_dot_grad_scl_alt(f,k+j-1,grad_f_tmp,uu,tmp,iadvec)
        if (ladd1) then
          ugradf(:,j)=ugradf(:,j)+tmp
        else
          ugradf(:,j)=tmp
        endif
!
      enddo
!
!  Adjustments for spherical coordinate system.
!  The following now works for general u.gradA.
!
      if (lspherical_coords) then
        ff=f(l1:l2,m,n,k:k+2)
        ugradf(:,1)=ugradf(:,1)-r1_mn*(uu(:,2)*ff(:,2)+uu(:,3)*ff(:,3))
        ugradf(:,2)=ugradf(:,2)+r1_mn*(uu(:,2)*ff(:,1)-uu(:,3)*ff(:,3)*cotth(m))
        ugradf(:,3)=ugradf(:,3)+r1_mn*(uu(:,3)*ff(:,1)+uu(:,3)*ff(:,2)*cotth(m))
      endif
!
!  The following now works for general u.gradA.
!
      if (lcylindrical_coords) then
        ff=f(l1:l2,m,n,k:k+2)
        ugradf(:,1)=ugradf(:,1)-rcyl_mn1*(uu(:,2)*ff(:,2))
        ugradf(:,2)=ugradf(:,2)+rcyl_mn1*(uu(:,2)*ff(:,1))
      endif
!
    endsubroutine u_dot_grad_vec_alt
!***********************************************************************
    subroutine u_dot_grad_mat(f,k,gradM,uu,ugradM,upwind)
!
!  Computes  u.grad(M)
!  where M is a second rank matrix.
!
!  07-aug-10/dhruba: coded
! 24-nov-11/dhruba: added upwinding
! 26-mar-12/MR: doupwind introduced
!
      intent(in) :: gradM,f,k
      intent(out) :: ugradM
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx,3,3,3) :: gradM
      real,dimension(nx,3) :: uu
      real, dimension (nx,3,3) :: ugradM
      integer :: k
      logical, optional :: upwind
      integer :: ipi,ipj,ipk
!
      if (k<1 .or. k>mfarray) then
        call fatal_error('u_dot_grad_mat','variable index is out of bounds')
        return
      endif
!
      call vec_dot_3tensor(uu,gradM,ugradM)
!
!  Test if Upwind is used.
!
      if (present(upwind)) then
        if (upwind) then
!
! The same operation needs to be done to each element
! of the matrix ugradM. We assume below that this matrix
! is symmetric. Otherwise the following should be rewritten.
!
          ipk=0
          do ipi=1,3
            do ipj=ipi,3
              call doupwind(f,k+ipk,uu,ugradM(1,ipi,ipj))
              ipk=ipk+1
            enddo
          enddo
          call symmetrise3x3_ut2lt(ugradM)
        endif
      endif
!
!  Spherical and cylindrical coordinates are not
!  implemented for this subroutine.
!
      if (lspherical_coords) then
        call fatal_error('u_dot_grad_mat','not implemented in sph-coordinates')
      endif
!
      if (lcylindrical_coords) then
        call fatal_error('u_dot_grad_mat','not implemented in cyl-coordinates')
      endif
!
    endsubroutine u_dot_grad_mat
!***********************************************************************
    subroutine u_dot_grad_scl(f,k,gradf,uu,ugradf,upwind,ladd)
!
!  Do advection-type term u.grad f_k.
!  Assumes gradf to be known, but takes f and k as arguments to be able
!  to calculate upwind correction.
!
!  28-Aug-2007/dintrans: attempt of upwinding in cylindrical coordinates
!  29-Aug-2007/dhruba: attempt of upwinding in spherical coordinates.
!  28-Sep-2009/MR: ladd added for incremental work
!  26-mar-12/MR: doupwind introduced
!
      use General, only: loptest
!
      intent(in) :: f,k,gradf,uu,upwind,ladd
      intent(out) :: ugradf
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx,3) :: uu,gradf
      real, dimension (nx) :: ugradf
      integer :: k
      logical, optional :: upwind, ladd
!
      if (k<1 .or. k>mfarray) then
        call fatal_error('u_dot_grad_scl','variable index is out of bounds')
        return
      endif
!
      call dot_mn(uu,gradf,ugradf,ladd)
!
!  Upwind correction
!
      if (loptest(upwind)) call doupwind(f,k,uu,ugradf)
!
    endsubroutine u_dot_grad_scl
!***********************************************************************
    subroutine u_dot_grad_scl_alt(f,k,gradf,uu,ugradf,iadvec,ladd)
!
!  Do advection-type term u.grad f_k.
!  Assumes gradf to be known, but takes f and k as arguments to be able
!  to calculate upwind correction.
!
!  28-Aug-2007/dintrans: attempt of upwinding in cylindrical coordinates
!  29-Aug-2007/dhruba: attempt of upwinding in spherical coordinates.
!  28-Sep-2009/MR: ladd added for incremental work
!  22-Jun-2011/dhruba: made this alternative version which also incorporated the
! kurganov-tadmore scheme.
!  26-mar-12/MR: doupwind introduced
!
      use General, only: loptest
!
      intent(in) :: f,k,gradf,uu,iadvec,ladd
      intent(out) :: ugradf
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx,3) :: uu,gradf
      real, dimension (nx) :: ugradf, udelf
      logical, optional :: ladd
      integer :: k,iadvec
      logical :: ladd1
!
! iadvec=0 normal scheme
! iadvec=1 upwinding
! iadvec=2 Kurganov-Tadmor (KT)
!
      ladd1=loptest(ladd)
!
      select case (iadvec)
      case (0)
        call dot_mn(uu,gradf,ugradf,ladd1)
      case (1)
        call dot_mn(uu,gradf,ugradf,ladd1)
        call doupwind(f,k,uu,ugradf)
      case (2)
!
! x, y and z directions respectively
!
        call u_grad_kurganov_tadmore(f,k,udelf,1)
        ugradf=udelf
        call u_grad_kurganov_tadmore(f,k,udelf,2)
        ugradf=ugradf+udelf
        call u_grad_kurganov_tadmore(f,k,udelf,3)
        ugradf=ugradf+udelf
!
      case default
        if (lroot) print*, 'sub:u_dot_grad: iadvec must be 0,1 or 2, iadvec=', &
              iadvec
          call fatal_error('u_dot_grad_scl','')
        endselect
!
    endsubroutine u_dot_grad_scl_alt
!***********************************************************************
    subroutine u_grad_kurganov_tadmore(f,k,udelf,j)
!
      use Deriv, only: der2_minmod
!
      intent(in) :: f,k,j
      intent(out) :: udelf
!
      real, dimension (mx,my,mz,mfarray) :: f
      real,dimension(nx) :: udelf
      real, dimension (nx) :: vel,velpj,velmj,amhalf,aphalf,delfj,delfjp1,delfjm1
      integer :: k,j
      integer :: ix,iix
!
      vel=f(l1:l2,m,n,iuu+j-1)
      call der2_minmod(f,k,delfj,delfjp1,delfjm1,j)
      select case(j)
      case(1)
        do ix=l1,l2
          iix=ix-nghost
          velpj(iix) = f(ix+1,m,n,iux)
          velmj(iix) = f(ix-1,m,n,iux)
        enddo
      case(2)
        velpj = f(l1:l2,m+1,n,iuy)
        velmj = f(l1:l2,m-1,n,iuy)
      case(3)
        velpj = f(l1:l2,m,n+1,iuz)
        velmj = f(l1:l2,m,n-1,iuz)
      case default
        call fatal_error('sub:u_grad_kurganov_tadmore:','wrong component')
      endselect
      aphalf = abs(vel+velpj)/2.
      amhalf = abs(vel+velmj)/2.
      udelf = (aphalf-vel)*(delfjp1+delfj)/2. + &
          (amhalf-vel)*(delfjm1+delfj)/2.
!
    endsubroutine u_grad_kurganov_tadmore
!***********************************************************************
    subroutine nou_dot_grad_scl(gradf,uu,ugradf,del6u,upwind,ladd)
!
!  Do advection-type term u.grad f_k.
!  Assumes gradf to be known, but takes f and k as arguments to be able
!  to calculate upwind correction
!
!  28-Aug-2007/dintrans: attempt of upwinding in cylindrical coordinates
!  29-Aug-2007/dhruba: attempt of upwinding in spherical coordinates.
!  28-Sep-2009/MR: ladd added for incremental work
!
      use General, only: loptest
!
      intent(in) :: gradf,uu,upwind,ladd
      intent(out) :: ugradf
!
      real, dimension (nx,3) :: uu,gradf,del6u
      real, dimension (nx) :: ugradf
      logical, optional :: upwind,ladd
!
      call dot_mn(uu,gradf,ugradf,loptest(ladd))
!
!  Upwind correction (currently just for z-direction).
!
      if (present(upwind)) then; if (upwind) then
!
!  x-direction.
!
        ugradf=ugradf-abs(uu(:,1))*del6u(:,1)
!
!  y-direction.
!
        if (lcartesian_coords) then
          ugradf=ugradf-abs(uu(:,2))*del6u(:,2)
        else
          if (lcylindrical_coords) &
             ugradf=ugradf-rcyl_mn1*abs(uu(:,2))*del6u(:,2)
          if (lspherical_coords) &
             ugradf=ugradf-r1_mn*abs(uu(:,2))*del6u(:,2)
        endif
!
!  z-direction.
!
        if ((lcartesian_coords).or.(lcylindrical_coords)) then
          ugradf=ugradf-abs(uu(:,3))*del6u(:,3)
        else
          if (lspherical_coords) &
             ugradf=ugradf-r1_mn*sin1th(m)*abs(uu(:,3))*del6u(:,3)
        endif
!
      endif; endif
!
    endsubroutine nou_dot_grad_scl
!***********************************************************************
    subroutine h_dot_grad_vec(hh,gradf,ff,hgradf)
!
!  h.gradf for vectors h and f.
!
!  23-mar-08/axel: adapted from u_dot_grad_vec
!
      use General, only: loptest
!
      intent(in) :: hh,gradf,ff
      intent(out) :: hgradf
!
      real, dimension (nx,3,3) :: gradf
      real, dimension (nx,3) :: hh,ff,hgradf
      real, dimension (nx) :: tmp
      integer :: j
!
!  Dot product for each of the three components of gradf .
!
      do j=1,3
        call h_dot_grad_scl(hh,gradf(:,j,:),tmp)
        hgradf(:,j)=tmp
      enddo
!
!  Adjustments for spherical coordinate system.
!  The following now works for general u.gradA .
!
      if (lspherical_coords) then
        hgradf(:,1)=hgradf(:,1)-r1_mn*(hh(:,2)*ff(:,2)+hh(:,3)*ff(:,3))
        hgradf(:,2)=hgradf(:,2)+r1_mn*(hh(:,2)*ff(:,1)-hh(:,3)*ff(:,3)*cotth(m))
        hgradf(:,3)=hgradf(:,3)+r1_mn*(hh(:,3)*ff(:,1)+hh(:,3)*ff(:,2)*cotth(m))
      endif
!
!  The following now works for general u.gradA .
!
      if (lcylindrical_coords) then
        hgradf(:,1)=hgradf(:,1)-rcyl_mn1*(hh(:,2)*ff(:,2))
        hgradf(:,2)=hgradf(:,2)+rcyl_mn1*(hh(:,2)*ff(:,1))
      endif
!
    endsubroutine h_dot_grad_vec
!***********************************************************************
    subroutine h_dot_grad_scl(hh,gradf,hgradf)
!
!  Do advection-type term h.grad f_k, but h is not taken from f array.
!
!  23-mar-08/axel: adapted from u_dot_grad_scl
!
      intent(in) :: hh,gradf
      intent(out) :: hgradf
!
      real, dimension (nx,3) :: hh,gradf
      real, dimension (nx) :: hgradf
!
      call dot_mn(hh,gradf,hgradf)
!
    endsubroutine h_dot_grad_scl
!***********************************************************************
    subroutine gradf_upw1st(f,uu,k,gradf)
!
!  Do advection-type term u.grad f_k for upwind 1st order der scheme.
!
      use Deriv, only: der_upwind1st
!
      intent(in) :: f,uu
      intent(out) :: gradf
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx,3) :: uu,gradf
      integer :: j,k
!
      do j=1,3
        call der_upwind1st(f,uu,k,gradf(:,j),j)
      enddo
!
    endsubroutine gradf_upw1st
!***********************************************************************
    subroutine rdim(file,mx_in,my_in,mz_in,mvar_in,maux_in,mglobal_in,&
        prec_in,nghost_in,ipx_in, ipy_in, ipz_in)
!
!  Read dimension from file.
!
!  15-sep-09/nils: adapted from rdim
!
      character (len=*) :: file
      character         :: prec_in
      integer           :: mx_in,my_in,mz_in
      integer           :: mvar_in,maux_in,mglobal_in,nghost_in
      integer           :: ipx_in, ipy_in, ipz_in
!
!  Every processor writes to their procN/dim.dat (with io_dist.f90).
!
      open(124,file=file,FORM='formatted')
      read(124,*) mx_in,my_in,mz_in,mvar_in,maux_in,mglobal_in
      read(124,*) prec_in
      read(124,*) nghost_in, nghost_in, nghost_in
      read(124,*) ipx_in, ipy_in, ipz_in
!
      close(124)
!
    endsubroutine rdim
!***********************************************************************
    subroutine read_snaptime(file,tout,nout,dtout,t)
!
!  Read in output time for next snapshot (or similar) from control file.
!
!  30-sep-97/axel: coded
!  24-aug-99/axel: allow for logarithmic spacing
!   9-sep-01/axel: adapted for MPI
!  10-sep-15/MR  : tout set to t if file is missing and dtout>0
!
      use Mpicomm, only: mpibcast_real, MPI_COMM_WORLD
!
      character (len=*), intent(in) :: file
      real, intent(out) :: tout
      integer, intent(out) :: nout
      real, intent(in) :: dtout
      double precision, intent(in) :: t
!
      integer, parameter :: lun = 31
      logical :: exist
      integer, parameter :: nbcast_array=2
      real, dimension(nbcast_array) :: bcast_array
      double precision :: t0
!
      if (lroot) then
!
!  Depending on whether or not file exists, we need to
!  either read or write tout and nout from or to the file.
!
        inquire(FILE=trim(file),EXIST=exist)
        open(lun,FILE=trim(file))
        if (exist) then
          read(lun,*) tout,nout
        else
!
!  Special treatment when dtout is negative.
!  Now tout and nout refer to the next snapshopt to be written.
!
          settout: if (dtout < 0.0) then
            tout = log10(t)
          elseif (dtout /= 0.0) then settout
            !  make sure the tout is a good time
            t0 = max(t - dt, 0.0D0)
            tout = t0 + (dble(dtout) - modulo(t0, dble(dtout)))
          else settout
            call warning("read_snaptime", "Writing snapshot every time step. ")
            tout = 0.0
          endif settout
          nout=1
          write(lun,*) tout,nout
        endif
        close(lun)
!
!  Broadcast tout and nout in one go.
!
        bcast_array(1) = tout
        bcast_array(2) = nout
      endif
!
      call mpibcast_real(bcast_array,nbcast_array,comm=MPI_COMM_WORLD)
      tout = bcast_array(1)
      nout = bcast_array(2)
!
    endsubroutine read_snaptime
!***********************************************************************
    subroutine update_snaptime(file,tout,nout,dtout,t,lout,ch,nowrite)
!
!  Check whether we need to write snapshot; if so, update the snapshot
!  file (e.g. tsnap.dat). Done by all processors.
!
!  30-sep-97/axel: coded
!  24-aug-99/axel: allow for logarithmic spacing
!  27-jul-15/MR  : try to fix a strange behavior with gfortran:
!                  when crashing, a big number of unmotivated snapshots
!                  is output -> test of NaN in t
!
      use General, only: itoa, notanumber_0d
!
      character (len=*), intent(in) :: file
      real, intent(inout) :: tout
      integer, intent(inout) :: nout
      real, intent(in) :: dtout
      double precision, intent(in) :: t
      logical, intent(inout) :: lout
      logical, intent(in), optional :: nowrite
      character (len=intlen), intent(out), optional :: ch
!
      integer, parameter :: lun = 31
      logical :: lwrite
      real :: t_sp   ! t in single precision for backwards compatibility
      logical, save :: lfirstcall=.true.
      real, save :: deltat_threshold
!
      if (notanumber_0d(t)) then
        lout=.false.
        return
      endif
!
!  Use t_sp as a shorthand for either t or lg(t).
!
      if (dtout<0.0) then
        t_sp=log10(t)
      else
        t_sp=t
      endif
!
!  Check if no writing tout is requested.
!
      lwrite = .true.
      if (present(nowrite)) lwrite = .not. nowrite
!
!  Generate a running file number, if requested.
!
      if (present (ch)) ch = itoa (nout)
!
!  Mark lout=.true. when time has exceeded the value of tout do while loop to
!  make sure tt is always larger than tout.
!  (otherwise slices are written just to catch up with tt.)
!
!  WL: Add possibility that there should be a small threshold in this
!      comparison. Needed for outputting at the exact tsnap, otherwise
!      a difference between tsp and tout to machine precision can be
!      interpreted as stating that the output is to be done at the next,
!      not the current, timestep.
!
      if (lfirstcall) then
        if (.not.loutput_varn_at_exact_tsnap) then
          deltat_threshold=0.0
        else
          deltat_threshold=dtmin
        endif
        lfirstcall=.false.
      endif
!
      if ((t_sp >= tout) .or. &
!      if (lout.or.t_sp    >= tout             .or. &
          (abs(t_sp-tout) <  deltat_threshold)) then
        tout=tout+abs(dtout)
!        if (.not.lout) tout=tout+abs(dtout)
        nout=nout+1
        lout=.true.
!
!  Write corresponding value of tout to file to make sure we have it, in case
!  the code craches. If the disk is full, however, we need to reset the values
!  manually.
!
        if (lroot .and. lwrite) then
          open(lun,FILE=trim(file))
          write(lun,*) tout,nout
          write(lun,*) 'This file is written automatically (routine'
          write(lun,*) 'update_snaptime in sub.f90). The values above give'
          write(lun,*) 'time and number of the *next* snapshot. These values'
          write(lun,*) 'are only read once in the beginning. You may adapt'
          write(lun,*) 'them by hand (eg after a crash).'
          close(lun)
        endif
      else
        lout=.false.
      endif
!
    endsubroutine update_snaptime
!***********************************************************************
    subroutine shift_dt(dt_)
!
!  Hack to make the code output the VARn files at EXACTLY the times
!  defined by dsnap, instead of slightly after it.
!
!  03-aug-11/wlad: coded
!
      use General, only: safe_character_assign
!
      real, intent(inout) :: dt_
      real, save :: tsnap
      integer, save :: nsnap
      character (len=fnlen) :: file
      logical, save :: lfirst_call=.true.
!
!  Read the output time defined by dsnap.
!
      if (lfirst_call) then
        call safe_character_assign(file,trim(datadir)//'/tsnap.dat')
        call read_snaptime(file,tsnap,nsnap,dsnap,t)
        lfirst_call=.false.
      endif
!
!  Adjust the time-step accordingly, so that the next timestepping
!  lands the simulation at the precise time defined by dsnap.
!
      if ((tsnap-t > dtmin).and.(t+dt_ > tsnap)) then
        dt_=tsnap-t
        lfirst_call=.true.
      endif
!
    endsubroutine shift_dt
!***********************************************************************
    subroutine set_dt(dt1_)

      use Mpicomm, only: mpiallreduce_max, MPI_COMM_WORLD

      real :: dt1_
      real :: dt1, dt1_local
      real, save :: dt1_last=0.0

      dt1_local=dt1_
      ! Timestep growth limiter
      if (ddt > 0.) dt1_local=max(dt1_local,dt1_last)
      call mpiallreduce_max(dt1_local,dt1,MPI_COMM_WORLD)
      dt=1.0/dt1
      if (loutput_varn_at_exact_tsnap) call shift_dt(dt)
      ! Timestep growth limiter
      if (ddt > 0.) dt1_last=dt1_local/ddt

    endsubroutine set_dt
!***********************************************************************
    subroutine vecout(lun,file,vv,thresh,nvec)
!
!  Write vectors to disc if their length exceeds thresh.
!
!  22-jul-03/axel: coded
!
      character (len=*) :: file
      real, dimension(nx,3) :: vv
      real, dimension(nx) :: v2
      real :: thresh,thresh2,dummy=0.
      integer :: l,lun,nvec
      real :: t_sp   ! t in single precision for backwards compatibility
!
      t_sp = t
!
!  Return if thresh=0 (default).
!
      if (thresh==0.) return
!
!  Open files when first data point.
!
      if (lfirstpoint) then
        open(lun,FILE=trim(file)//'.dat',form='unformatted',position='append')
        write(lun) 0,0,0,t_sp,dummy,dummy  !(marking first line)
        nvec=0
      endif
!
!  Write data.
!
      thresh2=thresh**2
      v2=vv(:,1)**2+vv(:,2)**2+vv(:,3)**2
      do l=1,nx
        if (v2(l)>=thresh2) then
          write(lun) l,m-nghost,n-nghost,vv(l,:)
          nvec=nvec+1
        endif
      enddo
!
!  Close file, and write number of vectors to a separate file.
!
      if (llastpoint) then
        close(lun)
        open(lun,FILE=trim(file)//'.num',position='append')
        write(lun,*) t_sp,nvec
        close(lun)
      endif
!
    endsubroutine vecout
!***********************************************************************
    subroutine despike(f,j,retval,factor)
!
!  Remove large spikes from
!
!  14-aug-06/tony: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension(nx) :: retval
      real, dimension (mx) :: tmp_penc
      real, dimension (mx) :: meanf
      real :: factor
      real, parameter :: t1 = 1./26.
      real, parameter :: t2 = 0.70710678/26.
      real, parameter :: t3 = 0.57735027/26.
      real, parameter :: t4 = 0.
      real, parameter, dimension (-1:1,-1:1,-1:1) :: interp3D = reshape(&
            (/ t3, t2, t3, &
               t2, t1, t2, &
               t3, t2, t3, &
               t2, t1, t2, &
               t1, t4, t1, &
               t2, t1, t2, &
               t3, t2, t3, &
               t2, t1, t2, &
               t3, t2, t3 /),&
            (/ 3,3,3 /))
      integer :: ii,jj,kk
      integer :: j
!
      meanf=0.
      if ((nxgrid/=1).and.(nygrid/=1).and.(nzgrid/=1)) then
        tmp_penc=f(:,m,n,j)
        do kk=-1,1
        do jj=-1,1
        do ii=-1,1
          if (ii/=0.or.jj/=0.or.kk/=0) &
          meanf(3:mx-2)=meanf(3:mx-2)+interp3D(ii,jj,kk)*(f(3+ii:mx-2+ii,m+jj,n+kk,j)-tmp_penc(3:mx-2))
        enddo
        enddo
        enddo
      else
        call fatal_error('shock_max3_pencil_interp', &
            'Tony got lazy and only implemented the 3D case')
      endif
!
!      factor1=1./factor
      retval=max(meanf(l1:l2)*factor,f(l1:l2,m,n,j))
!      retval=max(meanf(l1:l2)*factor,retval)
!
    endsubroutine despike
!***********************************************************************
    subroutine smooth_kernel(f,j,smth)
!
!  Smooth scalar field FF using predefined constant gaussian like kernel.
!
!  20-jul-06/tony: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension(nx) :: smth
      integer :: j,l
!
      do l=l1,l2
        smth(l-l1+1)=sum(smth_kernel*f(l-3:l+3,m-3:m+3,n-3:n+3,j))
      enddo
!
    endsubroutine smooth_kernel
!***********************************************************************
    subroutine smooth(f, ivar)
!
!  Smoothes the f-variable ivar with a polynomial kernel.  It assumes
!  that the boundary conditions for ivar have been applied, and the
!  ghost cells are not treated upon return.
!
!  23-jan-14/ccyang: coded.
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
      integer, intent(in) :: ivar
!
      real, dimension(-nghost:nghost,-nghost:nghost,-nghost:nghost) :: kernel = 0.0
      real, dimension(nx,ny,nz) :: work
      logical :: lfirstcall = .true.
      integer :: i1 = 0, i2 = 0, j1 = 0, j2 = 0, k1 = 0, k2 = 0
      integer :: i, j, k
!
!  Initialization at first call.
!
      init: if (lfirstcall) then
        call get_smooth_kernel(kernel)
        xdim: if (nxgrid > 1) then
          i1 = -nghost
          i2 = nghost
        endif xdim
        ydim: if (nygrid > 1) then
          j1 = -nghost
          j2 = nghost
        endif ydim
        zdim: if (nzgrid > 1) then
          k1 = -nghost
          k2 = nghost
        endif zdim
        lfirstcall = .false.
      endif init
!
!  Smooth the field (assuming boundary conditions for ivar has been applied).
!
      work = 0.0
      do k = k1, k2
        do j = j1, j2
          do i = i1, i2
            work = work + kernel(i,j,k) * f(l1+i:l2+i,m1+j:m2+j,n1+k:n2+k,ivar)
          enddo
        enddo
      enddo
      f(l1:l2,m1:m2,n1:n2,ivar) = work
!
    endsubroutine smooth
!***********************************************************************
    subroutine get_smooth_kernel(kernel)
!
!  Gets the smooth kernel used by subroutine smooth.
!
!  15-feb-14/ccyang: coded
!
      real, dimension(-nghost:nghost,-nghost:nghost,-nghost:nghost), intent(out) :: kernel
!
      integer, dimension(-nghost:nghost) :: b
      real, dimension(-nghost:nghost) :: weight
      integer :: k
!
      call binomial(b)
      weight = real(b)
      k = 2 * nghost + 1
      kernel = 0.0
      kernel(0,0,0) = 1.0
      if (nxgrid > 1) kernel(:,0,0) = kernel(0,0,0) * weight
      if (nygrid > 1) kernel(:,:,0) = spread(kernel(:,0,0), 2, k) * spread(weight, 1, k)
      if (nzgrid > 1) kernel = spread(kernel(:,:,0), 3, k) * spread(spread(weight, 1, k), 1, k)
      kernel = kernel / sum(kernel)
!
    endsubroutine get_smooth_kernel
!***********************************************************************
    subroutine binomial(b)
!
!  Finds the binomial coefficients.
!
!  23-jan-14/ccyang: coded.
!
      integer, dimension(:), intent(out) :: b
!
      integer :: n, i
!
      n = size(b) - 1
!
      b(1) = 1
      do i = 1, n - 1
        b(i+1) = b(i) * (n - i + 1) / i
      enddo
      b(n+1) = 1
!
    endsubroutine binomial
!***********************************************************************
    subroutine identify_bcs(varname_input,idx)
!
!  Print boundary conditions for scalar field.
!
!  19-jul-02/wolf: coded
!  29-may-04/axel: allowed variable name to be 8 chars long
!
      character (len=*) :: varname_input
      integer :: idx
!
      write(*,'(A,A10,",  x: <",A8,">, y: <",A8,">,  z: <",A8,">")') &
           'Bcs for ', varname_input, &
           trim(bcx(idx)), trim(bcy(idx)), trim(bcz(idx))
!
    endsubroutine identify_bcs
!***********************************************************************
    function noform(cname,lcomplex)
!
!  Given a string of the form `name(format)',
!  returns the name without format, fills empty space
!  of correct length (depending on format) with dashes.

!  For output as legend.dat and first line of time_series.dat.
!
!  22-jun-02/axel: coded
!  20-aug-13/MR: optional parameter lcomplex added
!  26-aug-13/MR: unnecessary p descriptors removed from cform
!
      use Cparam, only: max_col_width
      use General, only: loptest
!
      character (len=*) :: cname
      logical, optional :: lcomplex
      character (len=max_col_width) :: noform,cform,cnumber,dashes
      integer :: index_e,index_f,index_g,index_i,index_d,index_r,index1,index2
      integer :: iform0,iform1,iform2,length,number,number1,number2,io_code
!
      intent(in)  :: cname
!
!  Fill DASHES with, well, dashes.
!
      dashes = repeat('-', max_col_width)
!
!  Find position of left bracket to isolate format, cform.
!
      iform0=index(cname,' ')
      iform1=index(cname,'(')
      iform2=index(cname,')')
!
!  Set format; use default if not given.
!  Here we keep the parenthesis in cform.
!
      if (iform1>0) then
        cform=cname(iform1:iform2)
        length=iform1-1
      else
        cform='(e10.2)'
        length=iform0-1
      endif
!
!  Find length of formatted expression, examples: f10.2, e10.3, g12.1 .
!  index_1 is the position of the format type (f,e,g), and
!  index_d is the position of the dot.
!
      index_e=scan(cform,'eE')
      index_f=scan(cform,'fF')
      index_g=scan(cform,'gG')
      index_i=scan(cform,'iI')
      index_d=index(cform,'.')
      index_r=index(cform,')')
      index1=max(index_e,index_f,index_g,index_i)
      index2=index_d; if (index_d==0) index2=index_r
!
!  Calculate the length of the format and assemble expression for legend.
!
      cnumber=cform(index1+1:index2-1)
      read(cnumber,'(i4)',iostat=io_code) number
      if (io_code /= 0) then
!
!  Error while reading or end of file.
!
        print*,'noform: formatting problem'
        print*,'problematic cnumber= <',cnumber,'>'
        number=10
      endif
      if (loptest(lcomplex)) number = 2*number+4
      number1=max(0,(number-length)/2)
      number2=max(1,number-length-number1) ! at least one separating dash
!
!  Sanity check.
!
      if (number1+length+number2 > max_col_width) then
        call error("noform", &
                   "Increase max_col_width or sanitize print.in{,.double}")
      endif
!
      noform=dashes(1:number1)//cname(1:length)//dashes(1:number2)
!
    endfunction noform
!***********************************************************************
    function levi_civita(i,j,k)
!
!  Totally antisymmetric tensor.
!
!  20-jul-03/axel: coded
!
      real :: levi_civita
      integer :: i,j,k
!
      if ( &
        (i==1 .and. j==2 .and. k==3) .or. &
        (i==2 .and. j==3 .and. k==1) .or. &
        (i==3 .and. j==1 .and. k==2) ) then
        levi_civita=1.
      elseif ( &
        (i==3 .and. j==2 .and. k==1) .or. &
        (i==1 .and. j==3 .and. k==2) .or. &
        (i==2 .and. j==1 .and. k==3) ) then
        levi_civita=-1.
      else
        levi_civita=0.
      endif
!
    endfunction levi_civita
!***********************************************************************
    function kronecker_delta(i,j)
!
!  \delta_{ij} = 1 if i==j, 0 otherwise
!
!  28-oct-11/dhruba: coded
!
      real :: kronecker_delta
      integer :: i,j
!
      if (i==j)  then
        kronecker_delta = 1.
      else
        kronecker_delta=0.
      endif
!
    endfunction
!***********************************************************************
    function poly_1(coef, x)
!
!  Horner's scheme for polynomial evaluation.
!  Version for 1d array.
!
!  17-jan-02/wolf: coded
!
      real, dimension(:) :: coef
      real, dimension(:) :: x
      real, dimension(size(x,1)) :: poly_1
      integer :: Ncoef,i
!
      Ncoef = size(coef,1)
!
      poly_1 = coef(Ncoef)
      do i=Ncoef-1,1,-1
        poly_1 = poly_1*x+coef(i)
      enddo
!
    endfunction poly_1
!***********************************************************************
    function poly_0(coef, x)
!
!  Horner's scheme for polynomial evaluation.
!  Version for scalar.
!
!  17-jan-02/wolf: coded
!
      real, dimension(:) :: coef
      real :: x
      real :: poly_0
      integer :: Ncoef,i
!
      Ncoef = size(coef,1)
!
      poly_0 = coef(Ncoef)
      do i=Ncoef-1,1,-1
        poly_0 = poly_0*x+coef(i)
      enddo
!
    endfunction poly_0
!***********************************************************************
    function poly_3(coef, x)
!
!  Horner's scheme for polynomial evaluation.
!  Version for 3d array.
!
!  17-jan-02/wolf: coded
!
      real, dimension(:) :: coef
      real, dimension(:,:,:) :: x
      real, dimension(size(x,1),size(x,2),size(x,3)) :: poly_3
      integer :: Ncoef,i
!
      Ncoef = size(coef,1)
!
      poly_3 = coef(Ncoef)
      do i=Ncoef-1,1,-1
        poly_3 = poly_3*x+coef(i)
      enddo
!
    endfunction poly_3
!***********************************************************************
    subroutine lower_triangular_index(ij,i1,j1)
      integer,intent(out)::ij
      integer,intent(in) :: i1,j1
      integer :: ii,jj
      ii=i1;jj=j1
      if (i1.lt.j1) then
        ii=j1
        jj=i1
      endif
      ij=ii*(ii-1)/2 + jj
    endsubroutine lower_triangular_index
!***********************************************************************
    recursive function ylm(ell,emm,der) result (sph_har)
!
!  Spherical harmonic
!
!   24-nov-14/dhruba: copied from step
!   15-jun-17/MR: corrected  derivative
!
      real :: sph_har
      real :: theta,phi
      integer :: ell,emm
      real, optional :: der
!
      real :: cos2p,cost,sint,cosp,sinp
      integer :: aemm
      logical :: lother
!
! the one over pi, cosines and sines below may be pre-calculated
!
      lother=.false.
      cost=costh(m); sint=sinth(m); cosp=cosph(n); sinp=sinph(n)
      goto 1

      entry ylm_other(theta,phi,ell,emm,der) result (sph_har)

      lother=.true.
      cost=cos(theta); sint=sin(theta); cosp=cos(phi); sinp=sin(phi)
      
 1    aemm=iabs(emm)

      select case (ell)
        case (0)
          sph_har=(0.5)*sqrt(1./pi)
        case (1)
          select case(aemm)
            case (0)
              sph_har=(0.5)*sqrt(3./pi)*cost
            case (1) 
              sph_har=(0.5)*sqrt(3./(2*pi))*sint*cosp
              if (emm<0) sph_har = -sph_har       ! Condon-Shortley phase
            case default
              call fatal_error('sub:ylm','l=1 wrong m ')
            endselect
        case (2)
          if (aemm==2) cos2p=2*cosp*cosp-1
          select case(aemm)
            case (0)
              sph_har=(0.25)*sqrt(5./pi)*(3.*cost*cost-1.)
            case (1) 
              sph_har=-(0.5)*sqrt(15./(2*pi))*sint*cost*cosp
              if (emm<0) sph_har = -sph_har       ! Condon-Shortley phase
            case (2)
              sph_har=(0.25)*sqrt(15./(2*pi))*sint*sint*cos2p
            case default
              call fatal_error('sub:ylm','l=2 wrong m ')
            endselect
        case (9)
          select case(aemm)
            case (0)
              sph_har=sqrt(19./(4*pi))*(12155.*cost**9 - 25740.*cost**7 + 18018.*cost**5 - 4620.*cost**3 + 315.*cost)/128.
            case default
              call fatal_error('sub:ylm','l=9 wrong m ')
          endselect
        case (10)
          select case(aemm)
            case (0)
              sph_har=sqrt(21./(4*pi)) &
                      *(46189.*cost**10 - 109395.*cost**8 + 90090.*cost**6 - 30030.*cost**4 + 3465.*cost**2 - 63.)/256.
            case default
              call fatal_error('sub:ylm','l=10 wrong m ')
          endselect
        case default
          call fatal_error('sub:ylm','your ylm is not implemented')
      endselect

      if (present(der)) then
        der = ell*cost*sph_har/sint
        if (emm<ell) then
          if (lother) then
            der = der - sqrt((2.*ell+1.)*(ell-aemm)*(ell+aemm)/(2.*ell-1.))* &
                    ylm_other(theta,phi,ell-1,emm)/sint
          else
            der = der - sqrt((2.*ell+1.)*(ell-aemm)*(ell+aemm)/(2.*ell-1.))* &
                    ylm(ell-1,emm)/sint
          endif
        endif
      endif
!
    endfunction ylm
!***********************************************************************
    function step_scalar(x,x0,width)
!
!  Smooth unit step function centred at x0; implemented as tanh profile.
!
!   5-sep-08/dhruba: copied from step
!   9-nov-10/axel: no need to have the tini here
!
      real :: x
      real :: step_scalar
      real :: x0,width
!
!  check whether width is finite.
!
      if (width==0) call fatal_error('step_scalar','width must not be zero')
      step_scalar = 0.5*(1+tanh((x-x0)/width))
!
    endfunction step_scalar
!***********************************************************************
    function step_vector(x,x0,width)
!
!  Smooth unit step function centred at x0; implemented as tanh profile
!
!  23-jan-02/wolf: coded
!
      real, dimension(:) :: x
      real, dimension(size(x)) :: step_vector
      real :: x0,width
!
      step_vector = 0.5*(1+tanh((x-x0)/(width+tini)))
!
    endfunction step_vector
!***********************************************************************
    function der_step(x,x0,width)
!
!  Derivative of smooth unit STEP() function given above (i.e. a bump profile).
!  Adapt this if you change the STEP() profile, or you will run into
!  inconsistenies.
!  MR: perhaps to be merged with step
!
!  23-jan-02/wolf: coded
!
      real, dimension(:) :: x
      real, dimension(size(x,1)) :: der_step,arg
      real :: x0,width
!
!  Some argument gymnastics to avoid `floating overflow' for large
!  arguments.
!
      arg = abs((x-x0)/(width+tini))
      arg = min(arg,8.)         ! cosh^2(8) = 3e+27
      der_step = 0.5/(width*cosh(arg)**2)
!
      endfunction der_step
!***********************************************************************
    function der6_step(x,x0,width)
!
!  6th order derivative of smooth unit STEP() function given
!  above (i.e. a bump profile).
!  Adapt this if you change the STEP() profile, or you will run into
!  inconsistenies.
!
!  08-dec-09/dhruba: aped from der_step
!
      real, dimension(:) :: x
      real, dimension(size(x,1)) :: der6_step,arg,sechx,tanhx
      real :: x0,width
!
!  Some argument gymnastics to avoid `floating overflow' for large
!  arguments.
!
      arg = abs((x-x0)/(width+tini))
      tanhx=tanh(arg)
      arg = min(arg,8.)         ! cosh^2(8) = 3e+27
      sechx=1./cosh(arg)
      der6_step = (1./(2*width**6))*(&
             -272.0*(sechx**6)*tanhx+416.0*(sechx**4)*(tanhx**3) &
             -32.0*(sechx**2)*(tanhx**5) )
!
      endfunction der6_step
!***********************************************************************
    function stepdown(x,x0,width)
!
!  Smooth unit step function centred at x0; implemented as tanh profile
!
!  23-jan-02/wolf: coded
!
      real, dimension(:) :: x
      real, dimension(size(x,1)) :: stepdown
      real :: x0,width
!
      stepdown = -0.5*(1+tanh((x-x0)/(width+tini)))
!
    endfunction stepdown
!***********************************************************************
    function der_stepdown(x,x0,width)
!
!  Derivative of smooth unit STEPDOWN() function given above (i.e. a bump profile).
!  Adapt this if you change the STEP() profile, or you will run into
!  inconsistenies.
!
!  27-mar-10/dhruba: aped from der_step
!
      real, dimension(:) :: x
      real, dimension(size(x,1)) :: der_stepdown,arg
      real :: x0,width
!
!  Some argument gymnastics to avoid `floating overflow' for large
!  arguments.
!
      arg = abs((x-x0)/(width+tini))
      arg = min(arg,8.)         ! cosh^2(8) = 3e+27
      der_stepdown = -0.5/(width*cosh(arg)**2)
!
    endfunction der_stepdown
!***********************************************************************
    function cubic_step_pt(x,x0,width,shift)
!
!  Smooth unit step function with cubic (smooth) transition over [x0-w,x0+w].
!  Optional argument SHIFT shifts center:
!  for shift=1. the interval is [x0    ,x0+2*w],
!  for shift=-1. it is          [x0-2*w,x0    ].
!  This is to make sure the interior region is not affected.
!  Maximum slope is 3/2=1.5 times that of a linear profile.
!
!  This version is for scalar args.
!
!  18-apr-04/wolf: coded
!
      real :: x
      real :: cubic_step_pt,xi
      real :: x0,width
      real, optional :: shift
      real :: relshift
!
      if (present(shift)) then; relshift=shift; else; relshift=0.0; endif
      xi = (x-x0)/(width+tini) - relshift
      xi = max(xi,-1.0)
      xi = min(xi, 1.0)
      cubic_step_pt = 0.5 + xi*(0.75-xi**2*0.25)
!
    endfunction cubic_step_pt
!***********************************************************************
    function cubic_step_mn(x,x0,width,shift)
!
!  Smooth unit step function with cubic (smooth) transition over [x0-w,x0+w].
!  Version for 1d arg (in particular pencils).
!
!  18-apr-04/wolf: coded
!
      real, dimension(:) :: x
      real, dimension(size(x,1)) :: cubic_step_mn,xi
      real :: x0,width
      real, optional :: shift
      real :: relshift
!
      if (present(shift)) then; relshift=shift; else; relshift=0.0; endif
      xi = (x-x0)/(width+tini) - relshift
      xi = max(xi,-1.0)
      xi = min(xi, 1.0)
      cubic_step_mn = 0.5 + xi*(0.75-xi**2*0.25)
!
    endfunction cubic_step_mn
!***********************************************************************
    function cubic_der_step_pt(x,x0,width,shift)
!
!  Derivative of smooth unit step function, localized to [x0-w,x0+w].
!  This version is for scalar args.
!
!  12-jul-05/axel: adapted from cubic_step_pt
!
      real :: x
      real :: cubic_der_step_pt,xi
      real :: x0,width
      real, optional :: shift
      real :: relshift,width1
!
      if (present(shift)) then; relshift=shift; else; relshift=0.0; endif
      width1 = 1./(width+tini)
      xi = (x-x0)*width1 - relshift
      xi = max(xi,-1.0)
      xi = min(xi, 1.0)
      cubic_der_step_pt = (0.75-xi**2*0.75) * width1
!
    endfunction cubic_der_step_pt
!***********************************************************************
    function cubic_der_step_mn(x,x0,width,shift)
!
!  Derivative of smooth unit step function, localized to [x0-w,x0+w].
!  Version for 1d arg (in particular pencils).
!
!  12-jul-05/axel: adapted from cubic_step_mn
!
      real, dimension(:) :: x
      real, dimension(size(x,1)) :: cubic_der_step_mn,xi
      real :: x0,width
      real, optional :: shift
      real :: relshift,width1
!
      if (present(shift)) then; relshift=shift; else; relshift=0.0; endif
      width1 = 1./(width+tini)
      xi = (x-x0)*width1 - relshift
      xi = max(xi,-1.0)
      xi = min(xi, 1.0)
      cubic_der_step_mn = (0.75-xi**2*0.75) * width1
!
    endfunction cubic_der_step_mn
!***********************************************************************
    function quintic_step_pt(x,x0,width,shift)
!
!  Smooth unit step function with quintic (smooth) transition over [x0-w,x0+w].
!  Optional argument SHIFT shifts center:
!  for shift=1. the interval is [x0    ,x0+2*w],
!  for shift=-1. it is          [x0-2*w,x0    ].
!  Maximum slope is 15/8=1.875 times that of a linear profile.
!
!  This version is for scalar args.
!
!  09-aug-05/wolf: coded
!
      real :: x
      real :: quintic_step_pt,xi
      real :: x0,width
      real, optional :: shift
      real :: relshift
!
      if (present(shift)) then; relshift=shift; else; relshift=0.0; endif
      xi = (x-x0)/(width+tini) - relshift
      xi = max(xi,-1.0)
      xi = min(xi, 1.0)
      quintic_step_pt = 0.5 + xi*(0.9375 + xi**2*(-0.625 + xi**2*0.1875))
!
    endfunction quintic_step_pt
!***********************************************************************
    function quintic_step_mn(x,x0,width,shift)
!
!  Smooth unit step function with quintic (smooth) transition over [x0-w,x0+w].
!
!  Version for 1d arg (in particular pencils).
!
!  09-aug-05/wolf: coded
!
      real, dimension(:) :: x
      real, dimension(size(x,1)) :: quintic_step_mn,xi
      real :: x0,width
      real, optional :: shift
      real :: relshift
!
      if (present(shift)) then; relshift=shift; else; relshift=0.0; endif
      xi = (x-x0)/(width+tini) - relshift
      xi = max(xi,-1.0)
      xi = min(xi, 1.0)
      quintic_step_mn = 0.5 + xi*(0.9375 + xi**2*(-0.625 + xi**2*0.1875))
!
    endfunction quintic_step_mn
!***********************************************************************
    function quintic_der_step_pt(x,x0,width,shift)
!
!  Derivative of smooth unit step function, localized to [x0-w,x0+w].
!
!  This version is for scalar args.
!
!  09-aug-05/wolf: coded
!
      real :: x
      real :: quintic_der_step_pt,xi
      real :: x0,width
      real, optional :: shift
      real :: relshift,width1
!
      if (present(shift)) then; relshift=shift; else; relshift=0.0; endif
      width1 = 1./(width+tini)
      xi = (x-x0)*width1 - relshift
      xi = max(xi,-1.0)
      xi = min(xi, 1.0)
      quintic_der_step_pt = (0.9375 + xi**2*(-1.875 + xi**2*0.9375)) &
                              * width1
!
    endfunction quintic_der_step_pt
!***********************************************************************
    function quintic_der_step_mn(x,x0,width,shift)
!
!  Derivative of smooth unit step function, localized to [x0-w,x0+w].
!
!  Version for 1d arg (in particular pencils).
!
!  09-aug-05/wolf: coded
!
      real, dimension(:) :: x
      real, dimension(size(x,1)) :: quintic_der_step_mn,xi
      real :: x0,width
      real, optional :: shift
      real :: relshift,width1
!
      if (present(shift)) then; relshift=shift; else; relshift=0.0; endif
      width1 = 1./(width+tini)
      xi = (x-x0)*width1 - relshift
      xi = max(xi,-1.0)
      xi = min(xi, 1.0)
      quintic_der_step_mn = (0.9375 + xi**2*(-1.875 + xi**2*0.9375)) &
                              * width1
!
    endfunction quintic_der_step_mn
!***********************************************************************
    function sine_step_pt(x,x0,width,shift)
!
!  Smooth unit step function with sine (smooth) transition over [x0-w,x0+w].
!  Optional argument SHIFT shifts center:
!  for shift=1. the interval is [x0    ,x0+2*w],
!  for shift=-1. it is          [x0-2*w,x0    ].
!  Maximum slope is 15/8=1.875 times that of a linear profile.
!
!  This version is for scalar args.
!
!  13-jun-06/tobi: Adapted from cubic_step
!
      real :: x
      real :: sine_step_pt,xi
      real :: x0,width
      real, optional :: shift
      real :: relshift
!
      if (present(shift)) then; relshift=shift; else; relshift=0.0; endif
      xi = (x-x0)/(width+tini) - relshift
      xi = max(xi,-1.0)
      xi = min(xi, 1.0)
      sine_step_pt = 0.5*(1+sin(0.5*pi*xi))
!
    endfunction sine_step_pt
!***********************************************************************
    function sine_step_mn(x,x0,width,shift)
!
!  Smooth unit step function with sine (smooth) transition over [x0-w,x0+w].
!
!  Version for 1d arg (in particular pencils).
!
!  13-jun-06/tobi: Adapted from cubic_step
!
      real, dimension(:) :: x
      real, dimension(size(x,1)) :: sine_step_mn,xi
      real :: x0,width
      real, optional :: shift
      real :: relshift
!
      if (present(shift)) then; relshift=shift; else; relshift=0.0; endif
      xi = (x-x0)/(width+tini) - relshift
      xi = max(xi,-1.0)
      xi = min(xi, 1.0)
      sine_step_mn = 0.5*(1+sin(0.5*pi*xi))
!
    endfunction sine_step_mn
!***********************************************************************
    subroutine nan_inform(f,msg,region,int1,int2,int3,int4,lstop)
!
!  Check input array (f or df) for NaN, -Inf, Inf, and output location in
!  array.
!
!  30-apr-04/anders: coded
!  12-jun-04/anders: region or intervals supplied in call
!
      use General, only: notanumber
!
      real, dimension(:,:,:,:) :: f
      character (len=*) :: msg
      integer :: a,b,c,d,a1=1,a2=mx,b1=1,b2=my,c1=1,c2=mz,d1=1,d2=1
      integer, dimension(2), optional :: int1,int2,int3,int4
      character (len=*), optional :: region
      logical, optional :: lstop
!
!  Must set d2 according to whether f or df is considered.
!
      d2 = size(f,4)
!
!  Set intervals for different predescribed regions.
!
      if (present(region)) then
!
        select case (region)
        case ('f_array')
        case ('pencil')
          b1=m
          b2=m
          c1=n
          c2=n
        case ('default')
          call fatal_error('nan_inform','No such region')
        endselect
!
      endif
!
!  Overwrite with supplied intervals.
!
      if (present(int1)) then  ! x
        a1=int1(1)
        a2=int1(2)
      endif
!
      if (present(int2)) then  ! y
        b1=int2(1)
        b2=int2(2)
      endif
!
      if (present(int3)) then  ! z
        c1=int3(1)
        c2=int3(2)
      endif
!
      if (present(int4)) then  ! variable
        d1=int4(1)
        d2=int4(2)
      endif
!
!  Look for NaN and inf in resulting interval.
!
      do a=a1,a2; do b=b1,b2; do c=c1,c2; do d=d1,d2
        if (notanumber(f(a,b,c,d))) then
          print*,'nan_inform: NaN with message "', msg, &
              '" encountered in the variable ', varname(d)
          print*,'nan_inform: ', varname(d), ' = ', f(a,b,c,d)
          print*,'nan_inform: t, it, itsub   = ', t, it, itsub
          print*,'nan_inform: l, m, n, iproc = ', a, b, c, iproc_world
          print*,'----------------------------'
          if (present(lstop)) then
            if (lstop) call fatal_error('nan_stop','')
          endif
        endif
      enddo; enddo; enddo; enddo
!
    endsubroutine nan_inform
!***********************************************************************
    subroutine parse_bc(bc,bc12)
!
!  Parse boundary conditions, which may be in the form `a' (applies to
!  both `lower' and `upper' boundary) or `a:s' (use `a' for lower,
!  `s' for upper boundary.
!
!  24-jan-02/wolf: coded
!
      character (len=2*bclen+1), dimension(mcom) :: bc
      character (len=bclen), dimension(mcom,2) :: bc12
      integer :: j,isep
!
      intent(in) :: bc
      intent(out) :: bc12
!
      do j=1,mcom
        if (bc(j) == '') then ! will probably never happen due to default='p'
          if (lroot) print*, 'Empty boundary condition No. ', &
              j, 'in (x, y, or z)'
          !!call fatal_error('parse_bc','')
        endif
        isep = index(bc(j),':')
        if (isep > 0) then
          bc12(j,1) = bc(j)(1:isep-1)
          bc12(j,2) = bc(j)(isep+1:)
        else
          bc12(j,:) = bc(j)(1:bclen)
        endif
      enddo
!
    endsubroutine parse_bc
!***********************************************************************
    subroutine inverse_parse_bc(bc,bc12)
!
!  27-jun-12/joern+dhruba: coded
!
      character (len=2*bclen+1), dimension(mcom) :: bc
      character (len=bclen), dimension(mcom,2) :: bc12
      integer :: j
!
      intent(out) :: bc
      intent(in) :: bc12
!
      do j=1,mcom
        if ((bc12(j,1) == '') .or. (bc12(j,2) == '')) then
! will probably never happen due to default='p'
          if (lroot) print*, 'Empty boundary condition No. ', &
              j, 'in (x, y, or z)'
          call fatal_error('inverse_parse_bc','')
        endif
        bc(j)(1:bclen) = bc12(j,1)
        bc(j)(bclen+1:bclen+1) = ':'
        bc(j)(bclen+2:2*bclen+1)=bc12(j,2)
      enddo
!
    endsubroutine inverse_parse_bc
!***********************************************************************
    subroutine parse_bc_rad(bc,bc1,bc2)
!
!  Parse boundary conditions, which may be in the form `a' (applies to
!  both `lower' and `upper' boundary) or `a:s' (use `a' for lower,
!  `s' for upper boundary.
!
!   6-jul-03/axel: adapted from parse_bc
!
      character (len=2*bclen+1), dimension(3) :: bc
      character (len=bclen), dimension(3) :: bc1,bc2
      integer :: j,isep
!
      intent(in) :: bc
      intent(out) :: bc1,bc2
!
      do j=1,3
        if (bc(j) == '') then ! will probably never happen due to default='p'
          if (lroot) print*, 'Empty boundary condition No. ', &
                 j, 'in (x, y, or z)'
          call fatal_error('parse_bc','')
        endif
        isep = index(bc(j),':')
        if (isep > 0) then
          bc1(j) = bc(j)(1:isep-1)
          bc2(j) = bc(j)(isep+1:)
        else
          bc1(j) = bc(j)(1:bclen)
          bc2(j) = bc(j)(1:bclen)
        endif
      enddo
!
    endsubroutine parse_bc_rad
!***********************************************************************
    subroutine parse_bc_radg(bc,bc1,bc2)
!
!  Parse boundary conditions, which may be in the form `a' (applies to
!  both `lower' and `upper' boundary) or `a:s' (use `a' for lower,
!  `s' for upper boundary.
!
!   6-jul-03/axel: adapted from parse_bc
!
      character (len=2*bclen+1) :: bc
      character (len=bclen) :: bc1,bc2
      integer :: isep
!
      intent(in) :: bc
      intent(out) :: bc1,bc2
!
      if (bc == '') then
        if (lroot) print*, 'Empty boundary condition in (x, y, or z)'
        call fatal_error('parse_bc_radg','')
      endif
      isep = index(bc,':')
      if (isep > 0) then
        bc1 = bc(1:isep-1)
        bc2 = bc(isep+1:)
      else
        bc1 = bc(1:bclen)
        bc2 = bc(1:bclen)
      endif
!
    endsubroutine parse_bc_radg
!***********************************************************************
    subroutine parse_shell(strin,strout)
!
!  Parse string replacing all $XXXX sequences with appropriate
!  values from the environment. Return the parsed result in strout.
!
      use General, only: safe_character_assign
!
      character (len=*) :: strin, strout
      character (len=255) :: envname, chunk !, envvalue
      character (len=1) :: chr
      character (len=64), parameter :: envnamechars = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-'
      integer :: inptr, inlen, envstart, nameptr
!
      intent(in)    :: strin
      intent(inout)   :: strout
!
      inptr=1
      inlen=len(trim(strin))
      strout=''
!
dlrloop:do
        envstart =index(strin(inptr:inlen),'$')
        if (envstart <= 0) exit;
        chunk = trim(strin(inptr:envstart-1))
        if (envstart > inptr) call safe_character_assign(strout,trim(strout)//trim(chunk))
        inptr = envstart + 1;
        if (inptr > inlen) exit dlrloop
!
        nameptr = inptr
nameloop: do
          chr = trim(strin(nameptr:nameptr))
          if (index(envnamechars,chr) > 0) then
            nameptr=nameptr+1
          else
            exit nameloop
          endif
!
          if (nameptr > inlen) exit nameloop
        enddo nameloop
        if ((nameptr-1) >= inptr) then
         envname=trim(strin(inptr:nameptr-1))
! Commented pending a C replacement
!         call getenv(trim(envname),envvalue)
!         call safe_character_assign(strout,trim(strout)//trim(envvalue))
        endif
!
        inptr=nameptr
        if (inptr>inlen) exit dlrloop
!
      enddo dlrloop
!
      if (inptr <= inlen) then
        chunk = trim(strin(inptr:inlen))
        call safe_character_assign(strout,trim(strout)//trim(chunk))
      endif
!
    endsubroutine parse_shell
!***********************************************************************
    function control_file_exists(fname,delete)
!
!  Does the given control file exist in either ./ or ./runtime/ ?
!  If DELETE is true, delete the file after checking for existence.
!
!  26-jul-09/wolf: coded
!
      use File_io, only: parallel_file_exists
      use General, only: loptest
!
      logical :: control_file_exists
      character (len=*), intent(in) :: fname
      logical, optional, intent(in) :: delete
!
      control_file_exists = parallel_file_exists(trim(fname), loptest(delete))
      if (.not. control_file_exists) &
          control_file_exists = parallel_file_exists(trim("runtime/"//fname), loptest(delete))
!
    endfunction control_file_exists
!***********************************************************************
    function read_line_from_file(fname)
!
!  Read the first line from a file; return empty string if file is empty
!
!  4-oct-02/wolf: coded
!
      use File_io, only : file_exists
!
      character (len=linelen) :: read_line_from_file
      character (len=*) :: fname
!
      integer :: unit=1
      integer :: ierr=0
!
      read_line_from_file=char(0)
      if (.not. lroot) return
!
      if (file_exists(fname)) then
        open(unit,FILE=fname,IOSTAT=ierr)
        read(unit,'(A)',IOSTAT=ierr) read_line_from_file
        close(unit)
        if (ierr /= 0) read_line_from_file=char(0)
      endif
!
    endfunction read_line_from_file
!***********************************************************************
    subroutine get_nseed(nseed)
!
!  Get length of state of random number generator. The current seed can
!  be represented by nseed (4-byte) integers.
!  Different compilers have different lengths:
!    NAG: 1, Compaq: 2, Intel: 47, SGI: 64, NEC: 256
!
      use Mpicomm, only: lroot
      use General, only: random_seed_wrapper
!
      integer, intent(out) :: nseed
!
      call random_seed_wrapper(SIZE=nseed)
print*,'AXELsub nseed=',nseed
!
!  Test whether mseed is large enough for this machine.
!
      if (nseed > mseed) then
        if (lroot) print*, 'This machine requires mseed >= ', nseed, &
            ', but you have only ', mseed
        call fatal_error('get_nseed','Need to increase mseed')
      endif
!
    endsubroutine get_nseed
!***********************************************************************
    subroutine get_where(mask, indices, status)
!
!  Get the indices where mask is .true.
!
!  Note: The pointer argument indices will be reassociated with a newly
!      allocated array.  It is the user's responsibility to pass in a
!      disassociated pointer and deallocate it after use.
!
!  10-feb-15/ccyang: coded.
!
      logical, dimension(:), intent(in) :: mask
      integer, dimension(:), pointer :: indices
      integer, intent(out), optional :: status
!
      integer :: i, n, stat
!
!  Allocate the index array.
!
      allocate(indices(count(mask)), stat=stat)
      error: if (stat /= 0) then
        soft: if (present(status)) then
          status = stat
          return
        endif soft
        call fatal_error('get_where', 'unable to allocate the index array. ')
      endif error
!
!  Scan through the elements of mask.
!
      n = 0
      scan: do i = 1, size(mask)
        true: if (mask(i)) then
          n = n + 1
          indices(n) = i
        endif true
      enddo scan
      if (present(status)) status = 0
!
      endsubroutine get_where
!***********************************************************************
    subroutine write_dx_general(file,x00,y00,z00)
!
!  Write .general file for data explorer (aka DX).
!
!  04-oct-02/wolf: coded
!  08-oct-02/tony: use safe_character_assign() to detect string overflows
!  12-sep-13/mcnallcp: make the endianness detection reflect the platform
!
      use General, only: safe_character_append, date_time_string
!
      real :: x00,y00,z00
      character (len=*) :: file
      character (len=datelen) :: date
      character (len=linelen) :: field='',struct='',type='',dep=''
!
!  This is True for big endian, False of little endian
      logical, parameter :: bigendian = ichar(transfer(1,'a')) == 0
!
      call date_time_string(date)
!
!  Accumulate a few lines.
!
      if (lhydro    ) then
        call safe_character_append(field,  'uu, '       )
        call safe_character_append(struct, '3-vector, ' )
        call safe_character_append(type,   'float, '    )
        call safe_character_append(dep,    'positions, ')
      endif
      if (ldensity  ) then
        call safe_character_append(field,  'lnrho, '    )
        call safe_character_append(struct, 'scalar, '   )
        call safe_character_append(type,   'float, '    )
        call safe_character_append(dep,    'positions, ')
      endif
      if (lentropy  ) then
        call safe_character_append(field,  'ss, '       )
        call safe_character_append(struct, 'scalar, '   )
        call safe_character_append(type,   'float, '    )
        call safe_character_append(dep,    'positions, ')
      endif
      if (ltemperature .and. (.not. ltemperature_nolog) ) then
        call safe_character_append(field,  'lnTT, '       )
        call safe_character_append(struct, 'scalar, '   )
        call safe_character_append(type,   'float, '    )
        call safe_character_append(dep,    'positions, ')
      endif
      if (lmagnetic ) then
        call safe_character_append(field,  'aa, '       )
        call safe_character_append(struct, '3-vector, ' )
        call safe_character_append(type,   'float, '    )
        call safe_character_append(dep,    'positions, ')
      endif
      if (lheatflux ) then
        call safe_character_append(field,  'qq, '       )
        call safe_character_append(struct, '3-vector, ' )
        call safe_character_append(type,   'float, '    )
        call safe_character_append(dep,    'positions, ')
      endif
      if (lradiation) then
        call safe_character_append(field,  'e_rad, ff_rad, '       )
        call safe_character_append(struct, 'scalar, 3-vector, '    )
        call safe_character_append(type,   'float, float, '        )
        call safe_character_append(dep,    'positions, positions, ')
      endif
      if (lpscalar  ) then
        call safe_character_append(field,  'lncc, '     )
        call safe_character_append(struct, 'scalar, '   )
        call safe_character_append(type,   'float, '    )
        call safe_character_append(dep,    'positions, ')
      endif
!
!  Remove trailing comma.
!
      field  = field (1:len(trim(field ))-1)
      struct = struct(1:len(trim(struct))-1)
      type   = type  (1:len(trim(type  ))-1)
      dep    = dep   (1:len(trim(dep   ))-1)
!
!  Now write.
!
      open(1,FILE=file)
!
      write(1,'(A)'  ) '# Creator: The Pencil Code'
      write(1,'(A,A)') '# Date: ', trim(date)
      write(1,'(A,A)') 'file = ', trim(datadir)//'/proc0/var.dat'
      write(1,'(A,I4," x ",I4," x ",I4)') 'grid = ', mx, my, mz
      write(1,'(A,A)') '# NB: hardcoded assumption of ieee'
      if (bigendian) then
        write(1,'(A,A," ",A)') 'format = ', 'msb', 'ieee'
      else
        write(1,'(A,A," ",A)') 'format = ', 'lsb', 'ieee'
      endif
      write(1,'(A,A)') 'header = ', 'bytes 4'
      write(1,'(A,A)') 'interleaving = ', 'record'
      write(1,'(A,A)') 'majority = ', 'column'
      write(1,'(A,A)') 'field = ', trim(field)
      write(1,'(A,A)') 'structure = ', trim(struct)
      write(1,'(A,A)') 'type = ', trim(type)
      write(1,'(A,A)') 'dependency = ', trim(dep)
      write(1,'(A,A,6(", ",1PG12.4))') 'positions = ', &
           'regular, regular, regular', &
           x00, dx, y00, dy, z00, dz
      write(1,'(A)') ''
      write(1,'(A)') 'end'
!
      close(1)
!
    endsubroutine write_dx_general
!***********************************************************************
    subroutine write_zprof(fname,a)
!
!  Writes z-profile to a file.
!
!  10-jul-05/axel: coded
!
      use HDF5_IO, only: output_profile
!
      real, dimension(:), intent(in) :: a
      character (len=*), intent(in) :: fname

      if (size(a) == mz) then
        call output_profile(fname, z, a, 'z', lsave_name=.true., lhas_ghost=.true.)
      else
        call output_profile(fname, z(n1:n2), a, 'z', lsave_name=.true.)
      endif
 
    endsubroutine write_zprof
!***********************************************************************
    subroutine write_xprof(fname,a)
!
!  Writes x-profile to a file.
!
!  10-jul-05/axel: coded
!
      use HDF5_IO, only: output_profile
!
      real, dimension(:), intent(in) :: a
      character (len=*), intent(in) :: fname

      if (size(a) == mx) then
        call output_profile(fname, x, a, 'x', lsave_name=.true., lhas_ghost=.true.)
      else
        call output_profile(fname, x(l1:l2), a, 'x', lsave_name=.true.)
      endif
 
    endsubroutine write_xprof
!***********************************************************************
    subroutine remove_prof(type)
!
!  Remove profile file.
!
!  10-jul-05/axel: coded
!  05-Nov-2018/PABourdin: generalized to any direction
!
      use File_io, only: file_remove
      use General, only: safe_character_assign
!
      character, intent(in) :: type
!
      character (len=120) :: fname,listfile
      integer :: ierr, unit=2
!
      call file_remove(trim(directory)//'/'//'profile_'//type//'.h5')
!
      call safe_character_assign(listfile,trim(directory)//'/'//type//'/prof_list.dat')
!
!  Read list of file and remove them one by one.
!
      open(unit,file=listfile,status='old',iostat=ierr)
      if (ierr /= 0) return
      do while (it <= nt)
        read(unit,*,iostat=ierr) fname
        if (ierr /= 0) exit
        call file_remove(trim(directory)//'/'//type//'/prof_'//trim(fname)//'.dat')
      enddo
      close(unit)
!
!  Now delete this listfile altogether.
!
      call file_remove(listfile)
!
    endsubroutine remove_prof
!***********************************************************************
    subroutine blob(ampl,f,i,radius,xblob,yblob,zblob,radius_x)
!
!  Single blob.
!
!  27-jul-02/axel: coded
!
      integer :: i
      real, dimension (mx,my,mz,mfarray) :: f
      real, optional :: xblob,yblob,zblob
      real :: ampl,radius,x01=0.,y01=0.,z01=0.,fact, fact_x
      real, optional :: radius_x
!
!  Single  blob.
!
      if (present(xblob)) x01=xblob
      if (present(yblob)) y01=yblob
      if (present(zblob)) z01=zblob
      if (ampl==0) then
        if (lroot) print*,'ampl=0 in blob'
      else
        if (lroot.and.ip<14) print*,'blob: variable i,ampl=',i,ampl
        fact=1./radius**2
!
!  Possibility of elongated blob in the x direction
!
        if (present(radius_x)) then
          fact_x=1./radius_x**2
        else
          fact_x=fact
        endif
!
        f(:,:,:,i)=f(:,:,:,i)+ampl*( &
           spread(spread(exp(-fact_x*(x-x01)**2),2,my),3,mz) &
          *spread(spread(exp(-fact  *(y-y01)**2),1,mx),3,mz) &
          *spread(spread(exp(-fact  *(z-z01)**2),1,mx),2,my))
      endif
!
    endsubroutine blob
!***********************************************************************
    recursive function hypergeometric2F1(a,b,c,z,tol) result (hyp2F1)
!
!  DOCUMENT ME!!!!!!
!
      real, intent(in) :: a,b,c,z,tol
      real :: hyp2F1
      real :: fac
      integer :: n
!
      real :: aa,bb,cc
!
      aa=a; bb=b; cc=c
!
      fac=1
      hyp2F1=fac
      n=1
!
      if (z<=0.5) then
!
        do while (fac>tol)
          fac=fac*aa*bb*z/(cc*n)
          hyp2F1=hyp2F1+fac
          aa=aa+1
          bb=bb+1
          cc=cc+1
          n=n+1
        enddo
!
      else
!
        !!!!!!!! only valid for mu=-1 !!!!!!!!
        !hyp2F1=2*hypergeometric2F1(aa,bb,aa+bb-cc+1,1-z,tol)-sqrt(1-z)* &
               !2*hypergeometric2F1(cc-aa,cc-bb,cc-aa-bb+1,1-z,tol)
        hyp2F1=(gamma_function(cc)*gamma_function(cc-aa-bb))/ &
               (gamma_function(cc-aa)*gamma_function(cc-bb))* &
               hypergeometric2F1(aa,bb,aa+bb-cc+1,1-z,tol) &
              +(1-z)**(cc-aa-bb)* &
               (gamma_function(cc)*gamma_function(aa+bb-cc))/ &
               (gamma_function(aa)*gamma_function(bb))* &
               hypergeometric2F1(cc-aa,cc-bb,cc-aa-bb+1,1-z,tol)
!
      endif
!
    endfunction hypergeometric2F1
!***********************************************************************
    recursive function pi_function(x) result(pi_func)
!
!  Calculates the Pi-function using rational approximation.
!
!    Pi(x) = Gamma(x+1) = x!
!
!  Coefficients were determined using maple's minimax() function.
!
!  9-jun-04/tobi+wolf: coded
!
      real, intent(in) :: x
      real :: pi_func
      integer, parameter :: order=7
      real, dimension(order) :: coeff1,coeff2
      real :: enum,denom
      integer :: i
!
      coeff1=(/0.66761295020790986D00, &
               0.36946093910826145D00, &
               0.18669829780572704D00, &
               4.8801451277274492D-2, &
               1.36528684153155468D-2, &
               1.7488042503123817D-3, &
               3.6032044608268575D-4/)
!
      coeff2=(/0.66761295020791116D00, &
               0.754817592058897962D00, &
              -3.7915754844972276D-2, &
              -0.11379619871302534D00, &
               1.5035521280605477D-2, &
               3.1375176929984225D-3, &
              -5.5599617153443518D-4/)
!
      if (x>1) then
!
        pi_func=x*pi_function(x-1)
!
      elseif (x<0) then
!
        if (abs(x+1)<=epsilon(x)) then
          pi_func=pi_function(x+1)/epsilon(x)
        else
          pi_func=pi_function(x+1)/(x+1)
        endif
!
      else
!
        enum=coeff1(order)
        do i=order-1,1,-1
          enum=enum*x+coeff1(i)
        enddo
        denom=coeff2(order)
        do i=order-1,1,-1
          denom=denom*x+coeff2(i)
        enddo
        pi_func=enum/denom
!
      endif
!
    endfunction pi_function
!***********************************************************************
    function gamma_function(x)
!
!  Calculates the Gamma-function as
!
!    Gamma(x) = Pi(x-1)
!
!  9-jun-04/tobi+wolf: coded
!
      real, intent(in) :: x
      real :: gamma_function
!
      gamma_function=pi_function(x-1)
!
    endfunction gamma_function
!***********************************************************************
    subroutine tensor_diffusion_coef(gecr,ecr_ij,bij,bb,vKperp,vKpara,rhs,llog,gvKperp,gvKpara)
!
!  Calculates tensor diffusion with variable tensor (or constant tensor).
!  Calculates parts common to both variable and constant tensor first.
!
!  Note: ecr=lnecr in the below comment
!
!  Write diffusion tensor as K_ij = Kpara*ni*nj + (Kperp-Kpara)*del_ij.
!
!  vKperp*del2ecr + d_i(vKperp)d_i(ecr) + (vKpara-vKperp) d_i(n_i*n_j*d_j ecr)
!      + n_i*n_j*d_i(ecr)d_j(vKpara-vKperp)
!
!  = vKperp*del2ecr + gKperp.gecr + (vKpara-vKperp) (H.G + ni*nj*Gij)
!      + ni*nj*Gi*(vKpara_j - vKperp_j),
!  where H_i = (nj bij - 2 ni nj nk bk,j)/|b| and vKperp, vKpara are variable
!  diffusion coefficients.
!
!  Calculates (K.gecr).gecr
!  =  vKperp(gecr.gecr) + (vKpara-vKperp)*Gi(ni*nj*Gj)
!
!  Adds both parts into decr/dt.
!
!  10-oct-03/axel: adapted from pscalar
!  30-nov-03/snod: adapted from tensor_diff without variable diffusion
!  04-dec-03/snod: converted for evolution of lnecr (=ecr)
!   9-apr-04/axel: adapted for general purpose tensor diffusion
!  25-jun-05/bing:
!
      real, dimension (nx,3,3) :: ecr_ij,bij
      real, dimension (nx,3) :: gecr,bb,bunit,hhh,gvKperp1,gvKpara1,tmpv
      real, dimension (nx) :: abs_b,b1,del2ecr,gecr2,vKperp,vKpara
      real, dimension (nx) :: hhh2,quenchfactor,rhs,tmp,tmpi,tmpj,tmpk
      real :: limiter_tensordiff=3.
      integer :: i,j,k
      logical, optional :: llog
      real, optional, dimension (nx,3) :: gvKperp,gvKpara
!
      intent(in) :: bb,bij,gecr,ecr_ij
      intent(out) :: rhs
!
!  Calculate unit vector of bb.
!
!     call dot2_mn(bb,abs_b,PRECISE_SQRT=.true.)
      call dot2_mn(bb,abs_b,FAST_SQRT=.true.)
      b1=1./max(tini,abs_b)
      call multsv_mn(b1,bb,bunit)
!
!  Calculate first H_i.
!
      del2ecr=0.
      do i=1,3
        del2ecr=del2ecr+ecr_ij(:,i,i)
        hhh(:,i)=0.
        do j=1,3
          tmpj(:)=0.
          do k=1,3
            tmpj(:)=tmpj(:)-2.*bunit(:,k)*bij(:,k,j)
          enddo
          hhh(:,i)=hhh(:,i)+bunit(:,j)*(bij(:,i,j)+bunit(:,i)*tmpj(:))
        enddo
      enddo
      call multsv_mn(b1,hhh,tmpv)
!
!  Limit the length of H such that dxmin*H < 1, so we also multiply
!  by 1/sqrt(1.+dxmin^2*H^2).
!  And dot H with ecr gradient.
!
!     call dot2_mn(tmpv,hhh2,PRECISE_SQRT=.true.)
      call dot2_mn(tmpv,hhh2,FAST_SQRT=.true.)
      quenchfactor=1./max(1.,limiter_tensordiff*hhh2*dxmax)
      call multsv_mn(quenchfactor,tmpv,hhh)
      call dot_mn(hhh,gecr,tmp)
!
!  Dot Hessian matrix of ecr with bi*bj, and add into tmp.
!
      call multmv_mn(ecr_ij,bunit,hhh)
      call dot_mn(hhh,bunit,tmpj)
      tmp = tmp+tmpj
!
!  Calculate (Gi*ni)^2 needed for lnecr form; also add into tmp.
!
      gecr2=0.
      if (present(llog)) then
        call dot_mn(gecr,bunit,tmpi)
        tmp=tmp+tmpi**2
!
!  Calculate gecr2 - needed for lnecr form.
!
        call dot2_mn(gecr,gecr2)
      endif
!
!  If variable tensor, add extra terms and add result into decr/dt.
!
!  Set gvKpara, gvKperp.
!
     if (present(gvKpara)) then; gvKpara1=gvKpara; else; gvKpara1=0.; endif
     if (present(gvKperp)) then; gvKperp1=gvKperp; else; gvKperp1=0.; endif
!
!  Put d_i ecr d_i vKperp into tmpj.
!
      call dot_mn(gvKperp1,gecr,tmpj)
!
!  Nonuniform conductivities, add terms into tmpj.
!
      call dot(bunit,gvKpara1-gvKperp1,tmpi)
      call dot(bunit,gecr,tmpk)
      tmpj = tmpj+tmpi*tmpk
!
!  Calculate rhs.
!
      rhs=vKperp*(del2ecr+gecr2) + (vKpara-vKperp)*tmp + tmpj
!
    endsubroutine tensor_diffusion_coef
!***********************************************************************
    subroutine max_for_dt_nx_nx(f,maxf)
!
!  Like maxf = max(f,max), unless we have chosen to manipulate data
!  before taking the maximum value. Designed for calculation of time step,
!  where one may want to exclude certain regions, etc.
!
!  Would be nicer as an (assumed-size) array-valued function (as a plug-in
!  replacement for max), but this can be more than 2 times slower (NEC
!  SX-5, compared to about 15% slower with Intel F95) than a subroutine
!  call according to tests.
!
!  30-jan-04/wolf: coded
!
      real, dimension(nx) :: maxf,f
!
      intent(in)    :: f
      intent(inout) :: maxf
!
      maxf = max(f,maxf)
!
    endsubroutine max_for_dt_nx_nx
!***********************************************************************
    subroutine max_for_dt_1_nx(f,maxf)
!
!  Like max_for_dt_n_n, but with a different signature of argument shapes.
!
!  30-jan-04/wolf: coded
!
      real, dimension(nx) :: maxf
      real                :: f
!
      intent(in)    :: f
      intent(inout) :: maxf
!
      maxf = max(f,maxf)
!
    endsubroutine max_for_dt_1_nx
!***********************************************************************
    subroutine max_for_dt_1_1_1_nx(f1,f2,f3,maxf)
!
!  Like max_for_dt_n_n, but with a different signature of argument shapes.
!
!  30-jan-04/wolf: coded
!
      real, dimension(nx) :: maxf
      real                :: f1,f2,f3
!
      intent(in)    :: f1,f2,f3
      intent(inout) :: maxf
!
      maxf = max(f1,f2,f3,maxf)
!
    endsubroutine max_for_dt_1_1_1_nx
!***********************************************************************
    function pencil_multiply1(s,v)
!
!  The `*' operator may be extended through this function to allow
!  elementwise multiplication of a `pencil-scalar' with a `pencil-vector'.
!
!   6-Sep-05/tobi: coded
!
      real, dimension(nx), intent(in) :: s
      real, dimension(nx,3), intent(in) :: v
      real, dimension(nx,3) :: pencil_multiply1
!
      integer :: i
!
      do i=1,3; pencil_multiply1(:,i) = s(:) * v(:,i); enddo
!
    endfunction pencil_multiply1
!***********************************************************************
    function pencil_multiply2(v,s)
!
!  The `*' operator may be extended through this function to allow
!  elementwise multiplication of a `pencil-scalar' with a `pencil-vector'.
!
!   6-Sep-05/tobi: coded
!
      real, dimension(nx,3), intent(in) :: v
      real, dimension(nx), intent(in) :: s
      real, dimension(nx,3) :: pencil_multiply2
!
      integer :: i
!
      do i=1,3; pencil_multiply2(:,i) = v(:,i) * s(:); enddo
!
    endfunction pencil_multiply2
!***********************************************************************
    function pencil_add1(s,v)
!
!  The `+' operator may be extended through this function to allow
!  elementwise addition of a `pencil-scalar' to a `pencil-vector'.
!
!   6-Sep-05/tobi: coded
!
      real, dimension(nx), intent(in) :: s
      real, dimension(nx,3), intent(in) :: v
      real, dimension(nx,3) :: pencil_add1
!
      integer :: i
!
      do i=1,3; pencil_add1(:,i) = s(:) + v(:,i); enddo
!
    endfunction pencil_add1
!***********************************************************************
    function pencil_add2(v,s)
!
!  The `+' operator may be extended through this function to allow
!  elementwise addition of a `pencil-scalar' to a `pencil-vector'.
!
!   6-Sep-05/tobi: coded
!
      real, dimension(nx,3), intent(in) :: v
      real, dimension(nx), intent(in) :: s
      real, dimension(nx,3) :: pencil_add2
!
      integer :: i
!
      do i=1,3; pencil_add2(:,i) = v(:,i) + s(:); enddo
!
    endfunction pencil_add2
!***********************************************************************
    function pencil_divide1(s,v)
!
!  The `/' operator may be extended through this function to allow
!  elementwise division of a `pencil-scalar' by a `pencil-vector'.
!
!   6-Sep-05/tobi: coded
!
      real, dimension(nx), intent(in) :: s
      real, dimension(nx,3), intent(in) :: v
      real, dimension(nx,3) :: pencil_divide1
!
      integer :: i
!
      do i=1,3; pencil_divide1(:,i) = s(:) / v(:,i); enddo
!
    endfunction pencil_divide1
!***********************************************************************
    function pencil_divide2(v,s)
!
!  The `/' operator may be extended through this function to allow
!  elementwise division of a `pencil-vector' by a `pencil-scalar'.
!
!   6-Sep-05/tobi: coded
!
      real, dimension(nx,3), intent(in) :: v
      real, dimension(nx), intent(in) :: s
      real, dimension(nx,3) :: pencil_divide2
!
      integer :: i
!
      do i=1,3; pencil_divide2(:,i) = v(:,i) / s(:); enddo
!
    endfunction pencil_divide2
!***********************************************************************
    function pencil_subtract1(s,v)
!
!  The `-' operator may be extended through this function to allow
!  elementwise subtraction of a `pencil-vector' from a `pencil-scalar'.
!
!   6-Sep-05/tobi: coded
!
      real, dimension(nx), intent(in) :: s
      real, dimension(nx,3), intent(in) :: v
      real, dimension(nx,3) :: pencil_subtract1
!
      integer :: i
!
      do i=1,3; pencil_subtract1(:,i) = s(:) - v(:,i); enddo
!
    endfunction pencil_subtract1
!***********************************************************************
    function pencil_subtract2(v,s)
!
!  The `-' operator may be extended through this function to allow
!  elementwise subtraction of a `pencil-scalar' from a `pencil-vector'.
!
!   6-Sep-05/tobi: coded
!
      real, dimension(nx,3), intent(in) :: v
      real, dimension(nx), intent(in) :: s
      real, dimension(nx,3) :: pencil_subtract2
!
      integer :: i
!
      do i=1,3; pencil_subtract2(:,i) = v(:,i) - s(:); enddo
!
    endfunction pencil_subtract2
!***********************************************************************
    elemental real function one_minus_exp(x)
!
!  Returns 1 - exp(-x).
!
!  03-jan-15/ccyang: coded.
!
      real, intent(in) :: x
!
      if (x * x > epsilon(1.0)) then
        one_minus_exp = 1.0 - exp(-x)
      else
        one_minus_exp = x * (1.0 - 0.5 * x)
      endif
!
    endfunction one_minus_exp
!***********************************************************************
    elemental real function erfunc(x)
!
!  Error function from Numerical Recipes.
!
!  15-Jan-2007/dintrans: coded
!
      real, intent(in) :: x
!
      real :: dumerfc, t, z
!
      z = abs(x)
      t = 1.0 / ( 1.0 + 0.5 * z )
!
      dumerfc =  t * exp(-z * z - 1.26551223 + t *        &
          ( 1.00002368 + t * ( 0.37409196 + t *           &
          ( 0.09678418 + t * (-0.18628806 + t *           &
          ( 0.27886807 + t * (-1.13520398 + t *           &
          ( 1.48851587 + t * (-0.82215223 + t * 0.17087277 )))))))))
!
      if (x<0.0) dumerfc = 2.0 - dumerfc
      erfunc = 1.0 - dumerfc
!
    endfunction erfunc
!***********************************************************************
    subroutine power_law_mn(const,dist,power_law_index,output,xref)
!
!  General distance power law initial conditions.
!
!  24-feb-05/wlad: coded
!   4-jul-07/wlad: generalized for any power law case
!
      real, dimension(:) :: dist,output
      real :: const,power_law_index
      real, optional :: xref
!
      intent(in)  :: const,power_law_index
      intent(out) :: output
!
      if (present(xref)) dist=dist/xref
!
      if (rsmooth==0.) then
        output = const*dist**(-power_law_index)
      else
        output = const*(dist**2+rsmooth**2)**(-.5*power_law_index)
      endif
!
    endsubroutine power_law_mn
!***********************************************************************
    subroutine power_law_pt(const,dist,power_law_index,output,xref)
!
!  General distance power law initial conditions.
!
!  24-feb-05/wlad: coded
!   4-jul-07/wlad: generalized for any power law case
!
      real :: dist,output
      real :: const,power_law_index
      real, optional :: xref
!
      intent(in)  :: const,power_law_index
      intent(out) :: output
!
      if (present(xref)) dist=dist/xref
!
      if (rsmooth==0.) then
        output = const*dist**(-power_law_index)
      else
        output = const*(dist**2+rsmooth**2)**(-.5*power_law_index)
      endif
!
    endsubroutine power_law_pt
!***********************************************************************
    subroutine get_radial_distance(rrmn,rcylmn,e1_,e2_,e3_)
!
!  Calculate distance and its cylindrical projection for different
!  coordinate systems.
!
!  e1, e2, and e3 are the positions in the respective coordinate systems
!
!  15-mar-07/wlad : coded
!
      real, dimension(:),intent(out) :: rrmn,rcylmn
      real, dimension(size(rrmn,1)) :: xc
      real, intent(in), optional :: e1_,e2_,e3_
      real :: e1,e2,e3
      integer :: tmp
      logical :: lorigin
!
!  Check if we are dealing with distance from the origin.
!
      tmp=0 ; lorigin=.false.
      if (present(e1_)) then;e1=e1_;tmp=tmp+1;else;e1=0.;endif
      if (present(e2_)) then;e2=e2_;tmp=tmp+1;else;e2=0.;endif
      if (present(e3_)) then;e3=e3_;tmp=tmp+1;else;e3=0.;endif
      if (tmp==0) lorigin=.true.
!
!  Check if this array has size nx or mx.
!
!     select case (size(rrmn))
!     case (mx)
!       xc=x
!     case (nx)
!       xc=x(l1:l2)
!     case default
!       print*,'get_radial_distance: '//&
!            'the array has dimension=',size(rrmn),' is that correct?'
!       call fatal_error('get_radial_distance','')
!     endselect
!
!AB: the construct above doesn't work if there are no ghostzones,
!AB: i.e., if mx=nx. User therefore if statement.
!
      if (size(rrmn)==mx) then
        xc=x
      elseif (size(rrmn)==nx) then
        xc=x(l1:l2)
      else
        print*,'get_radial_distance: '//&
             'the array has dimension=',size(rrmn),' is that correct?'
        call fatal_error('get_radial_distance','')
      endif
!
!  Calculate the coordinate-free distance relative to the position (e1,e2,e3).
!
      if (lorigin) then
        if (coord_system=='cartesian') then
          rcylmn=sqrt(xc**2+y(m)**2)+tini
          rrmn  =sqrt(    rcylmn**2 +z(n)**2)
        elseif (coord_system=='cylindric') then
          rcylmn= xc                      +tini
          rrmn  =sqrt(  rcylmn**2+z(n)**2)
        elseif (coord_system=='spherical') then
          rcylmn=     xc*sinth(m)+tini
          rrmn  =     xc         +tini
        endif
      else
        if (coord_system=='cartesian') then
          rcylmn=sqrt((xc-e1)**2+(y(m)-e2)**2)+tini
          rrmn  =sqrt(       rcylmn**2+(z(n)-e3)**2)
        elseif (coord_system=='cylindric') then
          rcylmn=sqrt(xc**2+e1**2 - 2*xc*e1*cos(y(m)-e2))+tini
          rrmn  =sqrt(rcylmn**2+(z(n)-e3)**2)
        elseif (coord_system=='spherical') then
          rcylmn=sqrt((xc*sinth(m))**2 + (e1*sin(e2))**2 - &
               2*xc*e1*costh(m)*cos(e2))+tini
          rrmn  =sqrt(xc**2 + e1**2 - 2*xc*e1*&
               (costh(m)*cos(e2)+sinth(m)*sin(e2)*cos(z(n)-e3)))+tini
        endif
      endif
!
    endsubroutine get_radial_distance
!***********************************************************************
    function interp1(r,fr,nr,r0,ldescending)
!
!  20-dec-07/dintrans: coded
!  Note: if ldescending=T, then input arrays r and fr are in descending
!  order and we first reverse them.
!
      integer :: nr,istop,i,i1,i2
      real, dimension (nr) :: r,fr,r1,fr1
      real    :: r0,interp1
      logical, optional :: ldescending
!
      if (present(ldescending)) then
        if (ldescending) then
          r1=r(nr:1:-1)
          fr1=fr(nr:1:-1)
        else
          r1=r
          fr1=fr
        endif
      else
        r1=r
        fr1=fr
      endif
!
      if (r0 == r1(1)) then
        interp1=fr1(1)
        return
      elseif (r0 > r1(nr)) then
        interp1=fr1(nr)
        return
      else
        istop=0 ; i=1
        do while (istop /= 1)
          if (r1(i) >= r0) istop=1
          i=i+1
        enddo
        i1=i-2 ; i2=i-1
        interp1=(fr1(i1)*(r1(i2)-r0)+fr1(i2)*(r0-r1(i1)))/(r1(i2)-r1(i1))
      endif
!
    endfunction interp1
!***********************************************************************
    subroutine bspline_basis(k, x, b)
!
!  Computes the values of the non-zero B-spline basis functions
!  B_{i,k}(j+x) for i = j-k+1, j-k+2, ..., j.  The knot sequence {t_i)
!  is assumed to be infinite and be integers, i.e., t_i = i for all
!  integer i.
!
!  28-jul-15/ccyang: coded.
!
!  Input Arguments
!      k   Number of knot spans for each basis function, which has order
!          (k-1).
!      x   A number in [0,1).
!  Output Argument
!      b   An array of k elements, where b(i) = B_{j-k+i,k}(j+x).
!
      integer, intent(in) :: k
      real, intent(in) :: x
      real, dimension(k), intent(out) :: b
!
      integer :: i, j
!
!  Work up the order column by column.
!
      b = 0.0
      b(1) = 1.0
      order: do j = 2, k
        b(j) = x * b(j-1)
        do i = j - 1, 2, -1
          b(i) = (x - real(i-j)) * b(i-1) + (real(i) - x) * b(i)
        end do
        b(1) = (1.0 - x) * b(1)
        b(1:j) = b(1:j) / real(j-1)
      enddo order
!
    endsubroutine bspline_basis
!***********************************************************************
    subroutine bspline_interpolation(n, k, f, a, indx, shift)
!
!  Uses the B-spline interpolation to periodically shift a regular array
!  of data nodes.
!
!  31-jul-15/ccyang: coded.
!
!  Input/Output Argument
!      f   An array of node data; interpolated after shift on return.
!  Input Arguments
!      n   Number of nodes.
!      k   Number of knot spans for each basis function, which has order
!          (k-1).
!      a   Preconditioned by bspline_precondition() and then LU
!          decomposed by ludcmp().
!      indx
!          Index permutations returned by ludcmp().
!      shift
!          Shift in unit of array index.
!
      integer, intent(in) :: n, k
      real, dimension(n), intent(inout) :: f
      real, dimension(n,n), intent(in) :: a
      integer, dimension(n), intent(in) :: indx
      real, intent(in) :: shift
!
      real, dimension(:), allocatable, save :: bk
      integer :: k_old = -1
      real :: shift_old = 0.0
!
      real, dimension(n) :: b, c
      integer :: i, j
!
!  Solve for the coefficients for the B-spline basis functions.
!
      c = f
      call lubksb(a, indx, c)
!
!  Find the values of the basis functions at the interpolation points.
!
      j = ceiling(shift - 0.5)
      basis: if (k /= k_old .or. shift /= shift_old) then
        alloc: if (k /= k_old) then
          if (allocated(bk)) deallocate(bk)
          allocate(bk(k))
        endif alloc
        call bspline_basis(k, 0.5 - shift + real(j), bk)
        k_old = k
        shift_old = shift
      endif basis
      b = 0.0
      b(1:k) = bk
!
!  Make the interpolation.
!
      j = k + j
      do i = 1, n
        f(i) = sum(cshift(b, j - i) * c)
      enddo
!
    endsubroutine bspline_interpolation
!***********************************************************************
    subroutine bspline_precondition(n, k, a)
!
!  Sets up the linear system for the coefficients of the B-spline basis
!  functions, assuming periodic boundary conditions.  The linear system
!  reads
!
!      A x = f,
!
!  where
!          [ B_{1,k}(0.5)   B_{2,k}(0.5)   B_{3,k}(0.5)   ... B_{n,k}(0.5)   ]
!      A = [ B_{1,k}(1.5)   B_{2,k}(1.5)   B_{3,k}(1.5)   ... B_{n,k}(1.5)   ],
!          [ ...            ...            ...            ... ...            ]
!          [ B_{1,k}(n-0.5) B_{2,k}{n-0.5) B_{3,k}(n-0.5) ... B_{n,k}(n-0.5) ]
!
!      x = [ alpha_1, alpha_2, alpha_3, ..., alpha_n ]^T,
!
!  is the coefficients,
!
!      f = [ f_1, f_2, f_3, ..., f_n ],
!
!  is the node data.  The B-spline interpolation is then given by
!
!      f(x) = sum_i B_{i,k)(x)
!
!  with
!
!      f(i-0.5) = f_i, i = 1, 2, ..., n.
!
!  28-jul-15/ccyang: coded.
!
!  Input Arguments
!      n   Number of nodes.
!      k   Number of knot spans for each B-spline basis function, which
!          has order (k-1).
!  Output Argument
!      a   The square matrix for the linear system.
!
      integer, intent(in) :: n, k
      real, dimension(n,n), intent(out) :: a
!
      real, dimension(n) :: b
      integer :: i
!
!  Find the values of the B-spline basis functions.
!
      b = 0.0
      call bspline_basis(k, 0.5, b(1:k))
!
!  Cyclically assign the basis functions into the square matrix.
!
      do i = 1, n
        a(i,:) = cshift(b, k - i)
      enddo
!
    endsubroutine bspline_precondition
!***********************************************************************
    subroutine ludcmp(a,indx)
!
!  25-jun-09/rplasson: coded (adapted from numerical recipe)
!
!  Computes the LU decomposition of the matrix a.
!  The result is placed in the matrix a.
!  The row permutations are returned in indx.
!
      real, dimension(:,:), intent(INOUT) :: a
      integer, dimension(:), intent(OUT) :: indx
      real, dimension(size(a,1)) :: vv,swap
      integer :: j,n,imax
      integer, dimension(1) :: tmp
!
      n=size(a,1)
      if (n /= size(a,2)) call fatal_error('ludcmp','non square matrix')
      if (n /= size(indx)) call fatal_error('ludcmp','bad dimension for indx')
      vv=maxval(abs(a),dim=2)
      if (any(vv == 0.0)) call fatal_error('ludcmp','singular matrix')
      vv=1.0/vv
      do j=1,n
        tmp=maxloc(vv(j:n)*abs(a(j:n,j)))
        imax=(j-1)+tmp(1)
        if (j /= imax) then
          swap=a(imax,:)
          a(imax,:)=a(j,:)
          a(j,:)=swap
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if (a(j,j) == 0.0) a(j,j)=tiny(0.)
        a(j+1:n,j)=a(j+1:n,j)/a(j,j)
        a(j+1:n,j+1:n)=a(j+1:n,j+1:n)-spread(a(j+1:n,j),dim=2,ncopies=(n-j)) * &
            spread(a(j,j+1:n),dim=1,ncopies=(n-j))
      enddo
!
    endsubroutine ludcmp
!***********************************************************************
    subroutine lubksb(a,indx,b)
!
!  25-jun-09/rplasson: coded (adapted from numerical recipe)
!
!  Solves the equation A.X=B .
!  'a' must contain the LU decomposition of matrix A obtained by ludcmp.
!  'indx' is the permutation vector obtained by ludcmp.
!  'b' contains B, and returns the solution vector X.
!
      real, dimension(:,:), intent(IN) :: a
      integer, dimension(:), intent(IN) :: indx
      real, dimension(:), intent(INOUT) :: b
      integer :: i,n,ii,ll
      real :: summ
!
      n=size(a,1)
      if (n /= size(a,2)) call fatal_error('lubksb','non square matrix')
      if (n /= size(indx)) call fatal_error('lubksb','bad dimension for indx')
      ii=0
      do i=1,n
        ll=indx(i)
        summ=b(ll)
        b(ll)=b(i)
        if (ii /= 0) then
          summ=summ-dot_product(a(i,ii:i-1),b(ii:i-1))
        else if (summ /= 0.0) then
          ii=i
        endif
        b(i)=summ
      enddo
      do i=n,1,-1
        b(i) = (b(i)-dot_product(a(i,i+1:n),b(i+1:n)))/a(i,i)
      enddo
!
    endsubroutine lubksb
!***********************************************************************
    subroutine gij_psi(psif,ee,g)
!
!  Calculate gradient of a scalar field multiplied by a constant vector),
!  return matrix.
!
!  31-jul-07/dhruba: adapted from gij
!
      use Deriv, only: der
!
      real, dimension (mx,my,mz) :: psif
      real, dimension (3) :: ee
      real, dimension (nx,3,3) :: g
      real, dimension (nx) :: tmp
      integer :: i,j
!
      intent(in) :: psif
      intent(out) :: g
!
      do i=1,3; do j=1,3
          call der(psif*ee(i),tmp,j)
          g(:,i,j) = tmp
      enddo; enddo
!
    endsubroutine gij_psi
!***********************************************************************
    subroutine gij_psi_etc(psif,ee,aa,aij,Bij,del2,graddiv)
!
!  Calculate B_i,j = eps_ikl A_l,jk and A_l,kk .
!
!  1-aug-07/dhruba : adapted from gij_etc
!
      use Deriv, only: der2,derij
!
      real, dimension (mx,my,mz), intent (in) :: psif
      real, dimension(3), intent(in) :: ee
      real, dimension (nx,3,3), intent (out) :: bij
      real, dimension (nx,3,3), intent (in), optional :: aij
      real, dimension (nx,3), intent (out), optional :: del2,graddiv
      real, dimension (nx,3), intent (in), optional :: aa
!
!  Locally used variables.
!
      real, dimension (nx,3,3,3) :: d2A
      real, dimension (nx) :: tmp
      integer :: i,j
!
!  Calculate all mixed and non-mixed second derivatives
!  of the vector potential (A_k,ij).
!
!  Do not calculate both d^2 A/(dx dy) and d^2 A/(d dx)
!  (This wasn't spotted by me but by a guy from SGI...)
!  Note: for non-cartesian coordinates there are different correction terms,
!  see below.
!
      do i=1,3
        do j=1,3
          call der2(psif*ee(i),tmp,j); d2A(:,j,j,i)=tmp
        enddo
!!DHRUBA
        call derij(psif*ee(i),tmp,2,3); d2A(:,2,3,i)=tmp; d2A(:,3,2,i)=tmp
        call derij(psif*ee(i),tmp,3,1); d2A(:,3,1,i)=tmp; d2A(:,1,3,i)=tmp
        call derij(psif*ee(i),tmp,1,2); d2A(:,1,2,i)=tmp; d2A(:,2,1,i)=tmp
      enddo
!
!  Corrections for spherical polars from swapping mixed derivatives:
!  Psi_{,theta^ r^} = Psi_{,r^ theta^} - Psi_{,\theta^}/r
!  Psi_{,phi^ r^} = Psi_{,r^ phi^} - Psi_{,\phi^}/r
!  Psi_{,phi^ theta^} = Psi_{,theta^ phi^} - Psi_{,\phi^}*r^{-1}*cot(theta)
!
      if (lspherical_coords) then
        do i=1,3
          d2A(:,2,1,i)=d2A(:,2,1,i)-aij(:,i,2)*r1_mn
          d2A(:,3,1,i)=d2A(:,3,1,i)-aij(:,i,3)*r1_mn
          d2A(:,3,2,i)=d2A(:,3,2,i)-aij(:,i,3)*r1_mn*cotth(m)
        enddo
      endif
!
!  For cylindrical, only
!  Psi_{,phi^ pom^} = Psi_{,pom^ phi^} - Psi_{,\phi^}/pom .
!
      if (lcylindrical_coords) then
         do i=1,3
            d2A(:,2,1,i)=d2A(:,2,1,i)-aij(:,i,2)*rcyl_mn1
         enddo
      endif
!
!  Calculate b_i,j = eps_ikl A_l,kj, as well as optionally,
!  del2_i = A_i,jj and graddiv_i = A_j,ji .
!
      bij(:,1,:)=d2A(:,2,:,3)-d2A(:,3,:,2)
      bij(:,2,:)=d2A(:,3,:,1)-d2A(:,1,:,3)
      bij(:,3,:)=d2A(:,1,:,2)-d2A(:,2,:,1)
!
!  Corrections for spherical coordinates.
!
      if (lspherical_coords) then
        bij(:,3,2)=bij(:,3,2)+aij(:,2,2)*r1_mn
        bij(:,2,3)=bij(:,2,3)-aij(:,3,3)*r1_mn
        bij(:,1,3)=bij(:,1,3)+aij(:,3,3)*r1_mn*cotth(m)
        bij(:,3,1)=bij(:,3,1)+aij(:,2,1)*r1_mn         -aa(:,2)*r2_mn
        bij(:,2,1)=bij(:,2,1)-aij(:,3,1)*r1_mn         +aa(:,3)*r2_mn
        bij(:,1,2)=bij(:,1,2)+aij(:,3,2)*r1_mn*cotth(m)-aa(:,3)*r2_mn*sin2th(m)
      endif
!
!  Corrections for cylindrical coordinates.
!
      if (lcylindrical_coords) then
        bij(:,3,2)=bij(:,3,2)+ aij(:,2,2)*r1_mn
        bij(:,3,1)=bij(:,3,1)+(aij(:,2,1)+aij(:,1,2))*rcyl_mn1-aa(:,2)*rcyl_mn2
      endif
!
!  Calculate del2 and graddiv, if requested.
!
      if (present(graddiv)) then
!--     graddiv(:,:)=d2A(:,:,1,1)+d2A(:,:,2,2)+d2A(:,:,3,3)
        graddiv(:,:)=d2A(:,1,:,1)+d2A(:,2,:,2)+d2A(:,3,:,3)
        if (lspherical_coords) then
          graddiv(:,1)=graddiv(:,1)+aij(:,1,1)*r1_mn*2+ &
             aij(:,2,1)*r1_mn*cotth(m)-aa(:,2)*r2_mn*cotth(m)-aa(:,1)*r2_mn*2
          graddiv(:,2)=graddiv(:,2)+aij(:,1,2)*r1_mn*2+ &
             aij(:,2,2)*r1_mn*cotth(m)-aa(:,2)*r2_mn*sin2th(m)
          graddiv(:,3)=graddiv(:,3)+aij(:,1,3)*r1_mn*2+ &
             aij(:,2,3)*r1_mn*cotth(m)
        endif
      endif
!
      if (present(del2)) then
        del2(:,:)=d2A(:,1,1,:)+d2A(:,2,2,:)+d2A(:,3,3,:)
        if (lspherical_coords.and.present(aij).and.present(aa)) then
          del2(:,1)= del2(:,1)+&
            r1_mn*(2.*(aij(:,1,1)-aij(:,2,2)-aij(:,3,3)&
            -r1_mn*aa(:,1)-cotth(m)*r1_mn*aa(:,2) ) &
            +cotth(m)*aij(:,1,2) )
          del2(:,2)=del2(:,2)+&
            r1_mn*(2.*(aij(:,2,1)-cotth(m)*aij(:,3,3)&
            +aij(:,1,2) )&
            +cotth(m)*aij(:,2,2)-r1_mn*sin2th(m)*aa(:,2) )
          del2(:,3)=del2(:,3)+&
            r1_mn*(2.*(aij(:,3,1)+aij(:,1,3)&
            +cotth(m)*aij(:,2,3) ) &
            +cotth(m)*aij(:,3,2)-r1_mn*sin2th(m)*aa(:,3) )
        else
        endif
        if (lcylindrical_coords) call fatal_error('gij_etc', &
            'use del2=graddiv-curlcurl for cylindrical coords')
      endif
!
    endsubroutine gij_psi_etc
!***********************************************************************
    logical function location_in_proc(pos,lpos,mpos,npos)
!
! Finds out if a points defined by pos belongs to the subdomain of the processor.
! Note: pos contains only the relevant coordinates, that is, those of the *active* dimensions.
! If yes then also returns the nearest grid location of this point. (MR: "nearest" not exactly true)
!
!  -dec-10/dhruba: coded
!  28-dec-10/MR: changed into function
!  13-jan-11/MR: made dependent on dimensionality,
!                irrelevant position indices set to 1
!  17-jan-17/MR: made valid for all possible dimensionality cases;
!                set output for inactive dimensions to nghost+1
!
      real, dimension(dimensionality), intent(in)  :: pos
      integer,                         intent(out) :: lpos,mpos,npos
!
      logical :: linx,liny,linz
      integer :: i
!
      lpos=nghost+1 ; mpos=nghost+1 ; npos=nghost+1
      linx=.true. ; liny=.true. ; linz=.true.
!
      do i=1,dimensionality 
        select case (dim_mask(i))
          case (1); call xlocation(pos(i),lpos,linx); lpos=min(max(l1,lpos),l2)  !MR: minmax operation should be in *location
          case (2); call ylocation(pos(i),mpos,liny); mpos=min(max(m1,mpos),m2)
          case (3); call zlocation(pos(i),npos,linz); npos=min(max(n1,npos),n2)
        end select
      enddo
!
      location_in_proc = linx.and.liny.and.linz
!
    endfunction location_in_proc
!***********************************************************************
    subroutine xlocation(xpos,ixpos,lproc)
!
!  If xpos lies within this processor, then lproc=T and xpos=x(ixpos).
!  Otherwise lproc=F and ixpos=1.
!
!  18-nov-06/axel: coded
!  14-oct-08/ccyang: use half-closed interval and include the top-most plane
!  03-dec-09/MR: moved here from module slices
!  16-dec-10/joern: adapted from zlocation
!  18-jan-11/axel+joern: choose between l and l+1 depending on sign of x
!
      real :: xpos
      integer :: ixpos,l
      logical :: lproc
!
!  Use ixpos for testing if position is found.
!
      ixpos=-1
!
!  Run through all x positions until we hit the right interval.
!  If the right interval is found, jump out of the loop.
!  To make the output symmetric about zero, we must swap upper
!  and lower boundaries depending on the sign of x.
!
      do l=l1,l2
        if (x(l)<=xpos.and.x(l+1)>xpos) then
          if (x(l)<0) then
            ixpos=l+1
          else
            ixpos=l
          endif
          lproc=.true.
          exit
        endif
      enddo
!
!  If nothing is found, we set lproc=.false. and and put ixpos=1
!
      if (ixpos==-1) then
        ixpos=1
        lproc=.false.
      endif
!
    endsubroutine xlocation
!***********************************************************************
    subroutine ylocation(ypos,iypos,lproc)
!
!  If ypos lies within this processor, then lproc=T and ypos=y(iypos).
!  Otherwise lproc=F and iypos=1.
!
!  18-nov-06/axel: coded
!  14-oct-08/ccyang: use half-closed interval and include the top-most plane
!  03-dec-09/MR: adapted from module slices
!  16-dec-10/joern: took it from zlocation
!  18-jan-11/axel+joern: choose between m and m+1 depending on sign of y
!
      real :: ypos
      integer :: iypos,m
      logical :: lproc
!
!  Use iypos for testing if position is found.
!
      iypos=-1
!
!  Run through all y positions until we hit the right interval.
!  If the right interval is found, jump out of the loop.
!  To make the output symmetric about zero, we must swap upper
!  and lower boundaries depending on the sign of y.
!
      do m=m1,m2
        if (y(m)<=ypos.and.y(m+1)>ypos) then
          if (y(m)<0) then
            iypos=m+1
          else
            iypos=m
          endif
          lproc=.true.
          exit
        endif
      enddo
!
!  If nothing is found, we set lproc=.false. and and put iypos=1
!
      if (iypos==-1) then
        iypos=1
        lproc=.false.
      endif
!
    endsubroutine ylocation
!***********************************************************************
    subroutine zlocation(zpos,izpos,lproc)
!
!  If zpos lies within this processor, then lproc=T and zpos=z(izpos).
!  Otherwise lproc=F and izpos=1.
!
!  18-nov-06/axel: coded
!  14-oct-08/ccyang: use half-closed interval and include the top-most plane
!  03-dec-09/MR: moved here from module slices
!  18-jan-11/axel+joern: choose between n and n+1 depending on sign of z
!
      real :: zpos
      integer :: izpos,n
      logical :: lproc
!
!  Use izpos for testing if position is found.
!
      izpos=-1
!
!  Run through all z positions until we hit the right interval.
!  If the right interval is found, jump out of the loop.
!  To make the output symmetric about zero, we must swap upper
!  and lower boundaries depending on the sign of z.
!
      do n=n1,n2
        if (z(n)<=zpos.and.z(n+1)>zpos) then
          if (z(n)<0) then
            izpos=n+1
          else
            izpos=n
          endif
!AB: old  izpos=n
          lproc=.true.
          exit
        endif
      enddo
!
!  If nothing is found, we set lproc=.false. and and put izpos=1
!
      if (izpos==-1) then
        izpos=1
        lproc=.false.
      endif
!
    endsubroutine zlocation
!***********************************************************************
    subroutine position(ind,ip,ngrid,ind_loc,flag)
!
!  Determines local position ind_loc with respect to processor ip corresponding to global position ind
!  if grid has local extent ngrid. flag is set if ind_loc lies within the local range.
!  On return, ind_loc is corrected for number of ghost zones, thus can be used to index the f array.
!  ind_loc and flag are not altered if ind <= 0.
!
!  21-apr-15/MR: coded
!  27-may-18/MR: added ind_loc=-1 when ind outside proc range
!   6-may-19/MR: made +1 the "invalid" value to allow indexing with ind_loc
!
      integer, intent(inout):: ind
      integer, intent(in)   :: ip,ngrid
      integer, intent(inout):: ind_loc
      logical, intent(out)  :: flag
     
      if (ind>0) then
        ind_loc=ind-ip*ngrid
        if (ind_loc>=1.and.ind_loc<=ngrid) then
          flag=.true.
          ind_loc=ind_loc+nghost
        else
          flag=.false.
          ind_loc=1
        endif
      else
        if (ind_loc > nghost) ind = ind_loc + ip*ngrid
      endif

    endsubroutine position
!***********************************************************************
    subroutine fourier_single_mode(arr,idims,k,idir,amps,l2nd)
!
!  No parallelization in x allowed here
!
!  08-dec-10/MR: coded
!
      use mpicomm, only: mpireduce_sum
!
      implicit none
!
      integer, dimension(2)             , intent(in)  :: idims
      real, dimension(idims(1),idims(2)), intent(in)  :: arr
      real                              , intent(in)  :: k
      integer                           , intent(in)  :: idir
      real   , dimension(2,*)           , intent(out) :: amps
      logical                , optional , intent(in)  :: l2nd
!
      integer :: n,i,idim
      real, dimension(:)  , allocatable :: cg,sg
      real, dimension(:,:), allocatable :: buffer
      real :: fac
      logical :: l2ndl
!
      if (present(l2nd)) then
        l2ndl=l2nd
      else
        l2ndl=.false.
      endif
!
      if (l2ndl) then
        idim=idims(1)
      else
        idim=idims(2)
      endif
!
      select case (idir)
      case (1)    ; n=nxgrid
      case (2)    ; n=ny
      case (3)    ; n=nz
      case default; n=nxgrid
      end select
!
      if (idims(1)+idims(2)-idim/=n) then
        amps(:,1:idim)=0.
        return
      endif
!
      allocate(cg(n),sg(n))
!
      select case (idir)
      case (1)    ; cg=cos(k*xgrid)   ;sg=sin(k*xgrid)   ; fac=dx
      case (2)    ; cg=cos(k*y(m1:m2));sg=sin(k*y(m1:m2)); fac=dy
      case (3)    ; cg=cos(k*z(n1:n2));sg=sin(k*z(n1:n2)); fac=dz
      case default; cg=cos(k*xgrid)   ;sg=sin(k*xgrid)   ; fac=dx
      end select
!
      do i=1,idim
        if (l2ndl) then
          amps(:,i) = fac*(/ sum(arr(i,:)*cg), sum(arr(i,:)*sg) /)
        else
          amps(:,i) = fac*(/ sum(arr(:,i)*cg), sum(arr(:,i)*sg) /)
        endif
      enddo
!
      if (ncpus>1) then
        allocate(buffer(2,idim))
        select case (idir)
        case (2)
          if (nprocy>1) then
            call mpireduce_sum(amps,buffer,(/2,idim/),idir=2)
            if (ipy==0) amps(:,1:idim)=buffer      !result is in root of y-beams
          endif
        case (3)
          if (nprocz>1) then
            call mpireduce_sum(amps,buffer,(/2,idim/),idir=3)
            if (ipz==0) amps(:,1:idim)=buffer      !result is in root of z-beams
          endif
        end select
        deallocate(buffer)
      endif
!
      deallocate(cg,sg)
!
    endsubroutine fourier_single_mode
!***********************************************************************
    subroutine register_report_aux(name, index, ind_aux1, ind_aux2, ind_aux3, communicated)
!
!  Registers aux variable named 'name' if not already registered
!  (i.e. if index==0). Variable is scalar if ind_aux1,ind_aux2,
!  ind_aux3 are missing, vector with number of components equal to
!  number of present ind_aux* parameters. Index of variable and its
!  components (if any) are returned in index,ind_aux1,ind_aux2,ind_aux3
!  If already registered: outputs indices in index.pro
!
!  13-jan-11/MR: coded
!  29-may-14/ccyang: add optional argument communicated
!
      use FArrayManager, only: farray_register_auxiliary, farray_index_append
!
      implicit none
!
      integer,           intent(inout) :: index
      integer, optional, intent(inout) :: ind_aux1,ind_aux2,ind_aux3
      character (LEN=*), intent(in)    :: name
      logical, intent(in), optional :: communicated
!
      integer   :: vec
!
      vec=-1
!
      if ( present(ind_aux1) ) then
        vec=1
        if ( present(ind_aux2) ) then
          if ( present(ind_aux3) ) then
            vec=3
          else
            vec=2
          endif
        endif
      endif
!
      if (index == 0) then
        call farray_register_auxiliary(trim(name), index, vector=abs(vec), communicated=communicated)
        if (vec>=1) then
          ind_aux1=index
          if (vec>=2) then
            ind_aux2=index+1
            if (vec==3) ind_aux3=index+2
          endif
        endif
      else
        if (lroot) print*, 'register_report_aux: i'//trim(name)//' =', index
        call farray_index_append('i'//trim(name),index,vec)
      endif
!
    endsubroutine register_report_aux
!***********************************************************************
    subroutine unit_vector(bb,bb_hat)
!
!  Compute the unit vector for any given vector bb.
!  Tries to avoid division by zero.
!  Taken from http://nuclear.llnl.gov/CNP/apt/apt/aptvunb.html.
!
!  18-oct-11/bing: copied from bb_unitvec_shock in magnetic.f90
!
      real, dimension(nx,3) :: bb,bb_hat,bb2
      real, dimension(nx) :: a2,aerr2,bb_len
      integer :: j
      real :: tol
!
      intent(in) :: bb
      intent(out) :: bb_hat
!
!  Truncate small components to zero.
!
      tol=sqrt(tini)
!
      bb2 = bb**2
!
      aerr2 = tol * max(sum(bb2,2),1.)
!
      do j=1,3
        where (bb2(:,j) < aerr2)
          bb_hat(:,j) = 0.
        elsewhere
          bb_hat(:,j) = bb(:,j)
        endwhere
      enddo
!
!  Get unit vector.
!
      bb_len = sqrt(sum(bb_hat**2,2))
!
      do j=1,3; bb_hat(:,j) = bb_hat(:,j)/(bb_len+tini); enddo
!
!  Check if length is between 0. and 1.
!
      call dot2(bb_hat,a2)
!
      if (maxval(a2) > 1.+1e-6) &
          call fatal_error('unit_vector:','has not the length 1')
!
    endsubroutine unit_vector
!***********************************************************************
    subroutine calc_del6_for_upwind(f,k,hh,del6f_upwind,mask)
!
      use Deriv, only: der6, deri_3d_inds
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(nx,3)             :: del6f, hh
      real, dimension(nx)               :: del6f_upwind
      integer, dimension(nx)            :: indxs
      integer, intent(in), optional     :: mask
      integer                           :: j, k, msk
!
      intent(in) :: f,k,hh 
      intent(out) :: del6f_upwind
!
      msk=0
      if (present(mask)) msk=mask
!
      do j=1,3
!
        if (j==msk) then
          del6f(:,j) = 0.
        else
!
          if (lequidist(j) .or. lignore_nonequi) then
            call der6(f,k,del6f(1,j),j,UPWIND=.true.)
          else
            where(hh(:,j)>=0)
              indxs = 7
            elsewhere
              indxs = 8
            endwhere
            call deri_3d_inds(f(1,1,1,k),del6f(1,j),indxs,j,lnometric=.true.)
          endif
!
          del6f(:,j) = abs(hh(:,j))*del6f(:,j)
!
        endif
      enddo
!
      if (lcylindrical_coords) &
        del6f(:,2) = rcyl_mn1*del6f(:,2)
!
      if (lspherical_coords) then
        del6f(:,2) = r1_mn*del6f(:,2)
        del6f(:,3) = r1_mn*sin1th(m)*del6f(:,3)
      endif
!
      del6f_upwind = sum(del6f,2)
!
    endsubroutine calc_del6_for_upwind
!***********************************************************************    
    subroutine doupwind(f,k,uu,ugradf,mask)
!
!  Calculates upwind correction, works incrementally on ugradf
!
!  26-mar-12/MR: outsourced from routines u_dot_grad_mat, u_dot_grad_scl, u_dot_grad_scl_alt
!   9-apr-12/MR: optional parameter plus added
!  12-apr-12/MR: optional parameter modified
!   8-apr-17/wlyra: encapsulated the calculation of del6     
!
      real, dimension (mx,my,mz,mfarray), intent(IN)    :: f
      integer                                           :: k
      real, dimension (nx,3),             intent(IN)    :: uu
      real, dimension (nx),               intent(INOUT) :: ugradf
      integer,                            intent(IN), optional :: mask
      real, dimension (nx) :: del6f_upwind
      integer                :: msk
!
      msk=0
      if (present(mask)) then
        if ( mask>=1 .and. mask <=3 ) msk=mask
      endif
!
      call calc_del6_for_upwind(f,k,uu,del6f_upwind,msk)
!
      if (msk>0) then
        ugradf = ugradf+del6f_upwind
      else
        ugradf = ugradf-del6f_upwind
      endif
!
    endsubroutine doupwind
!***********************************************************************
    subroutine global_mean(f,inda,mean,indep,lexp)
!
!  Calculate global mean for a (several) field(s) selected by the index inda
!  (the index range inda - indep) in f.
!
!  15-oct-12/MR: adapted from remove_mean
!  22-aug-13/MR: rewritten into subroutine because of internal compiler error of gfortran
!  25-aug-13/MR: removed allocatable attribute from mean to adhere to f95
!
      use Mpicomm, only: mpiallreduce_sum
      use General, only: ioptest, loptest
!
      real, dimension (mx,my,mz,*), intent(in)           :: f
      integer,                      intent(in)           :: inda
      real, dimension(inda:),       intent(out)          :: mean
      integer,                      intent(in), optional :: indep
      logical,                      intent(in), optional :: lexp
!
      real, allocatable, dimension(:) :: mean_tmp
      integer :: j, inde
      logical :: lexpl
!
      inde = ioptest(indep,inda)
!
      allocate(mean_tmp(inda:inde))
!
!  initialize mean
!
      mean = 0.0
!
!  Compute mean for each field.
!
      lexpl = loptest(lexp)
      do j=inda,inde
        if (lexpl) then
          mean(j) = mean(j) + sum(exp(f(l1:l2,m1:m2,n1:n2,j)))
        else
          mean(j) = mean(j) + sum(f(l1:l2,m1:m2,n1:n2,j))
        endif
      enddo
!
!  Compute total mean for all processors
!
      call mpiallreduce_sum(mean,mean_tmp,inde-inda+1)
      mean = mean_tmp/nwgrid
!
    endsubroutine global_mean
!***********************************************************************
    subroutine remove_mean(f,inda,indep)
!
!  Substract mean from a (several) field(s) selected by the index inda
!  (the index range inda - indep) in f.
!
!  08-may-12/MR: adapted from remove_mean_flow
!
      use Mpicomm, only: mpiallreduce_sum
      use General, only: ioptest
!
      real, dimension (mx,my,mz,*), intent (inout)        :: f
      integer,                      intent (in)           :: inda
      integer,                      intent (in), optional :: indep
!
      real, allocatable, dimension(:) :: mean, mean_tmp
      integer :: m,n,j, inde
!
      inde = ioptest(indep,inda)
      allocate( mean(inda:inde), mean_tmp(inda:inde) )
!
!  initialize mean
!
      mean = 0.0
!
!  Go through all pencils.
!
      do n = n1,n2
      do m = m1,m2
!
!  Compute mean for each field.
!
        do j=inda,inde
          mean(j) = mean(j) + sum(f(l1:l2,m,n,j))
        enddo
      enddo
      enddo
!
!  Compute total mean for all processors
!
      call mpiallreduce_sum(mean/nwgrid,mean_tmp,inde-inda+1)
      mean = mean_tmp
!
!  Go through all pencils and subtract out the mean
!  separately for each field.
!
      do n = n1,n2
      do m = m1,m2
        do j=inda,inde
          f(l1:l2,m,n,j) = f(l1:l2,m,n,j) - mean(j)
        enddo
      enddo
      enddo
!
      if (lroot.and.ip<6) print*,'remove_mean: mean=',mean
!
      deallocate( mean, mean_tmp )
!
    endsubroutine remove_mean
!***********************************************************************
    subroutine insert_carray( array, insert, index, leng )
!
! inserts string vector insert in string vector array at position index; array is assumed
! to be settled until position leng; leng is updated.
! Oobs: routine does not check whether there is enough free space in
!       array for the insertion.
!
! 15-feb-2013/MR: parameter leng_insert removed (now derived from insert)
! 21-aug-2013/MR: moved from Testflow
!
      character(LEN=*), dimension(*) :: array
      character(LEN=*), dimension(:) :: insert
      integer                        :: index, leng, leng_insert, i
!
      intent(in)    :: index, insert
      intent(inout) :: leng, array
!
!  insert only if position index is meaningful
!
      if ( index>0.and.index<=leng+1 ) then
!
        leng_insert = size(insert)
        do i=leng,index,-1
          array(i+leng_insert) = array(i)
        enddo
!
        array(index:index+leng_insert-1) = insert
!
        leng = leng+leng_insert

      endif
!
    endsubroutine insert_carray
!***********************************************************************
    subroutine insert_rarray( array, insert, index, leng )
!
! 21-aug-2013/MR: derived from insert_carray for vectors of reals.
!
      integer            :: index, leng, leng_insert, i
      real, dimension(*) :: array
      real, dimension(:) :: insert
!
      intent(in)    :: index, insert
      intent(inout) :: leng, array
!
!  insert only if position index is meaningful
!
      if ( index>0.and.index<=leng+1 ) then
!
        leng_insert = size(insert)
        do i=leng,index,-1
          array(i+leng_insert) = array(i)
        enddo
!
        array(index:index+leng_insert-1) = insert
!
        leng = leng+leng_insert

      endif
!
    endsubroutine insert_rarray
!***********************************************************************
    subroutine insert_carray_mult( array, insert, mult, index, leng )
!
! inserts string insert mult times in string vector array at position index; array is assumed
! to be settled until position leng; leng is updated.
! Oobs: routine does not check whether there is enough free space in
!       array for the insertion.
!
! 15-feb-2013/MR: derived from insert_carray
!
      character(LEN=*), dimension(*) :: array
      character(LEN=*)               :: insert
      integer                        :: index, leng, mult, i
!
      intent(in)    :: index, insert, mult
      intent(inout) :: leng, array
!
!  insert only if position index is meaningful
!
      if ( index>0.and.index<=leng+1 ) then
!
        do i=leng,index,-1
          array(i+mult) = array(i)
        enddo
!
        array(index:index+mult-1) = insert
!
        leng = leng+mult
!
      endif
    endsubroutine insert_carray_mult
!***********************************************************************
    real function find_max_fvec(f,iv)
!
!  Find the maximum of the modulus of a vector field.
!
!  19-aug-2011/ccyang: coded
!  12-sep-2013/MR: outsourced from hydro
!
      use Mpicomm, only: mpiallreduce_max, MPI_COMM_WORLD
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      integer, intent(in) :: iv
!
      real :: umax1
!
!  Find the maximum.
!
      umax1 = sqrt(maxval(  f(l1:l2,m1:m2,n1:n2,iv  )**2 &
                          + f(l1:l2,m1:m2,n1:n2,iv+1)**2 &
                          + f(l1:l2,m1:m2,n1:n2,iv+2)**2  ))
      call mpiallreduce_max(umax1, find_max_fvec, MPI_COMM_WORLD)
!
    endfunction find_max_fvec
!***********************************************************************
    real function find_rms_fvec(f, iv)
!
!  Find the root-mean-square of the modulus of a vector field.
!
!  25-aug-2015/ccyang: coded
!
      use Mpicomm, only: mpiallreduce_sum
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      integer, intent(in) :: iv
!
      real :: s1, s
!
!  Find the rms.
!
      s1 = sum(f(l1:l2,m1:m2,n1:n2,iv:iv+2)**2)
      call mpiallreduce_sum(s1, s)
      find_rms_fvec = sqrt(s / real(nwgrid))
!
    endfunction find_rms_fvec
!***********************************************************************
    function find_xyrms_fvec(f, iv) result(rms)
!
!  Find the rms magnitude of a vector field in each z.
!
!  29-jun-14/ccyang: coded
!
      use Mpicomm, only: mpiallreduce_sum
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      integer, intent(in) :: iv
      real, dimension(nz) :: rms
!
      real, dimension(nz) :: rms_loc
      integer :: k
!
      forall(k = n1:n2) rms_loc(k-nghost) = sum(f(l1:l2,m1:m2,k,iv:iv+2)**2)
      call mpiallreduce_sum(rms_loc, rms, nz)
      rms = sqrt(rms / nxygrid)
!
    endfunction find_xyrms_fvec
!***********************************************************************
    subroutine finalize_aver_3D(nproc,idir,arrm)
!
!  Finalize calculation of an average of a 3D array arrm, by MPI communication
!  in direction idir, nproc - number of processors in this direction,
!  direction can be a beam or a layer
!
!  12-sep-2013/MR: coded
!
      use Mpicomm, only: mpiallreduce_sum
      use Messages, only: fatal_error
      use Yinyang_mpi, only: reduce_zsum
!
      integer,                intent(IN)   :: nproc,idir
      real, dimension(:,:,:), intent(INOUT):: arrm

      real, dimension(:,:,:), allocatable :: temp
      integer, dimension(3) :: sz
!
!  communicate in direction idir
!
        if (nproc>1) then
!
          sz=(/size(arrm,1),size(arrm,2),size(arrm,3)/)
          allocate(temp(sz(1),sz(2),sz(3)))
          if (lyinyang) then
            if (idir/=3) & 
              call fatal_error('finalize_aver_3D', &
                               'Not implemented for other than phi direction on Yin-Yang grid')
            !call reduce_zsum(arrm,temp)
            !call mpibcast_z_yy(temp)
          else
            call mpiallreduce_sum(arrm,temp,sz,idir=idir)
            arrm=temp
          endif
!
        endif
!
    endsubroutine finalize_aver_3D
!***********************************************************************
    subroutine finalize_aver_1D(nproc,idir,arrm)
!
!  Finalize calculation of an average of a 1D array arrm, by MPI communication
!  in direction idir, nproc - number of processors in this direction,
!  direction can be a beam or a layer
!
!  12-sep-2013/MR: coded
!
      use Mpicomm, only: mpiallreduce_sum
!
      integer,            intent(IN)   :: nproc,idir
      real, dimension(:), intent(INOUT):: arrm

      real, dimension(:), allocatable :: temp
      integer :: sz
!
!  communicate over direction idir
!
        if (nproc>1) then
!
          sz = size(arrm)
          allocate(temp(sz))
          if (lyinyang.and.idir==3) then
            !call zaverages_xy
          else
            call mpiallreduce_sum(arrm,temp,sz,idir=idir)
          endif
          arrm=temp
!
        endif
!
    endsubroutine finalize_aver_1D
!***********************************************************************
    subroutine finalize_aver_2D(nproc,idir,arrm)
!
!  Finalize calculation of an average of a 2D array arrm, by MPI communication
!  in direction idir, nproc - number of processors in this direction,
!  direction can be a beam or a layer
!
!  12-sep-2013/MR: coded
!
      use Mpicomm, only: mpiallreduce_sum
!
      integer,              intent(IN)   :: nproc,idir
      real, dimension(:,:), intent(INOUT):: arrm

      real, dimension(:,:), allocatable :: temp
      integer, dimension(2) :: sz
!
!  communicate over direction idir
!
        if (nproc>1) then
!
          sz=(/size(arrm,1),size(arrm,2)/)
          allocate(temp(sz(1),sz(2)))
          call mpiallreduce_sum(arrm,temp,sz,idir=idir)
          arrm=temp
!
        endif
!
    endsubroutine finalize_aver_2D
!***********************************************************************
    subroutine finalize_aver_4D(nproc,idir,arrm)
!
!  Finalize calculation of an average of a 4D array arrm, by MPI communication
!  in direction idir, nproc - number of processors in this direction,
!  direction can be a beam or a layer
!
!  12-sep-2013/MR: coded
!
      use Mpicomm, only: mpiallreduce_sum, mpireduce_sum   !!, mpicomm
!
!      include 'mpif.h'
!
      integer,                  intent(IN)   :: nproc,idir
      real, dimension(:,:,:,:), intent(INOUT):: arrm

      real, dimension(:,:,:,:), allocatable :: temp

      integer, dimension(4) :: sz
!
!  communicate over direction idir
!
        if (nproc>1) then
!
          sz=(/size(arrm,1),size(arrm,2),size(arrm,3),size(arrm,4)/)
          allocate(temp(sz(1),sz(2),sz(3),sz(4)))
          call mpiallreduce_sum(arrm,temp,sz,idir=idir)
!
          !!call MPI_BCAST(temp, product(sz), MPI_REAL, iprocx+nprocx*ipy, &
          !!mpicomm(idir),mpierr)
          arrm=temp
!
        endif
!
    endsubroutine finalize_aver_4D
!***********************************************************************
    subroutine calc_diffusive_flux(diffs,c_char,islope_limiter,h_slope_limited,flux)
!
!  Calculates diffusive flux from variable differences diffs acc. to Eqs. (6)-(10) in Rempel (2014).
!  c_char is vector of characteristic velocities.
!
!  23-sep-15/MR,joern,fred,petri: coded
!
      use General, only: notanumber

      real, dimension(:),intent(in) :: diffs,c_char
      real,              intent(in) :: h_slope_limited
      character(LEN=*),  intent(in) :: islope_limiter
      real, dimension(:),intent(out):: flux

      real, dimension(size(diffs)-1) :: slope
      real, dimension(size(diffs)-2) :: phi, right_left
      integer :: len

      len=size(diffs)

      call slope_limiter(diffs(2:),diffs(:len-1),slope,islope_limiter)
      right_left = diffs(2:len-1) - 0.5*(slope(2:) + slope(1:len-2))        ! = u_r - u_l

      call diff_flux(h_slope_limited, diffs(2:len-1), right_left, phi)
      flux = -0.5*c_char*phi*right_left
!!      flux = -0.5*0.3*phi*right_left
!!      print*, 'flux', flux(490:510)
if (notanumber(c_char)) then
   print*, 'CALC_DIFFUSIVE_FLUX: c_char=', len
   stop
endif 

    endsubroutine calc_diffusive_flux
!***********************************************************************
    subroutine slope_limiter(diff_right,diff_left,limited,type)
!
!  Returns limited slope in parameter limited, see Rempel (2014). 
!  Presently only limiter minmod is implemented. 
!
!  25-sep-15/MR,joern: coded
!  27-jan-16/MR: converted into non-elemental subroutine
!
      real, dimension(:), intent(IN) :: diff_left, diff_right
      real, dimension(:), intent(OUT):: limited
      character(LEN=*),   intent(IN) :: type
      integer :: len, ii

      len=size(diff_right)

      select case (type)
      case ('minmod-mine')
        where( sign(1.,diff_left)*sign(1.,diff_right)>0)
          limited = sign(min( 2.*abs(diff_left), 2.*abs(diff_right), &
                             0.5*abs(diff_right+diff_left) ), diff_left )
        elsewhere
          limited = 0.
        endwhere
      case ('minmod') 
         do ii=1,len
           if (diff_left(ii)*diff_right(ii)>0) then
             if (diff_left(ii)>0) limited(ii)=min(2.*diff_left(ii),2.*diff_right(ii), &
                                                  0.5*(diff_right(ii)+diff_left(ii)))
             if (diff_left(ii)<0) limited(ii)=max(2.*diff_left(ii),2.*diff_right(ii), &
                                                  0.5*(diff_right(ii)+diff_left(ii)))
           else
             limited(ii)=0
           endif
         enddo
      case ('superbee')
        limited=0.
        call fatal_error('slope_limiter','limiter not implemented')
      case ('')
        limited=0.
        call fatal_error('slope_limiter','limiter not implemented')
      case default
        limited=0.
        call fatal_error('slope_limiter','limiter not implemented')
      end select
    
    endsubroutine slope_limiter
!***********************************************************************
     subroutine diff_flux(h, diff_right, diff_lr, phi)
!
!  Calculates diffusive flux for one coordinate direction from u_i+1-u_i and u_r-u_l 
!  and returns it in parameter phi, see Rempel (2014). 
!
!  25-sep-15/MR,joern: coded
!  27-jan-16/MR: converted into non-elemental subroutine, because of malcompilation by gcc version 4.6.3
!
      real,               intent(IN) :: h
      real, dimension(:), intent(IN) :: diff_lr, diff_right
      real, dimension(:), intent(OUT):: phi

      where (diff_right*diff_lr>0.)
        phi=max(0.,1.+h*(diff_lr/diff_right-1.))
      elsewhere
        phi=0.
      endwhere
!
    endsubroutine diff_flux
!***********************************************************************
    subroutine periodic_fold_back(dd,Boxsize)
!
!  Comment on me!
!
      real,dimension(3),intent(in) :: Boxsize
      real,dimension(3),intent(inout) :: dd
!
      integer :: p,q,idim

      do idim=1,3
        q = floor(dd(idim)/(Boxsize(idim)/2))
        p = 0
        if (q.eq.1)  p = -1
        if (q.eq.-2) p = 1
        dd(idim) = dd(idim) + Boxsize(idim)*p
      enddo
!
    endsubroutine periodic_fold_back
!***********************************************************************
    subroutine calc_all_diff_fluxes(f,k,islope_limiter,h_slope_limited)
!
!  Calculates all <=3 components of the diffusive flux according to dimensionality. 
!
!  8-oct-15/MR: carved out from viscosity
!
      use General, only: notanumber
!
      real, dimension (:,:,:,:), intent(INOUT):: f
      integer,                   intent(IN)   :: k
      character(LEN=*),          intent(IN)   :: islope_limiter
      real,                      intent(IN)   :: h_slope_limited
             
      integer :: iff, nn, mm, ll
      real, dimension(mx-1) :: tmpx
      real, dimension(my-1) :: tmpy
      real, dimension(mz-1) :: tmpz

!      if (.not.lcartesian_coords) &
!        call fatal_error('calc_all_diff_fluxes', &
!                         'Slope-limited diffusion not implemented for curvilinear coordinates.')
! JW: comment out for testing

      f(:,:,:,iFF_diff1:iFF_diff2)=0.

      iff=iFF_diff

      if (nxgrid>1) then
        do nn=n1,n2; do mm=m1,m2
          tmpx = f(2:,mm,nn,k)-f(:mx-1,mm,nn,k)
if (notanumber(tmpx)) print*, 'TMPX:k,mm,nn=', k,mm,nn
!if (j==1) print*, 'TMPX:', tmpx
          call calc_diffusive_flux(tmpx,f(2:mx-2,mm,nn,iFF_char_c),islope_limiter,h_slope_limited,f(2:mx-2,mm,nn,iff))
if (notanumber(f(2:mx-2,mm,nn,iff))) print*, 'DIFFX:k,mm,nn=', k,mm,nn
        enddo; enddo
        iff=iff+1
      endif

      if (nygrid>1) then
        do nn=n1,n2; do ll=l1,l2
          tmpy = f(ll,2:,nn,k)-f(ll,:my-1,nn,k)
if (notanumber(tmpy)) print*, 'TMPY:k,mm,nn=', k,mm,nn
!if (j==2) print*, 'UY, iproc=:', iproc    !, '1-5,my-5:my-1'
!if (j==2) print'(16(e13.6,1x))', f(ll,:,nn,iuu+j-1)
          call calc_diffusive_flux(tmpy,f(ll,2:my-2,nn,iFF_char_c),islope_limiter,h_slope_limited,f(ll,2:my-2,nn,iff))
if (notanumber(f(ll,2:my-2,nn,iff))) print*, 'DIFFY:k,ll,mm=', k,ll,mm
        enddo; enddo
        iff=iff+1
      endif

      if (nzgrid>1) then
        do mm=m1,m2; do ll=l1,l2
if (notanumber(f(ll,mm,n1:n2,k))) print*, 'FARRAY(n1:n2):k,ll,mm=', k,ll,mm
if (notanumber(f(ll,mm,:,k))) print*, 'FARRAY(:):k,ll,mm=', k,ll,mm
          tmpz = f(ll,mm,2:,k)-f(ll,mm,:mz-1,k)
if (notanumber(tmpz)) print*, 'TMPZ:k,ll,mm=', k,ll,mm
          call calc_diffusive_flux(tmpz,f(ll,mm,2:mz-2,iFF_char_c),islope_limiter,h_slope_limited,f(ll,mm,2:mz-2,iff))
if (notanumber(f(ll,mm,2:mz-2,iff))) print*, 'DIFFZ:k,ll,mm=', k,ll,mm
        enddo; enddo
      endif
!
    endsubroutine calc_all_diff_fluxes
!***********************************************************************
    subroutine calc_sij2(f,sij2,lshear_rateofstrain)
!
!  Calculates square of traceless rate-of-strain tensor from f array.
!
!  16-dec-16/MR: Outsourced from hydro.
!
      real, dimension (mx,my,mz,mfarray), intent(IN) :: f
      real, dimension(nx),                intent(OUT):: sij2
      logical,                  optional, intent(IN) :: lshear_rateofstrain

      real, dimension(nx,3) :: uu
      real, dimension(nx,3,3) :: uij, sij

! uij from f
      call gij(f,iuu,uij,1)
      uu=f(l1:l2,m,n,iux:iuz)
! divu -> uij2
      call div_mn(uij,sij2,uu)
! sij
      call traceless_strain(uij,sij2,sij,uu,lshear_rateofstrain)
! sij^2
      call multm2_sym_mn(sij,sij2)

    endsubroutine calc_sij2
!***********************************************************************
    subroutine traceless_strain(uij,divu,sij,uu,lss)
!
!  Calculates traceless rate-of-strain tensor sij from derivative tensor uij
!  and divergence divu within each pencil;
!  curvilinear co-ordinates require optional velocity argument uu.
!  In-place operation is possible, i.e. uij and sij may refer to the same array.
!
!  16-oct-09/MR: carved out from calc_pencils_hydro
!  10-apr-11/MR: optional parameter lss added, replaces use of global lshear_rateofstrain
!
    use General, only: loptest

    real, dimension (nx,3,3)         :: uij, sij
    real, dimension (nx)             :: divu
    real, dimension (nx,3), optional :: uu
    logical,                optional :: lss
!
    integer :: i,j
    logical :: lshear_ROS
!
    intent(in)  :: uij, divu, lss
    intent(out) :: sij
!
    lshear_ROS=lshear.and.loptest(lss)
!
    do j=1,3
      sij(:,j,j)=uij(:,j,j)-(1./3.)*divu
      do i=j+1,3
        sij(:,i,j)=.5*(uij(:,i,j)+uij(:,j,i))
        sij(:,j,i)=sij(:,i,j)
      enddo
    enddo
!
    if (lspherical_coords.or.lcylindrical_coords) then
      if (.not.present(uu)) then
        call fatal_error('traceless_strain', &
            'Deformation matrix for curvilinear co-ordinates'// &
            'requires providing of the velocity itself')
        return
      endif
    endif
!
    if (lspherical_coords) then
! sij(:,1,1) remains unchanged in spherical coordinates
      sij(:,1,2)=sij(:,1,2)-.5*r1_mn*uu(:,2)
      sij(:,1,3)=sij(:,1,3)-.5*r1_mn*uu(:,3)
      sij(:,2,1)=sij(:,1,2)
      sij(:,2,2)=sij(:,2,2)+r1_mn*uu(:,1)
      sij(:,2,3)=sij(:,2,3)-.5*r1_mn*cotth(m)*uu(:,3)
      sij(:,3,1)=sij(:,1,3)
      sij(:,3,2)=sij(:,2,3)
      sij(:,3,3)=sij(:,3,3)+r1_mn*uu(:,1)+cotth(m)*r1_mn*uu(:,2)
    elseif (lcylindrical_coords) then
      sij(:,1,2)=sij(:,1,2)-.5*rcyl_mn1*uu(:,2)
      sij(:,2,2)=sij(:,2,2)+.5*rcyl_mn1*uu(:,1)
      sij(:,2,1)=sij(:,1,2)
    endif
!
    if (lshear_ROS) then
      sij(:,1,2)=sij(:,1,2)+Sshear
      sij(:,2,1)=sij(:,2,1)+Sshear
    endif
!
    endsubroutine traceless_strain
!***********************************************************************
    subroutine remove_mean_value(f,ind1,ind2)
!
!  Substract mean x-flow from the x-velocity field.
!  Useful to avoid unphysical winds in shearing box simulations.
!  Note: this is possibly not useful when there is rotation, because
!  then epicyclic motions don't usually grow catastrophically.
!
!  22-may-07/axel: adapted from remove_mean_momenta
!  15-dec-10/MR  : added parameters indux to make applicable to other ...
!  13-may-18/axel: moved to general, renamed, added end index
!
      use Mpicomm, only: mpiallreduce_sum
!
      real, dimension (mx,my,mz,mfarray), intent (inout) :: f
      integer,                            intent (in)    :: ind1,ind2
!
      real, dimension (nx) :: uu
      real, dimension (ind1:ind2) :: um, um_tmp
      integer :: m,n,j
      real    :: fac
!
!  initialize um and compute normalization factor fac
!
        um = 0.0
        fac = 1.0/nwgrid
!
!  Go through all pencils.
!
        do n = n1,n2
        do m = m1,m2
!
!  Compute mean flow in each of the 3 directions.
!
          do j=ind1,ind2
            uu = f(l1:l2,m,n,j)
            um(j) = um(j) + fac*sum(uu)
          enddo
        enddo
        enddo
!
!  Compute total sum for all processors
!
        call mpiallreduce_sum(um,um_tmp,ind2-ind1+1)
        um = um_tmp
!
!  Go through all pencils and subtract out the mean value
!  separately for each direction.
!
        do n = n1,n2
        do m = m1,m2
          do j=ind1,ind2
            f(l1:l2,m,n,j) = f(l1:l2,m,n,j) - um(j)
          enddo
        enddo
        enddo
!
    endsubroutine remove_mean_value
!***********************************************************************
    subroutine stagger_to_base_interp_1st(f,m,n,fint)
!
!  Performs first-order interpolation from the staggered grid to the original one.
!  The interpolation point is the center of the cuboid defined by the intervals 
!  (i1-1,i1), (j1-1,j1), (k1-1,k1) of the staggered grid.
!
!  10-jul-18/MR: coded
! 
      real, dimension(mx,my,mz), intent(IN) :: f      ! values on staggered grid
      integer,                   intent(IN) :: m,n    ! position in mn-loop at which interpolation is performed
      real, dimension(nx),       intent(OUT):: fint   ! interpoland
!
      real :: fac, facq, facc
      integer, dimension( 3,-1:0) :: ijk
      integer, dimension(-1:0) :: i2
      integer, dimension(-1:0), parameter :: shifts=(/-1,0/)
      integer :: k

      ijk(1,:)=l1; ijk(2,:)=m; ijk(3,:)=n

      fac=1./2.
      do k=1,3
        if (lactive_dimension(k)) ijk(k,:)=ijk(k,:)+shifts
!print*, 'k,ijk(k,:)=', k, ijk(k,:)
      enddo
      i2 = ijk(1,:)+nx-1

      if (dimensionality==1) then
        fint = (f(ijk(1,-1):i2(-1),ijk(2,-1),ijk(3,-1))+f(ijk(1,0):i2(0),ijk(2,0),ijk(3,0)))*fac
      else
        facq=fac*fac
        if (dimensionality==3) then
          facc=facq*fac
          fint = ( f(ijk(1,-1):i2(-1),ijk(2,-1),ijk(3,-1))+f(ijk(1, 0):i2( 0),ijk(2, 0),ijk(3, 0)) &
                  +f(ijk(1,-1):i2(-1),ijk(2,-1),ijk(3, 0))+f(ijk(1,-1):i2(-1),ijk(2, 0),ijk(3,-1)) &
                  +f(ijk(1, 0):i2( 0),ijk(2,-1),ijk(3,-1))+f(ijk(1,-1):i2(-1),ijk(2, 0),ijk(3, 0)) &
                  +f(ijk(1, 0):i2( 0),ijk(2,-1),ijk(3, 0))+f(ijk(1, 0):i2( 0),ijk(2, 0),ijk(3,-1)) )*facc
        elseif (lactive_dimension(3)) then
          fint = ( f(ijk(1,-1):i2(-1),ijk(2,-1),ijk(3,-1))+f(ijk(1, 0):i2( 0),ijk(2, 0),ijk(3,-1)) &
                  +f(ijk(1,-1):i2(-1),ijk(2,-1),ijk(3, 0))+f(ijk(1, 0):i2( 0),ijk(2, 0),ijk(3, 0)))*facq 
        else
!print*, ijk(1, 0),i2( 0),ijk(2,-1),n,f(ijk(1, 0):i2( 0),ijk(2,-1),n)
          fint = ( f(ijk(1,-1):i2(-1),ijk(2,-1),n)+f(ijk(1, 0):i2( 0),ijk(2, 0),n) &
                  +f(ijk(1,-1):i2(-1),ijk(2, 0),n)+f(ijk(1, 0):i2( 0),ijk(2,-1),n))*facq 
        endif
      endif

      endsubroutine stagger_to_base_interp_1st
!***********************************************************************
    subroutine stagger_to_base_interp_3rd(f,m,n,fint)
!
!  Performs third-order interpolation from the staggered grid to the original one.
!  The interpolation point is the center of the cuboid defined by the intervals 
!  (i1-1,i1), (j1-1,j1), (k1-1,k1) of the staggered grid.
!
!  10-jul-18/MR: coded
! 
      real, dimension(mx,my,mz), intent(IN) :: f      ! values on staggered grid
      integer,                   intent(IN) :: m,n    ! position in mn-loop at which interpolation is performed
      real, dimension(nx),       intent(OUT):: fint   ! interpoland
!
      real :: facm, facp, facmq, facpq, facpm, facpc, facpqm, facpmq, facmc
      integer, dimension( 3,-2:1) :: ijk
      integer, dimension(-2:1) :: i2
      integer, dimension(-2:1), parameter :: shifts=(/-2,-1,0,1/)
      integer :: k

      ijk(1,:)=l1; ijk(2,:)=m; ijk(3,:)=n

      facm=1./16.; facp=9./16.
      do k=1,3
        if (lactive_dimension(k)) ijk(k,:)=ijk(k,:)+shifts
!print*, 'k,ijk(k,:)=', k, ijk(k,:)
      enddo
      i2 = ijk(1,:)+nx-1

      if (dimensionality==1) then
        fint = (f(ijk(1,-1):i2(-1),ijk(2,-1),ijk(3,-1))+f(ijk(1,0):i2(0),ijk(2,0),ijk(3,0)))*facp &
              -(f(ijk(1,-2):i2(-2),ijk(2,-2),ijk(3,-2))+f(ijk(1,1):i2(1),ijk(2,1),ijk(3,1)))*facm
      else
        facmq=facm*facm; facpq=facp*facp; facpm=facp*facm
        if (dimensionality==2) then
          if (lactive_dimension(3)) then
            fint =  ( f(ijk(1,-1):i2(-1),ijk(2,-1),ijk(3,-1))+f(ijk(1, 0):i2( 0),ijk(2, 0),ijk(3, 0))         &
                     +f(ijk(1, 0):i2( 0),ijk(2, 0),ijk(3,-1))+f(ijk(1,-1):i2(-1),ijk(2,-1),ijk(3, 0)) )*facpq &
!
                   -( f(ijk(1,-1):i2(-1),ijk(2,-1),ijk(3,+1))+f(ijk(1,+1):i2(+1),ijk(2,+1),ijk(3,-1))         &
                     +f(ijk(1,-1):i2(-1),ijk(2,-1),ijk(3,-2))+f(ijk(1,-2):i2(-2),ijk(2,-2),ijk(3,-1))         &
                     +f(ijk(1, 0):i2( 0),ijk(2, 0),ijk(3,-2))+f(ijk(1,-2):i2(-2),ijk(2,-2),ijk(3, 0))         &
                     +f(ijk(1, 0):i2( 0),ijk(2, 0),ijk(3,+1))+f(ijk(1,+1):i2(+1),ijk(2,+1),ijk(3, 0)) )*facpm &
!
                   +( f(ijk(1,-2):i2(-2),ijk(2,-2),ijk(3,-2))+f(ijk(1,+1):i2(+1),ijk(2,+1),ijk(3,+1))         &
                     +f(ijk(1,-2):i2(-2),ijk(2,-2),ijk(3,+1))+f(ijk(1,+1):i2(+1),ijk(2,+1),ijk(3,-2)) )*facmq
          else
            fint =  ( f(ijk(1,-1):i2(-1),ijk(2,-1),n)+f(ijk(1, 0):i2( 0),ijk(2, 0),n)         &
                     +f(ijk(1,-1):i2(-1),ijk(2, 0),n)+f(ijk(1, 0):i2( 0),ijk(2,-1),n) )*facpq &
!
                   -( f(ijk(1,-1):i2(-1),ijk(2,+1),n)+f(ijk(1,+1):i2(+1),ijk(2,-1),n)         &
                     +f(ijk(1,-1):i2(-1),ijk(2,-2),n)+f(ijk(1,-2):i2(-2),ijk(2,-1),n)         &
                     +f(ijk(1, 0):i2( 0),ijk(2,-2),n)+f(ijk(1,-2):i2(-2),ijk(2, 0),n)         &
                     +f(ijk(1, 0):i2( 0),ijk(2,+1),n)+f(ijk(1,+1):i2(+1),ijk(2, 0),n) )*facpm &
!
                   +( f(ijk(1,-2):i2(-2),ijk(2,-2),n)+f(ijk(1,+1):i2(+1),ijk(2,+1),n)         &
                     +f(ijk(1,-2):i2(-2),ijk(2,+1),n)+f(ijk(1,+1):i2(+1),ijk(2,-2),n) )*facmq
          endif
        else
          facpc =facpq*facp
          facpqm=facpq*facm
          facpmq=facp*facmq
          facmc =facmq*facm
          fint = ( f(ijk(1,-1):i2(-1),ijk(2,-1),ijk(3,-1))+f(ijk(1, 0):i2( 0),ijk(2,-1),ijk(3,-1)) &
                  +f(ijk(1,-1):i2(-1),ijk(2, 0),ijk(3,-1))+f(ijk(1,-1):i2(-1),ijk(2,-1),ijk(3, 0)) &
                  +f(ijk(1, 0):i2( 0),ijk(2, 0),ijk(3,-1))+f(ijk(1, 0):i2( 0),ijk(2,-1),ijk(3, 0)) &
                  +f(ijk(1,-1):i2(-1),ijk(2, 0),ijk(3, 0))+f(ijk(1, 0):i2( 0),ijk(2, 0),ijk(3, 0)) )*facpc &  ! 8 terms
!
                -( f(ijk(1,-1):i2(-1),ijk(2,-1),ijk(3,-2))+f(ijk(1, 0):i2( 0),ijk(2,-1),ijk(3,-2))     &
                  +f(ijk(1,-1):i2(-1),ijk(2, 0),ijk(3,-2))+f(ijk(1, 0):i2( 0),ijk(2, 0),ijk(3,-2))     &
                  +f(ijk(1,-1):i2(-1),ijk(2,-1),ijk(3,+1))+f(ijk(1, 0):i2( 0),ijk(2,-1),ijk(3,+1))     &
                  +f(ijk(1,-1):i2(-1),ijk(2, 0),ijk(3,+1))+f(ijk(1, 0):i2( 0),ijk(2, 0),ijk(3,+1))     &
                  +f(ijk(1,-1):i2(-1),ijk(2,-2),ijk(3,-1))+f(ijk(1, 0):i2( 0),ijk(2,-2),ijk(3,-1))     &
                  +f(ijk(1,-1):i2(-1),ijk(2,-2),ijk(3, 0))+f(ijk(1, 0):i2( 0),ijk(2,-2),ijk(3, 0))     &
                  +f(ijk(1,-1):i2(-1),ijk(2,+1),ijk(3,-1))+f(ijk(1, 0):i2( 0),ijk(2,+1),ijk(3,-1))     &
                  +f(ijk(1,-1):i2(-1),ijk(2,+1),ijk(3, 0))+f(ijk(1, 0):i2( 0),ijk(2,+1),ijk(3, 0))     &
                  +f(ijk(1,-2):i2(-2),ijk(2,-1),ijk(3,-1))+f(ijk(1,-2):i2(-2),ijk(2, 0),ijk(3,-1))     &
                  +f(ijk(1,-2):i2(-2),ijk(2,-1),ijk(3, 0))+f(ijk(1,-2):i2(-2),ijk(2, 0),ijk(3, 0))     &
                  +f(ijk(1,+1):i2(+1),ijk(2,-1),ijk(3,-1))+f(ijk(1,+1):i2(+1),ijk(2, 0),ijk(3,-1))     &
                  +f(ijk(1,+1):i2(+1),ijk(2,-1),ijk(3, 0))+f(ijk(1,+1):i2(+1),ijk(2, 0),ijk(3, 0)) )*facpqm & ! 24 terms
!
                +( f(ijk(1,-2):i2(-2),ijk(2,-2),ijk(3,-1))+f(ijk(1,+1):i2(+1),ijk(2,-2),ijk(3,-1))     &
                  +f(ijk(1,-2):i2(-2),ijk(2,+1),ijk(3,-1))+f(ijk(1,+1):i2(+1),ijk(2,+1),ijk(3,-1))     &
                  +f(ijk(1,-2):i2(-2),ijk(2,-2),ijk(3, 0))+f(ijk(1,+1):i2(+1),ijk(2,-2),ijk(3, 0))     &
                  +f(ijk(1,-2):i2(-2),ijk(2,+1),ijk(3, 0))+f(ijk(1,+1):i2(+1),ijk(2,+1),ijk(3, 0))     &
                  +f(ijk(1,-2):i2(-2),ijk(2,-1),ijk(3,-2))+f(ijk(1,+1):i2(+1),ijk(2,-1),ijk(3,-2))     &
                  +f(ijk(1,-2):i2(-2),ijk(2,-1),ijk(3,+1))+f(ijk(1,+1):i2(+1),ijk(2,-1),ijk(3,+1))     &
                  +f(ijk(1,-2):i2(-2),ijk(2, 0),ijk(3,-2))+f(ijk(1,+1):i2(+1),ijk(2, 0),ijk(3,-2))     &
                  +f(ijk(1,-2):i2(-2),ijk(2, 0),ijk(3,+1))+f(ijk(1,+1):i2(+1),ijk(2, 0),ijk(3,+1))     &
                  +f(ijk(1,-1):i2(-1),ijk(2,-2),ijk(3,-2))+f(ijk(1,-1):i2(-1),ijk(2,+1),ijk(3,-2))     &
                  +f(ijk(1,-1):i2(-1),ijk(2,-2),ijk(3,+1))+f(ijk(1,-1):i2(-1),ijk(2,+1),ijk(3,+1))     &
                  +f(ijk(1, 0):i2( 0),ijk(2,-2),ijk(3,-2))+f(ijk(1, 0):i2( 0),ijk(2,+1),ijk(3,-2))     &
                  +f(ijk(1, 0):i2( 0),ijk(2,-2),ijk(3,+1))+f(ijk(1, 0):i2( 0),ijk(2,+1),ijk(3,+1)) )*facpmq & ! 24 terms
!
                -( f(ijk(1,-2):i2(-2),ijk(2,-2),ijk(3,-2))+f(ijk(1,+1):i2(+1),ijk(2,-2),ijk(3,-2))     &
                  +f(ijk(1,-2):i2(-2),ijk(2,+1),ijk(3,-2))+f(ijk(1,-2):i2(-2),ijk(2,-2),ijk(3,+1))     &
                  +f(ijk(1,+1):i2(+1),ijk(2,+1),ijk(3,-2))+f(ijk(1,+1):i2(+1),ijk(2,-2),ijk(3,+1))     &
                  +f(ijk(1,-2):i2(-2),ijk(2,+1),ijk(3,+1))+f(ijk(1,+1):i2(+1),ijk(2,+1),ijk(3,+1)) )*facmc    ! 8 terms
        endif
      endif

    endsubroutine stagger_to_base_interp_3rd
!***********************************************************************
endmodule Sub
