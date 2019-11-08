! $Id$
!
!  Module with general utility subroutines.
!
module General
!
  use Cparam
!
  implicit none
!
  private
!
  public :: gaunoise_number
  public :: safe_character_assign, safe_character_append, safe_character_prepend
  public :: lower_case
  public :: random_seed_wrapper
  public :: random_number_wrapper, random_gen, normal_deviate
  public :: parse_filename
!
  public :: keep_compiler_quiet
!
  public :: setup_mm_nn
  public :: find_index_range, find_index, find_index_range_hill, pos_in_array
  public :: find_proc
!
  public :: spline, tridag, pendag, complex_phase, erfcc
  public :: cspline
  public :: polynomial_interpolation
  public :: arcsinh
  public :: besselj_nu_int, calc_complete_ellints
  public :: bessj, cyclic
  public :: spline_derivative_double, spline_integral, linear_interpolate
  public :: itoa, count_bits, parser, write_full_columns
  public :: read_range, merge_ranges, add_merge_range, get_range_no, write_by_ranges, &
            write_by_ranges_1d_real, write_by_ranges_1d_cmplx, &
            write_by_ranges_2d_real, write_by_ranges_2d_cmplx
  public :: compress, fcompress
  public :: quick_sort, binary_search
  public :: date_time_string
  public :: backskip
  public :: lextend_vector
  public :: operator(.in.)
  public :: loptest, ioptest, roptest, doptest, coptest
  public :: indgen
  public :: ranges_dimensional
  public :: staggered_mean_scal, staggered_mean_vec
  public :: staggered_max_scal, staggered_max_vec
  public :: directory_names_std, numeric_precision
  public :: touch_file
  public :: var_is_vec
  public :: transform_cart_spher, transform_cart_spher_other, transform_spher_cart_yy
  public :: yy_transform_strip, yy_transform_strip_other, yin2yang_coors
  public :: transform_thph_yy, transform_thph_yy_other, merge_yin_yang
  public :: copy_kinked_strip_z, copy_kinked_strip_y, reset_triangle
  public :: transpose_mn
  public :: notanumber, notanumber_0d
  public :: reduce_grad_dim
  public :: meshgrid
  public :: linspace
  public :: linear_interpolate_2d
  public :: chk_time
  public :: get_species_nr
  public :: get_from_nml_str,get_from_nml_log,get_from_nml_real,convert_nml
!
  interface random_number_wrapper
    module procedure random_number_wrapper_0
    module procedure random_number_wrapper_1
    module procedure random_number_wrapper_3
  endinterface
!
  interface keep_compiler_quiet ! Overload `keep_compiler_quiet' function
    module procedure keep_compiler_quiet_r
    module procedure keep_compiler_quiet_r1d
    module procedure keep_compiler_quiet_r2d
    module procedure keep_compiler_quiet_r3d
    module procedure keep_compiler_quiet_r4d
    module procedure keep_compiler_quiet_p
    module procedure keep_compiler_quiet_bc
    module procedure keep_compiler_quiet_sl
    module procedure keep_compiler_quiet_i
    module procedure keep_compiler_quiet_i1d
    module procedure keep_compiler_quiet_i81d
    module procedure keep_compiler_quiet_i2d
    module procedure keep_compiler_quiet_i3d
    module procedure keep_compiler_quiet_l
    module procedure keep_compiler_quiet_l1d
    module procedure keep_compiler_quiet_c
  endinterface
!
  interface safe_character_append
    module procedure safe_character_append_2
    module procedure safe_character_append_3
  endinterface
!
  interface safe_character_prepend
    module procedure safe_character_prepend_2
  endinterface
!
  interface write_full_columns
    module procedure write_full_columns_real
    module procedure write_full_columns_cmplx
  endinterface
!
  interface read_range
    module procedure read_range_r
    module procedure read_range_i
  endinterface

  interface write_by_ranges
    module procedure write_by_ranges_1d_real
    module procedure write_by_ranges_1d_cmplx
    module procedure write_by_ranges_2d_real
    module procedure write_by_ranges_2d_cmplx
  endinterface

  interface lextend_vector
    module procedure lextend_vector_float
    module procedure lextend_vector_char
  endinterface
!
  interface polynomial_interpolation
    module procedure poly_interp_one
    module procedure poly_interp_fixorder
  endinterface
!
  interface pos_in_array
    module procedure pos_in_array_int
    module procedure pos_in_array_char
  endinterface
!
  interface in_array
    module procedure in_array_int
    module procedure in_array_char
  endinterface
!
  interface operator (.in.)
    module procedure in_array_int
    module procedure in_array_char
  endinterface
!
  interface notanumber          ! Overload the `notanumber' function
    module procedure notanumber_0
    module procedure notanumber_1
    module procedure notanumber_2
    module procedure notanumber_3
    module procedure notanumber_4
    module procedure notanumber_5
  endinterface
!
  interface meshgrid
    module procedure meshgrid_2d
    module procedure meshgrid_3d
  endinterface meshgrid

  interface compress
    module procedure compress_vec
    module procedure compress_arr
  endinterface
!
!  State and default generator of random numbers.
!
  integer, save, dimension(mseed) :: rstate=0
  character (len=labellen) :: random_gen='min_std'

  include 'general.h'
  !include 'general_f2003.h'
!
!***********************************************************************
    pure integer function find_proc(ipx, ipy, ipz)
!
!  Returns the rank of a process given its position in (ipx,ipy,ipz).
!
!  16-sep-15/ccyang: coded.
!
      use Cdata, only: lprocz_slowest
!
      integer, intent(in) :: ipx, ipy, ipz
!
      if (lprocz_slowest) then
        find_proc = ipz * nprocxy + ipy * nprocx + ipx
      else
        find_proc = ipy * nprocxz + ipz * nprocx + ipx
      endif
!
    endfunction find_proc
!***********************************************************************
    subroutine setup_mm_nn
!
!  Produce index-array for the sequence of points to be worked through:
!  Before the communication has been completed, the nghost=3 layers next
!  to the processor boundary (m1, m2, n1, or n2) cannot be used yet.
!  In the mean time we can calculate the interior points sufficiently far
!  away from the boundary points. Here we calculate the order in which
!  m and n are executed. At one point, necessary(imn)=.true., which is
!  the moment when all communication must be completed.
!
      use Cdata, only: mm,nn,imn_array,necessary,lroot,ip, &
                       lyinyang,lcutoff_corners,nycut,nzcut, &
                       lfirst_proc_y, lfirst_proc_z, llast_proc_y, llast_proc_z, iproc
!
      integer :: imn,m,n
      integer :: min_m1i_m2,max_m2i_m1
!
!  For non-parallel runs simply go through m and n from bottom left to to right.
!
      imn_array=0
      if (ncpus==1) then
        imn=1
        necessary(1)=.true.
        do n=n1,n2
          do m=m1,m2
            mm(imn)=m
            nn(imn)=n
            imn_array(m,n)=imn
            imn=imn+1
          enddo
        enddo
      else
!
!  For parallelized runs:
!  Do inner rectangle.
!
        imn=1
        do n=n1i+2,n2i-2
          do m=m1i+2,m2i-2
            if (imn_array(m,n) == 0) then
              mm(imn)=m
              nn(imn)=n
              imn_array(m,n)=imn
              imn=imn+1
            endif
          enddo
        enddo
        necessary(imn)=.true.
!
!  Do the upper stripe in the n-direction.
!
        do n=max(n2i-1,n1+1),n2
          do m=m1i+2,m2i-2
            if (imn_array(m,n) == 0) then
              mm(imn)=m
              nn(imn)=n
              imn_array(m,n)=imn
              imn=imn+1
            endif
          enddo
        enddo
!
!  Do the lower stripe in the n-direction.
!
        do n=n1,min(n1i+1,n2)
          do m=m1i+2,m2i-2
            if (imn_array(m,n) == 0) then
              mm(imn)=m
              nn(imn)=n
              imn_array(m,n)=imn
              imn=imn+1
            endif
          enddo
        enddo
!
!  Left and right hand boxes.
!  NOTE: need to have min(m1i,m2) instead of just m1i, and max(m2i,m1)
!  instead of just m2i, to make sure the case ny=1 works ok, and
!  also that the same m is not set in both loops.
!  ALSO: need to make sure the second loop starts not before the
!  first one ends; therefore max_m2i_m1+1=max(m2i,min_m1i_m2+1).
!
        min_m1i_m2=min(m1i+1,m2)
        max_m2i_m1=max(m2i-1,min_m1i_m2+1)
!
        do n=n1,n2
          do m=m1,min_m1i_m2
            if (imn_array(m,n) == 0) then
              mm(imn)=m
              nn(imn)=n
              imn_array(m,n)=imn
              imn=imn+1
            endif
          enddo
          do m=max_m2i_m1,m2
            if (imn_array(m,n) == 0) then
              mm(imn)=m
              nn(imn)=n
              imn_array(m,n)=imn
              imn=imn+1
            endif
          enddo
        enddo
      endif

      if (lyinyang.and.lcutoff_corners) then
        if (lfirst_proc_y) then
          if (lfirst_proc_z) call reset_triangle_inds(1,my-nycut+nghost, 1,mz-nzcut+nghost,imn_array)  ! lower left corner
          if (llast_proc_z ) call reset_triangle_inds(1,my-nycut+nghost,mz,nzcut+1-nghost ,imn_array)  ! upper left
        endif
        if (llast_proc_y) then
          if (lfirst_proc_z) call reset_triangle_inds(my,nycut+1-nghost,  1,mz-nzcut+nghost,imn_array) ! lower right
          if (llast_proc_z ) call reset_triangle_inds(my,nycut+1-nghost, mz,nzcut+1-nghost ,imn_array) ! upper right
        endif
      endif
!
!  Debugging output to be analysed with $PENCIL_HOME/utils/check-mm-nn.
!
      if (ip<=6) then
        if (lroot) then
          do imn=1,ny*nz
            if (necessary(imn)) write(*,'(A)') '==MM==NN==> Necessary'
            write(*,'(A,I5,I5)') '==MM==NN==> m,n= ', mm(imn), nn(imn)
          enddo
        endif
      endif
      !if (iproc==0) write(100,'(30(i3,1x))') imn_array
!
    endsubroutine setup_mm_nn
!***********************************************************************
    subroutine gaunoise_number(gn)
!
! Gaussian random number generator with unit variance and zero mean.
! Returns a pair of random numbers.
!
      real,dimension(2),intent(out) :: gn
      real :: r,p
      r = 0.0
      do while (r == 0.0)
        call random_number_wrapper(r)
      enddo
      call random_number_wrapper(p)
      gn(1)=sqrt(-2*log(r))*sin(twopi*p)
      gn(2)=sqrt(-2*log(r))*cos(twopi*p)
    endsubroutine gaunoise_number
!***********************************************************************
    subroutine random_number_wrapper_0(a,channel)
!
!  Fills a with a random number calculated with one of the generators
!  available with random_gen.
!
      use Cdata, only: lroot
!
      integer, optional :: channel
      real, intent(out) :: a
!
      select case (random_gen)
!
        case ('system')
          call random_number(a)
        case ('min_std')
          a=ran0(rstate(1))
        case ('nr_f90')
          if (present(channel)) then
            a=mars_ran2()
          else
            a=mars_ran()
          endif
 
        case default
          if (lroot) print*, 'No such random number generator: ', random_gen
          STOP 1                ! Return nonzero exit status
!
      endselect
!
    endsubroutine random_number_wrapper_0
!***********************************************************************
    subroutine random_number_wrapper_1(a,channel)
!
!  Fills a with an array of random numbers calculated with one of the
!  generators available with random_gen.
!
      use Cdata, only: lroot
!
      real, dimension(:), intent(out) :: a
      integer, optional :: channel
      integer :: i
!
      select case (random_gen)
!
        case ('system')
          call random_number(a)
        case ('min_std')
          do i=1,size(a,1)
            a(i)=ran0(rstate(1))
          enddo
        case ('nr_f90')
          if (present(channel)) then
            do i=1,size(a,1)
              a(i)=mars_ran2()
            enddo
          else
            do i=1,size(a,1)
              a(i)=mars_ran()
            enddo
          endif
        case default
          if (lroot) print*, 'No such random number generator: ', random_gen
          STOP 1                ! Return nonzero exit status
!
      endselect
!
    endsubroutine random_number_wrapper_1
!***********************************************************************
    subroutine random_number_wrapper_3(a,channel)
!
!  Fills a with a matrix of random numbers calculated with one of the
!  generators available with random_gen.
!
      use Cdata, only: lroot
!
      real, dimension(:,:,:), intent(out) :: a
      integer, optional :: channel
      integer :: i,j,k
!
      select case (random_gen)
!
        case ('system')
          call random_number(a)
        case ('min_std')
          do i=1,size(a,1); do j=1,size(a,2); do k=1,size(a,3)
            a(i,j,k)=ran0(rstate(1))
          enddo; enddo; enddo
        case ('nr_f90')
          do i=1,size(a,1); do j=1,size(a,2); do k=1,size(a,3)
            a(i,j,k)=mars_ran()
          enddo; enddo; enddo

          if (present(channel)) then
            do i=1,size(a,1); do j=1,size(a,2); do k=1,size(a,3)
              a(i,j,k)=mars_ran2()
            enddo; enddo; enddo
          else
            do i=1,size(a,1); do j=1,size(a,2); do k=1,size(a,3)
              a(i,j,k)=mars_ran()
            enddo; enddo; enddo
          endif
        case default
          if (lroot) print*, 'No such random number generator: ', random_gen
          STOP 1                ! Return nonzero exit status
!
      endselect
!
    endsubroutine random_number_wrapper_3
!***********************************************************************
    subroutine normal_deviate(a)
!
!  Return a normal deviate (Gaussian random number, mean=0, var=1) by means
!  of the "Ratio-of-uniforms" method.
!
!  28-jul-08/kapelrud: coded
!
      real,intent(out) :: a
!
      real :: u,v,x,y,q
      logical :: lloop
!
      lloop=.true.
      do while(lloop)
        call random_number_wrapper_0(u)
        call random_number_wrapper_0(v)
        v=1.7156*(v-0.5)
        x=u-0.449871
        y=abs(v)+0.386595
        q=x**2+y*(0.19600*y-0.25472*x)
        lloop=q>0.27597
        if (lloop) &
            lloop=(q>0.27846).or.(v**2>-4.0*log(u)*u**2)
      enddo
!
      a=v/u
!
    endsubroutine normal_deviate
!***********************************************************************
    subroutine random_seed_wrapper(size,put,get,channel)
!
!  Mimics the f90 random_seed routine.
!
      integer, optional :: size
      integer, optional :: channel
      integer, optional, dimension(:) :: put,get
!
      real :: dummy
      integer :: nseed
!
      select case (random_gen)
!
      case ('system')
        call random_seed(SIZE=nseed)
        if (present(size)) size=nseed
        if (present(get))  call random_seed(GET=get(1:nseed))
        if (present(put))  call random_seed(PUT=put(1:nseed))
      case ('min_std')
        nseed=1
        if (present(size)) size=nseed
        if (present(get))  get=rstate(1)
        if (present(put))  rstate(1)=put(1)
      case default ! 'nr_f90'
        nseed=2
        if (present(size)) size=nseed
        if (present(get)) get(1:nseed)=rstate(1:nseed)
        if (present(put)) then
          if (put(2)==0) then   ! state cannot be result from previous
                                ! call, so initialize

            if (present(channel)) then
              dummy = mars_ran2(put(1))
            else
              dummy = mars_ran(put(1))
            endif
          else
            rstate(1:nseed)=put(1:nseed)
          endif
        endif
!
      endselect
!
    endsubroutine random_seed_wrapper
!***********************************************************************
    function ran0(dummy)
!
!  The 'Minimal Standard' random number generator
!  by Lewis, Goodman and Miller.
!
!  28.08.02/nils: Adapted from Numerical Recipes
!
      integer, intent(inout) :: dummy
!
      integer :: k
      integer, parameter :: ia=16807,im=2147483647,iq=127773,ir=2836, &
           mask=123459876
      real, parameter :: am=1./im
      real :: ran0
!
      dummy=ieor(dummy,mask)
      k=dummy/iq
      dummy=ia*(dummy-k*iq)-ir*k
      if (dummy<0) dummy=dummy+im
      ran0=am*dummy
      dummy=ieor(dummy,mask)
!
    endfunction ran0
!***********************************************************************
    function mars_ran(init)
!
!  "Minimal" random number generator of Park and Miller, combined
!  with a Marsaglia shift sequence, with a period of supposedly
!  ~ 3.1x10^18.
!  Returns a uniform random number in (0, 1).
!  Call with (INIT=ival) to initialize.
!
!  26-sep-02/wolf: Implemented, following 'Numerical Recipes for F90'
!                  ran() routine
!
      implicit none
!
      integer, optional, intent(in) :: init
!
      real :: mars_ran
      real, save :: am   ! will be constant on a given platform
      integer, parameter :: ia=16807, im=2147483647, iq=127773, ir=2836
      integer :: k
      integer, save :: init1=1812
      logical, save :: first_call=.true.
!
!  Initialize.
!
      if (present(init) .or. (rstate(1) == 0) .or. (rstate(2) <= 0)) then
        if (present(init)) init1 = init
        am=nearest(1.0,-1.0)/im
        first_call=.false.
        rstate(1)=ieor(777755555,abs(init1))
        rstate(2)=ior(ieor(888889999,abs(init1)),1)
      elseif (first_call) then
        am=nearest(1.0,-1.0)/im
        first_call=.false.
      endif
!
!  Marsaglia shift sequence with period 2^32-1.
!
      rstate(1)=ieor(rstate(1),ishft(rstate(1),13))
      rstate(1)=ieor(rstate(1),ishft(rstate(1),-17))
      rstate(1)=ieor(rstate(1),ishft(rstate(1),5))
!
!  Park-Miller sequence by Schrage's method, period 2^31-2.
!
      k=rstate(2)/iq
      rstate(2)=ia*(rstate(2)-k*iq)-ir*k
      if (rstate(2) < 0) rstate(2)=rstate(2)+im
!
!  Combine the two generators with masking to ensure nonzero value.
!
      mars_ran=am*ior(iand(im,ieor(rstate(1),rstate(2))),1)
!
    endfunction mars_ran
!***********************************************************************
    function mars_ran2(init)
!
!  Same as mars_ran, but with another seed (for now)
!
!   8-nov-19/nils+axel: adapted from mars_ran
!
      implicit none
!
      integer, optional, intent(in) :: init
!
      real :: mars_ran2
      real, save :: am   ! will be constant on a given platform
      integer, parameter :: ia=16807, im=2147483647, iq=127773, ir=2836
      integer :: k
      integer, save :: init1=1812
      logical, save :: first_call=.true.
!
!  Initialize.
!
      if (present(init) .or. (rstate(1) == 0) .or. (rstate(2) <= 0)) then
        if (present(init)) init1 = init
        am=nearest(1.0,-1.0)/im
        first_call=.false.
        rstate(1)=ieor(777755555,abs(init1))
        rstate(2)=ior(ieor(888889999,abs(init1)),1)
      elseif (first_call) then
        am=nearest(1.0,-1.0)/im
        first_call=.false.
      endif
!
!  Marsaglia shift sequence with period 2^32-1.
!
      rstate(1)=ieor(rstate(1),ishft(rstate(1),13))
      rstate(1)=ieor(rstate(1),ishft(rstate(1),-17))
      rstate(1)=ieor(rstate(1),ishft(rstate(1),5))
!
!  Park-Miller sequence by Schrage's method, period 2^31-2.
!
      k=rstate(2)/iq
      rstate(2)=ia*(rstate(2)-k*iq)-ir*k
      if (rstate(2) < 0) rstate(2)=rstate(2)+im
!
!  Combine the two generators with masking to ensure nonzero value.
!
      mars_ran2=am*ior(iand(im,ieor(rstate(1),rstate(2))),1)
!
    endfunction mars_ran2
!***********************************************************************
    function nr_ran(iseed1)
!
!  (More or less) original routine from `Numerical Recipes in F90'. Not
!  sure we are allowed to distribute this.
!
!  28-aug-02/wolf: Adapted from Numerical Recipes
!
      implicit none
!
      integer, parameter :: ikind=kind(888889999)
      integer(ikind), intent(inout) :: iseed1
      real :: nr_ran
!
!  "Minimal" random number generator of Park and Miller combined
!  with a Marsaglia shift sequence. Returns a uniform random deviate
!  between 0.0 and 1.0 (exclusive of the endpoint values). This fully
!  portable, scalar generator has the "traditional" (not Fortran 90)
!  calling sequence with a random deviate as the returned function
!  value: call with iseed1 a negative integer to initialize;
!  thereafter, do not alter iseed1 except to reinitialize. The period
!  of this generator is about 3.1x10^18.
!
      real, save :: am
      integer(ikind), parameter :: ia=16807,im=2147483647,iq=127773,ir=2836
      integer(ikind), save      :: ix=-1,iy=-1,k
!
      if (iseed1 <= 0 .or. iy < 0) then   ! Initialize.
        am=nearest(1.0,-1.0)/im
        iy=ior(ieor(888889999,abs(iseed1)),1)
        ix=ieor(777755555,abs(iseed1))
        iseed1=abs(iseed1)+1   ! Set iseed1 positive.
      endif
      ix=ieor(ix,ishft(ix,13))   ! Marsaglia shift sequence with period 2^32-1.
      ix=ieor(ix,ishft(ix,-17))
      ix=ieor(ix,ishft(ix,5))
      k=iy/iq   ! Park-Miller sequence by Schrage's method,
      iy=ia*(iy-k*iq)-ir*k   ! period 231 - 2.
      if (iy < 0) iy=iy+im
      nr_ran=am*ior(iand(im,ieor(ix,iy)),1) ! Combine the two generators with
!                                           ! masking to ensure nonzero value.
    endfunction nr_ran
!***********************************************************************
    subroutine keep_compiler_quiet_r(v1,v2,v3,v4)
!
!  Call this to avoid compiler warnings about unused variables.
!  Optional arguments allow for more variables of the same shape+type.
!
!  04-aug-06/wolf: coded
!
      real     :: v1, v2, v3, v4
      optional ::     v2, v3, v4
!
      if (ALWAYS_FALSE) then
        write(0,*) 'keep_compiler_quiet_r: Never got here...'
        print*,                  v1
        if (present(v2)) print*, v2
        if (present(v3)) print*, v3
        if (present(v4)) print*, v4
        STOP 1
      endif
!
    endsubroutine keep_compiler_quiet_r
!***********************************************************************
    subroutine keep_compiler_quiet_r1d(v1,v2,v3,v4)
!
!  Call this to avoid compiler warnings about unused variables.
!  Optional arguments allow for more variables of the same shape+type.
!
!  04-aug-06/wolf: coded
!
      real, dimension(:) :: v1, v2, v3, v4
      optional           ::     v2, v3, v4
!
      if (ALWAYS_FALSE) then
        write(0,*) 'keep_compiler_quiet_r1d: Never got here...'
        print*,                  minval(v1)
        if (present(v2)) print*, minval(v2)
        if (present(v3)) print*, minval(v3)
        if (present(v4)) print*, minval(v4)
        STOP 1
      endif
!
    endsubroutine keep_compiler_quiet_r1d
!***********************************************************************
    subroutine keep_compiler_quiet_r2d(v1,v2,v3,v4)
!
!  Call this to avoid compiler warnings about unused variables.
!  Optional arguments allow for more variables of the same shape+type.
!
!  04-aug-06/wolf: coded
!
      real, dimension(:,:) :: v1, v2, v3, v4
      optional             ::     v2, v3, v4
!
      if (ALWAYS_FALSE) then
        write(0,*) 'keep_compiler_quiet_r2d: Never got here...'
        print*,                  minval(v1)
        if (present(v2)) print*, minval(v2)
        if (present(v3)) print*, minval(v3)
        if (present(v4)) print*, minval(v4)
        STOP 1
      endif
!
    endsubroutine keep_compiler_quiet_r2d
!***********************************************************************
    subroutine keep_compiler_quiet_r3d(v1,v2,v3,v4)
!
!  Call this to avoid compiler warnings about unused variables.
!  Optional arguments allow for more variables of the same shape+type.
!
!  04-aug-06/wolf: coded
!
      real, dimension(:,:,:) :: v1, v2, v3, v4
      optional               ::     v2, v3, v4
!
      if (ALWAYS_FALSE) then
        write(0,*) 'keep_compiler_quiet_r3d: Never got here...'
        print*,                  minval(v1)
        if (present(v2)) print*, minval(v2)
        if (present(v3)) print*, minval(v3)
        if (present(v4)) print*, minval(v4)
        STOP 1
      endif
!
    endsubroutine keep_compiler_quiet_r3d
!***********************************************************************
    subroutine keep_compiler_quiet_r4d(v1,v2,v3,v4)
!
!  Call this to avoid compiler warnings about unused variables.
!  Optional arguments allow for more variables of the same shape+type.
!
!  04-aug-06/wolf: coded
!
      real, dimension(:,:,:,:) :: v1, v2, v3, v4
      optional                 ::     v2, v3, v4
!
      if (ALWAYS_FALSE) then
        write(0,*) 'keep_compiler_quiet_r4d: never got here...'
        print*,                  minval(v1)
        if (present(v2)) print*, minval(v2)
        if (present(v3)) print*, minval(v3)
        if (present(v4)) print*, minval(v4)
        STOP 1
      endif
!
    endsubroutine keep_compiler_quiet_r4d
!***********************************************************************
    subroutine keep_compiler_quiet_p(v1,v2,v3,v4)
!
!  Call this to avoid compiler warnings about unused variables.
!  Optional arguments allow for more variables of the same shape+type.
!
!  04-aug-06/wolf: coded
!
      type (pencil_case) :: v1, v2, v3, v4
      optional           ::     v2, v3, v4
!
      if (ALWAYS_FALSE) then
        write(0,*) 'keep_compiler_quiet_p: Never got here...'
        print*,                  v1
        if (present(v2)) print*, v2
        if (present(v3)) print*, v3
        if (present(v4)) print*, v4
        STOP 1
      endif
!
    endsubroutine keep_compiler_quiet_p
!***********************************************************************
    subroutine keep_compiler_quiet_bc(v1,v2,v3,v4)
!
!  Call this to avoid compiler warnings about unused variables.
!  Optional arguments allow for more variables of the same shape+type.
!
!  04-aug-06/wolf: coded
!
      type (boundary_condition) :: v1, v2, v3, v4
      optional                  ::     v2, v3, v4
!
      if (ALWAYS_FALSE) then
        write(0,*) 'keep_compiler_quiet_bc: Never got here...'
        print*,                  v1
        if (present(v2)) print*, v2
        if (present(v3)) print*, v3
        if (present(v4)) print*, v4
        STOP 1
      endif
!
    endsubroutine keep_compiler_quiet_bc
!***********************************************************************
    subroutine keep_compiler_quiet_sl(v1,v2,v3,v4)
!
!  Call this to avoid compiler warnings about unused variables.
!  Optional arguments allow for more variables of the same shape+type.
!
!  04-aug-06/wolf: coded
!
      type (slice_data) :: v1, v2, v3, v4
      optional          ::     v2, v3, v4
!
      if (ALWAYS_FALSE) then
        write(0,*) 'keep_compiler_quiet_sl: Never got here...'
        print*,                  v1%index
        if (present(v2)) print*, v2%index
        if (present(v3)) print*, v3%index
        if (present(v4)) print*, v4%index
        STOP 1
      endif
!
    endsubroutine keep_compiler_quiet_sl
!***********************************************************************
    subroutine keep_compiler_quiet_i(v1,v2,v3,v4)
!
!  Call this to avoid compiler warnings about unused variables.
!  Optional arguments allow for more variables of the same shape+type.
!
!  04-aug-06/wolf: coded
!
      integer  :: v1, v2, v3, v4
      optional ::     v2, v3, v4
!
      if (ALWAYS_FALSE) then
        write(0,*) 'keep_compiler_quiet_i: Never got here...'
        print*,                  v1
        if (present(v2)) print*, v2
        if (present(v3)) print*, v3
        if (present(v4)) print*, v4
        STOP 1
      endif
!
    endsubroutine keep_compiler_quiet_i
!***********************************************************************
    subroutine keep_compiler_quiet_i81d(v1,v2,v3,v4)
!
!  Call this to avoid compiler warnings about unused variables.
!  Optional arguments allow for more variables of the same shape+type.
!
!  04-aug-06/wolf: coded
!
      integer(KIND=ikind8), dimension(:) :: v1, v2, v3, v4
      optional :: v2, v3, v4
!
      if (ALWAYS_FALSE) then
        write(0,*) 'keep_compiler_quiet_i: Never got here...'
        print*,                  v1(1)
        if (present(v2)) print*, v2(1)
        if (present(v3)) print*, v3(1)
        if (present(v4)) print*, v4(1)
        STOP 1
      endif
!
    endsubroutine keep_compiler_quiet_i81d
!***********************************************************************
    subroutine keep_compiler_quiet_i1d(v1,v2,v3,v4)
!
!  Call this to avoid compiler warnings about unused variables.
!  Optional arguments allow for more variables of the same shape+type.
!
!  04-aug-06/wolf: coded
!
      integer, dimension(:)  :: v1, v2, v3, v4
      optional               ::     v2, v3, v4
!
      if (ALWAYS_FALSE) then
        write(0,*) 'keep_compiler_quiet_i1d: Never got here...'
        print*,                  v1(1)
        if (present(v2)) print*, v2(1)
        if (present(v3)) print*, v3(1)
        if (present(v4)) print*, v4(1)
        STOP 1
      endif
!
    endsubroutine keep_compiler_quiet_i1d
!***********************************************************************
    subroutine keep_compiler_quiet_i2d(v1,v2,v3,v4)
!
!  Call this to avoid compiler warnings about unused variables.
!  Optional arguments allow for more variables of the same shape+type.
!
!  04-aug-06/wolf: coded
!
      integer, dimension(:,:)  :: v1, v2, v3, v4
      optional                 ::     v2, v3, v4
!
      if (ALWAYS_FALSE) then
        write(0,*) 'keep_compiler_quiet_i2d: Never got here...'
        print*,                  v1(1,1)
        if (present(v2)) print*, v2(1,1)
        if (present(v3)) print*, v3(1,1)
        if (present(v4)) print*, v4(1,1)
        STOP 1
      endif
!
    endsubroutine keep_compiler_quiet_i2d
!***********************************************************************
    subroutine keep_compiler_quiet_i3d(v1,v2,v3,v4)
!
!  Call this to avoid compiler warnings about unused variables.
!  Optional arguments allow for more variables of the same shape+type.
!
!  04-aug-06/wolf: coded
!
      integer, dimension(:,:,:)  :: v1, v2, v3, v4
      optional                   ::     v2, v3, v4
!
      if (ALWAYS_FALSE) then
        write(0,*) 'keep_compiler_quiet_i3d: Never got here...'
        print*,                  v1(1,1,1)
        if (present(v2)) print*, v2(1,1,1)
        if (present(v3)) print*, v3(1,1,1)
        if (present(v4)) print*, v4(1,1,1)
        STOP 1
      endif
!
    endsubroutine keep_compiler_quiet_i3d
!***********************************************************************
    subroutine keep_compiler_quiet_l1d(v1,v2,v3,v4)
!
!  Call this to avoid compiler warnings about unused variables.
!  Optional arguments allow for more variables of the same shape+type.
!
!  04-aug-06/wolf: coded
!
      logical, dimension(:)  :: v1, v2, v3, v4
      optional               ::     v2, v3, v4
!
      if (ALWAYS_FALSE) then
        write(0,*) 'keep_compiler_quiet_l1d: Never got here...'
        print*,                  v1(1)
        if (present(v2)) print*, v2(1)
        if (present(v3)) print*, v3(1)
        if (present(v4)) print*, v4(1)
        STOP 1
      endif
!
    endsubroutine keep_compiler_quiet_l1d
!***********************************************************************
    subroutine keep_compiler_quiet_l(v1,v2,v3,v4)
!
!  Call this to avoid compiler warnings about unused variables.
!  Optional arguments allow for more variables of the same shape+type.
!
!  04-aug-06/wolf: coded
!
      logical  :: v1, v2, v3, v4
      optional ::     v2, v3, v4
!
      if (ALWAYS_FALSE) then
        write(0,*) 'keep_compiler_quiet_l: Never got here...'
        print*,                  v1
        if (present(v2)) print*, v2
        if (present(v3)) print*, v3
        if (present(v4)) print*, v4
        STOP 1
      endif
!
    endsubroutine keep_compiler_quiet_l
!***********************************************************************
    subroutine keep_compiler_quiet_c(v1,v2,v3,v4)
!
!  Call this to avoid compiler warnings about unused variables.
!  Optional arguments allow for more variables of the same shape+type.
!
!  04-aug-06/wolf: coded
!
      character (len=*) :: v1, v2, v3, v4
      optional          ::     v2, v3, v4
!
      if (ALWAYS_FALSE) then
        write(0,*) 'keep_compiler_quiet_c: Never got here...'
        print*,                  v1
        if (present(v2)) print*, v2
        if (present(v3)) print*, v3
        if (present(v4)) print*, v4
        STOP 1
      endif
!
    endsubroutine keep_compiler_quiet_c
!***********************************************************************
    integer function count_bits(n)
!
!  Count non-zero bits. Returns the number of the left-most non-zero bit.
!
!  18-jun-2013/Bourdin.KIS: coded
!
      integer, intent(in) :: n
!
      integer :: num
!
      num = n
      count_bits = 0
      do while (num /= 0)
        num = ishft (num, -1)
        count_bits = count_bits + 1
      enddo
!
    endfunction count_bits
!***********************************************************************
    character (len=intlen) function itoa(n)
!
!  Convert integer to ASCII, similar to the C-stdlib 'itoa' function.
!  If the length of an integer changes, please adapt 'intlen' accordingly.
!
!  05-aug-2011/Bourdin.KIS: coded
!
      integer, intent(in) :: n
!
      write (itoa, '(I21)') n   ! (64 bit integer plus sign)
      itoa = adjustl(itoa)
!
    endfunction itoa
!***********************************************************************
  character function intochar(i)
!
  integer, intent(in) :: i
!
  write(intochar,'(i1)') i
!
  endfunction intochar
!***********************************************************************
    subroutine chk_time(label,time1,time2)
!
      integer :: time1,time2,count_rate
      character (len=*) :: label
!
!  Prints time in seconds.
!
      call system_clock(count=time2,count_rate=count_rate)
      print*,"chk_time: ",label,(time2-time1)/real(count_rate)
      time1=time2
!
    endsubroutine chk_time
!***********************************************************************
    subroutine parse_filename(filename,dirpart,filepart)
!
!  Split full pathname of a file into directory part and local filename part.
!
!  02-apr-04/wolf: coded
!
      character (len=*) :: filename,dirpart,filepart
      integer :: i
!
      intent(in)  :: filename
      intent(out) :: dirpart,filepart
!
      i = index(filename,'/',BACK=.true.) ! search last slash
      if (i>0) then
        call safe_character_assign(dirpart,filename(1:i-1))
        call safe_character_assign(filepart,trim(filename(i+1:)))
      else
        call safe_character_assign(dirpart,'.')
        call safe_character_assign(filepart,trim(filename))
      endif
!
    endsubroutine parse_filename
!***********************************************************************
    subroutine safe_character_assign(dest,src)
!
!  Do character string assignement with check against overflow.
!
!  08-oct-02/tony: coded
!  25-oct-02/axel: added tag in output to give name of routine
!
      character (len=*), intent(in):: src
      character (len=*), intent(out):: dest
      integer :: destLen, srcLen
!
      destLen = len(dest)
      srcLen = len(src)
!
      if (destLen<srcLen) then
        print*, "safe_character_assign: ", &
            "RUNTIME ERROR: FORCED STRING TRUNCATION WHEN ASSIGNING '" &
             //src//"' to ",destLen," characters"
        dest=src(1:destLen)
      else
        dest=src
      endif
!
    endsubroutine safe_character_assign
!***********************************************************************
    subroutine safe_character_append_2(str1,str2)
!
!  08-oct-02/wolf: coded
!
      character (len=*), intent(inout):: str1
      character (len=*), intent(in):: str2
!
      call safe_character_assign(str1, trim(str1) // str2)
!
    endsubroutine safe_character_append_2
!***********************************************************************
    subroutine safe_character_prepend_2(str1,str2)
!
!  23-feb-12/bing: adapted from safe_character_append_2
!
      character (len=*), intent(inout):: str1
      character (len=*), intent(in):: str2
!
      call safe_character_assign(str1, trim(str2) // trim(str1))
!
    endsubroutine safe_character_prepend_2
!***********************************************************************
    subroutine safe_character_append_3(str1,str2,str3)
!
!  08-oct-02/wolf: coded
!
      character (len=*), intent(inout):: str1
      character (len=*), intent(in):: str2,str3
!
      call safe_character_assign(str1, trim(str1) // trim(str2) // trim(str3))
!
    endsubroutine safe_character_append_3
!***********************************************************************
    function lower_case(input)
!
!  27-Sep-2015/PABourdin: coded
!
      character(*), intent(in) :: input
      character(len(input)) :: lower_case
!
      character(*), parameter :: lower_chars = 'abcdefghijklmnopqrstuvwxyz'
      character(*), parameter :: upper_chars = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
      integer pos, ind
!
      lower_case = input
      do pos = 1, len (input)
        ind = index (upper_chars, lower_case(pos:pos))
        if (ind /= 0) lower_case(pos:pos) = lower_chars(ind:ind)
      enddo
!
    endfunction lower_case
!***********************************************************************
    subroutine find_index_range(aa,naa,aa1,aa2,ii1,ii2)
!
!  Find index range (ii1,ii2) such that aa
!
!   9-mar-08/axel: coded
!
      use Cparam, only: nghost
!
      integer :: naa,ii,ii1,ii2
      real, dimension (naa) :: aa
      real :: aa1,aa2
!
      intent(in)  :: aa,naa,aa1,aa2
      intent(out) :: ii1,ii2
!
!  If not extent in this direction, set indices to interior values.
!
      if (naa==2*nghost+1) then
        ii1=nghost+1
        ii2=nghost+1
        return
      endif
!
!  Find lower index.
!
      ii1=naa
      do ii=1,naa
        if (aa(ii)>=aa1) then
          ii1=ii
          exit
        endif
      enddo
!
!  Find upper index.
!
      ii2=1
      do ii=naa,1,-1
        if (aa(ii)<=aa2) then
          ii2=ii
          return
        endif
      enddo
!
    endsubroutine find_index_range
!***********************************************************************
    function find_index_range_hill(aa,planes) result(range)
!
!  Find index range (ii1,ii2) if a "hill" between to "planes" with levels planes(1)
!  and planes(2) (possibly different). Inner "valleys" of the hill to the
!  planes' levels are not considered.
!
!   9-mar-16/MR: derived from find_index_range
!
      integer, dimension(:), intent(in) :: aa
      integer, dimension(2), intent(in) :: planes
      integer, dimension(2) :: range
!
      integer :: naa,ii
!
!  Find lower index.
!
      naa=size(aa)
      range(1)=naa
      do ii=1,naa
        if (aa(ii)/=planes(1)) then
          range(1)=ii
          exit
        endif
      enddo
!
!  Find upper index.
!
      range(2)=1
      do ii=naa,1,-1
        if (aa(ii)/=planes(2)) then
          range(2)=ii
          return
        endif
      enddo
!
    endfunction find_index_range_hill
!***********************************************************************
    pure integer function find_index(xa, x)
!
!  Returns the index of the element in array xa that is closest to x.
!
!  24-feb-13/ccyang: coded
!
      real, dimension(:), intent(in) :: xa
      real, intent(in) :: x
!
      integer, dimension(1) :: closest
!
      closest = minloc(abs(x - xa))
      find_index = closest(1)
!
    endfunction find_index
!***********************************************************************
    function spline_derivative(z,f)
!
!  Computes derivative of a given function using spline interpolation.
!
!  01-apr-03/tobi: originally coded by Aake Nordlund
!
      implicit none
      real, dimension(:) :: z
      real, dimension(:), intent(in) :: f
      real, dimension(size(z)) :: w1,w2,w3
      real, dimension(size(z)) :: d,spline_derivative
      real :: c
      integer :: mz,k
!
      mz=size(z)
!
      w1(1)=1./(z(2)-z(1))**2
      w3(1)=-1./(z(3)-z(2))**2
      w2(1)=w1(1)+w3(1)
      d(1)=2.*((f(2)-f(1))/(z(2)-z(1))**3 &
                  -(f(3)-f(2))/(z(3)-z(2))**3)
!
!  Interior points.
!
      w1(2:mz-1)=1./(z(2:mz-1)-z(1:mz-2))
      w3(2:mz-1)=1./(z(3:mz)-z(2:mz-1))
      w2(2:mz-1)=2.*(w1(2:mz-1)+w3(2:mz-1))
!
      d(2:mz-1)=3.*(w3(2:mz-1)**2*(f(3:mz)-f(2:mz-1)) &
           +w1(2:mz-1)**2*(f(2:mz-1)-f(1:mz-2)))
!
!  Last point.
!
      w1(mz)=1./(z(mz-1)-z(mz-2))**2
      w3(mz)=-1./(z(mz)-z(mz-1))**2
      w2(mz)=w1(mz)+w3(mz)
      d(mz)=2.*((f(mz-1)-f(mz-2))/(z(mz-1)-z(mz-2))**3 &
           -(f(mz)-f(mz-1))/(z(mz)-z(mz-1))**3)
!
!  Eliminate at first point.
!
      c=-w3(1)/w3(2)
      w1(1)=w1(1)+c*w1(2)
      w2(1)=w2(1)+c*w2(2)
      d(1)=d(1)+c*d(2)
      w3(1)=w2(1)
      w2(1)=w1(1)
!
!  Eliminate at last point.
!
      c=-w1(mz)/w1(mz-1)
      w2(mz)=w2(mz)+c*w2(mz-1)
      w3(mz)=w3(mz)+c*w3(mz-1)
      d(mz)=d(mz)+c*d(mz-1)
      w1(mz)=w2(mz)
      w2(mz)=w3(mz)
!
!  Eliminate subdiagonal.
!
      do k=2,mz
        c=-w1(k)/w2(k-1)
        w2(k)=w2(k)+c*w3(k-1)
        d(k)=d(k)+c*d(k-1)
      enddo
!
!  Backsubstitute.
!
      d(mz)=d(mz)/w2(mz)
      do k=mz-1,1,-1
        d(k)=(d(k)-w3(k)*d(k+1))/w2(k)
      enddo
!
      spline_derivative=d
!
    endfunction spline_derivative
!***********************************************************************
    function spline_derivative_double(z,f)
!
!  Computes derivative of a given function using spline interpolation.
!
!  01-apr-03/tobi: originally coded by Aake Nordlund
!  11-apr-03/axel: double precision version
!
      implicit none
      real, dimension(:) :: z
      double precision, dimension(:), intent(in) :: f
      double precision, dimension(size(z)) :: w1,w2,w3
      double precision, dimension(size(z)) :: d,spline_derivative_double
      double precision :: c
      integer :: mz,k
!
      mz=size(z)
!
      w1(1)=1./(z(2)-z(1))**2
      w3(1)=-1./(z(3)-z(2))**2
      w2(1)=w1(1)+w3(1)
      d(1)=2.*((f(2)-f(1))/(z(2)-z(1))**3 &
                  -(f(3)-f(2))/(z(3)-z(2))**3)
!
!  Interior points.
!
      w1(2:mz-1)=1./(z(2:mz-1)-z(1:mz-2))
      w3(2:mz-1)=1./(z(3:mz)-z(2:mz-1))
      w2(2:mz-1)=2.*(w1(2:mz-1)+w3(2:mz-1))
!
      d(2:mz-1)=3.*(w3(2:mz-1)**2*(f(3:mz)-f(2:mz-1)) &
           +w1(2:mz-1)**2*(f(2:mz-1)-f(1:mz-2)))
!
!  Last point.
!
      w1(mz)=1./(z(mz-1)-z(mz-2))**2
      w3(mz)=-1./(z(mz)-z(mz-1))**2
      w2(mz)=w1(mz)+w3(mz)
      d(mz)=2.*((f(mz-1)-f(mz-2))/(z(mz-1)-z(mz-2))**3 &
           -(f(mz)-f(mz-1))/(z(mz)-z(mz-1))**3)
!
!  Eliminate at first point.
!
      c=-w3(1)/w3(2)
      w1(1)=w1(1)+c*w1(2)
      w2(1)=w2(1)+c*w2(2)
      d(1)=d(1)+c*d(2)
      w3(1)=w2(1)
      w2(1)=w1(1)
!
!  Eliminate at last point.
!
      c=-w1(mz)/w1(mz-1)
      w2(mz)=w2(mz)+c*w2(mz-1)
      w3(mz)=w3(mz)+c*w3(mz-1)
      d(mz)=d(mz)+c*d(mz-1)
      w1(mz)=w2(mz)
      w2(mz)=w3(mz)
!
!  Eliminate subdiagonal.
!
      do k=2,mz
        c=-w1(k)/w2(k-1)
        w2(k)=w2(k)+c*w3(k-1)
        d(k)=d(k)+c*d(k-1)
      enddo
!
!  Backsubstitute.
!
      d(mz)=d(mz)/w2(mz)
      do k=mz-1,1,-1
        d(k)=(d(k)-w3(k)*d(k+1))/w2(k)
      enddo
!
      spline_derivative_double=d
!
    endfunction spline_derivative_double
!***********************************************************************
    function spline_integral(z,f,q0) result (q)
!
!  Computes integral of a given function using spline interpolation.
!
!  01-apr-03/tobi: originally coded by Aake Nordlund
!  02-jun-15/MR: fix for size of f equal to 1.
!
      implicit none
      real, dimension(:) :: z
      real, dimension(:) :: f
      real, dimension(size(z)) :: df,dz,q
      real, optional :: q0
      integer :: mz,k
!
      mz=size(z)
      if (mz==1) then
        q = f(1)
        return
      endif
!
      q(1)=0.
      if (present(q0)) q(1)=q0
      df=spline_derivative(z,f)
      dz(2:mz)=z(2:mz)-z(1:mz-1)
!
      q(2:mz)=.5*dz(2:mz)*(f(1:mz-1)+f(2:mz)) &
              +(1./12.)*dz(2:mz)**2*(df(1:mz-1)-df(2:mz))
!
      do k=2,mz
        q(k)=q(k)+q(k-1)
      enddo
!
    endfunction spline_integral
!***********************************************************************
    function spline_integral_double(z,f,q0)
!
!  Computes integral of a given function using spline interpolation.
!
!  01-apr-03/tobi: originally coded by Aake Nordlund
!  11-apr-03/axel: double precision version
!
      implicit none
      real, dimension(:) :: z
      double precision, dimension(:) :: f
      real, dimension(size(z)) :: dz
      double precision, dimension(size(z)) :: df
      double precision, dimension(size(z)) :: q,spline_integral_double
      double precision, optional :: q0
      integer :: mz,k
!
      mz=size(z)
!
      q(1)=0.
      if (present(q0)) q(1)=q0
      df=spline_derivative_double(z,f)
      dz(2:mz)=z(2:mz)-z(1:mz-1)
!
      q(2:mz)=.5*dz(2:mz)*(f(1:mz-1)+f(2:mz)) &
              +(1./12.)*dz(2:mz)**2*(df(1:mz-1)-df(2:mz))
!
      do k=2,mz
        q(k)=q(k)+q(k-1)
      enddo
!
      spline_integral_double=q
!
    endfunction spline_integral_double
!***********************************************************************
    pure subroutine tridag(a,b,c,r,u,err,msg)
!
!  Solves a tridiagonal system.
!
!  01-apr-03/tobi: from Numerical Recipes (p42-43).
!
!  | b1 c1 0  ...            | | u1   |   | r1   |
!  | a2 b2 c2 ...            | | u2   |   | r2   |
!  | 0  a3 b3 c3             | | u3   | = | r3   |
!  |          ...            | | ...  |   | ...  |
!  |          an-1 bn-1 cn-1 | | un-1 |   | rn-1 |
!  |          0    a_n  b_n  | | un   |   | rn   |
!
      implicit none
      real, dimension(:), intent(in) :: a,b,c,r
      real, dimension(:), intent(out) :: u
      real, dimension(size(b)) :: gam
      logical, intent(out), optional :: err
      character(len=*), intent(out), optional :: msg
      integer :: n,j
      real :: bet
!
      if (present(err)) err=.false.
      n=size(b)
      bet=b(1)
      if (bet==0.0) then
        if (present(msg)) msg = 'tridag: Error at code stage 1'
        if (present(err)) err = .true.
      endif
!
      u(1)=r(1)/bet
      do j=2,n
        gam(j)=c(j-1)/bet
        bet=b(j)-a(j)*gam(j)
        if (bet==0.0) then
          if (present(msg)) msg = 'tridag: Error at code stage 2'
          if (present(err)) err = .true.
          return
        endif
        u(j)=(r(j)-a(j)*u(j-1))/bet
      enddo
!
      do j=n-1,1,-1
        u(j)=u(j)-gam(j+1)*u(j+1)
      enddo
!
    endsubroutine tridag
!***********************************************************************
    subroutine tridag_double(a,b,c,r,u,err)
!
!  Solves tridiagonal system.
!
!  01-apr-03/tobi: from Numerical Recipes (p42-43).
!  11-apr-03/axel: double precision version.
!
!  | b1 c1 0  ...            | | u1   |   | r1   |
!  | a2 b2 c2 ...            | | u2   |   | r2   |
!  | 0  a3 b3 c3             | | u3   | = | r3   |
!  |          ...            | | ...  |   | ...  |
!  |          an-1 bn-1 cn-1 | | un-1 |   | rn-1 |
!  |          0    a_n  b_n  | | un   |   | rn   |
!
      implicit none
      double precision, dimension(:), intent(in) :: a,b,c,r
      double precision, dimension(:), intent(out) :: u
      double precision, dimension(size(b)) :: gam
      logical, intent(out), optional :: err
      integer :: n,j
      double precision :: bet
!
      if (present(err)) err=.false.
      n=size(b)
      bet=b(1)
      if (bet==0.0) then
        print*,'tridag_double: Error at code stage 1'
        if (present(err)) err=.true.
        return
      endif
!
      u(1)=r(1)/bet
      do j=2,n
        gam(j)=c(j-1)/bet
        bet=b(j)-a(j)*gam(j)
        if (bet==0.0) then
          print*,'tridag_double: Error at code stage 2'
          if (present(err)) err=.true.
          return
        endif
        u(j)=(r(j)-a(j)*u(j-1))/bet
      enddo
      do j=n-1,1,-1
        u(j)=u(j)-gam(j+1)*u(j+1)
      enddo
!
    endsubroutine tridag_double
!***********************************************************************
    subroutine pendag(n,a,b,c,d,e,r,u)
!
!  01-apr-00/John Crowe (Newcastle): written
!  10-avr-12/dintrans: recoded entirely because the old version did
!  not work (I realized that by inverting pentadiagonal systems with known
!  analytical solutions, e.g. laplacian of cos/sin functions)
!
      real, dimension(:), intent(in)  :: a,b,c,d,e,r
      real, dimension(:), intent(out) :: u
      real, dimension(size(r)+1) :: w,beta,alpha,cg,h
      integer :: k,n
!
      w(1)=c(1)
      beta(1)=0.0
      beta(2)=d(1)/w(1)
      alpha(1)=0.0
      alpha(2)=e(1)/w(1)
      alpha(n)=0.0
      alpha(n+1)=0.0
!
      do k=2,n
        cg(k)=b(k)-a(k)*beta(k-1)
        w(k)=c(k)-a(k)*alpha(k-1)-cg(k)*beta(k)
        if (w(k)==0.0) then
          print*,"w(k)=0.0 in pendag"
          stop
        endif
        beta(k+1)=(d(k)-cg(k)*alpha(k))/w(k)
        alpha(k+1)=e(k)/w(k)
      enddo
!
      h(1)=0.0
      h(2)=r(1)/w(1)
      do k=2,n
        h(k+1)=(r(k)-a(k)*h(k-1)-cg(k)*h(k))/w(k)
      enddo
!
      u(n)=h(n+1)
      u(n-1)=h(n)-beta(n)*u(n)
      do k=n-2,1,-1
        u(k)=h(k+1)-beta(k+1)*u(k+1)-alpha(k+1)*u(k+2)
      enddo
!
    endsubroutine pendag
!***********************************************************************
    pure subroutine spline(arrx,arry,x2,S,psize1,psize2,err,msg)
!
!  Interpolates in x2 a natural cubic spline with knots defined by the 1d
!  arrays arrx and arry.
!
!  25-mar-05/wlad : coded
!
      integer, intent(in) :: psize1,psize2
      integer :: i,j,ct1,ct2
      real, dimension (psize1) :: arrx,arry,h,h1,a,b,d,sol
      real, dimension (psize2) :: x2,S
      real, parameter :: fac=0.1666666
      logical, intent(out), optional :: err
      character(len=*), intent(out), optional :: msg
!
      intent(in)  :: arrx,arry,x2
      intent(out) :: S
!
      if (present(err)) err=.false.
      ct1 = psize1
      ct2 = psize2
!
!  Short-circuit if there is only 1 knot.
!
      if (ct1 == 1) then
        S = arry(1)
        return
      endif
!
!  Breaks if x is not monotonically increasing.
!
      do i=1,ct1-1
        if (arrx(i+1)<=arrx(i)) then
          if (present(msg)) msg = 'spline x:y in x2:y2 : vector x is not monotonically increasing'
          if (present(err)) err = .true.
          return
        endif
      enddo
!
!  Step h.
!
      h(1:ct1-1) = arrx(2:ct1) - arrx(1:ct1-1)
      h(ct1) = h(ct1-1)
      h1=1./h
!
!  Coefficients for tridiagonal system.
!
      a(2:ct1) = h(1:ct1-1)
      a(1) = a(2)
!
      b(2:ct1) = 2*(h(1:ct1-1) + h(2:ct1))
      b(1) = b(2)
!
      !c = h
!
      d(2:ct1-1) = 6*((arry(3:ct1) - arry(2:ct1-1))*h1(2:ct1-1) - (arry(2:ct1-1) - arry(1:ct1-2))*h1(1:ct1-2))
      d(1) = 0. ; d(ct1) = 0.
!
      if (present(msg)) then
        call tridag(a,b,h,d,sol,err,msg)
      else
        call tridag(a,b,h,d,sol,err)
      endif
      if (err) return
!
!  Interpolation formula.
!
      do j=1,ct2
        do i=1,ct1-1
!
          if ((x2(j)>=arrx(i)).and.(x2(j)<=arrx(i+1))) then
!
!  Substitute 1/6. by 0.1666666 to avoid divisions.
!
            S(j) = (fac*h1(i)) * (sol(i+1)*(x2(j)-arrx(i))**3 + sol(i)*(arrx(i+1) - x2(j))**3)  + &
                    (x2(j) - arrx(i))*(arry(i+1)*h1(i) - h(i)*sol(i+1)*fac)                          + &
                    (arrx(i+1) - x2(j))*(arry(i)*h1(i) - h(i)*sol(i)*fac)
           endif
!
        enddo
!
!  Use border values beyond this interval - should perhaps allow for linear
!  interpolation.
!
        if (x2(j)<=arrx(1)) then
          S(j) = arry(1)
        elseif (x2(j)>=arrx(ct1)) then
          S(j) = arry(ct1)
        endif
      enddo
!
    endsubroutine spline
!***********************************************************************
    subroutine cspline(y1, x2, y2, nonperiodic, tvd, posdef)
!
!  Use cubic spline interpolants to interpolate y1 into y2 at x2,
!  assuming y1(i) is at x = i-1 for all i.
!  If parameter nonperiodic is present and is .true., the natural spline
!  is used; periodic boundary conditions are assumed, otherwise.
!  If parameter tvd is present and is .true., those interpolated points
!  that violate the TVD condition are downgraded to linear interpolation.
!  If parameter posdef is present and is .true., those interpolated
!  points that dip below zero are downgraded to linear interpolation.
!
!  16-dec-14/ccyang: coded.
!
      real, dimension(:), intent(in) :: y1, x2
      real, dimension(size(x2)), intent(out) :: y2
      logical, intent(in), optional :: nonperiodic, tvd, posdef
!
      integer, dimension(size(x2)) :: inode
      real, dimension(size(y1)) :: b, c, d
      logical :: lperiodic, ltvd, lposdef
      integer :: n1, i, n
      real :: a, s
!
!  Check the switches.
!
      lperiodic = .true.
      if (present(nonperiodic)) lperiodic = .not. nonperiodic
      ltvd = .false.
      if (present(tvd)) ltvd = tvd
      lposdef = .false.
      if (present(posdef)) lposdef = posdef
!
!  Find the interpolants.
!
      if (lperiodic) then
        call spline_coeff_periodic(y1, b, c, d)
      else
        call spline_coeff(y1, b, c, d)
      endif
!
!  Interpolate.
!
      n1 = size(y1)
      inode = floor(x2)
      y2 = x2 - real(inode)
      interp: if (lperiodic) then
        inode = modulo(inode, n1) + 1
        y2 = y1(inode) + (b(inode) + (c(inode) + d(inode) * y2) * y2) * y2
      else interp
        inode = inode + 1
        n = n1 - 1
        s = b(n) + 2.0 * c(n) + 3.0 * d(n)
        where (inode < 1)
          y2 = y1(1) + b(1) * x2
        elsewhere (inode > n)
          y2 = y1(n1) + s * (x2 - real(n))
        elsewhere
          y2 = y1(inode) + (b(inode) + (c(inode) + d(inode) * y2) * y2) * y2
        endwhere
      endif interp
!
!  Check TVD and/or positive definiteness.
!
      tvdpd: if (ltvd .or. lposdef) then
        scan: do i = 1, size(x2)
          n = inode(i)
          if (.not. lperiodic .and. (n < 1 .or. n >= n1)) then
            cycle scan  ! Allow for extrapolation.
          elseif (lperiodic .and. n >= n1) then
            a = y1(1)
          else
            a = y1(n + 1)
          endif
          downgrade: if (ltvd .and. (y1(n) - y2(i)) * (y2(i) - a) < 0.0 .or. &
                         lposdef .and. y2(i) < 0.0) then
!           Downgrade to linear interpolation
            s = x2(i) - floor(x2(i))
            y2(i) = y1(n) * (1.0 - s) + a * s
          endif downgrade
        enddo scan
      endif tvdpd
!
    endsubroutine cspline
!***********************************************************************
    pure subroutine spline_coeff(y, b, c, d, yp1, yp2, err, msg)
!
!  Calculates the coefficients of the cubic spline interpolants
!
!    S_i(x) = y(i) + b(i)*(x-i) + c(i)*(x-i)^2 + d(i)*(x-i)^3,
!
!  for x in [i, i+1] and with S_i(i) = y(i) and S_i(i+1) = y(i+1), where
!  i = 1, 2, ..., n-1.
!
!  If yp1 is present, S'(1) = yp1; S''(1) = 0, otherwise.
!  If yp2 is present, S'(n) = yp2; S''(n) = 0, otherwise.
!
!  15-dec-14/ccyang: coded.
!
      real, dimension(:), intent(in) :: y
      real, dimension(size(y)), intent(out) :: b, c, d
      real, intent(in), optional :: yp1, yp2
      character(len=*), intent(out), optional :: msg
      logical, intent(out), optional :: err
!
      real, parameter :: onethird = 1.0 / 3.0
      real, dimension(size(y)) :: p, q, r
      character(len=linelen) :: msg1
      logical :: err1
      integer :: n
!
!  Solve the tridiagonal system for coefficient c's.
!
      n = size(y)
      p = 4.0
      q = 1.0
      r(2:n-1) = y(3:n) - 2.0 * y(2:n-1) + y(1:n-2)
!
!     Left boundary condition.
!
      left: if (present(yp1)) then
        p(1) = 2.0
        r(1) = y(2) - y(1) - yp1
      else left
        q(1) = 0.0
        r(1) = 0.0
      endif left
!
!     Right boundary condition
!
      right: if (present(yp2)) then
        p(n) = 2.0
        r(n) = yp2 - y(n) + y(n-1)
      else right
        q(n) = 0.0
        r(n) = 0.0
      endif right
!
      r = 3.0 * r
      call tridag(q, p, q, r, c, err1, msg1)
      if (present(err)) err = err1
      if (present(msg)) msg = msg1
!
!  Find the rest of the coefficients.
!
      b(1:n-1) = y(2:n) - y(1:n-1) - onethird * (c(2:n) + 2.0 * c(1:n-1))
      d(1:n-1) = onethird * (c(2:n) - c(1:n-1))
!
    endsubroutine spline_coeff
!***********************************************************************
    subroutine spline_coeff_periodic(y, b, c, d)
!
!  Calculates the coefficients of the cubic spline interpolants with
!  periodic boundary conditions.  The nodes start from one and the step
!  size is assumed to be unity.  The interpolants are
!
!    S_i(x) = y_i + b_i (x - i) + c_i (x - i)^2 + d_i (x - i)^3,
!
!  for x in [i, i+1] and i = 1, 2, ..., n.
!
!  13-dec-14/ccyang: coded.
!
      real, dimension(:), intent(in) :: y
      real, dimension(size(y)), intent(out) :: b, c, d
!
      real, parameter :: onethird = 1.0 / 3.0
      real, dimension(size(y)) :: r
      integer :: n
!
!  Solve the cyclic system for coefficient c's.
!
      n = size(y)
      r(2:n-1) = y(3:n) - 2.0 * y(2:n-1) + y(1:n-2)
      r(1) = y(2) - 2.0 * y(1) + y(n)
      r(n) = y(1) - 2.0 * y(n) + y(n-1)
      r = 3.0 * r
      call cyclic(spread(1.0, 1, n), spread(4.0, 1, n), spread(1.0, 1, n), 1.0, 1.0, r, c, n)
!
!  Find the rest of the coefficients.
!
      b(1:n-1) = y(2:n) - y(1:n-1) - onethird * (c(2:n) + 2.0 * c(1:n-1))
      d(1:n-1) = onethird * (c(2:n) - c(1:n-1))
      b(n) = y(1) - y(n) - onethird * (c(1) + 2.0 * c(n))
      d(n) = onethird * (c(1) - c(n))
!
    endsubroutine spline_coeff_periodic
!***********************************************************************
    subroutine poly_interp_one(xa, ya, x, y, istatus, message)
!
!  Uses polynomial interpolation to interpolate (xa, ya) to (x, y),
!  where y(n) is the n-th order interpolation, 0 <= n <= size(xa)-1.
!
!  01-oct-14/ccyang: adapted from Numerical Recipes.
!
      real, dimension(:), intent(in) :: xa, ya
      real, intent(in) :: x
      real, dimension(0:size(xa)-1), intent(out) :: y
      integer, intent(out), optional :: istatus
      character(len=*), intent(out), optional :: message
!
      real, dimension(size(xa)) :: c, d, den, ho
      integer :: m, n, ns
      real :: dy
!
!  Check the sizes of the input arrays.
!
      n = size(xa)
      incompatible: if (size(ya) /= n) then
        if (present(istatus)) istatus = -1
        if (present(message)) message = 'Input arrays xa and ya are incompatible. '
        return
      endif incompatible
!
!  Initialize the tableau of c's and d's.
!
      c = ya
      d = ya
      ho = xa - x
!
!  Find index ns of closest table entry.
!
      ns = find_index(xa, x)
!
!  Initial approximation to y.
!
      y(0) = ya(ns)
      ns = ns - 1
!
!  For each column of the tableau, loop over the current c's and d's and update them.
!
      update: do m = 1, n - 1
        den(1:n-m) = ho(1:n-m) - ho(1+m:n)
        failure: if (any(den(1:n-m) == 0.0)) then
          if (present(istatus)) istatus = -2
          if (present(message)) message = 'calculation failure'
          return
        endif failure
        den(1:n-m) = (c(2:n-m+1) - d(1:n-m)) / den(1:n-m)
        d(1:n-m) = ho(1+m:n) * den(1:n-m)
        c(1:n-m) = ho(1:n-m) * den(1:n-m)
        take_side: if (2 * ns < n - m) then
          dy = c(ns + 1)
        else take_side
          dy = d(ns)
          ns = ns - 1
        endif take_side
        y(m) = y(m-1) + dy
      enddo update
!
!  Clean exit
!
      if (present(istatus)) istatus = 0
!
    endsubroutine poly_interp_one
!***********************************************************************
    subroutine poly_interp_fixorder(xa, ya, x, y, norder, tvd, posdef, istatus, message)
!
!  Uses polynomials of norder to interpolate (xa, ya) to each of (x, y).
!  If tvd or posdef is present and set true, the order will be reduced
!  at places where the total variation diminishing or positive
!  definiteness is violated, respectively.
!
!  01-oct-14/ccyang: coded
!
      real, dimension(:), intent(in) :: xa, ya
      real, dimension(:), intent(in) :: x
      real, dimension(:), intent(out) :: y
      integer, intent(in) :: norder
      logical, intent(in), optional :: tvd, posdef
      integer, intent(out), optional :: istatus
      character(len=*), intent(out), optional :: message
!
      real, dimension(0:norder) :: yi
      character(len=256) :: msg
      logical :: tvd1, posdef1, left, lodd, ok
      integer :: ord, noh, nxa, ix, ix1, ix2
      integer :: i, n, istat
!
!  Check the dimension of the output arrays.
!
      n = size(x)
      incompatible: if (size(y) /= n) then
        if (present(istatus)) istatus = -3
        if (present(message)) message = 'Arrays x and y are incompatible.'
        return
      endif incompatible
!
!  Check if total variation diminishing or positive definiteness is turned on.
!
      tvd_on: if (present(tvd)) then
        tvd1 = tvd
      else tvd_on
        tvd1 = .false.
      endif tvd_on
!
      pos_on: if (present(posdef)) then
        posdef1 = posdef
      else pos_on
        posdef1 = .false.
      endif pos_on
!
!  Interpolate each point.
!
      istat = 0
      nxa = size(xa)
      noh = norder / 2
      lodd = mod(norder, 2) /= 0
!
      loop: do i = 1, n
!
!  Find the index range to construct the interpolant.
!
        ix = find_index(xa, x(i))
        left = x(i) < xa(ix)
        ix1 = ix - noh
        ix2 = ix + noh
        odd: if (lodd) then
          side: if (left) then
            ix1 = ix1 - 1
          else side
            ix2 = ix2 + 1
          endif side
        endif odd
        ix1 = max(ix1, 1)
        ix2 = min(ix2, nxa)
        ord = ix2 - ix1
!
!  Send for polynomial interpolation.
!
        call poly_interp_one(xa(ix1:ix2), ya(ix1:ix2), x(i), yi, istat, msg)
        if (istat /= 0) exit loop
!
!  Check the total variation and/or positive definiteness.
!
        order: if (tvd1 .or. posdef1) then
          bracket: if (left) then
            ix1 = max(ix - 1, 1)
            ix2 = ix
          else bracket
            ix1 = ix
            ix2 = min(ix + 1, nxa)
          endif bracket
!
          reduce: do while (ord > 0)
            ok = .true.
            if (tvd1 .and. (yi(ord) - ya(ix1)) * (ya(ix2) - yi(ord)) < 0.0) ok = .false.
            if (posdef1 .and. yi(ord) < 0.0) ok = .false.
            if (ok) exit reduce
            ord = ord - 1
          enddo reduce
        endif order
!
        y(i) = yi(ord)
!
      enddo loop
!
!  Error handling
!
      if (present(istatus)) istatus = istat
      if (present(message) .and. istat /= 0) message = msg
!
    endsubroutine poly_interp_fixorder
!***********************************************************************
    function complex_phase(z)
!
!  Takes complex number and returns Theta where
!  z = A*exp(i*theta).
!
!  17-may-06/anders+jeff: coded
!
      real :: c,re,im,complex_phase
      complex :: z
!
      complex_phase=0.0
!
      c=abs(z)
      re=real(z)
      im=aimag(z)
! I
      if ((re>=0.0).and.(im>=0.0)) complex_phase=     asin(im/c)
! II
      if ((re< 0.0).and.(im>=0.0)) complex_phase=  pi-asin(im/c)
! III
      if ((re< 0.0).and.(im< 0.0)) complex_phase=  pi-asin(im/c)
! IV
      if ((re>=0.0).and.(im< 0.0)) complex_phase=twopi+asin(im/c)
!
    endfunction complex_phase
!***********************************************************************
    function erfcc(x)
!
!  Numerical Recipes routine.
!
!  12-jul-2005/joishi: added, translated syntax to f90, and pencilized
!  21-jul-2006/joishi: generalized.
!
       real :: x(:)
!
       real, dimension(size(x)) :: erfcc,t,z
!
      z=abs(x)
      t=1/(1.0+0.5*z)
      erfcc=t*exp(-z*z-1.26551223+t*(1.00002368+t*(.37409196+t* &
            (.09678418+t*(-.18628806+t*(.27886807+t*(-1.13520398+t*  &
            (1.48851587+t*(-.82215223+t*.17087277)))))))))
      where (x<0.0) erfcc=2.0-erfcc
!
      return
!
    endfunction erfcc
!***********************************************************************
    elemental real function arcsinh(x)
!
!  Returns the inverse hyperbolic sine of x.
!
!  24-dec-14/ccyang: coded
!
      real, intent(in) :: x
!
      arcsinh = log(x + sqrt(x * x + 1.0))
!
    endfunction arcsinh
!***********************************************************************
    subroutine besselj_nu_int(res,nu,arg,loversample)
!
!  Calculate the cylindrical bessel function
!  with integer index. The function in gsl_wrapper.c
!  only calculates the cylindrical Bessel functions
!  with real index. The amount of factorials in the
!  real index Bessel function leads to over and underflows
!  as the index goes only moderately high.
!
!                 _
!             1  /  pi
!  J_m(z) = ____ |     (cos(z*sin(theta)-m*theta)) dtheta
!                |
!            pi _/  0
!
!  The function defines its own theta from 0 to pi for the
!  integration, with the same number of points as the
!  azimuthal direction.
!
!  06-03-08/wlad: coded
!
      real, dimension(:),allocatable :: angle,a
      real :: arg,res,d_angle
      integer :: i,nu,nnt
      logical, optional :: loversample
!
      intent(in)  :: nu,arg
      intent(out) :: res
!
!  Possibility of very high resolution.
!  Useful in start time, for instance.
!
      nnt=max(100,nygrid)
      if (present(loversample)) then
        if (loversample) nnt=30000
      endif
!
      allocate(angle(nnt))
      allocate(a(nnt))
!
      d_angle=pi/(nnt-1)
!
      do i=1,nygrid
        angle(i)=(i-1)*d_angle
      enddo
      a=cos(arg*sin(angle)-nu*angle)
!
      res=pi_1*d_angle*(sum(a(2:nnt-1))+.5*(a(1)+a(nnt)))
!
    endsubroutine besselj_nu_int
!***********************************************************************
    subroutine calc_complete_ellints(mu,Kappa_mu,E_mu,loversample)
!
!  Calculate the complete elliptic integrals of first (K)
!  and second kind (E)
!
!              _
!             /  pi/2
!  K(mu)  =   |       1/sqrt(1-mu*sin(x)) dx
!             |
!            _/  0
!
!  The function defines its own theta from 0 to pi for the
!  integration, with the same number of points as the
!  azimuthal direction, or 100 points if nygrid<100. The
!  integration is performed with the trapezoidal rule.  As
!  K(mu) is not defined at the point mu=1, we set K(1)=0
!
!  The complete elliptic integral of second kind
!
!              _
!             /  pi/2
!  E(mu)  =   |       sqrt(mu*sin(x)) dx
!             |
!            _/  0
!
!  is defined everywhere and does not need this fix.
!
!
      real, dimension(:),allocatable :: angle,a_K,a_E
      real :: mu,d_angle,Kappa_mu
      real, optional :: E_mu
      integer :: i,nnt
      logical, optional :: loversample
!
!  Possibility of very high resolution.
!  Useful in start time, for instance.
!
      nnt=max(100,nygrid)
      if (present(loversample)) then
        if (loversample) nnt=30000
      endif
!
      allocate(angle(nnt))
      allocate(a_K(nnt))
!
      d_angle=.5*pi/(nnt-1)
      do i=1,nnt
        angle(i)=(i-1)*d_angle
      enddo
!
!  Elliptic integral of first kind.
!
      a_K=d_angle/sqrt(1-(mu*sin(angle))**2)
      if (mu == 1 ) then
        Kappa_mu=0
      else
        Kappa_mu=sum(a_K(2:nnt-1)) + .5*(a_K(1)+a_K(nnt))
      endif
!
! Elliptic integral of second kind.
!
      if (present(E_mu)) then
        allocate(a_E(nnt))
        a_E=d_angle*sqrt(1-(mu*sin(angle))**2)
        E_mu=sum(a_E(2:nnt-1)) + .5*(a_E(1)+a_E(nnt))
      endif
!
    endsubroutine calc_complete_ellints
!***********************************************************************
!
!***********************************************************************
!*                                                                     *
!*    Program to calculate the first kind Bessel function of integer   *
!*    order N, for any REAL X, using the function BESSJ(N,X).          *
!*                                                                     *
!* ------------------------------------------------------------------- *
!*                                                                     *
!*    SAMPLE RUN:                                                      *
!*                                                                     *
!*    (Calculate Bessel function for N=2, X=0.75).                     *
!*                                                                     *
!*    Bessel function of order  2 for X =  0.7500:                     *
!*                                                                     *
!*         Y =  0.67073997E-01                                         *
!*                                                                     *
!* ------------------------------------------------------------------- *
!*   Reference: From Numath Library By Tuan Dang Trong in Fortran 77.  *
!*                                                                     *
!*                               F90 Release 1.0 By J-P Moreau, Paris. *
!***********************************************************************
!PROGRAM TBESSJ
!
!  REAL*8  BESSI, X, Y
!  INTEGER N
!
!  N=2
!  X=0.75
!
!  Y = BESSJ(N,X)
!
!  write(*,10) N, X
!  write(*,20) Y
!
!  stop
!
!10 format (/' Bessel Function of order ',I2,' for X=',F8.4,':')
!20 format(/'      Y = ',E15.8/)
!
!STOP
!END
!
     FUNCTION BESSJ (N,X)
!
!     This subroutine calculates the first kind modified Bessel function
!     of integer order N, for any REAL X. We use here the classical
!     recursion formula, when X > N. For X < N, the Miller's algorithm
!     is used to avoid overflows.
!     REFERENCE:
!     C.W.CLENSHAW, CHEBYSHEV SERIES FOR MATHEMATICAL FUNCTIONS,
!     MATHEMATICAL TABLES, VOL.5, 1962.
!
      integer, parameter :: IACC = 40
      real, parameter :: BIGNO = 1.E10, BIGNI = 1.E-10
      real :: X,BESSJ,TOX,BJM,BJ,BJP,SUM1
      integer :: N,J,M,JSUM
!
      IF (N==0) THEN
      BESSJ = BESSJ0(X)
      RETURN
      ENDIF
      IF (N==1) THEN
      BESSJ = BESSJ1(X)
      RETURN
      ENDIF
      IF (X==0.) THEN
      BESSJ = 0.
      RETURN
      ENDIF
      TOX = 2./X
      IF (X>FLOAT(N)) THEN
      BJM = BESSJ0(X)
      BJ  = BESSJ1(X)
      DO 11 J = 1,N-1
      BJP = J*TOX*BJ-BJM
      BJM = BJ
      BJ  = BJP
   11 CONTINUE
      BESSJ = BJ
      ELSE
      M = 2*((N+INT(SQRT(FLOAT(IACC*N))))/2)
      BESSJ = 0.
      JSUM = 0
      SUM1 = 0.
      BJP = 0.
      BJ  = 1.
      DO 12 J = M,1,-1
      BJM = J*TOX*BJ-BJP
      BJP = BJ
      BJ  = BJM
      IF (ABS(BJ)>BIGNO) THEN
      BJ  = BJ*BIGNI
      BJP = BJP*BIGNI
      BESSJ = BESSJ*BIGNI
      SUM1 = SUM1*BIGNI
      ENDIF
      IF (JSUM/=0) SUM1 = SUM1+BJ
      JSUM = 1-JSUM
      IF (J==N) BESSJ = BJP
   12 CONTINUE
      SUM1 = 2.*SUM1-BJ
      BESSJ = BESSJ/SUM1
      ENDIF
      RETURN
    endfunction BESSJ
!***********************************************************************
    FUNCTION BESSJ0 (X)
      REAL :: X,BESSJ0,AX,FR,FS,Z,FP,FQ,XX
!
!     This subroutine calculates the First Kind Bessel Function of
!     order 0, for any real number X. The polynomial approximation by
!     series of Chebyshev polynomials is used for 0<X<8 and 0<8/X<1.
!     REFERENCES:
!     M.ABRAMOWITZ,I.A.STEGUN, HANDBOOK OF MATHEMATICAL FUNCTIONS, 1965.
!     C.W.CLENSHAW, NATIONAL PHYSICAL LABORATORY MATHEMATICAL TABLES,
!     VOL.5, 1962.
!
      REAL :: Y,P1,P2,P3,P4,P5,R1,R2,R3,R4,R5,R6  &
               ,Q1,Q2,Q3,Q4,Q5,S1,S2,S3,S4,S5,S6
      DATA P1,P2,P3,P4,P5 /1.D0,-.1098628627D-2,.2734510407D-4, &
      -.2073370639D-5,.2093887211D-6 /
      DATA Q1,Q2,Q3,Q4,Q5 /-.1562499995D-1,.1430488765D-3, &
      -.6911147651D-5,.7621095161D-6,-.9349451520D-7 /
      DATA R1,R2,R3,R4,R5,R6 /57568490574.,-13362590354., &
      651619640.7,-11214424.18,77392.33017,-184.9052456 /
      DATA S1,S2,S3,S4,S5,S6 /57568490411.,1029532985., &
      9494680.718,59272.64853,267.8532712,1. /
      IF(X==0.) GO TO 1
      AX = ABS (X)
      IF (AX<8.) THEN
      Y = X*X
      FR = R1+Y*(R2+Y*(R3+Y*(R4+Y*(R5+Y*R6))))
      FS = S1+Y*(S2+Y*(S3+Y*(S4+Y*(S5+Y*S6))))
      BESSJ0 = FR/FS
      ELSE
      Z = 8./AX
      Y = Z*Z
      XX = AX-.785398164
      FP = P1+Y*(P2+Y*(P3+Y*(P4+Y*P5)))
      FQ = Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*Q5)))
      BESSJ0 = SQRT(.636619772/AX)*(FP*COS(XX)-Z*FQ*SIN(XX))
      ENDIF
      RETURN
    1 BESSJ0 = 1.
      RETURN
    endfunction BESSJ0
!***********************************************************************
    FUNCTION BESSJ1 (X)
      REAL :: X,BESSJ1,AX,FR,FS,Z,FP,FQ,XX
!     This subroutine calculates the First Kind Bessel Function of
!     order 1, for any real number X. The polynomial approximation by
!     series of Chebyshev polynomials is used for 0<X<8 and 0<8/X<1.
!     REFERENCES:
!     M.ABRAMOWITZ,I.A.STEGUN, HANDBOOK OF MATHEMATICAL FUNCTIONS, 1965.
!     C.W.CLENSHAW, NATIONAL PHYSICAL LABORATORY MATHEMATICAL TABLES,
!     VOL.5, 1962.
      REAL :: Y,P1,P2,P3,P4,P5,P6,R1,R2,R3,R4,R5,R6  &
               ,Q1,Q2,Q3,Q4,Q5,S1,S2,S3,S4,S5,S6
      DATA P1,P2,P3,P4,P5 /1.D0,.183105D-2,-.3516396496D-4,  &
      .2457520174D-5,-.240337019D-6 /,P6 /.636619772D0 /
      DATA Q1,Q2,Q3,Q4,Q5 /.04687499995D0,-.2002690873D-3,   &
      .8449199096D-5,-.88228987D-6,.105787412D-6 /
      DATA R1,R2,R3,R4,R5,R6 /72362614232.,-7895059235., &
      242396853.1,-2972611.439,15704.48260,-30.16036606 /
      DATA S1,S2,S3,S4,S5,S6 /144725228442.,2300535178., &
      18583304.74,99447.43394,376.9991397,1. /
!
      AX = ABS(X)
      IF (AX<8.) THEN
      Y = X*X
      FR = R1+Y*(R2+Y*(R3+Y*(R4+Y*(R5+Y*R6))))
      FS = S1+Y*(S2+Y*(S3+Y*(S4+Y*(S5+Y*S6))))
      BESSJ1 = X*(FR/FS)
      ELSE
      Z = 8./AX
      Y = Z*Z
      XX = AX-2.35619491
      FP = P1+Y*(P2+Y*(P3+Y*(P4+Y*P5)))
      FQ = Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*Q5)))
      BESSJ1 = SQRT(P6/AX)*(COS(XX)*FP-Z*SIN(XX)*FQ)*SIGN(S6,X)
      ENDIF
      RETURN
    endfunction BESSJ1
!
!***********************************************************************
    subroutine cyclic(a,b,c,alpha,beta,r,x,n)
!
!  Inversion of a tridiagonal system with periodic BC (alpha and beta
!  coefficients in the left and right corners). Used in the ADI scheme of the
!  implicit_physics module.
!  Note: this subroutine is using twice the tridag one written above by tobi.
!
!  | b1 c1 0  ...       beta | | x1   |   | r1   |
!  | a2 b2 c2 ...            | | x2   |   | r2   |
!  | 0  a3 b3 c3             | | x3   | = | r3   |
!  |          ...            | | ...  |   | ...  |
!  |          an-1 bn-1 cn-1 | | xn-1 |   | rn-1 |
!  | alpha    0    a_n  b_n  | | xn   |   | rn   |
!
! 08-Sep-07/gastine+dintrans: coded from Numerical Recipes (p67-68).
!
      implicit none
!
      integer :: i,n
      real, dimension(n) :: a,b,c,r,x,bb,u,z
      real    :: alpha,beta,gamma,fact
      logical :: err
      character(len=255) :: msg
!
      if (n<=2) stop "cyclic in the general module: n too small"
      gamma=-b(1)
      bb(1)=b(1)-gamma
      bb(n)=b(n)-alpha*beta/gamma
      do i=2,n-1
        bb(i)=b(i)
      enddo
      call tridag(a,bb,c,r,x,err,msg)
      if (err) print *, 'cyclic: ', trim(msg)
      u(1)=gamma
      u(n)=alpha
      do i=2,n-1
        u(i)=0.
      enddo
      call tridag(a,bb,c,u,z,err,msg)
      if (err) print *, 'cyclic: ', trim(msg)
      fact=(x(1)+beta*x(n)/gamma)/(1.+z(1)+beta*z(n)/gamma)
      do i=1,n
        x(i)=x(i)-fact*z(i)
      enddo
!
      return
!
    endsubroutine cyclic
!***********************************************************************
   function linear_interpolate_2d(f,xx,yy,xxp,lcheck) result(gp)
!
!  Interpolate the value of g to arbitrary (xp, yp, zp) coordinate
!  using the linear interpolation formula
!
!    g(x,y,z) = B*x*y + E*x + F*y + G*z + H .
!
!  The coefficients are determined by the 8 grid points surrounding the
!  interpolation point.
!
!  30-dec-04/anders: coded
!  04-nov-10/nils: moved from particles_map to general
!  22-apr-11/MR: changed to logical function to get rid of dependence on
!                module Messages
!
      use Cdata
!
      real, dimension(:,:) :: f
      real, dimension(:) :: xx,yy
      real, dimension(2) :: xxp

      real :: gp
!
      real, parameter :: eps=1.e-5
      real :: g1, g2, g3, g4
      real :: xp0, yp0
      real, save :: dx1, dy1
      integer :: ix0, iy0
      logical :: lcheck
!
      intent(in) :: f, xx, yy, xxp, lcheck
!
!  Determine index value of lowest lying corner point of grid box surrounding
!  the interpolation point.
!
      do ix0=1,size(xx)
        if ( xx(ix0)>=xxp(1) ) exit
      enddo
      if (ix0>1) ix0=ix0-1

      do iy0=1,size(yy)
        if ( yy(iy0)>=xxp(2) ) exit
      enddo
      if (iy0>1) iy0=iy0-1
!
!  Check if the grid point interval is really correct.
!
      if ( xx(ix0)-eps<=xxp(1) .and. xx(ix0+1)+eps>=xxp(1) .and. &
           yy(iy0)-eps<=xxp(2) .and. yy(iy0+1)+eps>=xxp(2) )  then
        ! Everything okay
      else
        if (lcheck) then
          print*, 'linear_interpolate: Interpolation point does not ' // &
                  'lie within the calculated grid point interval.'
          print'(a,f8.3,1x,f8.3,1x,f8.3)', 'xp, xp0, xp1 = ', xxp(1), xx(ix0), xx(ix0+1)
          print'(a,f8.3,1x,f8.3,1x,f8.3)', 'yp, yp0, yp1 = ', xxp(2), yy(iy0), yy(iy0+1)
        endif
        return
      endif
!
!  Redefine the interpolation point in coordinates relative to lowest corner.
!
      xp0=xxp(1)-xx(ix0)
      yp0=xxp(2)-yy(iy0)
!
!  Calculate derived grid spacing parameters needed for interpolation.
!
      dx1=1./(xx(2)-xx(1))
      dy1=1./(yy(2)-yy(1))
!
!  Function values at all corners.
!
      g1=f(ix0  ,iy0  )
      g2=f(ix0+1,iy0  )
      g3=f(ix0  ,iy0+1)
      g4=f(ix0+1,iy0+1)
!
!  Interpolation formula.
!
      gp = g1 + xp0*dx1*(-g1+g2) + yp0*dy1*(-g1+g3) + xp0*yp0*dx1*dy1*(g1-g2-g3+g4)

    endfunction linear_interpolate_2d
!***********************************************************************
   logical function linear_interpolate(f,ivar1,ivar2,xxp,gp,inear,lcheck)
!
!  Interpolate the value of g to arbitrary (xp, yp, zp) coordinate
!  using the linear interpolation formula
!
!    g(x,y,z) = A*x*y*z + B*x*y + C*x*z + D*y*z + E*x + F*y + G*z + H .
!
!  The coefficients are determined by the 8 grid points surrounding the
!  interpolation point.
!
!  30-dec-04/anders: coded
!  04-nov-10/nils: moved from particles_map to general
!  22-apr-11/MR: changed to logical function to get rid of dependence on
!                module Messages
!
      use Cdata
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: ivar1, ivar2
      real, dimension (3) :: xxp
      real, dimension (ivar2-ivar1+1) :: gp
      integer, dimension (3) :: inear
!
      real, dimension (ivar2-ivar1+1) :: g1, g2, g3, g4, g5, g6, g7, g8
      real :: xp0, yp0, zp0
      real, save :: dxdydz1, dxdy1, dxdz1, dydz1, dx1, dy1, dz1
      integer :: i, ix0, iy0, iz0
      logical :: lfirstcall=.true.,lcheck
!
      intent(in)  :: f, xxp, ivar1, lcheck
      intent(out) :: gp
!
!  Determine index value of lowest lying corner point of grid box surrounding
!  the interpolation point.
!
      linear_interpolate = .true.
!
      ix0=inear(1); iy0=inear(2); iz0=inear(3)
      if ( (x(ix0)>xxp(1)) .and. nxgrid/=1) ix0=ix0-1
      if ( (y(iy0)>xxp(2)) .and. nygrid/=1) iy0=iy0-1
      if ( (z(iz0)>xxp(3)) .and. nzgrid/=1) iz0=iz0-1
!
!  Check if the grid point interval is really correct.
!
      if ((x(ix0)<=xxp(1) .and. x(ix0+1)>=xxp(1) .or. nxgrid==1) .and. &
          (y(iy0)<=xxp(2) .and. y(iy0+1)>=xxp(2) .or. nygrid==1) .and. &
          (z(iz0)<=xxp(3) .and. z(iz0+1)>=xxp(3) .or. nzgrid==1)) then
        ! Everything okay
      else
        print*, 'linear_interpolate: Interpolation point does not ' // &
            'lie within the calculated grid point interval.'
        print*, 'iproc = ', iproc_world
        print*, 'mx, x(1), x(mx) = ', mx, x(1), x(mx)
        print*, 'my, y(1), y(my) = ', my, y(1), y(my)
        print*, 'mz, z(1), z(mz) = ', mz, z(1), z(mz)
        print*, 'ix0, iy0, iz0 = ', ix0, iy0, iz0
        print*, 'xp, xp0, xp1 = ', xxp(1), x(ix0), x(ix0+1)
        print*, 'yp, yp0, yp1 = ', xxp(2), y(iy0), y(iy0+1)
        print*, 'zp, zp0, zp1 = ', xxp(3), z(iz0), z(iz0+1)
        linear_interpolate = .false.
        return
      endif
!
!  Redefine the interpolation point in coordinates relative to lowest corner.
!  Set it equal to 0 for dimensions having 1 grid points; this will make sure
!  that the interpolation is bilinear for 2D grids.
!
      xp0=0; yp0=0; zp0=0
      if (nxgrid/=1) xp0=xxp(1)-x(ix0)
      if (nygrid/=1) yp0=xxp(2)-y(iy0)
      if (nzgrid/=1) zp0=xxp(3)-z(iz0)
!
!  Calculate derived grid spacing parameters needed for interpolation.
!  For an equidistant grid we only need to do this at the first call.
!
      if (lequidist(1)) then
        if (lfirstcall) dx1=dx_1(ix0) !1/dx
      else
        dx1=1/(x(ix0+1)-x(ix0))
      endif
!
      if (lequidist(2)) then
        if (lfirstcall) dy1=dy_1(iy0)
      else
        dy1=1/(y(iy0+1)-y(iy0))
      endif
!
      if (lequidist(3)) then
        if (lfirstcall) dz1=dz_1(iz0)
      else
        dz1=1/(z(iz0+1)-z(iz0))
      endif
!
      if ( (.not. all(lequidist)) .or. lfirstcall) then
        dxdy1=dx1*dy1; dxdz1=dx1*dz1; dydz1=dy1*dz1
        dxdydz1=dx1*dy1*dz1
      endif
!
!  Function values at all corners.
!
      g1=f(ix0  ,iy0  ,iz0  ,ivar1:ivar2)
      g2=f(ix0+1,iy0  ,iz0  ,ivar1:ivar2)
      g3=f(ix0  ,iy0+1,iz0  ,ivar1:ivar2)
      g4=f(ix0+1,iy0+1,iz0  ,ivar1:ivar2)
      g5=f(ix0  ,iy0  ,iz0+1,ivar1:ivar2)
      g6=f(ix0+1,iy0  ,iz0+1,ivar1:ivar2)
      g7=f(ix0  ,iy0+1,iz0+1,ivar1:ivar2)
      g8=f(ix0+1,iy0+1,iz0+1,ivar1:ivar2)
!
!  Interpolation formula.
!
      gp = g1 + xp0*dx1*(-g1+g2) + yp0*dy1*(-g1+g3) + zp0*dz1*(-g1+g5) + &
          xp0*yp0*dxdy1*(g1-g2-g3+g4) + xp0*zp0*dxdz1*(g1-g2-g5+g6) + &
          yp0*zp0*dydz1*(g1-g3-g5+g7) + &
          xp0*yp0*zp0*dxdydz1*(-g1+g2+g3-g4+g5-g6-g7+g8)
!
!  Do a reality check on the interpolation scheme.
!
      if (lcheck) then
        do i=1,ivar2-ivar1+1
          if (gp(i)>max(g1(i),g2(i),g3(i),g4(i),g5(i),g6(i),g7(i),g8(i))) then
            print*, 'linear_interpolate: interpolated value is LARGER than'
            print*, 'linear_interpolate: a values at the corner points!'
            print*, 'linear_interpolate: x0, y0, z0=', &
                x(ix0), y(iy0), z(iz0)
            print*, 'linear_interpolate: i, gp(i)=', i, gp(i)
            print*, 'linear_interpolate: g1...g8=', &
                g1(i), g2(i), g3(i), g4(i), g5(i), g6(i), g7(i), g8(i)
            print*, '------------------'
          endif
          if (gp(i)<min(g1(i),g2(i),g3(i),g4(i),g5(i),g6(i),g7(i),g8(i))) then
            print*, 'linear_interpolate: interpolated value is smaller than'
            print*, 'linear_interpolate: a values at the corner points!'
            print*, 'linear_interpolate: xxp=', xxp
            print*, 'linear_interpolate: x0, y0, z0=', &
                x(ix0), y(iy0), z(iz0)
            print*, 'linear_interpolate: i, gp(i)=', i, gp(i)
            print*, 'linear_interpolate: g1...g8=', &
                g1(i), g2(i), g3(i), g4(i), g5(i), g6(i), g7(i), g8(i)
            print*, '------------------'
          endif
        enddo
      endif
!
      if (lfirstcall) lfirstcall=.false.
!
    endfunction linear_interpolate
!***********************************************************************
    integer function parser( zeile, feld, tz )
!
! Parses string zeile according to separator character tz;
! Puts up substrings consecutively into feld.
! Returns number of found substrings.
!
!  10-apr-11/MR: coded
!  18-nov-13/MR: removed dummy parameter lenf, can be inferred from feld
!
      character (len=*),               intent(in)  :: zeile
      character (len=*), dimension(:), intent(out) :: feld
      character,                       intent(in)  :: tz
!
      integer :: ind, inda, lenz, lenf
!
      parser = 0
      inda = 1
      lenz = len_trim(zeile)
      lenf = size(feld)
!
      do while ( inda <= lenz )
!
        ind = index( zeile(inda:), tz )
!
        if ( inda+ind-1 < lenz .and. parser == lenf ) then
          print*, 'Parser - Warning: too many substrings!'
          return
        endif
!
        parser = parser+1
!
        if ( ind == 0 ) then
!
          feld(parser) = trim(zeile(inda:))
          return
!
        else
!
          if ( ind>1 ) then
            feld(parser) = trim(zeile(inda:inda+ind-2))
          else
            feld(parser) = ''
          endif
!
          inda = inda+ind
!
        endif
!
      enddo
!
    endfunction parser
!***********************************************************************
  subroutine write_full_columns_real(unit,buffer,range,unfilled,ncol,fmt)
!
! range-wise output of a real vector in ncol columns
! unfilled (inout) - number of unfilled slots in last written line
!
!  20-apr-11/MR: coded
!  29-jan-14/MR: introduced is for range step; inserted some trim calls
!  05-feb-14/MR: corrected wrong placement of is definition
!
    integer,                        intent(in)    :: unit
    real,    dimension(*),          intent(in)    :: buffer
    integer, dimension(3),          intent(in)    :: range
    integer,                        intent(inout) :: unfilled
    integer,              optional, intent(in)    :: ncol
    character(len=*),     optional, intent(in)    :: fmt
!
    integer              :: ncoll, nd, ia, ie, is, rest
    character(len=intlen):: str
    character(len=intlen):: fmtl, fmth
!
    if ( present(fmt) ) then
      fmtl = fmt
    else
      fmtl = 'e10.2'
    endif
!
    nd = get_range_no(range,1)
    if ( nd==0 ) return
!
    ncoll = ioptest(ncol,8)
    is = range(3)
!
    if ( unfilled > 0 ) then
!
      fmth = fmtl
      if ( nd>unfilled ) then
        ie = range(1)+(unfilled-1)*is
        str=itoa(unfilled)
      else
        ie = range(2)
        if (nd<unfilled) fmth = trim(fmtl)//'$'
        str=itoa(nd)
      endif
      nd = nd-unfilled
!
      write(unit,'(1p,'//trim(str)//trim(fmth)//')') buffer(range(1):ie:is)
!
      if ( nd>0 ) then
        ia = ie + is
      else
        unfilled = -nd
        return
      endif
    else
      ia = range(1)
    endif
!
    rest = mod(nd,ncoll)
    ie = (nd-rest-1)*is + ia
!
    if ( rest<nd ) then
      str=itoa(ncoll)
      write(unit,'(1p,'//trim(str)//trim(fmtl)//')') buffer(ia:ie:is)
    endif
!
    if ( rest > 0 ) then
      str=itoa(rest)
      write(unit,'(1p,'//trim(str)//trim(fmtl)//'$)') buffer(ie+is:range(2):is)
      unfilled = ncoll-rest
    else
      unfilled = 0
    endif
!
  endsubroutine write_full_columns_real
!***********************************************************************
  subroutine write_full_columns_cmplx(unit,buffer,range,unfilled,ncol,fmt)
!
! range-wise output of a complex vector in ncol columns
! unfilled (inout) - number of unfilled slots in last written line
!
!  20-apr-11/MR: coded
!  29-jan-14/MR: introduced is for range step; inserted some trim calls
!  05-feb-14/MR: corrected wrong placement of is definition
!
    integer,                        intent(in)    :: unit
    complex, dimension(*),          intent(in)    :: buffer
    integer, dimension(3),          intent(in)    :: range
    integer,                        intent(inout) :: unfilled
    integer,              optional, intent(in)    :: ncol
    character(len=*),     optional, intent(in)    :: fmt
!
    integer              :: ncoll, nd, ia, ie, is, rest
    character(len=intlen):: str
    character(len=intlen):: fmtl, fmth
!
    if ( present(fmt) ) then
      fmtl = '('//fmt//',1x,'//fmt//')'
    else
      fmtl = '(e10.2,1x,e10.2)'
    endif
!
    nd = get_range_no(range,1)
    if ( nd==0 ) return
!
    ncoll = ioptest(ncol,8)
    is = range(3)
!
    if ( unfilled > 0 ) then
!
      fmth = fmtl
      if ( nd>unfilled ) then
        ie = range(1)+(unfilled-1)*is
        str=itoa(unfilled)
      else
        ie = range(2)
        if (nd<unfilled) fmth = trim(fmtl)//'$'
        str=itoa(nd)
      endif
      nd = nd-unfilled
!
      write(unit,'(1p,'//trim(str)//trim(fmth)//')') buffer(range(1):ie:is)
!
      if ( nd>0 ) then
        ia = ie + is
      else
        unfilled = -nd
        return
      endif
    else
      ia = range(1)
    endif
!
    rest = mod(nd,ncoll)
    ie = (nd-rest-1)*is + ia
!
    if ( rest<nd ) then
      str=itoa(ncoll)
      write(unit,'(1p,'//trim(str)//trim(fmtl)//')') buffer(ia:ie:is)
    endif
!
    if ( rest > 0 ) then
      str=itoa(rest)
      write(unit,'(1p,'//trim(str)//trim(fmtl)//'$)') buffer(ie+is:range(2):is)
      unfilled = ncoll-rest
    else
      unfilled = 0
    endif
!
  endsubroutine write_full_columns_cmplx
!***********************************************************************
    function add_merge_range( ranges, ie, range ) result (iel)
!  
!  Checks whether range is contained (maybe partly) in any of ranges(:,1:ie) and stores it or its possible
!  fragments in ranges(:,ie+1:). Returns new total number of ranges.
!
!  20-apr-11/MR: coded
!
    integer, dimension(:,:),intent(INOUT):: ranges
    integer,                intent(IN)   :: ie
    integer, dimension(2),  intent(IN)   :: range

    integer :: iel

    integer, dimension(2,2*size(ranges,2)) :: ranges_loc
    integer :: i,ieloc,iloc
 
    if (ie==0) then
      iel=1
      ranges(:,1) = range
    else
      ranges_loc(:,1) = range
      ieloc=1
      do i=1,ie
        iloc=1
        do while (iloc<=ieloc)
          call analyze_range(ranges(:,i),ranges_loc,iloc,ieloc)
        enddo
      enddo
      iel=ie
      do i=1,ieloc
        if (ranges_loc(1,ieloc)>0) then
          iel=iel+1
          ranges(:,iel)=ranges_loc(:,ieloc)
        endif  
      enddo
    endif

    endfunction add_merge_range
!***********************************************************************
    subroutine analyze_range(testrange,store,istore,iestore)
!
!  Analyzes range in store(:,istore) with respect to testrange.
!  If fully contained it is marked deleted, if partly contained
!  leftover(s) is(are) stored in store(:,istore) (and store(:,iestore+1).
!
!  20-apr-11/MR: coded
!
    integer, dimension(2),  intent(IN)   :: testrange
    integer, dimension(:,:),intent(INOUT):: store
    integer,                intent(INOUT):: istore,iestore

    integer, dimension(2) :: range

      range=store(:,istore)

      if (range(1)>=testrange(1) .and. range(1)<=testrange(2)) then

        if (range(2)<=testrange(2)) then
          store(:,istore) = (/0,0/)                      ! range is completely absorbed
        else
          store(:,istore) = (/testrange(2)+1,range(2)/)  ! left leftover
        endif

      elseif (range(1)<testrange(1)) then
        if (range(2)>=testrange(1)) then
          store(:,istore) = (/range(1),testrange(1)-1/)  ! right leftover
        elseif (range(2)<testrange(1)) then
          continue                                       ! left unchanged
        else
          store(:,istore) = (/range(1),testrange(1)-1/)  ! two leftovers: left
          iestore=iestore+1
          store(:,iestore) = (/testrange(2)+1,range(2)/) ! right
          istore=istore+1
        endif
      else
        continue                                         ! right unchanged
      endif
      istore=istore+1

    endsubroutine analyze_range
!***********************************************************************
    function merge_ranges( ranges, ie, range, ia, istore )
!
! merges ranges ia through ie in vector ranges with range;
! return value indicates whether range has been modified
!
! 20-apr-11/MR: coded
! 10-feb-14/MR: case of different stepsizes implemented
!
      logical                                         :: merge_ranges
      integer, dimension(:,:),          intent(inout) :: ranges
      integer,                          intent(inout) :: ie
      integer, dimension(3),            intent(in)    :: range
      integer,                optional, intent(in)    :: ia
      integer,                optional, intent(inout) :: istore

      integer :: i, j, step, stepi, i1, ial, istl, iar, ier, i2, iai, iei, nor, ne, nstore, deleted
      integer, dimension(:), allocatable :: store
!
      merge_ranges=.false.
      ial =ioptest(ia,1)
      istl=ioptest(istore,ie)

      iar=range(1); ier=range(2); step=range(3)
!
      do i=ial,ie
!
        iai=ranges(1,i); iei=ranges(2,i); stepi=ranges(3,i)
!
! has ranges(:,i) same stepsize as range?
!
        if ( stepi == step ) then
!
          if ( iar <= iei+stepi .and. &
               ier >= iai-stepi .and. &
               mod(iar-iai,stepi) == 0 ) then
!
! join range with ranges(:,i)
!
            ranges(1,i) = min(iai, iar)
            ranges(2,i) = max(iei, ier)
!
! range absorbed
!
            merge_ranges=.true.
            return
!
          endif
        endif
      enddo
!
      nor  = get_range_no( range, 1 )
      allocate( store(nor) )
!
! expand range into store
!
      j=1
      do i=iar,ier,step
        store(j)=i
        j=j+1
      enddo
!
! marker for deleted elements in store: safely not belonging to range
!
      deleted = iar-step
      nstore = nor

      do i=ial,ie
!
        iai=ranges(1,i); iei=ranges(2,i); stepi=ranges(3,i)

        if (  iar <= iei .and. ier >= iai ) then
!
          do i2=1,nor
            do i1=iai,iei,stepi
              if (store(i2)==i1) then
!
! if element von range found in ranges(:,i) mark as deleted
!
                store(i2)=deleted
                nstore=nstore-1
                exit
              endif
            enddo
          enddo
!
! if nstore==0: range absorbed
!
          if (nstore==0)  then
            merge_ranges=.true.
            return
          endif
!
        endif
      enddo
!
      if ( nstore<nor ) then
!
! compress store by dropping deleted elements
!
        i=1; ne=nor
        do while (i<=ne)
          if ( store(i)==deleted ) then
            if ( i<ne ) store(i:ne-1) = store(i+1:ne)
            ne = ne-1
          else
            i=i+1
          endif
        enddo
        call find_ranges(store(1:ne),ranges,istl)
        merge_ranges=.true.
!
      elseif (istl==size(ranges,2)) then
        print*, 'merge_ranges: Warning - cannot store further ranges!'
        return
      else
!
! range is stored without modification
!
        if (.not.present(istore)) then
          istl=istl+1
          ranges(:,istl) = range
        endif
      endif
!
      if (present(istore)) then
        istore=istl
      else
        ie=istl
      endif
!
    endfunction merge_ranges
!***********************************************************************
    subroutine compress_vec(vec,mask,len)
!
!  Compresses vector vec guided by mask: vec(i) is removed where mask(i)=T.
!  Length of compressed vec is returned in len.
!  If len>0 on input it is taken as the actual length of vec.
!
!  12-jan-17/MR: coded
!
      real,    dimension(:), intent(inout):: vec
      logical, dimension(:), intent(in)   :: mask
      integer,               intent(inout):: len

      integer :: i

      if (len<=0) len=size(vec)

      do i=len,1,-1
        if ( mask(i) ) then
          if ( i<len ) vec(i:len-1) = vec(i+1:len)
          len = len-1
        endif
      enddo

    endsubroutine compress_vec
!***********************************************************************
    subroutine compress_arr(arr,mask,len)
!
!  Compresses array arr in 2nd dimension guided by mask: arr(:,i) is removed where mask(i)=T.
!  2nd dimension of compressed arr is returned in len.
!  If len>0 on input it is taken as the actual 2nd dimension of arr.
!
!  12-jan-17/MR: coded
!
      real,    dimension(:,:),intent(inout):: arr
      logical, dimension(:),  intent(in)   :: mask
      integer,                intent(inout):: len

      integer :: i

      if (len<=0) len=size(arr,2)

      do i=len,1,-1
        if ( mask(i) ) then
          if ( i<len ) arr(:,i:len-1) = arr(:,i+1:len)
          len = len-1
        endif
      enddo

    endsubroutine compress_arr
!***********************************************************************
    function fcompress(arr,mask,len) result(comparr)
!
!  Returns array arr, compressed in 2nd dimension guided by mask: arr(:,i) is removed where mask(i)=T.
!  2nd dimension of compressed arr is returned in len.
!  If len>0 on input it is taken as the actual 2nd dimension of arr.
!
!  12-jan-17/MR: coded
!
      real,    dimension(:,:),intent(inout):: arr
      logical, dimension(:),  intent(in)   :: mask
      integer,                intent(inout):: len
!
      real,    dimension(size(arr,1),size(arr,2))  :: comparr

      integer :: i,icomp

      if (len<=0) len=size(arr,2)
      icomp=0

      do i=1,len
        if ( .not.mask(i) ) then
          comparr(:,i) = arr(:,i)
          icomp = icomp+1
        endif
      enddo

      len=icomp

    endfunction fcompress
!***********************************************************************
    subroutine find_ranges(list,ranges,ie)
!
! extracts ranges start:stop:stop from a list and adds them to array of
! ie ranges, updates ie accordingly
!
! 4-feb-14/MR: coded
!
      integer, dimension(:),  intent(in)   :: list
      integer, dimension(:,:),intent(OUT)  :: ranges
      integer,                intent(INOUT):: ie

      integer :: i, j, ia, is, ja, len, sum

      len=size(list); i=1

      do while (i<=len)

        ia=i
        if (i<len) then
          is = list(i+1)-list(i)
          sum=list(i+1); ja=i+2

          do j=ja,len
            sum=sum+is
            if (list(j)/=sum) exit
          enddo
        else
          j=len+1; is=1
        endif
!
        if (ie==size(ranges,2)) then
          print*, 'find_ranges: Warning - cannot generate further ranges!'
          return
        else
          ie=ie+1
          ranges(:,ie) = (/list(ia),list(j-1),is/)
          i=j
        endif
!
      enddo
!
    endsubroutine find_ranges
!***********************************************************************
    logical function read_range_r( crange, range, defrange )
!
! reads a range (start,stop,step) of reals from string crange of the shape  <start>:<stop>[:<step>]
! if <step> missing it is set to 1.
! adjusts range not to exceed limits of default range 'defrange'
! crange is interpreted as far as possible, missing data are taken from defrange
!
! 20-apr-11/MR: coded
! 20-sep-12/MR: made defrange optional, changed order of dummy arguments
! 10-feb-14/MR: use of defrange debugged
!
      character (len=20), intent(in)           :: crange
      real, dimension(3), intent(out)          :: range
      real, dimension(3), intent(in), optional :: defrange
!
      integer :: isep, isep1, ios, ios1, ios2, lenrng
      real    :: tmp
      logical :: ldefr
!
      ldefr = present(defrange)

      ios=0; ios1=0; ios2=0
!
      if ( crange /= '' ) then
!
        if (ldefr) range = defrange
!
        isep = index(crange,':')
        lenrng = len_trim(crange)
!
        if ( isep > 0 ) then
!
          if ( isep > 1 ) then
            read( crange(1:isep-1),*,IOSTAT=ios ) range(1)
            if ( ldefr .and. ios == 0 ) &
              range(1) = min(max(defrange(1),range(1)),defrange(2))
          endif
!
          if ( isep < lenrng ) then
!
            isep1 = index(crange(isep+1:lenrng),':')+isep
!
            if ( isep1 == isep ) isep1 = lenrng+1
!
            if ( isep1 > isep+1 ) then
              read( crange(isep+1:isep1-1),*,IOSTAT=ios1 ) range(2)
              if ( ldefr .and. ios1 == 0 ) &
                range(2) = min(max(defrange(1),range(2)),defrange(2))
            endif
!
            if ( isep1 < lenrng ) then
!
              range(3)=1.
              read( crange(isep1+1:lenrng),*,IOSTAT=ios2 ) range(3)
!
              if ( ios2 == 0 ) range(3) = abs(range(3))
!
            endif
          endif
!
          if ( range(1) > range(2) ) then
            tmp = range(1); range(1) = range(2); range(2) = tmp
          endif
!
        else
!
          read( crange, *,IOSTAT=ios ) range(1)
          if ( ldefr .and. ios == 0 ) &
            range(1) = min(max(defrange(1),range(1)),defrange(2))
          range(2) = range(1)
!
        endif
!
        if ( ios/=0 .or. ios1/=0 .or. ios2/=0 ) &
          print*, 'read_range_r - Warning: invalid data in range!'
!
        read_range_r = .true.
      else
        read_range_r = .false.
      endif
!
    endfunction read_range_r
!***********************************************************************
    logical function read_range_i( crange, range, defrange )
!
! reads a range (start,stop,step) of integers from string crange of the shape <start>:<stop>[:<step>]
! if <step> missing it is set to 1
! adjusts range not to exceed limits of default range defrange
! crange is interpreted as far as possible, missing data are taken from defrange
!
! 20-sep-12/MR: coded
! 10-feb-14/MR: use of defrange debugged
!
      character (len=20)   , intent(in)           :: crange
      integer, dimension(3), intent(out)          :: range
      integer, dimension(3), intent(in), optional :: defrange
!
      integer :: isep, isep1, ios, ios1, ios2, lenrng, tmp
      logical :: ldefr
!
      ldefr = present(defrange)
!
      ios=0; ios1=0; ios2=0
!
      if ( crange /= '' ) then
!
        if (ldefr) range = defrange
!
        isep = index(crange,':')
        lenrng = len_trim(crange)
!
        if ( isep > 0 ) then
!
          if ( isep > 1 ) then
            read( crange(1:isep-1),*,IOSTAT=ios ) range(1)
            if ( ldefr.and.ios == 0 ) &
              range(1) = min(max(defrange(1),range(1)),defrange(2))
          endif
!
          if ( isep < lenrng ) then
!
            isep1 = index(crange(isep+1:lenrng),':')+isep
!
            if ( isep1 == isep ) isep1 = lenrng+1
!
            if ( isep1 > isep+1 ) then
              read( crange(isep+1:isep1-1),*,IOSTAT=ios1 ) range(2)
              if ( ldefr.and.ios1 == 0 ) &
                range(2) = min(max(defrange(1),range(2)),defrange(2))
            endif
!
            if ( isep1 < lenrng ) then
!
              range(3)=1
              read( crange(isep1+1:lenrng),*,IOSTAT=ios2 ) range(3)
!
              if ( ios2 == 0 ) range(3) = abs(range(3))
!
            endif
          endif
!
          if ( range(1) > range(2) ) then
            tmp = range(1); range(1) = range(2); range(2) = tmp
          endif
!
        else
!
          read( crange, *,IOSTAT=ios ) range(1)
          if ( ldefr.and.ios == 0 ) &
            range(1) = min(max(defrange(1),range(1)),defrange(2))
          range(2) = range(1)
!
        endif
!
        if ( ios/=0 .or. ios1/=0 .or. ios2/=0 ) &
          print*, 'read_range_i - Warning: invalid data in range!'
!
        read_range_i = .true.
      else
        read_range_i = .false.
      endif
!
    endfunction read_range_i
!***********************************************************************
    integer function get_range_no( ranges, nr )
!
! determines total number of elements selected by all ranges in ranges
!
! 20-apr-11/MR: coded
! 18-nov-13/MR: made nr mandatory
!
      integer, dimension(3,*), intent(in) :: ranges
      integer,                 intent(in) :: nr
!
      integer :: i
!
      get_range_no = 0
!
      do i=1,nr
        if ( ranges(1,i) > 0 ) then
          get_range_no = get_range_no + &
                         ceiling( (ranges(2,i)-ranges(1,i)+1.)/ranges(3,i))
        else
          exit
        endif
      enddo
!
    endfunction get_range_no
!******************************************************************************
  subroutine write_by_ranges_2d_real( unit, buffer, xranges, yranges, trans )
!
! writes a real 2D array controlled by lists of ranges
! xranges and yranges for each of the two dimensions; output optionally transposed
!
! 10-may-11/MR: coded
! 27-jan-14/MR: loops truncated if all valid ranges processed
! 29-jan-14/MR: changed declaration of buffer*
!
    use Cdata, only: nk_max
!
    integer,                 intent(in)           :: unit
    integer, dimension(3,*), intent(in)           :: xranges, yranges
    real,    dimension(:,:), intent(in)           :: buffer
    logical,                 intent(in), optional :: trans
!
    integer :: i,j,jl,unfilled
    logical :: transl
!
    unfilled = 0
!
    if ( present(trans) ) then
      transl = trans
    else
      transl = .false.
    endif
!
    do j=1,nk_max
      if ( yranges(1,j) == 0 ) exit
      do jl=yranges(1,j),yranges(2,j),yranges(3,j)
        do i=1,nk_max
          if ( xranges(1,i) == 0 ) exit
          if ( transl ) then
            call write_full_columns_real( unit, buffer(jl,:), xranges(1,i), unfilled )
          else
            call write_full_columns_real( unit, buffer(:,jl), xranges(1,i), unfilled )
          endif
        enddo
      enddo
    enddo
!
    if ( unfilled > 0 ) write( unit,'(a)')
!
  endsubroutine write_by_ranges_2d_real
!***********************************************************************
  subroutine write_by_ranges_2d_cmplx( unit, buffer, xranges, yranges, trans )
!
! writes a real or complex 2D array controlled by lists of ranges
! xranges and yranges for each of the two dimensions; output optionally transposed
!
! 10-may-11/MR: coded
! 27-jan-14/MR: loops truncated if all valid ranges processed
! 29-jan-14/MR: changed declaration of buffer*
!
    use Cdata, only: nk_max
!
    integer,                 intent(in)           :: unit
    integer, dimension(3,*), intent(in)           :: xranges, yranges
    complex, dimension(:,:), intent(in)           :: buffer
    logical,                 intent(in), optional :: trans
!
    integer :: i,j,jl,unfilled
    logical :: transl
!
    unfilled = 0
!
    if ( present(trans) ) then
      transl = trans
    else
      transl = .false.
    endif
!
    do j=1,nk_max
      if ( yranges(1,j) == 0 ) exit
      do jl=yranges(1,j),yranges(2,j),yranges(3,j)
        do i=1,nk_max
          if ( xranges(1,i) == 0 ) exit
          if ( transl ) then
            call write_full_columns_cmplx( unit, buffer(jl,:), xranges(1,i), unfilled )
          else
            call write_full_columns_cmplx( unit, buffer(:,jl), xranges(1,i), unfilled )
          endif
        enddo
      enddo
    enddo
!
    if ( unfilled > 0 ) write( unit,'(a)')
!
  endsubroutine write_by_ranges_2d_cmplx
!***********************************************************************
  subroutine write_by_ranges_1d_real(unit,buffer,ranges)
!
! writes a real vector controlled by lists of ranges
! output optionally transposed
!
! 10-may-11/MR: coded
! 27-jan-14/MR: loop truncated if all valid ranges processed
! 29-jan-14/MR: changed declaration of buffer*
!
    use Cdata, only: nk_max
!
    integer,                 intent(in) :: unit
    real,    dimension(:)  , intent(in) :: buffer
    integer, dimension(3,*), intent(in) :: ranges
!
    integer :: unfilled, i
!
    unfilled = 0
    do i=1,nk_max
      if ( ranges(1,i) == 0 ) exit
      call write_full_columns_real( unit, buffer, ranges(1,i), unfilled )
    enddo
!
    if ( unfilled > 0 ) write(unit,'(a)')
!
  endsubroutine write_by_ranges_1d_real
!***********************************************************************
  subroutine write_by_ranges_1d_cmplx(unit,buffer,ranges)
!
! writes a complex vector controlled by lists of ranges
! output optionally transposed
!
! 10-may-11/MR: coded
! 27-jan-14/MR: loop truncated if all valid ranges processed
! 29-jan-14/MR: changed declaration of buffer*
!
    use Cdata, only: nk_max
!
    integer,                 intent(in) :: unit
    complex, dimension(:)  , intent(in) :: buffer
    integer, dimension(3,*), intent(in) :: ranges
!
    integer :: unfilled, i
!
    unfilled = 0
    do i=1,nk_max
      if ( ranges(1,i) == 0 ) exit
      call write_full_columns_cmplx( unit, buffer, ranges(1,i), unfilled )
    enddo
!
    if ( unfilled > 0 ) write(unit,'(a)')
!
  endsubroutine write_by_ranges_1d_cmplx
!***********************************************************************
    subroutine date_time_string(date)
!
!  Return current date and time as a string.
!  Subroutine, because nested writes don't work on some machines, so
!  calling a function like
!    print*, date_time_string()
!  may crash mysteriously.
!
!  4-oct-02/wolf: coded
!  4-nov-11/MR: moved from Sub to avoid circular dep's; crash due to too short parameter date avoided
!
      intent (out) :: date
!
      character (len=*) :: date
      integer, dimension(8) :: values
      character (len=3), dimension(12) :: month = &
           (/ 'jan', 'feb', 'mar', 'apr', 'may', 'jun', &
              'jul', 'aug', 'sep', 'oct', 'nov', 'dec' /)
      character (len=datelen) datestr
!
      if (len(date) < 20) &
          print*, 'date_time_string: WARNING -- string arg "date" too short'
!
      call date_and_time(VALUES=values)
      write(datestr,'(I2.2,"-",A3,"-",I4.2," ",I2.2,":",I2.2,":",I2.2)') &
           values(3), month(values(2)), values(1), &
           values(5), values(6), values(7)
      date=trim(datestr)
!
! TEMPORARY DEBUGGING STUFF
! SOMETIMES THIS ROUTINE PRINTS '***' WHEN IT SHOULDN'T
!
      if (index(date,'*')>0) then
        open(11,FILE='date_time_string.debug')
        write(11,*) 'This file was generated because sub$date_time_string()'
        write(11,*) 'produced a strange result.'
        write(11,*)
        write(11,*) 'date = <', date, '>'
        write(11,*) 'values = ', values
        write(11,*) 'i.e.'
        write(11,*) 'values(1) = ', values(1)
        write(11,*) 'values(2) = ', values(2)
        write(11,*) 'values(3) = ', values(3)
        write(11,*) 'values(4) = ', values(4)
        write(11,*) 'values(5) = ', values(5)
        write(11,*) 'values(6) = ', values(6)
        write(11,*) 'values(7) = ', values(7)
        close(11)
      endif
!
!  END OF TEMPORARY DEBUGGING STUFF
!
    endsubroutine date_time_string
!***********************************************************************
    logical function backskip(unit,count)
!
! sets record pointer back by count positions
!
!  3-nov-11/MR: coded
! 16-nov-11/MR: changed into logical function to signal I/O errors
!
    integer,           intent(in) :: unit
    integer, optional, intent(in) :: count
!
    integer :: i,n,iostat
!
    if (present(count)) then
      n=count
    else
      n=1
    endif
!
    backskip = .true.
!
    do i=1,n
      backspace(unit,IOSTAT=iostat)
      if (iostat/=0) return
    enddo
!
    backskip = .false.
!
    endfunction backskip
!***********************************************************************
    logical function lextend_vector_float(vector,newlen)
!
!  Checks whether a vector of floats can be used up to length newlen.
!  Surrogate for a real extension routine possible with FORTRAN 2003.
!
!  16-may-12/MR: coded
!
      real, dimension(:), intent(in) :: vector
      integer           , intent(in) :: newlen
!
      lextend_vector_float = newlen<=size(vector)
!
    endfunction lextend_vector_float
!***********************************************************************
    logical function lextend_vector_char(vector,newlen)
!
!  Checks whether a vector of chars can be used up to length newlen.
!  Surrogate for a real extension routine possible with FORTRAN 2003.
!
!  16-may-12/MR: coded
!
      character (len=*), dimension(:), intent(in) :: vector
      integer          , intent(in) :: newlen
!
      lextend_vector_char = newlen<=size(vector)
!
    endfunction lextend_vector_char
!***********************************************************************
  integer function pos_in_array_int(needle, haystack)
!
!  finds the position of a number in an array
!  returns 0 if string is not found
!
!  28-May-2015/Bourdin.KIS: reworked
!
    integer, intent(in) :: needle
    integer, dimension(:), intent(in) :: haystack
    integer :: pos

    do pos = 1, size(haystack)
      if (needle == haystack(pos)) then
        pos_in_array_int = pos
        return
      endif
    enddo

    pos_in_array_int = 0

  endfunction pos_in_array_int
!***********************************************************************
  integer function pos_in_array_char(needle, haystack)
!
!  finds the position of a string in an array
!  returns 0 if string is not found
!
!  28-May-2015/Bourdin.KIS: reworked
!
    character (len=*), intent(in) :: needle
    character (len=*), dimension(:), intent(in) :: haystack
    integer :: pos

    pos_in_array_char = -1

    do pos = 1, size(haystack)
      if (needle == haystack(pos)) then
         pos_in_array_char = pos
         return
      endif
    enddo

    pos_in_array_char = 0

  endfunction pos_in_array_char
!***********************************************************************
  logical function in_array_int(needle, haystack)
!
!  tests if a number is contained in a given array
!
!  28-May-2015/Bourdin.KIS: coded
!
    integer, intent(in) :: needle
    integer, dimension(:), intent(in) :: haystack

    in_array_int = .false.

    if (pos_in_array (needle, haystack) > 0) in_array_int = .true.

  endfunction in_array_int
!***********************************************************************
  logical function in_array_char(needle, haystack)
!
!  tests if a string str is contained in a vector of strings cvec
!
!  28-May-2015/Bourdin.KIS: coded, inspired by MR's string_in_array
!
    character (len=*), intent(in) :: needle
    character (len=*), dimension(:), intent(in) :: haystack

    in_array_char = .false.

    if (pos_in_array (needle, haystack) > 0) in_array_char = .true.

  endfunction in_array_char
!***********************************************************************
    pure logical function loptest(lopt,ldef)
!
!  returns value of optional logical parameter opt if present,
!  otherwise the default value ldef, if present, .false. if not
!
!  20-aug-13/MR: coded
!  26-aug-13/MR: optional default value ldef added
!
      logical, optional, intent(in) :: lopt, ldef

      if (present(lopt)) then
        loptest=lopt
      else if (present(ldef)) then
        loptest=ldef
      else
        loptest=.false.
      endif

    endfunction loptest
!***********************************************************************
    pure integer function ioptest(iopt,idef)
!
!  returns value of optional integer parameter iopt if present,
!  otherwise the default value idef, if present, zero, if not.
!
!  20-aug-13/MR: coded
!
      integer, optional, intent(in) :: iopt, idef

      if (present(iopt)) then
        ioptest=iopt
      elseif (present(idef)) then
        ioptest=idef
      else
        ioptest=0
      endif

    endfunction ioptest
!***********************************************************************
    pure real function roptest(ropt,rdef)
!
!  returns value of optional real parameter ropt if present,
!  otherwise the default value rdef, if present, zero, if not.
!
!  20-aug-13/MR: coded
!
      real, optional, intent(in) :: ropt, rdef

      if (present(ropt)) then
        roptest=ropt
      elseif (present(rdef)) then
        roptest=rdef
      else
        roptest=0.
      endif

    endfunction roptest
!***********************************************************************
    pure real(KIND=rkind8) function doptest(dopt,ddef)
!
!  returns value of optional real*8 parameter dopt if present,
!  otherwise the default value ddef, if present, zero, if not.
!
!  20-aug-13/MR: coded
!
      real(KIND=rkind8), optional, intent(in) :: dopt, ddef

      if (present(dopt)) then
        doptest=dopt
      elseif (present(ddef)) then
        doptest=ddef
      else
        doptest=0.
      endif

    endfunction doptest
!***********************************************************************
    pure function coptest(copt,cdef)
!
!  returns value of optional character parameter copt if present,
!  otherwise the default value cdef, if present, '', if not.
!
!  27-jan-15/MR: coded
!
      character(len=2*labellen) :: coptest
      character(len=*), optional, intent(in) :: copt, cdef

      if (present(copt)) then
        coptest=copt
      elseif (present(cdef)) then
        coptest=cdef
      else
        coptest=''
      endif

    endfunction coptest
!***********************************************************************
    pure integer function binary_search(x, v)
!
!  Uses binary search to find the index of the element in array v that
!  matches x; return 0 if no match is found.
!
!  v must already be in ascending order.
!
!  29-feb-16/ccyang: coded.
!
      integer, dimension(:), intent(in) :: v
      integer, intent(in) :: x
!
      integer :: low, high, mid
!
      binary_search = 0
      low = 1
      high = size(v)
      binary: do while (low <= high)
        mid = (low + high) / 2
        if (x < v(mid)) then
          high = mid - 1
        elseif (x > v(mid)) then
          low = mid + 1
        else
          binary_search = mid
          exit binary
        endif
      enddo binary
!
    endfunction binary_search
!***********************************************************************
    SUBROUTINE quick_sort(list, order)
!
! Quick sort routine from:
! Brainerd, W.S., Goldberg, C.H. & Adams, J.C. (1990) "Programmer's Guide to
! Fortran 90", McGraw-Hill  ISBN 0-07-000248-7, pages 149-150.
! Modified by Alan Miller to include an associated integer array which gives
! the positions of the elements in the original order.
!
! 5-feb-14/MR: added
!
      IMPLICIT NONE
      INTEGER, DIMENSION(:), INTENT(INOUT):: list
      INTEGER, DIMENSION(*), INTENT(OUT)  :: order
!
! Local variable
!
      INTEGER :: i

      DO i = 1, SIZE(list)
        order(i) = i
      ENDDO
!
      CALL quick_sort_1(1, SIZE(list))
!
      CONTAINS
!----------------------------------------------------------------------------
        RECURSIVE SUBROUTINE quick_sort_1(left_end, right_end)
!
          INTEGER, INTENT(IN) :: left_end, right_end
!
! Local variables
!
          INTEGER             :: i, j, itemp, reference, temp
          INTEGER, PARAMETER  :: max_simple_sort_size = 6
!
          IF (right_end < left_end + max_simple_sort_size) THEN
! Use interchange sort for small lists
            CALL interchange_sort(left_end, right_end)
!
          ELSE
! Use partition ("quick") sort
            reference = list((left_end + right_end)/2)
            i = left_end - 1; j = right_end + 1
!
            DO
! Scan list from left end until element >= reference is found
              DO
                i = i + 1
                IF (list(i) >= reference) EXIT
              ENDDO
! Scan list from right end until element <= reference is found
              DO
                j = j - 1
                IF (list(j) <= reference) EXIT
              ENDDO
!
              IF (i < j) THEN
! Swap two out-of-order elements
                temp = list(i); list(i) = list(j); list(j) = temp
                itemp = order(i); order(i) = order(j); order(j) = itemp
              ELSE IF (i == j) THEN
                i = i + 1
                EXIT
              ELSE
                EXIT
              ENDIF
            ENDDO
!
            IF (left_end < j)  CALL quick_sort_1(left_end, j)
            IF (i < right_end) CALL quick_sort_1(i, right_end)
          ENDIF
!
        END SUBROUTINE quick_sort_1
!----------------------------------------------------------------------------
        SUBROUTINE interchange_sort(left_end, right_end)

          INTEGER, INTENT(IN) :: left_end, right_end
!
! Local variables
!
          INTEGER :: i, j, itemp, temp
!
          DO i = left_end, right_end - 1
            DO j = i+1, right_end
              IF (list(i) > list(j)) THEN
                temp = list(i); list(i) = list(j); list(j) = temp
                itemp = order(i); order(i) = order(j); order(j) = itemp
              ENDIF
            ENDDO
          ENDDO
!
        ENDSUBROUTINE interchange_sort
!----------------------------------------------------------------------------
    ENDSUBROUTINE quick_sort
!****************************************************************************
    function indgen(n)
!
! Generates vector of integers  1,...,n, analogous to IDL-indgen, but starting
! at 1.
!
! 5-feb-14/MR: coded
!
      integer, intent(in) :: n
      integer, dimension(n) :: indgen

      integer :: i

      do i=1,n
        indgen(i)=i
      enddo

    endfunction indgen
!****************************************************************************
    subroutine ranges_dimensional(jrange)
!
!   Determines indices of the non-degenerate dimensions (out of {1,2,3} for {x,y,z})
!   and returns them in vector jrange of length dimensionality. Order of indices
!   is always ascending.
!
!   9-oct-15/MR: coded
! 
      integer, dimension(:), intent(OUT) :: jrange

      if (dimensionality==3) then
        jrange=(/1,2,3/)
      else if (dimensionality==2) then
        if (nxgrid==1) then
          jrange=(/2,3/)
        else if (nygrid==1) then
          jrange=(/1,3/)
        else if(nzgrid==1) then
          jrange=(/1,2/)
        endif
      else
        if (nxgrid/=1) then
          jrange=(/1/)
        else if(nygrid/=1) then
          jrange=(/2/)
        else if(nzgrid/=1) then
          jrange=(/3/)
        endif
      endif

    endsubroutine ranges_dimensional
!***********************************************************************
    subroutine staggered_mean_vec(f,k,jmean,weight)
!
!   Calculates mean squared modulus of a vector quantity in the centre of a grid cell.
!
!   9-oct-15/MR: coded
!  25-oct-16/JW: rewritten
!
      real, dimension (mx,my,mz,mfarray), intent(inout):: f
      integer,                            intent(in)   :: k,jmean
      real,                               intent(in)   :: weight

      real, parameter :: i8_1=1/8., i4_1=1/4., i2_1=1/2.
      integer         :: ll,mm,nn
!
      if (dimensionality==3) then
!     
        do ll=2,mx-2; do mm=2,my-2; do nn=2,mz-2
          f(ll,mm,nn,jmean) = f(ll,mm,nn,jmean) &
          +weight*i8_1*(maxval(abs(f(ll,mm,nn,      k:k+2))) &
                      + maxval(abs(f(ll,mm,nn+1,    k:k+2))) &
                      + maxval(abs(f(ll,mm+1,nn,    k:k+2))) &
                      + maxval(abs(f(ll,mm+1,nn+1,  k:k+2))) &
                      + maxval(abs(f(ll+1,mm,nn,    k:k+2))) &
                      + maxval(abs(f(ll+1,mm,nn+1,  k:k+2))) &
                      + maxval(abs(f(ll+1,mm+1,nn,  k:k+2))) &
                      + maxval(abs(f(ll+1,mm+1,nn+1,k:k+2))))
        enddo; enddo; enddo
!
      elseif (dimensionality==1) then
!
        if (nxgrid/=1) then
          do ll=2,mx-2; do mm=m1,m2; do nn=n1,n2
            f(ll,mm,nn,jmean) = f(ll,mm,nn,jmean) &
            +weight*i2_1*(maxval(abs(f(ll  ,mm,nn,k:k+2))) &
                        + maxval(abs(f(ll+1,mm,nn,k:k+2))))
          enddo; enddo; enddo
        elseif (nygrid/=1) then
          do ll=l1,l2; do mm=2,my-2; do nn=n1,n2
            f(ll,mm,nn,jmean) = f(ll,mm,nn,jmean) &
            +weight*i2_1*(maxval(abs(f(ll,mm  ,nn,k:k+2))) &
                        + maxval(abs(f(ll,mm+1,nn,k:k+2))))
          enddo; enddo; enddo
        else
          do ll=l1,l2; do mm=m1,m2; do nn=2,mz-2
            f(ll,mm,nn,jmean) = f(ll,mm,nn,jmean) &
            +weight*i2_1*(maxval(abs(f(ll,mm,nn,  k:k+2))) &
                        + maxval(abs(f(ll,mm,nn+1,k:k+2))))
          enddo; enddo; enddo      
        endif
!
!       dimensionality==2
!
      elseif (nzgrid==1) then   !  x-y
        do ll=2,mx-2; do mm=2,my-2; do nn=n1,n2
          f(ll,mm,nn,jmean) = f(ll,mm,nn,jmean) &
          +weight*i4_1*(maxval(abs(f(ll  ,mm  ,nn,k:k+2))) &
                      + maxval(abs(f(ll  ,mm+1,nn,k:k+2))) &
                      + maxval(abs(f(ll+1,mm  ,nn,k:k+2))) &                 
                      + maxval(abs(f(ll+1,mm+1,nn,k:k+2))))
        enddo; enddo; enddo
      elseif (nygrid==1) then   !  x-z
         do ll=2,mx-2; do mm=m1,m2; do nn=2,mz-2
          f(ll,mm,nn,jmean) = f(ll,mm,nn,jmean) &
          +weight*i4_1*(maxval(abs(f(ll  ,mm,nn,  k:k+2))) &
                      + maxval(abs(f(ll  ,mm,nn+1,k:k+2))) &
                      + maxval(abs(f(ll+1,mm,nn,  k:k+2))) &                 
                      + maxval(abs(f(ll+1,mm,nn+1,k:k+2))))
        enddo; enddo; enddo
      else                      !  y-z
        do ll=l1,l2; do mm=2,my-2; do nn=2,mz-2
          f(ll,mm,nn,jmean) = f(ll,mm,nn,jmean) &
          +weight*i4_1*(maxval(abs(f(ll,mm  ,nn,  k:k+2))) &
                      + maxval(abs(f(ll,mm  ,nn+1,k:k+2))) &
                      + maxval(abs(f(ll,mm+1,nn,  k:k+2))) &                 
                      + maxval(abs(f(ll,mm+1,nn+1,k:k+2))))
        enddo; enddo; enddo
      endif
    endsubroutine staggered_mean_vec
!***********************************************************************
    subroutine staggered_mean_scal(f,k,jmean,weight)
!
!   Calculates squared mean of a scalar quantity in the centre of a grid cell.
!
!   9-oct-15/MR: coded
!  25-oct-16/JW: rewritten
!
      real, dimension (mx,my,mz,mfarray), intent(inout):: f
      integer,                            intent(in)   :: k,jmean
      real,                               intent(in)   :: weight
!
      real, parameter :: i8_1=1/8., i4_1=1/4., i2_1=1/2.
      integer         :: ll,mm,nn
!
      if (dimensionality==3) then
!     
        do ll=2,mx-2; do mm=2,my-2; do nn=2,mz-2
          f(ll,mm,nn,jmean) = f(ll,mm,nn,jmean) &
          +weight*i8_1*(abs(f(ll,mm,nn,      k)) &
                      + abs(f(ll,mm,nn+1,    k)) &
                      + abs(f(ll,mm+1,nn,    k)) &
                      + abs(f(ll,mm+1,nn+1,  k)) &
                      + abs(f(ll+1,mm,nn,    k)) &
                      + abs(f(ll+1,mm,nn+1,  k)) &
                      + abs(f(ll+1,mm+1,nn,  k)) &
                      + abs(f(ll+1,mm+1,nn+1,k)))
        enddo; enddo; enddo
!
      elseif (dimensionality==1) then
!
        if (nxgrid/=1) then
          do ll=2,mx-2; do mm=m1,m2; do nn=n1,n2
            f(ll,mm,nn,jmean) = f(ll,mm,nn,jmean) &
            +weight*i2_1*(abs(f(ll  ,mm,nn,k)) &
                        + abs(f(ll+1,mm,nn,k)))
          enddo; enddo; enddo
        elseif (nygrid/=1) then
          do ll=l1,l2; do mm=2,my-2; do nn=n1,n2
            f(ll,mm,nn,jmean) = f(ll,mm,nn,jmean) &
            +weight*i2_1*(abs(f(ll,mm  ,nn,k)) &
                        + abs(f(ll,mm+1,nn,k)))
          enddo; enddo; enddo
        else
          do ll=l1,l2; do mm=m1,m2; do nn=2,mz-2
            f(ll,mm,nn,jmean) = f(ll,mm,nn,jmean) &
            +weight*i2_1*(abs(f(ll,mm,nn,  k)) &
                        + abs(f(ll,mm,nn+1,k)))
          enddo; enddo; enddo      
        endif
!
!       dimensionality==2
!
      elseif (nzgrid==1) then   !  x-y
        do ll=2,mx-2; do mm=2,my-2; do nn=n1,n2
          f(ll,mm,nn,jmean) = f(ll,mm,nn,jmean) &
          +weight*i4_1*(abs(f(ll  ,mm  ,nn,k)) &
                      + abs(f(ll  ,mm+1,nn,k)) &
                      + abs(f(ll+1,mm  ,nn,k)) &                 
                      + abs(f(ll+1,mm+1,nn,k)))
        enddo; enddo; enddo
      elseif (nygrid==1) then   !  x-z
         do ll=2,mx-2; do mm=m1,m2; do nn=2,mz-2
          f(ll,mm,nn,jmean) = f(ll,mm,nn,jmean) &
          +weight*i4_1*(abs(f(ll  ,mm,nn,  k)) &
                      + abs(f(ll  ,mm,nn+1,k)) &
                      + abs(f(ll+1,mm,nn,  k)) &                 
                      + abs(f(ll+1,mm,nn+1,k)))
        enddo; enddo; enddo
      else                      !  y-z
        do ll=l1,l2; do mm=2,my-2; do nn=2,mz-2
          f(ll,mm,nn,jmean) = f(ll,mm,nn,jmean) &
          +weight*i4_1*(abs(f(ll,mm  ,nn,  k)) &
                      + abs(f(ll,mm  ,nn+1,k)) &
                      + abs(f(ll,mm+1,nn,  k)) &                 
                      + abs(f(ll,mm+1,nn+1,k)))
        enddo; enddo; enddo
      endif

    endsubroutine staggered_mean_scal
!***********************************************************************
    subroutine staggered_max_vec(f,k,jmax,weight)
!
!   Calculates max values of modulus of a vector quantity in the centre of a grid cell.
!
!   22-dec-15/JW: coded, adapted from staggered_mean_vec
!
      real, dimension (mx,my,mz,mfarray), intent(inout):: f
      integer,                            intent(in)   :: k,jmax
      real,                               intent(in)   :: weight
      integer                                          :: ll,mm,nn
!
      if (dimensionality==3) then
!    
        do ll=2,mx-2; do mm=2,my-2; do nn=2,mz-2
          f(ll,mm,nn,jmax) = max(f(ll,mm,nn,jmax), &
          weight*max(maxval(abs(f(ll,mm,nn,      k:k+2))) &
                     ,maxval(abs(f(ll,mm,nn+1,    k:k+2))) &
                     ,maxval(abs(f(ll,mm+1,nn,    k:k+2))) &
                     ,maxval(abs(f(ll,mm+1,nn+1,  k:k+2))) &
                     ,maxval(abs(f(ll+1,mm,nn,    k:k+2))) &
                     ,maxval(abs(f(ll+1,mm,nn+1,  k:k+2))) &
                     ,maxval(abs(f(ll+1,mm+1,nn,  k:k+2))) &
                     ,maxval(abs(f(ll+1,mm+1,nn+1,k:k+2)))))
        enddo; enddo; enddo

      elseif (dimensionality==1) then
        if (nxgrid/=1) then
          do ll=2,mx-2; do mm=m1,m2; do nn=n1,n2
            f(ll,mm,nn,jmax) = max(f(ll,mm,nn,jmax), &
            weight*max(maxval(abs(f(ll,mm,nn,  k:k+2))) &
                       ,maxval(abs(f(ll+1,mm,nn,k:k+2)))))
          enddo; enddo; enddo
        elseif (nygrid/=1) then
          do ll=l1,l2; do mm=2,my-2; do nn=n1,n2
            f(ll,mm,nn,jmax) = max(f(ll,mm,nn,jmax), &
            weight*max(maxval(abs(f(ll,mm,nn,  k:k+2))) &
                       ,maxval(abs(f(ll,mm+1,nn,k:k+2)))))
          enddo; enddo; enddo
        else
          do ll=l1,l2; do mm=m1,m2; do nn=2,mz-2
            f(ll,mm,nn,jmax) = max(f(ll,mm,nn,jmax), &
            weight*max(maxval(abs(f(ll,mm,nn,  k:k+2))) &
                       ,maxval(abs(f(ll,mm,nn+1,k:k+2)))))
          enddo; enddo; enddo      
        endif
      elseif (nzgrid==1) then   !  x-y
        do ll=2,mx-2; do mm=2,my-2; do nn=n1,n2
          f(ll,mm,nn,jmax) = max(f(ll,mm,nn,jmax), &
          weight*max(maxval(abs(f(ll,mm,nn,    k:k+2))) &
                     ,maxval(abs(f(ll,mm+1,nn,  k:k+2))) &
                     ,maxval(abs(f(ll+1,mm,nn,  k:k+2))) &                 
                     ,maxval(abs(f(ll+1,mm+1,nn,k:k+2)))))
        enddo; enddo; enddo
      elseif (nygrid==1) then   !  x-z
         do ll=2,mx-2; do mm=m1,m2; do nn=2,mz-2
          f(ll,mm,nn,jmax) = max(f(ll,mm,nn,jmax), &
          weight*max(maxval(abs(f(ll,mm,nn,    k:k+2))) &
                     ,maxval(abs(f(ll,mm,nn+1,  k:k+2))) &
                     ,maxval(abs(f(ll+1,mm,nn,  k:k+2))) &                 
                     ,maxval(abs(f(ll+1,mm,nn+1,k:k+2)))))
        enddo; enddo; enddo
      else                      !  y-z
        do ll=l1,l2; do mm=2,my-2; do nn=2,mz-2
          f(ll,mm,nn,jmax) = max(f(ll,mm,nn,jmax), &
           weight*max(maxval(abs(f(ll,mm,nn,    k:k+2))) &
                     ,maxval(abs(f(ll,mm,nn+1,  k:k+2))) &
                     ,maxval(abs(f(ll,mm+1,nn,  k:k+2))) &                 
                     ,maxval(abs(f(ll,mm+1,nn+1,k:k+2)))))
        enddo; enddo; enddo
      endif

    endsubroutine staggered_max_vec
!***********************************************************************
    subroutine staggered_max_scal(f,k,jmax,weight)
!
!   Calculates max values of modulus of a scalar quantity in the centre of a grid cell.
!
!   22-dec-15/JW: coded, adapted from staggered_mean_scal
!
      real, dimension (mx,my,mz,mfarray), intent(inout):: f
      integer,                            intent(in)   :: k,jmax
      real,                               intent(in)   :: weight
      integer                                          :: ll,mm,nn
!
      if (dimensionality==3) then
!    
        do ll=2,mx-2; do mm=2,my-2; do nn=2,mz-2
          f(ll,mm,nn,jmax) = max(f(ll,mm,nn,jmax), &
          weight*max(abs(f(ll,mm,nn,      k)) &
                     ,abs(f(ll,mm,nn+1,    k)) &
                     ,abs(f(ll,mm+1,nn,    k)) &
                     ,abs(f(ll,mm+1,nn+1,  k)) &
                     ,abs(f(ll+1,mm,nn,    k)) &
                     ,abs(f(ll+1,mm,nn+1,  k)) &
                     ,abs(f(ll+1,mm+1,nn,  k)) &
                     ,abs(f(ll+1,mm+1,nn+1,k))))
        enddo; enddo; enddo

      elseif (dimensionality==1) then
        if (nxgrid/=1) then
          do ll=2,mx-2; do mm=m1,m2; do nn=n1,n2
            f(ll,mm,nn,jmax) = max(f(ll,mm,nn,jmax), &
            weight*max(abs(f(ll,mm,nn,  k)) &
                       ,abs(f(ll+1,mm,nn,k))))
          enddo; enddo; enddo
        elseif (nygrid/=1) then
          do ll=l1,l2; do mm=2,my-2; do nn=n1,n2
            f(ll,mm,nn,jmax) = max(f(ll,mm,nn,jmax), &
            weight*max(abs(f(ll,mm,nn,  k)) &
                       ,abs(f(ll,mm+1,nn,k))))
          enddo; enddo; enddo
        else
          do ll=l1,l2; do mm=m1,m2; do nn=2,mz-2
            f(ll,mm,nn,jmax) = max(f(ll,mm,nn,jmax), &
            weight*max(abs(f(ll,mm,nn,  k)) &
                       ,abs(f(ll,mm,nn+1,k))))
          enddo; enddo; enddo
        endif
      elseif (nzgrid==1) then   !  x-y
        do ll=2,mx-2; do mm=2,my-2; do nn=n1,n2
          f(ll,mm,nn,jmax) = max(f(ll,mm,nn,jmax), &
          weight*max(abs(f(ll,mm,nn,    k)) &
                     ,abs(f(ll,mm+1,nn,  k)) &
                     ,abs(f(ll+1,mm,nn,  k)) &                 
                     ,abs(f(ll+1,mm+1,nn,k))))
        enddo; enddo; enddo
      elseif (nygrid==1) then   !  x-z
         do ll=2,mx-2; do mm=m1,m2; do nn=2,mz-2
          f(ll,mm,nn,jmax) = max(f(ll,mm,nn,jmax), &
          weight*max(abs(f(ll,mm,nn,    k)) &
                     ,abs(f(ll,mm,nn+1,  k)) &
                     ,abs(f(ll+1,mm,nn,  k)) &                 
                     ,abs(f(ll+1,mm,nn+1,k))))
        enddo; enddo; enddo
      else                      !  y-z
        do ll=l1,l2; do mm=2,my-2; do nn=2,mz-2
          f(ll,mm,nn,jmax) = max(f(ll,mm,nn,jmax), &
          weight*max(abs(f(ll,mm,nn,    k)) &
                     ,abs(f(ll,mm,nn+1,  k)) &
                     ,abs(f(ll,mm+1,nn,  k)) &                 
                     ,abs(f(ll,mm+1,nn+1,k))))
        enddo; enddo; enddo
      endif

    endsubroutine staggered_max_scal
!***********************************************************************
    subroutine directory_names_std(lproc)
!
!  Set up the directory names:
!  set directory name for the output (one subdirectory for each processor)
!  if datadir_snap (where var.dat, VAR# go) is empty, initialize to datadir
!
!  02-oct-2002/wolf: coded
!  23-may-2017/axel: added directory_prestart, not to be erased after start
!
      use Cdata, only: iproc_world, directory, datadir, directory_dist, &
                       datadir_snap, directory_snap, directory_collect, &
                       datadir_prestart, directory_prestart

      logical, optional :: lproc

      character (len=intlen) :: chproc
!
      chproc=itoa(iproc_world)
      call safe_character_assign(directory, trim(datadir)//'/proc'//chproc)
      call safe_character_assign(directory_dist, &
                                            trim(datadir_snap)//'/proc'//chproc)
      call safe_character_assign(directory_prestart, &
                                            trim(datadir_prestart)//'/proc'//chproc)
!
      if (loptest(lproc)) then
        call safe_character_assign(directory_snap, &
                                            trim(datadir_snap)//'/proc'//chproc)
      else
        call safe_character_assign(directory_snap, &
                                            trim(datadir_snap)//'/allprocs')
      endif
      call safe_character_assign(directory_collect, &
                                            trim (datadir_snap)//'/allprocs')
!
    endsubroutine directory_names_std
!****************************************************************************
    character function numeric_precision()
!
!  Return 'S' if running in single, 'D' if running in double precision.
!
!  12-jul-06/wolf: extracted from wdim()
!
      integer :: real_prec
!
      real_prec = precision(1.)
      if (real_prec==6 .or. real_prec==7) then
        numeric_precision = 'S'
      elseif (real_prec == 15) then
        numeric_precision = 'D'
      else
        print*, 'WARNING: encountered unknown precision ', real_prec
        numeric_precision = '?'
      endif
!
    endfunction numeric_precision
!***********************************************************************
    subroutine touch_file(file)
!
!  Touches a given file (used for code locking).
!
!  25-may-03/axel: coded
!  24-mar-10/Bourdin.KIS: moved here from sub.f90
!
      character(len=*) :: file
!
      integer :: unit = 1
!
      open (unit, FILE=file)
      close (unit)
!
    endsubroutine touch_file
!***********************************************************************
    logical function var_is_vec(j)
!
!  Checks whether variable (starting) at slot j in f array is vector field.
!  At the moment magnetic field, including testfields and velocity including
!  testflows are taken into account.
!
! 4-dec-2015/MR: coded
!
      use Cdata, only: iuu,iux,iuz,iaa,iax,iaz,iaatest,ntestfield,iuutest,ntestflow

      integer :: j

      var_is_vec = iuu>0.and.j>=iux.and.j<=iuz .or.  &
                   iaa>0.and.j>=iax.and.j<=iaz .or.  &
                   iaatest>0.and.j>=iaatest.and.j<=iaatest+ntestfield-1 .or. &
                   iuutest>0.and.j>=iuutest.and.j<=iuutest+ntestflow-1

    endfunction var_is_vec
!***********************************************************************
    subroutine transform_cart_spher(f,ith1,ith2,iph1,iph2,j)
!
!  Transforms a vector given in f array in slots j to j+2 on the rectangle
!  ith1:ith2 x iph1:iph2 from Cartesian to spherical basis. Works in-place.
!
! 4-dec-2015/MR: coded
!
      use Cdata, only: cosph, sinph, costh, sinth    !, iproc, lyang

      real, dimension(:,:,:,:) :: f
      integer :: ith1, ith2, iph1, iph2, j

      real, dimension(size(f,1)) :: tmp12,tmp3
      integer :: ith,iph
!      real, dimension(size(f,1)) :: tmp

      do ith=ith1,ith2; do iph=iph1,iph2
!tmp=sum(f(:,ith,iph,j:j+2)**2,2)
        tmp12=cosph(iph)*f(:,ith,iph,j)+sinph(iph)*f(:,ith,iph,j+1)
        tmp3 =f(:,ith,iph,j+2)
        f(:,ith,iph,j+2) = -sinph(iph)*f(:,ith,iph,j)+cosph(iph)*f(:,ith,iph,j+1)
        f(:,ith,iph,j  ) =  sinth(ith)*tmp12 + costh(ith)*tmp3
        f(:,ith,iph,j+1) =  costh(ith)*tmp12 - sinth(ith)*tmp3

!if (any(abs(tmp-sum(f(:,ith,iph,j:j+2)**2,2)) > 1.e-7)) print*, 'iproc,ith,iph=', iproc,ith,iph,maxval(tmp)

      enddo; enddo

    endsubroutine transform_cart_spher
!***********************************************************************
    subroutine transform_cart_spher_other(arr,th,ph)
!
!  Transforms a vector given in f array in slots j to j+2 on the rectangle
!  ith1:ith2 x iph1:iph2 from Cartesian to spherical basis. Works in-place.
!
! 4-dec-2015/MR: coded
!
      real, dimension(nx,3), intent(INOUT) :: arr
      real,                  intent(IN)    :: th,ph

      real, dimension(nx) :: tmp12, tmp3

        tmp12=cos(ph)*arr(:,1)+sin(ph)*arr(:,2)
        tmp3 =arr(:,3)
        arr(:,3) = -sin(ph)*arr(:,1)+ cos(ph)*arr(:,2)
        arr(:,1) =  sin(th)*tmp12   + cos(th)*tmp3
        arr(:,2) =  cos(th)*tmp12   - sin(th)*tmp3

    endsubroutine transform_cart_spher_other
!***********************************************************************
    subroutine transform_spher_cart_yy(f,ith1,ith2,iph1,iph2,dest,lyy)
!
!  Transforms a vector given on the rectangle ith1:ith2 x iph1:ip2 from spherical
!  to Cartesian basis. For lyy=T in addition the Yin-Yang transform is performed.
!  Not in-place capable!
!
! 4-dec-2015/MR: coded
!
      use Cdata, only: cosph, sinph, costh, sinth

      real, dimension(:,:,:,:), intent(IN) :: f
      integer                 , intent(IN) :: ith1, ith2, iph1, iph2
      real, dimension(size(f,1),ith2-ith1+1,iph2-iph1+1,3), intent(OUT) :: dest
      logical, optional, intent(IN) :: lyy

      real, dimension(mx) :: tmp12
      integer :: i,j,itd,ipd,ju,jo

      ju=max(n1,iph1); jo=min(n2,iph2)
      do i=max(m1,ith1),min(m2,ith2); do j=ju,jo

        tmp12=sinth(i)*f(:,i,j,1)+costh(i)*f(:,i,j,2)
        itd=i-ith1+1; ipd=j-iph1+1
        if (loptest(lyy)) then
          dest(:,itd,ipd,1) = -(cosph(j)*tmp12 - sinph(j)*f(:,i,j,3))
          dest(:,itd,ipd,2) = -(costh(i)*f(:,i,j,1)-sinth(i)*f(:,i,j,2))
          dest(:,itd,ipd,3) = -(sinph(j)*tmp12 + cosph(j)*f(:,i,j,3))
        else
          dest(:,itd,ipd,1) = cosph(j)*tmp12 - sinph(j)*f(:,i,j,3)
          dest(:,itd,ipd,2) = sinph(j)*tmp12 + cosph(j)*f(:,i,j,3)
          dest(:,itd,ipd,3) = costh(i)*f(:,i,j,1)-sinth(i)*f(:,i,j,2)
        endif
!if (any(abs(sum(dest(:,itd,ipd,:)**2,2)-sum(f(:,i,j,:)**2,2)) > 2.e-14 )) print*, 'i,j=', i,j,maxval(abs(sum(dest(:,itd,ipd,:)**2,2)-sum(f(:,i,j,:)**2,2)))
      enddo; enddo

    endsubroutine transform_spher_cart_yy
!***********************************************************************
    subroutine yy_transform_strip(ith1,ith2,iph1,iph2,thphprime,ith_shift_,iph_shift_)
!
!  Transform coordinates of ghost zones of Yin or Yang grid to other grid.
!  Strip is defined by index ranges (ith1,ith2), (iph1,iph2) with respect
!  to local grid, in particular sinth, costh etc.
!
!  4-dec-15/MR: coded
! 12-mar-16/MR: entry yy_transform_strip_other added
! 10-sep-18/MR: added parameters ith_shift_,iph_shift_ for creating slanted strips
! 25-oct-18/PABourdin: entry yy_transform_strip_other removed and simplified code
!
      use Cdata, only: costh,sinth,cosph,sinph, y, z

      integer,               intent(IN) :: ith1,ith2,iph1,iph2
      real, dimension(:,:,:),intent(OUT):: thphprime
      integer, optional,     intent(IN) :: iph_shift_, ith_shift_

      integer :: i,j,itp,jtp,iph_shift,ith_shift,ii,jj,ishift,jshift
      real :: sth, cth, xprime, yprime, zprime, sprime
      logical :: ltransp

      if (ith2<ith1.or.iph2<iph1) return

      ltransp = ith2-ith1+1 /= size(thphprime,2)

      iph_shift=ioptest(iph_shift_)        ! default is zero
      ith_shift=ioptest(ith_shift_)        ! default is zero

      thphprime(1,:,:)=0.                  ! to indicate undefined values
!
      jshift=iph_shift
      do i=ith1,ith2
        ishift=ith_shift
        do j=iph1,iph2
!
!  Rotate by Pi about z axis, then by Pi/2 about x axis.
!  No distinction between Yin and Yang as transformation matrix is self-inverse.
!
          jj=j+jshift
          if (jj<1.or.jj>mz) cycle

          ii=i+ishift
          if (ishift/=0) ishift=ishift+ith_shift

          if (ii<1.or.ii>my) cycle

          sth=sinth(ii); cth=costh(ii)
          xprime = -cosph(jj)*sth
          zprime = -sinph(jj)*sth
          yprime = -cth

          sprime = sqrt(xprime**2 + yprime**2)

          if (ltransp) then
            jtp = i-ith1+1; itp = j-iph1+1
          else
            itp = i-ith1+1; jtp = j-iph1+1
          endif

          thphprime(1,itp,jtp) = datan2(dble(sprime),dble(zprime))
          thphprime(2,itp,jtp) = datan2(dble(yprime),dble(xprime))
          if (thphprime(2,itp,jtp)<0.) thphprime(2,itp,jtp) = thphprime(2,itp,jtp) + twopi
!
          !thphprime(1,itp,jtp) = ii   !!!
          !thphprime(2,itp,jtp) = jj   !!!
!
        enddo
        if (jshift/=0) jshift=jshift+iph_shift
      enddo

    endsubroutine yy_transform_strip
!***********************************************************************
    subroutine yy_transform_strip_other(th,ph,thphprime,ith_shift_,iph_shift_)
!
!  Transform coordinates of ghost zones of Yin or Yang grid to other grid.
!  Strip is defined by by vectors th and ph not related to local grid.
!
!  4-dec-15/MR: coded
! 12-mar-16/MR: entry yy_transform_strip_other added
! 10-sep-18/MR: added parameters ith_shift_,iph_shift_ for creating slanted strips
! 25-oct-18/PABourdin: entry yy_transform_strip_other removed and simplified code
!
      use Cdata, only: costh,sinth,cosph,sinph, y, z

      real, dimension(:),    intent(IN) :: th,ph
      real, dimension(:,:,:),intent(OUT):: thphprime
      integer, optional,     intent(IN) :: iph_shift_, ith_shift_

      integer :: i,j,itp,jtp,iph_shift,ith_shift,ii,jj,ishift,jshift
      real :: sth, cth, xprime, yprime, zprime, sprime
      logical :: ltransp

      ltransp = size(th) /= size(thphprime,2)

      iph_shift=ioptest(iph_shift_)        ! default is zero
      ith_shift=ioptest(ith_shift_)        ! default is zero

      thphprime(1,:,:)=0.                  ! to indicate undefined values
!
      jshift=iph_shift
      do i=1,size(th)
        ishift=ith_shift; 
        do j=1,size(ph)
!
!  Rotate by Pi about z axis, then by Pi/2 about x axis.
!  No distinction between Yin and Yang as transformation matrix is self-inverse.
!
          ii=i+ishift
          sth=sin(th(ii)); cth=cos(th(ii))

          jj=j+jshift
          xprime = -cos(ph(jj))*sth
          zprime = -sin(ph(jj))*sth
          yprime = -cth

          sprime = sqrt(xprime**2 + yprime**2)

          if (ltransp) then
            jtp = i; itp = j
          else
            itp = i; jtp = j
          endif

          thphprime(1,itp,jtp) = atan2(sprime,zprime)
          thphprime(2,itp,jtp) = atan2(yprime,xprime)
          if (thphprime(2,itp,jtp)<0.) thphprime(2,itp,jtp) = thphprime(2,itp,jtp) + twopi
!
          !thphprime(1,itp,jtp) = ii   !!!
          !thphprime(2,itp,jtp) = jj   !!!
!
          if (ishift/=0) ishift=ishift+ith_shift
        enddo
        if (jshift/=0) jshift=jshift+iph_shift
      enddo

    endsubroutine yy_transform_strip_other
!***********************************************************************
    subroutine copy_kinked_strip_z(len_cut,jstart_,source,dest,j,shift,leftright,ladd_)
!
! Copies (optionally cumulatively) variable j from rectangular strip source into variable j of kinked strip
! in dest, consisting of a horizontal (y-aligned) part with 
! length len_cut and width nghost and a slanted part with width 2*nghost.
!
! 20-jan-19/MR: coded
!
      use Cdata, only: iproc, iproc_world
      use Cdata, only: lyang,lroot

      real, dimension(:,:,:,:), intent(IN)   :: source
      real, dimension(:,:,:,:), intent(INOUT):: dest
      integer,                  intent(IN)   :: len_cut,jstart_,j,shift
      logical,                  intent(IN)   :: leftright
      logical, optional,        intent(IN)   :: ladd_

      integer :: istart,iend,jstart,jend
      logical :: ladd

      jstart=jstart_
      jend=jstart+nghost-1
      ladd=loptest(ladd_)

      if (leftright) then   !upper and lower left
        istart=1; iend=my-len_cut
        jstart=jstart_+shift*nghost; jend=jstart+nghost-1
        call copy_with_shift(istart,iend,jstart,jend,source(:,my-len_cut+1:,:,j),dest(:,:,:,j),0,shift,ladd) ! outer skew band
        jstart=jstart_; jend=jstart+nghost-1
        call copy_with_shift(istart,iend,jstart,jend,source(:,:,:,j),dest(:,:,:,j),0,shift,ladd) ! main skew band

        if (ladd) then
          dest(:,my-len_cut+1:my,jstart:jend,j) = dest  (:,my-len_cut+1:my,jstart:jend,j) &      ! horizonzal band 
                                                 +source(:,2*(my-len_cut)+1:my,1:nghost,j)
        else
          dest(:,my-len_cut+1:my,jstart:jend,j) = source(:,2*(my-len_cut)+1:my,1:nghost,j)
        endif
!if (notanumber(source(:,my-len_cut+1:my,1:nghost,j))) print*, 'source(:,my-len_cut+1:my,1:nghost,j): iproc,j=', iproc, iproc_world, j
      else                  !upper and lower right
        if (ladd) then
          dest(:,:len_cut,jstart:jend,j) = dest  (:,:len_cut,jstart:jend,j) &                    ! horizonzal band 
                                          +source(:,:len_cut,1:nghost,j)
        else
          dest(:,:len_cut,jstart:jend,j) = source(:,:len_cut,1:nghost,j)
        endif
if (notanumber(source(:,1:len_cut,1:nghost,j))) print*, 'source(:,1:len_cut,1:nghost,j): iproc,j=', iproc, iproc_world, j
        istart=len_cut+1; iend=my
        call copy_with_shift(istart,iend,jstart,jend,source(:,:,:,j),dest(:,:,:,j),0,shift,ladd) ! main skew band
        jstart=jstart_-shift*nghost; jend=jstart+nghost-1
        call copy_with_shift(istart,iend,jstart,jend,source(:,my-len_cut+1:,:,j),dest(:,:,:,j),0,shift,ladd) ! outer skew band
      endif

    endsubroutine copy_kinked_strip_z
!***********************************************************************
    subroutine copy_kinked_strip_y(len_cut,istart_,source,dest,j,shift,leftright,ladd_)
!
! Copies (optionally cumulatively) variable j from rectangular strip source into variable j of kinked strip
! in dest, consisting of a vertical (z-aligned) part with 
! length len_cut and width nghost and a slanted part with width 2*nghost.
!
! 20-jan-19/MR: coded
!
      use Cdata, only: iproc, iproc_world

      real, dimension(:,:,:,:), intent(IN)   :: source
      real, dimension(:,:,:,:), intent(INOUT):: dest
      integer,                  intent(IN)   :: len_cut,istart_,j,shift
      logical,                  intent(IN)   :: leftright
      logical, optional,        intent(IN)   :: ladd_

      integer :: istart,iend,jstart,jend
      logical :: ladd

      istart=istart_
      iend=istart+nghost-1
      ladd=loptest(ladd_)
      
      if (leftright) then   !upper and lower left
        jstart=1; jend=mz-len_cut
        istart=istart_+shift*nghost; iend=istart+nghost-1
        call copy_with_shift(istart,iend,jstart,jend,source(:,:,mz-len_cut+1:,j),dest(:,:,:,j),shift,0,ladd)   ! outer skew band
        istart=istart_; iend=istart+nghost-1
        call copy_with_shift(istart,iend,jstart,jend,source(:,:,:,j),dest(:,:,:,j),shift,0,ladd)   ! main skew band
        if (ladd) then                                                                             ! vertical band
          dest(:,istart:iend,mz-len_cut+1:mz,j) = dest  (:,istart:iend,mz-len_cut+1:mz,j) &
                                                 +source(:,1:nghost,2*(mz-len_cut)+1:mz,j)
        else
          dest(:,istart:iend,mz-len_cut+1:mz,j) = source(:,1:nghost,2*(mz-len_cut)+1:mz,j)
        endif
!if (notanumber(source(:,1:nghost,mz-len_cut+1:mz,j))) print*, 'source(:,1:nghost,mz-len_cut+1:mz,j): iproc,j=', iproc, iproc_world, j
      else                  !upper and lower right
        jstart=len_cut+1; jend=mz
        if (ladd) then
          dest(:,istart:iend,1:len_cut,j) = dest  (:,istart:iend,1:len_cut,j) &                    ! vertical band
                                           +source(:,1:nghost,1:len_cut,j)
        else
          dest(:,istart:iend,1:len_cut,j) = source(:,1:nghost,1:len_cut,j)
        endif
if (notanumber(source(:,1:nghost,1:len_cut,j))) print*, 'source(:,1:nghost,1:len_cut,j): iproc,j=', iproc, iproc_world,j
        jstart=len_cut+1; jend=mz
        call copy_with_shift(istart,iend,jstart,jend,source(:,:,:,j),dest(:,:,:,j),shift,0,ladd)   ! main skew band
        istart=istart_-shift*nghost; iend=istart+nghost-1
        call copy_with_shift(istart,iend,jstart,jend,source(:,:,mz-len_cut+1:,j),dest(:,:,:,j),shift,0,ladd)   ! outer skew band
      endif

    endsubroutine copy_kinked_strip_y
!***********************************************************************
    subroutine copy_with_shift(i1, i2, j1, j2, source, dest, ishift_, jshift_,ladd)
!
! Copies (optionally cumulatively) into slanted strip.
! Sets i1,i2 (or j1,j2) to last interval with shift.
!
! 20-sep-18/MR: coded
!
      use Cdata, only: iproc, iproc_world

      real, dimension(:,:,:), intent(IN)   :: source
      real, dimension(:,:,:), intent(INOUT):: dest
      integer,                intent(INOUT):: i1, i2, j1, j2
      integer,                intent(IN)   :: ishift_, jshift_
      logical, optional,      intent(IN)   :: ladd

      integer :: ishift, jshift, i, j, ii, jj, id, jd, is, js, iimax, jjmax

      iimax=size(dest,2); jjmax=size(dest,3)

      jshift=jshift_
      if (ishift_/=0) then; is=1; else; is=i1; endif

      do i=i1,i2
        ishift=ishift_

        if (jshift_/=0) then; js=1; else; js=j1; endif

        do j=j1,j2
!
          ii=i+ishift
          jj=j+jshift

          if (ii>=1.and.jj>=1.and.ii<=iimax.and.jj<=jjmax) then
            if (loptest(ladd)) then
              dest(:,ii,jj) = dest(:,ii,jj) + source(:,is,js)
            else
              dest(:,ii,jj) = source(:,is,js)
            endif
          endif
if (notanumber(source(:,is,js))) print*, 'source(:,is,js): iproc,j=', iproc, iproc_world, j
          js=js+1
!
          if (ishift/=0) ishift=ishift+ishift_
        enddo
        is=is+1
        if (jshift/=0) jshift=jshift+jshift_
      enddo

      if (ishift_/=0) then; id=i2-i1; i2=ii+ishift_; i1=i2-id; endif 
      if (jshift_/=0) then; jd=j2-j1; j2=jj+jshift_; j1=j2-jd; endif 
!
    endsubroutine copy_with_shift
!***********************************************************************
    subroutine reset_triangle(i1,i2,j1,j2,f)
!
!  Sets triangular area of float array f, defined by i1,i2,j1,j2, to zero.
!  "Right angle" of triangle is at (i1,j1).
! 
!  28-jan-2019/MR: coded
!
    integer,                intent(IN) :: i1,i2,j1,j2
    real, dimension(:,:,:), intent(OUT):: f

    integer :: istep, jstep, i, j

      istep = sign(1,i2-i1); jstep = sign(1,j2-j1)

      do i=i1,i2,istep
        do j=j1,j2,jstep
          if ( jstep*(j-j2)+istep*(i-i1) <= 0 ) then
            f(:,i,j)=0.
          endif
        enddo
      enddo

    endsubroutine reset_triangle
!***********************************************************************
    subroutine reset_triangle_inds(i1,i2,j1,j2,inds_arr)
!
!  Sets triangular area of integer array inds_arr, defined by i1,i2,j1,j2, to zero.
!  "Right angle" of triangle is at (i1,j1).
! 
!  28-jan-2019/MR: coded
!
    integer,                 intent(IN) :: i1,i2,j1,j2
    integer, dimension(:,:), intent(OUT):: inds_arr

    integer :: istep, jstep, i, j

      istep = sign(1,i2-i1); jstep = sign(1,j2-j1)

      do i=i1,i2,istep
        do j=j1,j2,jstep
          if ( jstep*(j-j2)+istep*(i-i1) <= 0 ) then
            inds_arr(i,j)=0
          endif
        enddo
      enddo

    endsubroutine reset_triangle_inds
!***********************************************************************
    subroutine transform_thph_yy( vec, powers, transformed, theta, phi )
!
!  Transforms theta and phi components of vector vec defined with the Yang
!  grid basis to the Yin grid basis using theta and phi coordinates of the Yang grid.
!  Without optional parameters theta, phi: for use on pencils within mn-loop.
!  Note that components of transformed are undefined if corresponding power
!  in mask powers is 0.
!  If theta, phi are present, it is for use outside mn-loop or for other
!  coordinates than y(m), z(n).
!
! 30-mar-2016/MR: coded
! 28-jun-2016/MR: added optional parameters
! 28-aug-2017/MR: cared for singularity at pole
!
      use Cdata, only: costh, sinth, cosph, sinph, m, n

      real,    dimension(:,:),intent(IN) :: vec
      integer, dimension(3),  intent(IN) :: powers
      real,    dimension(:,:),intent(OUT):: transformed
      real, optional,         intent(IN) :: theta, phi

      real :: sinth1, a, b, sinp, sisisq

      if (any(powers(2:3)/=0)) then
        if (present(theta)) then
          sinp=sin(phi)
          sisisq=sqrt(1.-(sinp*sin(theta))**2)
          if (sisisq==0.) then                 ! i.e. at pole of other grid -> theta and phi components indefined
            a=0.; b=0.
          else
            sinth1=1./sisisq
            a=cos(phi)*sinth1; b=sinp*cos(theta)*sinth1
          endif
        else
          sisisq=sqrt(1.-(sinth(m)*sinph(n))**2)
          if (sisisq==0.) then
            a=0.; b=0.
          else
            sinth1=1./sisisq
            a=cosph(n)*sinth1; b=sinph(n)*costh(m)*sinth1
          endif
        endif
      endif

      if (powers(1)/=0) &
        transformed(:,1) = vec(:,1)
      if (powers(2)/=0) &
        transformed(:,2) = b*vec(:,2) + a*vec(:,3)
      if (powers(3)/=0) &
        transformed(:,3) =-a*vec(:,2) + b*vec(:,3)

    endsubroutine transform_thph_yy
!***********************************************************************
    subroutine transform_thph_yy_other( vec,indthl,indthu,indphl,indphu,transformed )
!
!  Transforms theta and phi components of a vector field vec defined with the Yang grid basis
!  to the Yin grid basis using theta and phi coordinates of the Yang grid.
!  The r component is transfered unaltered.
!  Both vec and vec_transformed must have shape (dimx,ny,nz,3).
!  For use outside mn-loop.
!
! 30-mar-2016/MR: coded
! 28-aug-2017/MR: cared for singularity at pole
! 30-jan-2019/MR: index bounds now parameters; transferes r component unaltered
!
      use Cdata, only: costh, sinth, cosph, sinph

      real, dimension(:,:,:,:), intent(IN) :: vec
      integer,                  intent(IN) :: indthl,indthu,indphl,indphu
      real, dimension(:,:,:,:), intent(OUT):: transformed

      integer :: ig,jg,jl,ju
      real :: sisisq,sinth1,a,b

      jl=max(m1,indthl); ju=min(m2,indthu)

      do ig=max(n1,indphl),min(n2,indphu)
        do jg=jl,ju

          sisisq=sqrt(1.-(sinth(jg)*sinph(ig))**2)
          if (sisisq==0.) then                     ! i.e. at pole of other grid -> theta and phi components indefined
            a=0.; b=0.
          else
            sinth1=1./sisisq
            a=cosph(ig)*sinth1; b=sinph(ig)*costh(jg)*sinth1
          endif

          transformed(:,jg-indthl+1,ig-indphl+1,1) =   vec(:,jg,ig,1)
          transformed(:,jg-indthl+1,ig-indphl+1,2) = b*vec(:,jg,ig,2) + a*vec(:,jg,ig,3)
          transformed(:,jg-indthl+1,ig-indphl+1,3) =-a*vec(:,jg,ig,2) + b*vec(:,jg,ig,3)

        enddo
      enddo

    endsubroutine transform_thph_yy_other
!***********************************************************************
    subroutine transpose_mn(a,b,ladd)
!
!  Transpose mxn pencil array, b=transpose(a), on pencil arrays.
!
!   7-aug-10/dhruba: coded
!
      use Cdata, only: lroot
!
      real, dimension(:,:,:) :: a,b
      logical, optional :: ladd
!
      intent(in) :: a
      intent(out) :: b
      integer :: i,j
      logical :: ladd_
!
      if (size(a,2)/=size(b,3) .or. size(a,3)/=size(b,2)) then
        if (lroot) print*, 'ERROR -- transpose_mn: inconsistent dimensions of a and b!'
        stop
      endif

      ladd_=loptest(ladd)
      do i=1,size(b,2)
        do j=1,size(b,3)
          if (ladd_) then
            b(:,i,j)=b(:,i,j)+a(:,j,i)
          else
            b(:,i,j)=a(:,j,i)
          endif
        enddo
      enddo
!
    endsubroutine transpose_mn
!***********************************************************************
    function notanumber_0(f)
!
!  Check for NaN or Inf values.
!  Not well tested with all compilers and options, but avoids false
!  positives in a case where the previous implementation had problems
!  Version for scalars
!
!  22-Jul-11/sven+philippe: coded
!
      logical :: notanumber_0
      real :: f
!
      notanumber_0 = .not. ((f <= huge(f)) .or. (f > huge_real))
!
    endfunction notanumber_0
!***********************************************************************
   function notanumber_0d(f)
!
!  Check for NaN or Inf values.
!  Not well tested with all compilers and options, but avoids false
!  positives in a case where the previous implementation had problems
!  Version for scalars
!
!  27-Jul-15/MR: adapted
!
     logical :: notanumber_0d
     double precision :: f
!
     notanumber_0d = .not. ((f <= huge(f)) .or. (f > huge_double))
!
    endfunction notanumber_0d
!***********************************************************************
    function notanumber_1(f)
!
!  Check for NaN or Inf values.
!  Not well tested with all compilers and options, but avoids false
!  positives in a case where the previous implementation had problems
!  Version for 1-d arrays
!
!  22-Jul-11/sven+philippe: coded
!
      logical :: notanumber_1
      real, dimension(:) :: f
!
      notanumber_1 = any(.not. ((f <= huge(f)) .or. (f > huge_real)))
!
    endfunction notanumber_1
!***********************************************************************
    function notanumber_2(f)
!
!  Check for NaN or Inf values.
!  Not well tested with all compilers and options, but avoids false
!  positives in a case where the previous implementation had problems
!  Version for 2-d arrays
!
!  22-Jul-11/sven+philippe: coded
!
      logical :: notanumber_2
      real, dimension(:,:) :: f
!
      notanumber_2 = any(.not. ((f <= huge(f)) .or. (f > huge_real)))
!
    endfunction notanumber_2
!***********************************************************************
    function notanumber_3(f)
!
!  Check for NaN or Inf values.
!  Not well tested with all compilers and options, but avoids false
!  positives in a case where the previous implementation had problems
!  Version for 3-d arrays
!
!  22-Jul-11/sven+philippe: coded
!
      logical :: notanumber_3
      real, dimension(:,:,:) :: f
!
      notanumber_3 = any(.not. ((f <= huge(f)) .or. (f > huge_real)))
!
    endfunction notanumber_3
!***********************************************************************
    function notanumber_4(f)
!
!  Check for NaN or Inf values.
!  Not well tested with all compilers and options, but avoids false
!  positives in a case where the previous implementation had problems
!  Version for 4-d arrays
!
!  22-Jul-11/sven+philippe: coded
!
      logical :: notanumber_4
      real, dimension(:,:,:,:) :: f
!
      notanumber_4 = any(.not. ((f <= huge(f)) .or. (f > huge_real)))
!
    endfunction notanumber_4
!***********************************************************************
    function notanumber_5(f)
!
!  Check for NaN or Inf values.
!  Not well tested with all compilers and options, but avoids false
!  positives in a case where the previous implementation had problems
!  Version for 4-d arrays
!
!  22-Jul-11/sven+philippe: coded
!
      logical :: notanumber_5
      real, dimension(:,:,:,:,:) :: f
!
      notanumber_5 = any(.not. ((f <= huge(f)) .or. (f > huge_real)))
!
    endfunction notanumber_5
!***********************************************************************
    subroutine reduce_grad_dim(g)
!
!  Compresses a gradient vector according to dimensionality using precalculated
!  dim_mask.
!
!  28-Jan-16/MR: coded
!
      use Cdata, only: dim_mask

      real, dimension(:,:) :: g

      if (dimensionality==3) return

      g(:,1:dimensionality)=g(:,dim_mask(1:dimensionality))

    endsubroutine reduce_grad_dim
!****************************************************************************
    function merge_yin_yang(y,z,dy,dz,yzyang,yz,inds) result (nok)
!
!  Merges Yin and Yang grids: Yin, given by y,z,dy,dz, remains unchanged while Yang,
!  given by yzyang of dimension 2 x size(y)*size(z) which is assumed to be
!  transformed into Yin basis, is clipped by removing points which lie within Yin.
!  Output is merged coordinate array yz[2,*] and index vector inds into
!  yzyang selecting the unclipped points. inds can be used to merge data arrays accordingly.
!  Returns number of unclipped points in Yang.
!
!  30-mar-16/MR: cloned from corresp. IDL routine
!
    real, dimension(:),    intent(IN) :: y,z
    real,                  intent(IN) :: dy,dz
    real, dimension(:,:) , intent(IN) :: yzyang
    real, dimension(:,:) , intent(OUT):: yz
    integer, dimension(:), intent(OUT):: inds
    integer :: nok

    integer :: ind,i,ny,nz,nyz

    ny=size(y); nz=size(z); nyz=ny*nz
!
!  Determine indices of unclipped points of Yang in yzyang.
!
    nok=0
    do i=1,nyz
      if (yzyang(1,i) < minval(y)-dy .or. yzyang(1,i) > maxval(y)+dy .or. &
          yzyang(2,i) < minval(z)-dz .or. yzyang(2,i) > maxval(z)+dz) then
        nok=nok+1
        inds(nok)=i
      endif
    enddo
!
!  Put Yin coordinates in 2D array.
!
    ind=1
    do i=1,ny
      yz(1,ind:ind+nz-1) = y(i)
      yz(2,ind:ind+nz-1) = z
      ind=ind+nz
    enddo
!
!  Add Yang coordinates.
!
    yz(:,nyz+1:nyz+nok)=yzyang(:,inds(:nok))

  endfunction merge_yin_yang
!****************************************************************************
  subroutine yin2yang_coors(costh,sinth,cosph,sinph,yz)
!
!  Transforms (theta,phi) coordinates of Yin or Yang grid, given by cos(theta)
!  etc., into (theta',phi') coordinates of the other grid, stored in array yz
!  of dimension 2 x size(theta)*size(phi).
!
!  30-mar-16/MR: cloned from corresp. IDL routine
!
      real, dimension(:)  , intent(IN) :: sinth, costh, sinph, cosph
      real, dimension(:,:), intent(OUT):: yz

      integer :: ind, i, j
      real :: sth, cth, xprime, yprime, zprime, sprime

      ind=1
      do i=1,size(sinth)

        sth=sinth(i); cth=costh(i)

        do j=1,size(sinph)
!
!  Rotate by Pi about z axis, then by Pi/2 about x axis.
!  No distinction between Yin and Yang as transformation matrix is self-inverse.
!
          xprime = -cosph(j)*sth
          yprime = -cth
          zprime = -sinph(j)*sth

          sprime = sqrt(xprime**2 + yprime**2)

          yz(1,ind) = atan2(sprime,zprime)
          yz(2,ind) = atan2(yprime,xprime)
          if (yz(2,ind) < 0.) yz(2,ind) = yz(2,ind)+twopi
          ind = ind+1

        enddo
      enddo

  endsubroutine yin2yang_coors
!****************************************************************************
    subroutine meshgrid_2d(array1,array2,output_array1,output_array2)
!
! extends 1d arrays array1 and array2 into 2d arrays of size
! length(array1)xlength(array2). Analagous to numpy.meshgrid
!
! 19-aug-16/vince: adapted
! 16-jun-17/MR: double -> real
!
      real, intent(in), dimension(:) :: array1
      real, intent(in), dimension(:) :: array2
      real, intent(out), dimension(:,:) :: output_array1
      real, intent(out), dimension(:,:) :: output_array2
!    
      output_array1=spread(array1,2,size(array2))
      output_array2=spread(array2,1,size(array1))
!
    endsubroutine meshgrid_2d
!***********************************************************************
    subroutine meshgrid_3d(array1,array2,array3,output_array1,output_array2,output_array3)
!
! extends 1d arrays array1, array2, and array3 into 3d arrays
! of shape [size(array1),size(array2),size(array3)].
!
      real, intent(in), dimension(:) :: array1
      real, intent(in), dimension(:) :: array2
      real, intent(in), dimension(:) :: array3
      real, intent(out), dimension(:,:,:) :: output_array1
      real, intent(out), dimension(:,:,:) :: output_array2
      real, intent(out), dimension(:,:,:) :: output_array3
!
      output_array1(:,:,1)=spread(array1,2,size(array2))
      output_array1=spread(output_array1(:,:,1),3,size(array3))
!
      output_array2(:,:,1)=spread(array2,1,size(array1))
      output_array2=spread(output_array2(:,:,1),3,size(array3))
!
      output_array3(:,1,:)=spread(array3,1,size(array1))
      output_array3=spread(output_array3(:,1,:),2,size(array2))
!
    endsubroutine meshgrid_3d
!***********************************************************************
    function linspace(start,end,n,step_size,endpoint)
!
! Returns a vector of n linearly spaced values on the interval [start,end]
! if endpoint=.true. (default) or [start,end) if endpoint=.false.
! If a step_size argument is provided, returns the distance between adjacent
! values.
! Analagous to numpy.linspace
!
! 19-aug-16/vince: adapted
!
      real, intent(in) :: start,end
      logical, intent(in), optional :: endpoint !whether or not to include the endpoint
      integer, intent(in) :: n
      real, intent(out), optional :: step_size !can choose to get the step size in the array
      real, dimension(n) :: linspace
      integer :: i
      real :: step
      logical :: include_endpoint
!
        !by default, include the endpoint
      if(present(endpoint)) then
         include_endpoint=endpoint
      else
         include_endpoint=.true.
      end if
!
! the number of elements in the array is fixed, so the step size determines whether or not
! the endpoint is included as the last element of the array
!
      if(include_endpoint) then
         step=(end-start)/(n-1)
      else
         step=(end-start)/n
      end if
!
! if a step size variable is provided, return the step size used
!
      if (present(step_size)) then
        step_size=step
      end if
!
! fill the array
!
      do i=1,n
        linspace(i)=(start+step*(i-1))
      end do
!
    endfunction linspace
!***********************************************************************
    subroutine newton_arr(x,func,add,fmax,dxmax,itmax)
!
!  Newton iteration for a 2d-array of decoupled nonlinear equations.
!  Not yet tested.
!
!  12-sep-16/MR: coded
!
      real, dimension(:,:),           intent(INOUT) :: x
      !external :: func
      real, dimension(:,:), optional, intent(IN) :: add
      real,                 optional, intent(IN) :: fmax, dxmax
      integer,              optional, intent(IN) :: itmax

      interface
        pure subroutine func(x,f,df)
          real, intent(in) :: x
          real, intent(out):: f, df
        endsubroutine
      endinterface

      real, parameter :: damp=0.1
      integer :: itmax_, it, i, j
      real :: fmax_, dxmax_, dXi, fXi, dfdXi
      logical :: ladd

      fmax_ =roptest(fmax,1e-5)
      dxmax_=roptest(dxmax,1e-4)
      itmax_=ioptest(itmax,1000)
      ladd  =present(add)

      do i=1,size(x,1); do j=1,size(x,2)

        it=0;  fXi=fmax_
        do while (abs(fXi)>=fmax_)
!  
          call func(x(i,j),fXi,dfdXi)
          if (ladd) fXi=fXi+add(i,j)
          dXi=damp*fXi/dfdXi
          if (dXi<dxmax) exit
          x(i,j)=x(i,j)-dXi
!
          it=it+1
          if (it>=itmax_) exit
!
        enddo

      enddo; enddo
      
    endsubroutine newton_arr
!***********************************************************************
    integer function ld(ix)
!
!  Simple dyadic logarithm for integer argument ix without rounding.
!  Return value is negative if ix is not a power of 2,
!  and equal to impossible_int if ix<=0.
!
!  10-dec-17/MR: coded
!
      integer :: ix
      real :: q

      if (ix<=0) then
        ld=impossible_int
        return
      endif

      ld=0; q=ix
      do 

        q=ix/2
        if (q<1.) exit
        ld=ld+1

      enddo
 
      if (2**ld /= ix) ld=-ld 

    endfunction ld
!****************************************************************************** 
    function get_species_nr(name,stem,maxnr,caller) result (ispec)
!
! Extracts species number of a multi-species quantity from name if base name of quantity is
! stem and maximum legal number is maxnr. caller is the calling routine.
! If no number conatined, 1 is returned.
!
! 10-may-18/MR: coded
!
    use Cdata, only: lroot

    character (LEN=*) name, stem, caller
    integer :: maxnr
    integer :: ispec

    integer :: lenstem, ios

    lenstem=len(stem)

    if (name(lenstem+1:)==' ') then
      ispec=1
    else
      ispec=0
      read(name(lenstem+1:),'(i3)',iostat=ios) ispec
      if (ios/=0) then
        if (lroot) &
          print*, trim(caller)//': Warning - unreadable species number in slice name "'// &                        
                  trim(name)//'"'
        return
      endif
    endif

    if (ispec>maxnr) then
      if (lroot) &
        print*, trim(caller)//': Warning - species number in slice name "'// &
                trim(name)//'" bigger than maximum '//trim(itoa(maxnr))
      ispec=0
    endif

    endfunction get_species_nr
!***********************************************************************
    subroutine convert_nml(str,lvec)
!
!  Converts a string of the form *<integer number>[TtFf] into a vector of
!  logicals.
!   
!  12-jun-19/MR: coded

      character(LEN=*) :: str
      logical, dimension(:) :: lvec

      character :: ch
      integer :: i,j,datapos,mult

      j=0; i=1
      do
        ch=str(i:i)
        if (ch=='*') then
          datapos=scan(str(i+1:),'TtFf')
          read(str(i+1:i+datapos-1),*) mult
          lvec(j+1:j+mult) = str(datapos:datapos)=='T' .or. str(datapos:datapos)=='t'
          j=j+mult
          i=datapos
        elseif (ch=='T'.or.ch=='t') then
          j=j+1
          lvec(j)=.true.
        elseif (ch=='F'.or.ch=='f') then
          j=j+1
          lvec(j)=.false.
        else !invalid
        endif
        i=i+1
        if (i>len(str)) exit 
      enddo

    endsubroutine convert_nml
!***********************************************************************
    subroutine extract_from_nml(name,nml,lvec)
!
! Extracts (greps) data item with name "name" from a namelist file (default: data/param2.nml)
! and stores result in file "tmp"
! 
! 12-jun-19/MR: coded
!
      use Syscalls, only: system_cmd

      character(LEN=*)           :: name
      character(LEN=*), optional :: nml
      logical,          optional :: lvec

      character(LEN=256) :: cmd

      cmd="grep '"//trim(name)//" *=' "//trim(coptest(nml,'data/param2.nml'))//" | sed -e's/^.*= *//'"
      if (loptest(lvec)) cmd=trim(cmd)//" -e's/\([1-9][0-9]*\)\*/*\1/g'"
      cmd=trim(cmd)//" > tmp"
      call system_cmd(cmd)

    endsubroutine extract_from_nml
!***********************************************************************
    function get_from_nml_str(name,lfound,nml,lvec) result (res)
!
! Returns the value of a string variable with name "name", previously extracted from a namelist file
! and stored in file "tmp"
! 
! 12-jun-19/MR: coded
!
      character(LEN=*)           :: name
      logical                    :: lfound
      character(LEN=*), optional :: nml
      logical,          optional :: lvec
      character(LEN=128) :: res
      
      call extract_from_nml(name,nml,lvec)
 
      res=''
      inquire(file='tmp',exist=lfound)
      if (lfound) then
        open(1,file='tmp')
        read(1,*,err=1,end=1) res
        goto 2
1       lfound=.false.
2       close(1,status='delete')
      endif

    endfunction get_from_nml_str
!***********************************************************************
    function get_from_nml_log(name,lfound,nml,lvec) result (res)
!
! Returns the value of a logical variable with name "name", previously extracted from a namelist file
! and stored in file "tmp"
! 
! 12-jun-19/MR: coded
!
      character(LEN=*)           :: name
      logical                    :: lfound
      character(LEN=*), optional :: nml
      logical,          optional :: lvec
      logical :: res

      call extract_from_nml(name,nml,lvec)

      res=.false.
      inquire(file='tmp',exist=lfound)
      if (lfound) then
        open(1,file='tmp')
        read(1,*,err=1,end=1) res
        goto 2
1       lfound=.false.
2       close(1,status='delete')
      endif

    endfunction get_from_nml_log
!***********************************************************************
    function get_from_nml_real(name,lfound,nml,lvec) result (res)
!
! Returns the value of a real variable with name "name", previously extracted from a namelist file
! and stored in file "tmp"
! 
! 12-jun-19/MR: coded
!
      character(LEN=*)           :: name
      logical                    :: lfound
      character(LEN=*), optional :: nml
      logical,          optional :: lvec
      real :: res

      call extract_from_nml(name,nml,lvec)

      res=0.
      inquire(file='tmp',exist=lfound)
      if (lfound) then
        open(1,file='tmp')
        read(1,*,err=1,end=1) res
        goto 2
1       lfound=.false.
2       close(1,status='delete')
      endif

    endfunction get_from_nml_real
!***********************************************************************
  endmodule General
