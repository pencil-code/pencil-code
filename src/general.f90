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
  public :: safe_character_assign,safe_character_append,safe_character_prepend
  public :: random_seed_wrapper
  public :: random_number_wrapper, random_gen, normal_deviate
  public :: parse_filename
!
  public :: setup_mm_nn
  public :: input_persistent_general, output_persistent_general
  public :: find_index_range
!
  public :: spline,tridag,pendag,complex_phase,erfcc
  public :: besselj_nu_int,calc_complete_ellints
  public :: bessj,cyclic
  public :: spline_integral,linear_interpolate
  public :: chn, parser, write_full_columns
  public :: read_range, merge_ranges, get_range_no, write_by_ranges, &
            write_by_ranges_1d_real, write_by_ranges_1d_cmplx, &
            write_by_ranges_2d_real, write_by_ranges_2d_cmplx
!
  include 'record_types.h'
!
  interface random_number_wrapper
    module procedure random_number_wrapper_0
    module procedure random_number_wrapper_1
    module procedure random_number_wrapper_3
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
  interface write_by_ranges
    module procedure write_by_ranges_1d_real
    module procedure write_by_ranges_1d_cmplx
    module procedure write_by_ranges_2d_real
    module procedure write_by_ranges_2d_cmplx
  endinterface
!
!  State and default generator of random numbers.
!
  integer, save, dimension(mseed) :: rstate=0
  character (len=labellen) :: random_gen='min_std'
!
  contains
!***********************************************************************
    subroutine setup_mm_nn()
!
!  Produce index-array for the sequence of points to be worked through:
!  Before the communication has been completed, the nghost=3 layers next
!  to the processor boundary (m1, m2, n1, or n2) cannot be used yet.
!  In the mean time we can calculate the interior points sufficiently far
!  away from the boundary points. Here we calculate the order in which
!  m and n are executed. At one point, necessary(imn)=.true., which is
!  the moment when all communication must be completed.
!
      use Cdata, only: mm,nn,imn_array,necessary,lroot
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
!
!  Debugging output to be analysed with $PENCIL_HOME/utils/check-mm-nn.
!  Uncommenting manually, since we can't use ip here (as it is not yet
!  read from run.in).
!
      if (.false.) then
        if (lroot) then
          do imn=1,ny*nz
            if (necessary(imn)) write(*,'(A)') '==MM==NN==> Necessary'
            write(*,'(A,I5,I5)') '==MM==NN==> m,n= ', mm(imn), nn(imn)
          enddo
        endif
      endif
!
    endsubroutine setup_mm_nn
!***********************************************************************
    subroutine input_persistent_general(id,lun,done)
!
!  Reads seed from a snapshot.
!
      use Cdata, only: seed,nseed
!
      integer :: id,lun
      logical :: done
!
      call random_seed_wrapper(GET=seed)
      if (id==id_record_RANDOM_SEEDS) then
        read (lun) seed(1:nseed)
        call random_seed_wrapper(PUT=seed)
        done=.true.
      endif
!
    endsubroutine input_persistent_general
!***********************************************************************
    subroutine output_persistent_general(lun)
!
!  Writes seed to a snapshot.
!
      use Cdata, only: seed,nseed
!
      integer :: lun
!
      call random_seed_wrapper(GET=seed)
      write (lun) id_record_RANDOM_SEEDS
      write (lun) seed(1:nseed)
!
    endsubroutine output_persistent_general
!***********************************************************************
    subroutine random_number_wrapper_0(a)
!
!  Fills a with a random number calculated with one of the generators
!  available with random_gen.
!
      real :: a
      real, dimension(1) :: b
!
      intent(out) :: a
!
!     b = a                     ! not needed unless numbers are non-Markovian
!
      call random_number_wrapper(b)
      a = b(1)
!
    endsubroutine random_number_wrapper_0
!***********************************************************************
    subroutine random_number_wrapper_1(a)
!
!  Fills a with an array of random numbers calculated with one of the
!  generators available with random_gen.
!
      use Cdata, only: lroot
!
      real, dimension(:) :: a
      integer :: i
!
      intent(out) :: a
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
        do i=1,size(a,1)
          a(i)=mars_ran()
        enddo
      case default
        if (lroot) print*, 'No such random number generator: ', random_gen
        STOP 1                ! Return nonzero exit status
!
     endselect
!
    endsubroutine random_number_wrapper_1
!***********************************************************************
    subroutine random_number_wrapper_3(a)
!
!  Fills a with a matrix of random numbers calculated with one of the
!  generators available with random_gen.
!
      use Cdata, only: lroot
!
      real, dimension(:,:,:) :: a
      integer :: i,j,k
!
      intent(out) :: a
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
    subroutine random_seed_wrapper(size,put,get)
!
!  Mimics the f90 random_seed routine.
!
      integer, optional :: size
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
        if (present(get)) then
          get(1:nseed)=rstate(1:nseed)
        endif
        if (present(put)) then
          if (put(2)==0) then   ! state cannot be result from previous
                                ! call, so initialize
            dummy = mars_ran(put(1))
          else
            rstate(1:nseed)=put(1:nseed)
          endif
        endif
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
      integer, parameter :: ia=16807,im=2147483647,iq=127773,ir=2836, &
           mask=123459876
      real, parameter :: am=1./im
      real :: ran0
      integer :: dummy,k
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
!  26-sep-02/wolf: Implemented, following `Numerical Recipes for F90'
!                  ran() routine
!
      implicit none
!
      real :: mars_ran
      real, save :: am=impossible    ! will be constant on a given platform
      integer, optional, intent(in) :: init
      integer, parameter :: ia=16807,im=2147483647,iq=127773,ir=2836
      integer :: k,init1=1812   ! default value
      logical, save :: first_call=.true.
!
!ajwm This doesn't appear to always get set!
      if (first_call) then
        am=nearest(1.0,-1.0)/im
        first_call=.false.
      endif
      if (present(init) .or. rstate(1)==0 .or. rstate(2)<=0) then
!
!  Initialize.
!
        if (present(init)) init1 = init
        am=nearest(1.0,-1.0)/im
        rstate(1)=ieor(777755555,abs(init1))
        rstate(2)=ior(ieor(888889999,abs(init1)),1)
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
    subroutine chn(n,ch,label)
!
!  Make a character out of a number.
!  Take care of numbers that have less than 4 digits.
!
!  30-sep-97/axel: coded
!
      character (len=5) :: ch
      character (len=*), optional :: label
      integer :: n
!
      intent(in) :: n,label
      intent(out) :: ch
!
      ch='     '
      if (n<0) stop 'chn: lt1'
      if (n<10) then
        write(ch(1:1),'(i1)') n
      elseif (n<100) then
        write(ch(1:2),'(i2)') n
      elseif (n<1000) then
        write(ch(1:3),'(i3)') n
      elseif (n<10000) then
        write(ch(1:4),'(i4)') n
      elseif (n<100000) then
        write(ch(1:5),'(i5)') n
      else
        if (present(label)) print*, 'CHN: <', label, '>'
        print*,'CHN: n=',n
        stop "CHN: n too large"
      endif
!
    endsubroutine chn
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
      character (len=*), intent(inout):: dest
      integer :: destLen, srcLen
!
      destLen = len(dest)
      srcLen = len(src)
!
      if (destLen<srcLen) then
        print*, "safe_character_assign: ", &
            "RUNTIME ERROR: FORCED STRING TRUNCATION WHEN ASSIGNING '" &
             //src//"' to '"//dest//"'"
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
      return
!
      entry safe_character_prepend_2(str1,str2)
!
      call safe_character_assign(str1, trim(str2) // trim(str1))
!
    endsubroutine safe_character_append_2
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
    subroutine input_array(file,a,dimx,dimy,dimz,dimv)
!
!  Generalized form of input, allows specifying dimension.
!
!  27-sep-03/axel: coded
!
      character (len=*) :: file
      integer :: dimx,dimy,dimz,dimv
      real, dimension (dimx,dimy,dimz,dimv) :: a
!
      open(1,FILE=file,FORM='unformatted')
      read(1) a
      close(1)
!
    endsubroutine input_array
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
        goto 99
      endif
!
!  Find lower index.
!
      ii1=naa
      do ii=1,naa
        if (aa(ii)>=aa1) then
          ii1=ii
          goto 10
        endif
      enddo
!
!  Find upper index.
!
10    ii2=1
      do ii=naa,1,-1
        if (aa(ii)<=aa2) then
          ii2=ii
          goto 99
        endif
      enddo
!
99    continue
!
    endsubroutine find_index_range
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
    function spline_integral(z,f,q0)
!
!  Computes integral of a given function using spline interpolation.
!
!  01-apr-03/tobi: originally coded by Aake Nordlund
!
      implicit none
      real, dimension(:) :: z
      real, dimension(:) :: f
      real, dimension(size(z)) :: df,dz
      real, dimension(size(z)) :: q,spline_integral
      real, optional :: q0
      integer :: mz,k
!
      mz=size(z)
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
      spline_integral=q
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
    subroutine tridag(a,b,c,r,u,err)
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
      integer :: n,j
      real :: bet
!
      if (present(err)) err=.false.
      n=size(b)
      bet=b(1)
      if (bet==0.0) then
        print*,'tridag: Error at code stage 1'
        if (present(err)) err=.true.
        return
      endif
!
      u(1)=r(1)/bet
      do j=2,n
        gam(j)=c(j-1)/bet
        bet=b(j)-a(j)*gam(j)
        if (bet==0.0) then
          print*,'tridag: Error at code stage 2'
          if (present(err)) err=.true.
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
    subroutine pendag (m,a,b,c,d,e,r)
!
!  Solve pentadiagonal system of M linear equations.
!  A,B,C,D,E are the diagonals (A:subsub, B:sub, C:main, etc.).
!  R is the rhs on input and contains the solution on output.
!
!  01-apr-00/John Crowe (Newcastle): written
!
      implicit none
!
      integer :: i,m
      real :: x
      real, dimension (m) :: a,b,c,d,e,r
!
!  Eliminate sub-diagonals.
!
      x    = b(2)/c(1)
      c(2) = c(2)-d(1)*x
      d(2) = d(2)-e(1)*x
      r(2) = r(2)-r(1)*x
!
      do i=3,m-1,1
!
        x    = a(i)/c(i-2)
        b(i) = b(i)-d(i-2)*x
        c(i) = c(i)-e(i-2)*x
        r(i) = r(i)-r(i-2)*x
!
        x    = b(i)/c(i-1)
        c(i) = c(i)-d(i-1)*x
        d(i) = d(i)-e(i-1)*x
        r(i) = r(i)-r(i-1)*x
!
      enddo
!
      x    = a(m)/c(m-2)
      b(m) = b(m)-d(m-2)*x
      c(m) = c(m)-e(m-2)*x
      r(m) = r(m)-r(m-2)*x
!
      x    = b(m)/c(m-1)
      c(m) = c(m)-d(m-1)*x
      r(m) = r(m)-r(m-1)*x
!
!  Eliminate super-diagonals.
!
      r(m-1) = r(m-1) - d(m-1)*r(m)/c(m)
!
      do i=m-2,1 ,-1
        r(i)   = r(i) - d(i)*r(i+1)/c(i+1) - e(i)*r(i+2)/c(i+2)
      enddo
!
!  Reduce c's to unity, leaving the answers in r.
!
      do i=1,m,1
        r(i) = r(i) / c(i)
      enddo
!
    endsubroutine pendag
!***********************************************************************
    subroutine spline(arrx,arry,x2,S,psize1,psize2,err)
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
      real :: fac=0.1666666
      logical, intent(out), optional :: err
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
          print*,'spline x:y in x2:y2 : vector x is not monotonically increasing'
          if (present(err)) err=.true.
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
      call tridag(a,b,h,d,sol)
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
    function complex_phase(z)
!
!  Takes complex number and returns Theta where
!  z = A*exp(i*theta).
!
!  17-may-06/anders+jeff: coded
!
      use Cdata, only: pi
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
      if ((re>=0.0).and.(im< 0.0)) complex_phase=2*pi+asin(im/c)
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
    subroutine besselj_nu_int(res,nu,arg,loversample)
!
      use Cdata, only: pi,pi_1
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
      use Cdata, only : pi
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
      if (present(E_mu)) allocate(a_E(nnt))
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
!  X=0.75D0
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
      real, parameter :: BIGNO = 1.D10, BIGNI = 1.D-10
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
      DATA R1,R2,R3,R4,R5,R6 /57568490574.D0,-13362590354.D0, &
      651619640.7D0,-11214424.18D0,77392.33017D0,-184.9052456D0 /
      DATA S1,S2,S3,S4,S5,S6 /57568490411.D0,1029532985.D0, &
      9494680.718D0,59272.64853D0,267.8532712D0,1.D0 /
      IF(X==0.D0) GO TO 1
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
    1 BESSJ0 = 1.D0
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
      DATA R1,R2,R3,R4,R5,R6 /72362614232.D0,-7895059235.D0, &
      242396853.1D0,-2972611.439D0,15704.48260D0,-30.16036606D0 /
      DATA S1,S2,S3,S4,S5,S6 /144725228442.D0,2300535178.D0, &
      18583304.74D0,99447.43394D0,376.9991397D0,1.D0 /
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
!
      if (n<=2) stop "cyclic in the general module: n too small"
      gamma=-b(1)
      bb(1)=b(1)-gamma
      bb(n)=b(n)-alpha*beta/gamma
      do i=2,n-1
        bb(i)=b(i)
      enddo
      call tridag(a,bb,c,r,x)
      u(1)=gamma
      u(n)=alpha
      do i=2,n-1
        u(i)=0.
      enddo
      call tridag(a,bb,c,u,z)
      fact=(x(1)+beta*x(n)/gamma)/(1.+z(1)+beta*z(n)/gamma)
      do i=1,n
        x(i)=x(i)-fact*z(i)
      enddo
!
      return
!
    endsubroutine cyclic
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
!  module Messages
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
        print*, 'iproc = ', iproc
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
    integer function parser( zeile, feld, lenf, tz )
!
! Parses string zeile according to separator character tz;
! Puts up to lenf substrings consecutively into feld.
! Returns number of found substrings.
!
!  10-apr-11/MR: coded
!
      character (LEN=*),               intent(in)  :: zeile
      character (LEN=*), dimension(*), intent(out) :: feld
      integer,                         intent(in)  :: lenf
      character,                       intent(in)  :: tz
!
      integer :: ind, inda, lenz
!
      parser = 0
      inda = 1
      lenz = len_trim(zeile)
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
! range-wise output of a real or complex vector in ncol columns
! unfilled (inout) - number of unfilled slots in last written line
!
!  20-apr-11/MR: coded
!
    integer,                        intent(in)    :: unit
    real,    dimension(*),          intent(in)    :: buffer
    complex, dimension(*),          intent(in)    :: buffer_cmplx
    integer, dimension(3),          intent(in)    :: range
    integer,                        intent(inout) :: unfilled
    integer,              optional, intent(in)    :: ncol
    character(LEN=*),     optional, intent(in)    :: fmt
!
    integer          :: ncoll, nd, ia, ie, rest
    character(LEN=5) :: str
    character(LEN=20):: fmtl, fmth
    logical          :: lcomplex
!
    lcomplex = .false.
!
    if ( present(fmt) ) then
      fmtl = fmt
    else
      fmtl = 'e10.2'
    endif
!
    goto 1
!
  entry write_full_columns_cmplx(unit,buffer_cmplx,range,unfilled,ncol,fmt)
!
    lcomplex = .true.
!
    if ( present(fmt) ) then
      fmtl = '('//fmt//',1x,'//fmt//')'
    else
      fmtl = '(e10.2,1x,e10.2)'
    endif
!
 1  nd = get_range_no(range,1)
    if ( nd==0 ) return
!
    if ( present(ncol) ) then
      ncoll = ncol
    else
      ncoll = 8
    endif
!
    if ( unfilled > 0 ) then
!
      fmth = fmtl
      if ( nd>unfilled ) then
        ie = range(1)+(unfilled-1)*range(3)
        call chn(unfilled,str)
      else
        ie = range(2)
        if (nd<unfilled) fmth = trim(fmtl)//'$'
        call chn(nd,str)
      endif
!
      nd = nd-unfilled
!
      if (lcomplex) then
        write(1,'(1p,'//str//trim(fmth)//')') buffer_cmplx(range(1):ie:range(3))
      else
        write(1,'(1p,'//str//trim(fmth)//')') buffer(range(1):ie:range(3))
!        print*, 'str,buffer, nd=', str, buffer(range(1):ie:range(3)), nd
      endif
!
      if ( nd>0 ) then
        ia = ie + range(3)
      else
        unfilled = -nd
        return
      endif
    else
      ia = range(1)
    endif
!
    rest = mod(nd,ncoll)
    ie = (nd-rest-1)*range(3) + ia
!
    if ( rest<nd ) then
      call chn(ncoll,str)
      if (lcomplex) then
        write(1,'(1p,'//str//trim(fmtl)//')') buffer_cmplx(ia:ie:range(3))
      else
        write(1,'(1p,'//str//trim(fmtl)//')') buffer(ia:ie:range(3))
      endif
    endif
!
    if ( rest > 0 ) then
!
      call chn(rest,str)
      if (lcomplex) then
        write(1,'(1p,'//str//trim(fmtl)//'$)') buffer_cmplx(ie+range(3):range(2):range(3))
      else
        write(1,'(1p,'//str//trim(fmtl)//'$)') buffer(ie+range(3):range(2):range(3))
      endif
!
      unfilled = ncoll-rest
    else
      unfilled = 0
    endif
!
  endsubroutine write_full_columns_real
!***********************************************************************
    subroutine merge_ranges( ranges, ie, range, ia )
!
! merges ranges ia through ie in vector ranges with range
!
! 20-apr-11/MR: coded
!
      integer, dimension(3,*),          intent(inout) :: ranges
      integer,                          intent(in)    :: ie
      integer, dimension(3),            intent(inout) :: range
      integer,                optional, intent(in)    :: ia
      integer :: i, step, ial
!
      if ( present(ia) ) then
        ial = ia
      else
        ial = 1
      endif
!
      do i=ial,ie
!
        if ( ranges(1,i) > 0) then
!
          step = range(3)
          if ( ranges(3,i) == step ) then
!
            if ( range(1) <= ranges(2,i)+step .and. &
                 range(2) >= ranges(1,i)-step .and. &
                 mod(range(1)-ranges(1,i),step) == 0 ) then
!
              ranges(1,i) = min( ranges(1,i), range(1) )
              ranges(2,i) = max( ranges(2,i), range(2) )
              range = 0
!
            endif
!
          else
!           no implementation yet
          endif
        endif
!
      enddo
!
    endsubroutine merge_ranges
!***********************************************************************
    logical function read_range( crange, defrange, range )
!
! reads a range (start,stop,step) from string crange of the shape [0...9][:[0...9][:[0...9]]]
! adjusts range not to exceed limits of default range defrange
!
! 20-apr-11/MR: coded
!
      character (LEN=20)   , intent(in)  :: crange
      integer, dimension(3), intent(in)  :: defrange
      integer, dimension(3), intent(out) :: range
!
      integer :: isep, isep1, ios, ios1, ios2, lenrng, tmp
!
      ios=0; ios1=0; ios2=0
!
      if ( crange /= '' ) then
!
        range = defrange
!
        isep = index(crange,':')
        lenrng = len_trim(crange)
!
        if ( isep > 0 ) then
!
          if ( isep > 1 ) then
            read( crange(1:isep-1),*,IOSTAT=ios ) range(1)
            if ( ios == 0 ) &
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
              if ( ios1 == 0 ) &
                range(2) = min(max(defrange(1),range(2)),defrange(2))
            endif
!
            if ( isep1 < lenrng ) then
!
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
          if ( ios == 0 ) &
            range(1) = min(max(defrange(1),range(1)),defrange(2))
          range(2) = range(1)
!
        endif
!
        if ( ios/=0 .or. ios1/=0 .or. ios2/=0 ) &
          print*, 'read_range - Warning: invalid data in range!'
!
        read_range = .true.
      else
        read_range = .false.
      endif
!
    endfunction read_range
!***********************************************************************
    integer function get_range_no( ranges, nr )
!
! determines total number of elements selected by all ranges in ranges
!
! 20-apr-11/MR: coded
!
      integer, dimension(3,*), intent(in) :: ranges
      integer,       optional, intent(in) :: nr
!
      integer :: i, nrl
!
      get_range_no = 0
!
      if ( present(nr) ) then
        nrl = nr
      else
        nrl = 10
      endif
!
      do i=1,nrl
        if ( ranges(1,i) > 0 ) &
          get_range_no = get_range_no + &
                      ceiling( (ranges(2,i)-ranges(1,i)+1.)/ranges(3,i))
      enddo
!
    endfunction get_range_no
!******************************************************************************
  subroutine write_by_ranges_2d_real( unit, buffer, xranges, yranges, trans )
!
! writes a real or complex 2D array controlled by lists of ranges
! xranges and yranges for each of the two dimensions; output optionally transposed
!
! 10-may-11/MR: coded
!
    integer,                           intent(in)           :: unit
    integer, dimension(3,*)          , intent(in)           :: xranges, yranges
    real,    dimension(nxgrid,nygrid), intent(in)           :: buffer
    complex, dimension(nxgrid,nygrid), intent(in)           :: buffer_cmplx
    logical,                           intent(in), optional :: trans
!
    integer :: i,j,il,jl,unfilled
    logical :: transl, lcomplex
!
    lcomplex = .false.
    goto 1
!
  entry write_by_ranges_2d_cmplx( unit, buffer_cmplx, xranges, yranges, trans )
!
    lcomplex = .true.
!
 1  unfilled = 0
!
    if ( present(trans) ) then
      transl = trans
    else
      transl = .false.
    endif
!
    do j=1,10
      if ( yranges(1,j) > 0 ) then
        do jl=yranges(1,j),yranges(2,j),yranges(3,j)
          do i=1,10
            if ( xranges(1,i) > 0 ) then
              if ( transl ) then
                if (lcomplex) then
                  call write_full_columns_cmplx( unit, buffer_cmplx(jl,:), xranges(1,i), unfilled )
                else
                  call write_full_columns_real( unit, buffer(jl,:), xranges(1,i), unfilled )
                endif
              else
                if (lcomplex) then
                  call write_full_columns_cmplx( unit, buffer_cmplx(1,jl), xranges(1,i), unfilled )
                else
                  call write_full_columns_real( unit, buffer(1,jl), xranges(1,i), unfilled )
                endif
              endif
            endif
          enddo
        enddo
      endif
    enddo
!
    if ( unfilled > 0 ) write( unit,'(a)')
!
  endsubroutine write_by_ranges_2d_real
!***********************************************************************
  subroutine write_by_ranges_1d_real(unit,buffer,ranges)
!
! writes a real or complex vector controlled by lists of ranges
! output optionally transposed
!
! 10-may-11/MR: coded
!
    integer,                 intent(in) :: unit
    real,    dimension(*)  , intent(in) :: buffer
    complex, dimension(*)  , intent(in) :: buffer_cmplx
    integer, dimension(3,*), intent(in) :: ranges
!
    integer :: unfilled, i
    logical :: lcomplex
!
    lcomplex = .false.
    goto 1
!
  entry  write_by_ranges_1d_cmplx(unit,buffer_cmplx,ranges)
    lcomplex = .true.
!
 1  unfilled = 0
    do i=1,10
      if ( ranges(1,i) > 0 ) then
        if (lcomplex) then
          call write_full_columns_cmplx( unit, buffer_cmplx, ranges(1,i), unfilled )
        else
          call write_full_columns_real( unit, buffer, ranges(1,i), unfilled )
        endif
      endif
    enddo
!
    if ( unfilled > 0 ) write(unit,'(a)')
!
  endsubroutine write_by_ranges_1d_real
!***********************************************************************
endmodule General
