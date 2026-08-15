module Quiet

  use Cparam

  public :: keep_compiler_quiet
  public :: keep_compiler_quiet_dble
!
  interface keep_compiler_quiet ! Overload `keep_compiler_quiet' function
    module procedure keep_compiler_quiet_r
    module procedure keep_compiler_quiet_r1d
    module procedure keep_compiler_quiet_r2d
    module procedure keep_compiler_quiet_r3d
    module procedure keep_compiler_quiet_r4d
    module procedure keep_compiler_quiet_r5d
    module procedure keep_compiler_quiet_p
    module procedure keep_compiler_quiet_bc
    module procedure keep_compiler_quiet_sl
    module procedure keep_compiler_quiet_i
    module procedure keep_compiler_quiet_i1
    module procedure keep_compiler_quiet_i1d
    module procedure keep_compiler_quiet_i1d1
    module procedure keep_compiler_quiet_i81d
    module procedure keep_compiler_quiet_i2d
    module procedure keep_compiler_quiet_i3d
    module procedure keep_compiler_quiet_l
    module procedure keep_compiler_quiet_l1d
    module procedure keep_compiler_quiet_c
    module procedure keep_compiler_quiet_c1d
  endinterface

contains
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
    subroutine keep_compiler_quiet_dble(v1,v2,v3,v4)
!
      real(KIND=rkind8)     :: v1,v2,v3,v4
      optional ::                 v2,v3,v4
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
    endsubroutine keep_compiler_quiet_dble
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
    subroutine keep_compiler_quiet_r5d(v1,v2,v3,v4)
!
!  Call this to avoid compiler warnings about unused variables.
!  Optional arguments allow for more variables of the same shape+type.
!
!  4-apr-25/TP: coded
!
      real, dimension(:,:,:,:,:) :: v1, v2, v3, v4
      optional                 ::     v2, v3, v4
!
      if (ALWAYS_FALSE) then
        write(0,*) 'keep_compiler_quiet_r5d: never got here...'
        print*,                  minval(v1)
        if (present(v2)) print*, minval(v2)
        if (present(v3)) print*, minval(v3)
        if (present(v4)) print*, minval(v4)
        STOP 1
      endif
!
    endsubroutine keep_compiler_quiet_r5d
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
    subroutine keep_compiler_quiet_i1(v1,v2,v3,v4)
!
!  Call this to avoid compiler warnings about unused variables.
!  Optional arguments allow for more variables of the same shape+type.
!
!  04-aug-06/wolf: coded
!
      integer(KIND=ikind1)  :: v1, v2, v3, v4
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
    endsubroutine keep_compiler_quiet_i1
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
    subroutine keep_compiler_quiet_i1d1(v1,v2,v3,v4)
!
!  Call this to avoid compiler warnings about unused variables.
!  Optional arguments allow for more variables of the same shape+type.
!
!  04-aug-06/wolf: coded
!
      integer(KIND=ikind1), dimension(:)  :: v1, v2, v3, v4
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
    endsubroutine keep_compiler_quiet_i1d1
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
    subroutine keep_compiler_quiet_c1d(v1, v2, v3, v4)
!
!  Call this to avoid compiler warnings about unused variables.
!  Optional arguments allow for more variables of the same shape+type.
!
!  10-nov-20/ccyang: coded
!
      character(len=*), dimension(:), intent(in) :: v1
      character(len=*), dimension(:), intent(in), optional :: v2, v3, v4
!
      if (ALWAYS_FALSE) then
        write(0,*) 'keep_compiler_quiet_c1d: Never got here...'
        print *, v1
        if (present(v2)) print *, v2
        if (present(v3)) print *, v3
        if (present(v4)) print *, v4
        STOP 1
      endif
!
    endsubroutine keep_compiler_quiet_c1d
!***********************************************************************
endmodule Quiet
