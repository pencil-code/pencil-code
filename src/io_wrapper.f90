module IO_wrapper
!
  use Cparam, only: mx,my,mz,mparray,ikind8,max_int
  use Cdata, only: luse_alt_io
!
  implicit none

  private
  public :: initialize_io_wrapper
!
  external caller1, caller2_str1, caller3_str1, caller4_str1, caller5_str12, caller7_str67, caller0, dlclose_c
  integer(KIND=ikind8), external :: dlopen_c, dlsym_c
  character(LEN=128), external :: dlerror_c

  interface read_persist_wrap
    module procedure read_persist_logical_0D_wrap
    module procedure read_persist_logical_1D_wrap
    module procedure read_persist_int_0D_wrap
    module procedure read_persist_int_1D_wrap
    module procedure read_persist_real_0D_wrap
    module procedure read_persist_real_1D_wrap
  endinterface

  integer, parameter :: N_ROUTINES=13, &
                        I_INPUT_SNAP=1, I_INPUT_SNAP_FINALIZE=2, &
                        I_INPUT_PART_SNAP=3, I_INPUT_POINTMASS=4, &
                        I_INPUT_GLOBALS=5, I_READ_PERSIST_ID=6, &
                        I_READ_PERSIST_LOGICAL_0D=7, &
                        I_READ_PERSIST_LOGICAL_1D=8, &
                        I_READ_PERSIST_INT_0D=9, &
                        I_READ_PERSIST_INT_1D=10, &
                        I_READ_PERSIST_REAL_0D=11, &
                        I_READ_PERSIST_REAL_1D=12, &
                        I_RGRID=13 

  character(LEN=29), dimension(N_ROUTINES) :: routines=(/'input_snap                   ', &
                                                         'input_snap_finalize          ', &
                                                         'input_part_snap              ', &
                                                         'input_pointmass              ', &
                                                         'input_globals                ', &
                                                         'read_persist_id              ', &
                                                         'read_persist_logical_0d      ', &
                                                         'read_persist_logical_1d      ', &
                                                         'read_persist_int_0d          ', &
                                                         'read_persist_int_1d          ', &
                                                         'read_persist_real_0d         ', &
                                                         'read_persist_real_1d         ', &
                                                         'rgrid                        ' /)
  integer(KIND=ikind8) :: libhandle
  integer(KIND=ikind8), dimension(n_routines) :: sub_handles

  contains
!****************************************************************************
  subroutine initialize_io_wrapper

    use Messages, only: warning,fatal_error

    integer, parameter :: unit=11
    integer, parameter :: RTLD_LAZY=0, RTLD_NOW=0  !!!

    integer :: j
    integer(KIND=ikind8) :: sub_handle
    character(LEN=8) :: mod_prefix, mod_infix, mod_suffix
    character(LEN=128) :: error

    if (luse_alt_io) then

      libhandle=dlopen_c('src/io_alt.so'//char(0),RTLD_NOW)

      if (libhandle==0) then
        call warning('initialize_io_wrapper','Could not open src/io_alt.so. Continue w/o alternative I/O')
        luse_alt_io=.false.
        return
      endif

      call getenv("MODULE_PREFIX", mod_prefix)
      call getenv("MODULE_INFIX", mod_infix)
      call getenv("MODULE_SUFFIX", mod_suffix)

      do j=1,N_ROUTINES
        sub_handle=dlsym_c(libhandle,trim(mod_prefix)//'io'// &
            trim(mod_infix)//trim(routines(j))//trim(mod_suffix)//char(0))
        if (sub_handle==0) &
          call fatal_error('initialize_io_wrapper','Error for symbol '// &
                           trim(routines(j))//' in module IO: '//trim(dlerror_c()))
        sub_handles(j) = sub_handle
      enddo

    endif

  endsubroutine initialize_io_wrapper
!***********************************************************************
  subroutine finalize_io_wrapper

    if (luse_alt_io) call dlclose_c(libhandle)

  endsubroutine finalize_io_wrapper
!***********************************************************************
    subroutine input_snap_wrap(file,a,nv,mode)
!
!  manages reading of snapshot from different precision
!
!  24-oct-13/MR: coded
!
      use IO, only: input_snap

      character (len=*), intent(in) :: file
      integer, intent(in) :: nv, mode
      real, dimension (mx,my,mz,nv), intent(out) :: a

      if (luse_alt_io) then
        call caller4_str1(sub_handles(I_INPUT_SNAP),file,a,nv,mode)
      else
        call input_snap(file,a,nv,mode)
      endif

    endsubroutine input_snap_wrap
!***********************************************************************
    subroutine input_snap_finalize_wrap
!
!  Close snapshot file.
!
      use IO, only: input_snap_finalize

      if (luse_alt_io) then
        call caller0(sub_handles(I_INPUT_SNAP_FINALIZE))
      else
        call input_snap_finalize
      endif

    endsubroutine input_snap_finalize_wrap
!***********************************************************************
    subroutine input_part_snap_wrap(ipar, ap, mv, nv, npar_total, file, label)
!
!  Read particle snapshot file, mesh and time are read in 'input_snap'.
!
      use IO, only: input_part_snap

      integer, intent(in) :: mv
      integer, dimension (mv), intent(out) :: ipar
      real, dimension (mv,mparray), intent(out) :: ap
      integer, intent(out) :: nv, npar_total
      character (len=*), intent(in) :: file
      character (len=*), optional, intent(in) :: label

      if (luse_alt_io) then
        call caller7_str67(sub_handles(I_INPUT_PART_SNAP),ipar,ap,mv,nv,npar_total,file,label)
      else
        call input_part_snap(ipar, ap, mv, nv, npar_total, file, label)
      endif

    endsubroutine input_part_snap_wrap
!***********************************************************************
    subroutine input_pointmass_wrap(file, labels, fq, mv, nc)
!
!  Read pointmass snapshot file.
!
      use IO, only: input_pointmass

      character (len=*), intent(in) :: file
      integer, intent(in) :: mv, nc
      character (len=*), dimension (nc), intent(in) :: labels
      real, dimension (mv,nc), intent(out) :: fq

      if (luse_alt_io) then
        !!!call caller5_str12(sub_handles(I_INPUT_POINTMASS),file,labels,fq,mv,nc)
      else
        call input_pointmass(file,labels,fq,mv,nc)
      endif

    endsubroutine input_pointmass_wrap
!***********************************************************************
    subroutine input_globals_wrap(file, a, nv)
!
!  Read globals snapshot file, ignoring mesh.
!
      use IO, only: input_globals

      character (len=*) :: file
      integer :: nv
      real, dimension (mx,my,mz,nv) :: a

      if (luse_alt_io) then
        call caller3_str1(sub_handles(I_INPUT_GLOBALS),file,a,nv)
      else
        call input_globals(file,a,nv)
      endif

    endsubroutine input_globals_wrap    
!***********************************************************************
    logical function read_persist_id_wrap(label, id, lerror_prone)
!
!  Read persistent block ID from snapshot file.
!
      use IO, only: read_persist_id

      character (len=*), intent(in) :: label
      integer, intent(out) :: id
      logical, intent(in), optional :: lerror_prone

      if (luse_alt_io) then
        call caller3_str1(sub_handles(I_READ_PERSIST_ID),label,id,lerror_prone)
      else
        read_persist_id_wrap=read_persist_id(label,id,lerror_prone)
      endif
!
    endfunction read_persist_id_wrap
!***********************************************************************
    logical function read_persist_logical_0D_wrap(label, value)
!
!  Read persistent data from snapshot file.
!
      use IO, only: read_persist

      character (len=*), intent(in) :: label
      logical, intent(out) :: value
!
      if (luse_alt_io) then
        call caller2_str1(sub_handles(I_READ_PERSIST_LOGICAL_0D),label,value)
      else
        read_persist_logical_0D_wrap=read_persist(label,value)
      endif
!
    endfunction read_persist_logical_0D_wrap
!***********************************************************************
    logical function read_persist_logical_1D_wrap(label, value)
!
!  Read persistent data from snapshot file.
!
      use IO, only: read_persist

      character (len=*), intent(in) :: label
      logical, dimension(:), intent(out) :: value
!
      if (luse_alt_io) then
        call caller2_str1(sub_handles(I_READ_PERSIST_LOGICAL_1D),label,value)
      else
        read_persist_logical_1D_wrap=read_persist(label,value)
      endif
!
    endfunction read_persist_logical_1D_wrap
!***********************************************************************
    logical function read_persist_int_0D_wrap(label, value)
!
!  Read persistent data from snapshot file.
!
      use IO, only: read_persist

      character (len=*), intent(in) :: label
      integer, intent(out) :: value
!
      if (luse_alt_io) then
        call caller2_str1(sub_handles(I_READ_PERSIST_INT_0D),label,value)
      else
        read_persist_int_0D_wrap=read_persist(label,value)
      endif
!
    endfunction read_persist_int_0D_wrap
!***********************************************************************
    logical function read_persist_int_1D_wrap(label, value)
!
!  Read persistent data from snapshot file.
!
      use IO, only: read_persist

      character (len=*), intent(in) :: label
      integer, dimension(:), intent(out) :: value
!
      if (luse_alt_io) then
        call caller2_str1(sub_handles(I_READ_PERSIST_INT_1D),label,value)
      else
        read_persist_int_1D_wrap=read_persist(label,value)
      endif
!
    endfunction read_persist_int_1D_wrap
!***********************************************************************
    logical function read_persist_real_0D_wrap(label, value)
!
!  Read persistent data from snapshot file.
!
!  13-Dec-2011/Bourdin.KIS: coded
!  23-oct-2013/MR: modified for reading of different precision
!
      use IO, only: read_persist

      character (len=*), intent(in) :: label
      real, intent(out) :: value
!
      if (luse_alt_io) then
        call caller2_str1(sub_handles(I_READ_PERSIST_REAL_0D),label,value)
      else
        read_persist_real_0D_wrap=read_persist(label,value)
      endif
!
    endfunction read_persist_real_0D_wrap
!***********************************************************************
    logical function read_persist_real_1D_wrap(label, value)
!
!  Read persistent data from snapshot file.
!
      use IO, only: read_persist

      character (len=*), intent(in) :: label
      real, dimension(:), intent(out) :: value
!
      if (luse_alt_io) then
        call caller2_str1(sub_handles(I_READ_PERSIST_REAL_1D),label,value)
      else
        read_persist_real_1D_wrap=read_persist(label,value)
      endif
!
    endfunction read_persist_real_1D_wrap
!***********************************************************************
    subroutine rgrid_wrap(file)
!
!  Read processor-local part of grid coordinates.
!
      use IO, only: rgrid

      character (len=*) :: file
!
      if (luse_alt_io) then
        call caller1(sub_handles(I_RGRID),file)
      else
        call rgrid(file)
      endif

    endsubroutine rgrid_wrap
!***********************************************************************
end module 
