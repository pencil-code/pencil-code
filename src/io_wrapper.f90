!
module Io
!
  use Cparam, only: mx,my,mz,mparray,labellen,fnlen
  use Cdata, only: lstart, lroot
!
  implicit none

  include 'io.h'
!
! Indicates if IO is done distributed (each proc writes into a procdir)
! or collectively (eg. by specialized IO-nodes or by MPI-IO).
!
  logical :: lcollective_IO
  character (len=labellen) :: IO_strategy

  logical :: lswitched_to_out=.false., lwr_grid=.false.

  contains
!****************************************************************************
  subroutine register_io

    use Io_in, only: register_io_ => register_io, &
                     IO_strategy_ => IO_strategy, lcollective_IO_ => lcollective_IO
    use Io_out, only: IO_strategy_out => IO_strategy

    use Messages, only: warning
    use File_io, only: file_exists
    use Syscalls, only: system_cmd

    if (lstart.or.file_exists('IO_LOCK')) then
      if (lstart) then
        call warning('register_io','Alternative IO strategy for input will be ignored in start.x'// &
                     ' AND in subsequent executions of run.x')
        if (lroot) call system_cmd("touch IO_LOCK")
      else
        call warning('register_io','IO_LOCK file found - alternative IO strategy for input will be ignored.')
      endif
      call switch_io
    else 
      IO_strategy = IO_strategy_
      lcollective_IO = lcollective_IO_
      call register_io_
      lwr_grid = .true.
    endif

  endsubroutine register_io
!***********************************************************************
  subroutine switch_io

     use Io_out, only: register_io_ => register_io, &
                       IO_strategy_ => IO_strategy, &
                       lcollective_IO_ => lcollective_IO, &
                       directory_names_ => directory_names
     use Syscalls, only: system_cmd

    if (.not.lswitched_to_out) then

      lswitched_to_out=.true.
      IO_strategy = IO_strategy_
      lcollective_IO = lcollective_IO_

      call directory_names_ 
      call register_io_

     if (lroot) call system_cmd("touch IO_LOCK")

    endif

  endsubroutine switch_io
!***********************************************************************
  subroutine finalize_io
!
     use Io_in, only: finalize_io_in => finalize_io
     use Io_out, only: finalize_io_out => finalize_io

     call finalize_io_in
     call finalize_io_out
     !call system_cmd("sed -i -e's/^\( *IO_IN\)/#\1/' src/Makefile.local")

  endsubroutine finalize_io
!***********************************************************************
    subroutine log_filename_to_file(filename, flist)
!
     use Io_in, only: log_filename_to_file_in => log_filename_to_file
     use Io_out, only: log_filename_to_file_out => log_filename_to_file
!
      character (len=*) :: filename, flist

      if (lswitched_to_out) then
        call log_filename_to_file_out(filename, flist)
      else
        call log_filename_to_file_in(filename, flist)
      endif

    endsubroutine log_filename_to_file
!***********************************************************************
    subroutine directory_names
!
     use Io_in, only: directory_names_in => directory_names
     use Io_out, only: directory_names_out => directory_names

      if (lswitched_to_out) then
        call directory_names_out
      else
        call directory_names_in
      endif

    endsubroutine directory_names
!***********************************************************************
    subroutine input_snap(file,a,nv,mode)
!
      use Io_in, only: input_snap_in => input_snap
      use Io_out, only: input_snap_out => input_snap
!
      character (len=*), intent(in) :: file
      integer, intent(in) :: nv, mode
      real, dimension (mx,my,mz,nv), intent(out) :: a

      if (lswitched_to_out) then
        call input_snap_out(file,a,nv,mode)
      else
        call input_snap_in(file,a,nv,mode)
      endif

    endsubroutine input_snap
!***********************************************************************
    subroutine input_snap_finalize
!
      use Io_in, only: input_snap_finalize_in => input_snap_finalize
      use Io_out, only: input_snap_finalize_out => input_snap_finalize

      if (lswitched_to_out) then
        call input_snap_finalize_out
      else
        call input_snap_finalize_in
      endif

    endsubroutine input_snap_finalize
!***********************************************************************
    subroutine input_slice_real_arr(file, time, pos, data)
!
!  Read a slice file.
!
!  13-nov-20/ccyang: wrapper
!
      use HDF5_IO, only: hdf5_input_slice
!
      character(len=*), intent(in) :: file
      real, intent(out):: time, pos
      real, dimension(:,:,:), intent(out):: data
!
      call hdf5_input_slice(file, time, pos, data)
!
    endsubroutine input_slice_real_arr
!***********************************************************************
    subroutine input_slice_scat_arr(file, pos, data, ivar, nt)
!
!  Read a slice file.
!
!  13-nov-20/ccyang: wrapper
!
      use General, only: scattered_array
      use HDF5_IO, only: hdf5_input_slice
!
      character(len=*), intent(in) :: file
      real, intent(out):: pos
      type(scattered_array), pointer :: data   !intent(inout)
      integer, intent(in) :: ivar
      integer, intent(out):: nt
!
      call hdf5_input_slice(file, pos, data, ivar, nt)
!
    endsubroutine input_slice_scat_arr
!***********************************************************************
    subroutine input_part_snap(ipar, ap, mv, nv, npar_total, file, label)
!
      use Io_in, only: input_part_snap_in => input_part_snap
      use Io_out, only: input_part_snap_out => input_part_snap

      integer, intent(in) :: mv
      integer, dimension (mv), intent(out) :: ipar
      real, dimension (mv,mparray), intent(out) :: ap
      integer, intent(out) :: nv, npar_total
      character (len=*), intent(in) :: file
      character (len=*), optional, intent(in) :: label

      if (lswitched_to_out) then
        call input_part_snap_out(ipar, ap, mv, nv, npar_total, file, label)
      else
        call input_part_snap_in(ipar, ap, mv, nv, npar_total, file, label)
      endif

    endsubroutine input_part_snap
!***********************************************************************
    subroutine input_pointmass(file, labels, fq, mv, nc)
!
      use Io_in, only: input_pointmass_in => input_pointmass
      use Io_out, only: input_pointmass_out => input_pointmass

      character (len=*), intent(in) :: file
      integer, intent(in) :: mv, nc
      character (len=*), dimension (nc), intent(in) :: labels
      real, dimension (mv,nc), intent(out) :: fq

      if (lswitched_to_out) then
        call input_pointmass_out(file, labels, fq, mv, nc)
      else
        call input_pointmass_in(file, labels, fq, mv, nc)
      endif

    endsubroutine input_pointmass
!***********************************************************************
    subroutine input_globals(file, a, nv)
!
      use Io_in, only: input_globals_in => input_globals
      use Io_out, only: input_globals_out => input_globals

      character (len=*) :: file
      integer :: nv
      real, dimension (mx,my,mz,nv) :: a

      if (lswitched_to_out) then
        call input_globals_out(file, a, nv)
      else
        call input_globals_in(file, a, nv)
      endif

    endsubroutine input_globals
!***********************************************************************
    logical function init_read_persist(file)
!
      use Io_in, only: init_read_persist_in => init_read_persist
      use Io_out, only: init_read_persist_out => init_read_persist
!
      character (len=*), intent(in), optional :: file
!
      if (lswitched_to_out) then
        init_read_persist=init_read_persist_out(file)
      else
        init_read_persist=init_read_persist_in(file)
      endif

    endfunction init_read_persist
!***********************************************************************
    logical function read_persist_id(label, id, lerror_prone) 

      use Io_in, only: read_persist_id_in => read_persist_id
      use Io_out, only: read_persist_id_out => read_persist_id 
!
      character (len=*), intent(in) :: label
      integer, intent(out) :: id
      logical, intent(in), optional :: lerror_prone

      if (lswitched_to_out) then
        read_persist_id = read_persist_id_out(label, id, lerror_prone)
      else
        read_persist_id = read_persist_id_in(label, id, lerror_prone)
      endif
!
    endfunction read_persist_id
!***********************************************************************
    subroutine rgrid(file)

      use Io_in, only: rgrid_in => rgrid
      use Io_out, only: rgrid_out => rgrid
!
      character (len=*) :: file
!
      if (lswitched_to_out) then
        call rgrid_out(file)
      else
        call rgrid_in(file)
      endif

    endsubroutine rgrid
!***********************************************************************
    subroutine rproc_bounds(file)

      use Io_in, only: rproc_bounds_in => rproc_bounds
      use Io_out, only: rproc_bounds_out => rproc_bounds
!
      character (len=*) :: file
!
     if (lswitched_to_out) then
        call rproc_bounds_out(file)
      else
        call rproc_bounds_in(file)
      endif
!
    endsubroutine rproc_bounds
!***********************************************************************
    subroutine wgrid(file,mxout,myout,mzout,lwrite)

      use Io_out, only: wgrid_ => wgrid
      use General, only: loptest
!
      character (len=*) :: file
      integer, optional :: mxout,myout,mzout
      logical, optional :: lwrite
!
      call switch_io
      call wgrid_(file,mxout,myout,mzout,loptest(lwrite,lwr_grid))

    endsubroutine wgrid
!***********************************************************************
    subroutine wproc_bounds(file)

      use Io_out, only: wproc_bounds_ => wproc_bounds
!
      character (len=*) :: file
!
      call switch_io
      call wproc_bounds_(file)

    endsubroutine wproc_bounds
!***********************************************************************
    subroutine output_snap(a,nv1,nv2,file) 

      use Io_out, only: output_snap_ => output_snap 
!
      real, dimension (:,:,:,:),  intent(IN) :: a
      integer,           optional,intent(IN) :: nv1,nv2
      character (len=*), optional,intent(IN) :: file
!
      call switch_io
      call output_snap_(a,nv1,nv2,file)
!
    endsubroutine output_snap
!***********************************************************************
    subroutine output_snap_finalize

      use Io_out, only: output_snap_finalize_ => output_snap_finalize

      call switch_io
      call output_snap_finalize_ 
!
    endsubroutine output_snap_finalize
!***********************************************************************
    subroutine output_average_2D(label, nc, name, data, time, lwrite, header)
!
!  Output 2D average to a file.
!
!  01-dec-2022/ccyang: stub
!
      use General, only: keep_compiler_quiet
!
      integer, intent(in) :: nc
      character(len=fmtlen), dimension(nc), intent(in) :: name
      character(len=*), intent(in) :: label
      real, dimension(:,:,:), intent(in) :: data
      real, intent(in) :: time
      logical, intent(in) :: lwrite
      real, dimension(:), optional, intent(in) :: header
!
      call keep_compiler_quiet(label)
      call keep_compiler_quiet(nc)
      call keep_compiler_quiet(name)
      call keep_compiler_quiet(data)
      call keep_compiler_quiet(time)
      call keep_compiler_quiet(lwrite)
      if (present(header)) call keep_compiler_quiet(header)
!
      call fatal_error("output_average_2D", "not implemented")
!
    endsubroutine output_average_2D
!***********************************************************************
    subroutine output_slice_position()
!
!  Record slice positions.
!
!  13-nov-20/ccyang: wrapper
!
      use HDF5_IO, only: hdf5_output_slice_position
!
      call hdf5_output_slice_position()
!
    endsubroutine output_slice_position
!***********************************************************************
    subroutine output_slice(lwrite, time, label, suffix, pos, grid_pos, data)
!
!  Append to a slice file
!
!  13-nov-20/ccyang: wrapper
!
      use HDF5_IO, only: hdf5_output_slice
!
      logical, intent(in) :: lwrite
      real, intent(in) :: time
      character(len=*), intent(in) :: label, suffix
      real, intent(in) :: pos
      integer, intent(in) :: grid_pos
      real, dimension(:,:), pointer :: data
!
      call hdf5_output_slice(lwrite, time, label, suffix, pos, grid_pos, data)
!
    endsubroutine output_slice
!***********************************************************************
    subroutine output_part_snap(ipar, ap, mv, nv, file, label, ltruncate) 

      use Io_out, only: output_part_snap_ => output_part_snap 
!
      integer, intent(in) :: mv, nv
      integer, dimension (mv), intent(in) :: ipar
      real, dimension (mv,mparray), intent(in) :: ap
      character (len=*), intent(in) :: file
      character (len=*), optional, intent(in) :: label
      logical, optional, intent(in) :: ltruncate
!
      call switch_io
      call output_part_snap_(ipar, ap, mv, nv, file, label, ltruncate)
!
    endsubroutine output_part_snap
!***********************************************************************
    subroutine output_stalker_init(num, nv, snap, ID) 

      use Io_out, only: output_stalker_init_ => output_stalker_init 
!
      integer, intent(in) :: num, nv, snap
      integer, dimension(nv), intent(in) :: ID
!
      call switch_io
      call output_stalker_init_(num, nv, snap, ID)
!
    endsubroutine output_stalker_init
!***********************************************************************
    subroutine output_stalker(label, mv, nv, data, nvar, lfinalize) 

      use Io_out, only: output_stalker_ => output_stalker 
!
      character (len=*), intent(in) :: label
      integer, intent(in) :: mv, nv
      real, dimension (mv), intent(in) :: data
      integer, intent(in), optional :: nvar
      logical, intent(in), optional :: lfinalize
!
      call switch_io
      call output_stalker_(label, mv, nv, data, nvar, lfinalize)
!
    endsubroutine output_stalker
!***********************************************************************
    subroutine output_part_finalize

      use Io_out, only: output_part_finalize_ => output_part_finalize

      call switch_io
      call output_part_finalize_ 
!
    endsubroutine output_part_finalize
!***********************************************************************
    subroutine output_pointmass(file, labels, fq, mv, nc) 

      use Io_out, only: output_pointmass_ => output_pointmass 
!
      character (len=*), intent(in) :: file
      integer, intent(in) :: mv, nc
      character (len=*), dimension (nc), intent(in) :: labels
      real, dimension (mv,nc), intent(in) :: fq
!
      call switch_io
      call output_pointmass_(file, labels, fq, mv, nc)
!
    endsubroutine output_pointmass
!***********************************************************************
    subroutine output_globals(file, a, nv, label) 

      use Io_out, only: output_globals_ => output_globals 

      integer :: nv
      real, dimension (mx,my,mz,nv) :: a
      character (len=*) :: file
      character (len=*), intent(in), optional :: label
!
      call switch_io
      call output_globals_(file, a, nv, label)

    endsubroutine output_globals
!***********************************************************************
    logical function persist_exists(label)

      use Io_in, only: persist_exists_ => persist_exists 
!
      character (len=*), intent(in) :: label
!
      persist_exists=persist_exists_(label)
!
    endfunction persist_exists
!***********************************************************************
    logical function init_write_persist(file)

      use Io_out, only: init_write_persist_ => init_write_persist
!
      character (len=*), intent(in), optional :: file
!
      call switch_io
      init_write_persist=init_write_persist_(file)

    endfunction init_write_persist
!***********************************************************************
    logical function write_persist_id(label, id) 

      use Io_out, only: write_persist_id_ => write_persist_id 
!
      character (len=*), intent(in) :: label
      integer, intent(in) :: id
!
      call switch_io
      write_persist_id = write_persist_id_(label, id)
!
    endfunction write_persist_id
!***********************************************************************
    logical function read_persist_logical_0D(label, value)

      use Io_in, only: read_persist_logical_0D_in => read_persist_logical_0D
      use Io_out, only: read_persist_logical_0D_out => read_persist_logical_0D
!
      character (len=*), intent(in) :: label
      logical, intent(out) :: value
!
      if (lswitched_to_out) then
        read_persist_logical_0D = read_persist_logical_0D_out(label, value)
      else  
        read_persist_logical_0D = read_persist_logical_0D_in(label, value)
      endif
!
    endfunction read_persist_logical_0D
!***********************************************************************
    logical function read_persist_logical_1D(label, value)

      use Io_in, only: read_persist_logical_1D_in => read_persist_logical_1D
      use Io_out, only: read_persist_logical_1D_out => read_persist_logical_1D
!
      character (len=*), intent(in) :: label
      logical, dimension(:), intent(out) :: value
!
      if (lswitched_to_out) then
        read_persist_logical_1D = read_persist_logical_1D_out(label, value)
      else  
        read_persist_logical_1D = read_persist_logical_1D_in(label, value)
      endif
!
    endfunction read_persist_logical_1D
!***********************************************************************
    logical function read_persist_int_0D(label, value)

      use Io_in, only: read_persist_int_0D_in => read_persist_int_0D
      use Io_out, only: read_persist_int_0D_out => read_persist_int_0D
!
      character (len=*), intent(in) :: label
      integer, intent(out) :: value
!
      if (lswitched_to_out) then
        read_persist_int_0D = read_persist_int_0D_out(label, value)
      else
        read_persist_int_0D = read_persist_int_0D_in(label, value)
      endif
!
    endfunction read_persist_int_0D
!***********************************************************************
    logical function read_persist_int_1D(label, value)

      use Io_in, only: read_persist_int_1D_in => read_persist_int_1D
      use Io_out, only: read_persist_int_1D_out => read_persist_int_1D
!
      character (len=*), intent(in) :: label
      integer, dimension(:), intent(out) :: value

      if (lswitched_to_out) then
        read_persist_int_1D = read_persist_int_1D_out(label, value)
      else
        read_persist_int_1D = read_persist_int_1D_in(label, value)
      endif
!
    endfunction read_persist_int_1D
!***********************************************************************
    logical function read_persist_real_0D(label, value)

      use Io_in, only: read_persist_real_0D_in => read_persist_real_0D
      use Io_out, only: read_persist_real_0D_out => read_persist_real_0D
!
      character (len=*), intent(in) :: label
      real, intent(out) :: value
!
      if (lswitched_to_out) then
        read_persist_real_0D = read_persist_real_0D_out(label, value)
      else
        read_persist_real_0D = read_persist_real_0D_in(label, value)
      endif
!
    endfunction read_persist_real_0D
!***********************************************************************
    logical function read_persist_real_1D(label, value)

      use Io_in, only: read_persist_real_1D_in => read_persist_real_1D
      use Io_out, only: read_persist_real_1D_out => read_persist_real_1D
!
      character (len=*), intent(in) :: label
      real, dimension(:), intent(out) :: value
!
      if (lswitched_to_out) then
        read_persist_real_1D = read_persist_real_1D_out(label, value)
      else
        read_persist_real_1D = read_persist_real_1D_in(label, value)
      endif
!
    endfunction read_persist_real_1D
!***********************************************************************
    logical function read_persist_torus_rect(label,value)
!
!  Read persistent data from snapshot file.
!
!  16-May-2020/MR: coded
!
      use Geometrical_types
      use Io_in, only: read_persist_torus_rect_in => read_persist_torus_rect
      use Io_out, only: read_persist_torus_rect_out => read_persist_torus_rect

      character (len=*), intent(in) :: label
      type(torus_rect), intent(out) :: value
!
      if (lswitched_to_out) then
        read_persist_torus_rect = read_persist_torus_rect_out(label, value)
      else
        read_persist_torus_rect = read_persist_torus_rect_in(label, value)
      endif
!
    endfunction read_persist_torus_rect
!***********************************************************************
    logical function write_persist_logical_0D(label, id, value)

      use Io_out, only: write_persist_logical_0D_ => write_persist_logical_0D
!
      character (len=*), intent(in) :: label
      integer, intent(in) :: id
      logical, intent(in) :: value
!
      call switch_io
      write_persist_logical_0D = write_persist_logical_0D_(label, id, value)
!
    endfunction write_persist_logical_0D
!***********************************************************************
    logical function write_persist_logical_1D(label, id, value)

      use Io_out, only: write_persist_logical_1D_ => write_persist_logical_1D
!
      character (len=*), intent(in) :: label
      integer, intent(in) :: id
      logical, dimension(:), intent(in) :: value
!
      call switch_io
      write_persist_logical_1D = write_persist_logical_1D_(label, id, value)
!
    endfunction write_persist_logical_1D
!***********************************************************************
    logical function write_persist_int_0D(label, id, value)

      use Io_out, only: write_persist_int_0D_ => write_persist_int_0D
!
      character (len=*), intent(in) :: label
      integer, intent(in) :: id
      integer, intent(in) :: value
!
      call switch_io
      write_persist_int_0D = write_persist_int_0D_(label, id, value)
!
    endfunction write_persist_int_0D
!***********************************************************************
    logical function write_persist_int_1D(label, id, value)

      use Io_out, only: write_persist_int_1D_ => write_persist_int_1D
!
      character (len=*), intent(in) :: label
      integer, intent(in) :: id
      integer, dimension(:), intent(in) :: value
!
      call switch_io
      write_persist_int_1D = write_persist_int_1D_(label, id, value)
!
    endfunction write_persist_int_1D
!***********************************************************************
    logical function write_persist_real_0D(label, id, value)

      use Io_out, only: write_persist_real_0D_ => write_persist_real_0D
!
      character (len=*), intent(in) :: label
      integer, intent(in) :: id
      real, intent(in) :: value
!
      call switch_io
      write_persist_real_0D = write_persist_real_0D_(label, id, value)
!
    endfunction write_persist_real_0D
!***********************************************************************
    logical function write_persist_real_1D(label, id, value)

      use Io_out, only: write_persist_real_1D_ => write_persist_real_1D
!
      character (len=*), intent(in) :: label
      integer, intent(in) :: id
      real, dimension(:), intent(in) :: value
!
      call switch_io
      write_persist_real_1D = write_persist_real_1D_(label, id, value)
!
    endfunction write_persist_real_1D
!***********************************************************************
    logical function write_persist_torus_rect(label, id, value)
!
!  Write persistent data to snapshot file.
!
!  16-May-2020/MR: coded
!
      use Geometrical_types
      use Io_out, only: write_persist_torus_rect_ => write_persist_torus_rect

      character (len=*), intent(in) :: label
      integer, intent(in) :: id
      type(torus_rect), intent(in) :: value
!
      call switch_io
      write_persist_torus_rect = write_persist_torus_rect_(label, id, value)
!
    endfunction write_persist_torus_rect
!***********************************************************************
end module Io
