! $Id$
!***********************************************************************
!
! This program combines one set of PVAR files from individual processor
! directories into one file under data/allprocs/.
!
! 07-nov-20/ccyang: coded
!
!***********************************************************************
program pc_collect_pvar
!
  use Cparam
!
  implicit none
!
  integer, parameter :: lun_in = 1, lun_out = 2
!
  integer, dimension(ncpus) :: npar_proc = 0
  character(len=fnlen) :: pvar = "pvar.dat"
  integer :: ierr
!
! Sanity check.
!
  if (sanity_check()) stop "failed sanity check. "
!
! Prompt for the name of the data file.
!
  print *, "Enter the name of the particle data file (e.g., PVAR0): "
  read *, pvar
!
! Get number of particles in each process.
!
  call get_npar(pvar, npar_proc)
!
! Read and combine particle data files.
!
  call read_and_combine(pvar, npar_proc, ierr)
  if (ierr /= 0) stop "unable to combine the data files. "
!
contains
!***********************************************************************
  subroutine get_npar(pvar, npar_proc)
!
! Reads number of particles in each process from the corresponding data
! file.
!
! 07-nov-20/ccyang: coded
!
! Input Arguments
!     pvar
!         Name of the particle data file
!
! Output Arguments
!     npar_proc
!         Number of particles in each process
!
    character(len=*), intent(in) :: pvar
    integer, dimension(ncpus), intent(out) :: npar_proc
!
    character(len=fnlen) :: fpath
    integer :: i, n
!
! Name the data file under each process directory.
!
    proc: do i = 1, ncpus
!
! Open the data file.
!
      call open_pvar(i - 1, pvar, lun_in, fpath)
!
! Read the number of particles.
!
      read(unit=lun_in) n
      print *, trim(fpath), ":", n, "particles"
      npar_proc(i) = n
!
! Close the file.
!
      close(unit=lun_in)
    enddo proc
!
  endsubroutine get_npar
!***********************************************************************
  subroutine open_pvar(iproc, pvar, lun, fpath)
!
! Open one particle data file under one process directory.
!
! 07-nov-20/ccyang: coded
!
! Input Arguments
!     iproc
!         Process ID
!     pvar
!         Name of the particle data file
!     lun
!         I/O unit assigned to the file
!
! Output Arguments
!     fpath
!         Path to the data file
!
    use Cdata, only: datadir, datadir_snap, directory, iproc_world
    use General, only: directory_names_std
!
    character(len=*), intent(in) :: pvar
    integer, intent(in) :: iproc, lun
    character(len=*), intent(out) :: fpath
!
! Compose the path to the data file.
!
    datadir_snap = datadir
    iproc_world = iproc
    call directory_names_std()
    fpath = trim(directory) // "/" // pvar
!
! Open the file.
!
    open(unit=lun, file=fpath, status="old", action="read", form="unformatted")
!
  endsubroutine open_pvar
!***********************************************************************
  subroutine read_and_combine(pvar, npar_proc, ierr)
!
! Reads particle data from individual process directories, combines
! them, and then writes one single combined file under data/allprocs/.
!
! 07-nov-20/ccyang: coded
!
! Input Arguments
!     pvar
!         Name of the particle data file
!     npar_proc
!         Number of particles in each process
!
! Output Arguments
!     ierr
!         Error code; zero for success
!
    use Cdata, only: directory_snap
!
    integer, dimension(ncpus), intent(in) :: npar_proc
    character(len=*), intent(in) :: pvar
    integer, intent(out) :: ierr
!
    real, allocatable, dimension(:,:) :: fp
    integer, allocatable, dimension(:) :: ipar, ibuf
!
    character(len=fnlen) :: fout, fin
    integer :: intlen, reallen, rlen, npar_tot, offset, nbuf, irec1, irec2, irec
    integer :: i, j, np
    real :: t = -1.0, tread
!
! Inquire IO lengths.
!
    inquire(iolength=intlen) 0
    inquire(iolength=reallen) 0.0
    print 10, "integer length = ", intlen, " bytes"
    print 10, "real length = ", reallen, " bytes"
    10 format (1x, a, i2, a)
    errlen: if (mod(reallen, intlen) /= 0) then
      print *, "read_and_combine: integer and real are not aligned. "
      ierr = -1
      return
    endif errlen
    rlen = reallen / intlen
!
! Open target file.
!
    fout = trim(directory_snap) // "/" // trim(pvar)
    open(unit=lun_out, file=fout, status="replace", action="write", form="unformatted", access="direct", recl=intlen)
!
! Find and write total number of particles.
!
    npar_tot = sum(npar_proc)
    write(unit=lun_out, rec=1) npar_tot
    offset = rlen * npar_tot
!
! Initialize record pointers.
!
    irec1 = 1
    irec2 = irec1 + npar_tot
!
    proc: do i = 1, ncpus
!
! Open each particle data file.
!
      call open_pvar(i - 1, pvar, lun_in, fin)
      print *, "Processing ", trim(fin), "......"
!
! Read and check number of particles.
!
      read(unit=lun_in) np
      errnr1: if (np /= npar_proc(i)) then
        print *, "read_and_combine: inconsistent number of particles. "
        print *, "read_and_combine: i, np, npar_proc = ", i, np, npar_proc(i)
        ierr = -2
        return
      endif errnr1
!
! Allocate working arrays for read.
!
      nbuf = rlen * np
      allocate(ipar(np), fp(np,mparray), ibuf(nbuf), stat=ierr)
      errall: if (ierr /= 0) then
        print *, "read_and_combine: unable to allocate working arrays"
        return
      endif errall
!
! Read particle data.
!
      read(unit=lun_in) ipar
      read(unit=lun_in) fp
      read(unit=lun_in) tread
!
! Check time consistency.
!
      errt1: if (t < 0.0) then
        t = tread
      elseif (t /= tread) then errt1
        print *, "read_and_combine: inconsistent time"
        print *, "read_and_combine: t, tread = ", t, tread
        ierr = -3
        return
      endif errt1
!
! Close data file.
!
      close(unit=lun_in)
!
! Write particle tags.
!
      call write_buffer(lun_out, irec1, ipar)
      irec1 = irec1 + np
!
! Write particle data.
!
      irec = irec2
      wfp: do j = 1, mparray
        ibuf = transfer(fp(:,j), ibuf)
        call write_buffer(lun_out, irec, ibuf)
        irec = irec + offset
      enddo wfp
      irec2 = irec2 + nbuf
!
! Deallocate working arrays.
!
      deallocate(ipar, fp, ibuf)
!
    enddo proc
!
! Write time.
!
    allocate(ibuf(rlen))
    ibuf = transfer(t, ibuf)
    irec = irec1 + mparray * offset
    call write_buffer(lun_out, irec, ibuf)
    deallocate(ibuf)
!
! Close the target file.
!
    close(unit=lun_out)
    ierr = 0
!
  endsubroutine read_and_combine
!***********************************************************************
  logical function sanity_check() result(error)
!
! Returns .false. if the preconditions of the program are met; .true.,
! otherwise.  An error message is also printed when the latter occurs.
!
! 07-nov-20/ccyang: coded
!
    use Cdata, only: datadir
!
    logical :: exists
!
    error = .true.
!
! Check if the file datadir.in exists.
!
    inquire(file='datadir.in', exist=exists)
    err1: if (exists) then
      print *, "sanity_check: not implemented for file datadir.in. "
      return
    endif err1
!
! Check if the file data/directory_snap exists.
!
    inquire(file=trim(datadir)//'/directory_snap', exist=exists)
    err2: if (exists) then
      print *, "sanity_check: not implemented for file ", trim(datadir), "/directory_snap. "
      return
    endif err2
!
    error = .false.
!
  endfunction sanity_check
!***********************************************************************
  subroutine write_buffer(lun, pos, ibuf)
!
! Writes an array of integers to a direct-access file.
!
! 07-nov-20/ccyang: coded
!
! Input Arguments
!     lun
!         I/O unit associated with a direct-access file for write, with a record length of one integer
!     pos
!         Record position right before where the first integer to be written
!     ibuf
!         The array of integers
!
    integer, dimension(:), intent(in) :: ibuf
    integer, intent(in) :: lun, pos
!
    integer :: i
!
! Write the integers.
!
    do i = 1, size(ibuf)
      write(unit=lun, rec=pos+i) ibuf(i)
    enddo
!
  endsubroutine write_buffer
!***********************************************************************
endprogram pc_collect_pvar
