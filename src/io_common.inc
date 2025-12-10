!***********************************************************************
    subroutine rproc_bounds(file)
!
!   Import processor boundaries from file.
!
!   10-jul-08/kapelrud: coded
!   16-Feb-2012/Bourdin.KIS: rewritten
!
      use Mpicomm, only: mpibcast

      character (len=*) :: file
      integer :: ierr
!
      real(KIND=rkind4), dimension(0:nprocx):: procx_boundssg
      real(KIND=rkind4), dimension(0:nprocy):: procy_boundssg
      real(KIND=rkind4), dimension(0:nprocz):: procz_boundssg
!
      real(KIND=rkind8), dimension(0:nprocx):: procx_boundsdb
      real(KIND=rkind8), dimension(0:nprocy):: procy_boundsdb
      real(KIND=rkind8), dimension(0:nprocz):: procz_boundsdb
!
      if (lroot) then

        open(lun_input,FILE=file,FORM='unformatted',status='old',action='read',iostat=ierr)
        if (ierr /= 0) call fatal_error("rproc_bounds", "Cannot open "//trim(file)//" for reading")

        if (lread_from_other_prec) then
          if (kind(x)==rkind4) then
            call input_proc_bounds(procx_boundsdb,procy_boundsdb,procz_boundsdb)
            procx_bounds=procx_boundsdb; procy_bounds=procy_boundsdb; procz_bounds=procz_boundsdb
          elseif (kind(x)==rkind8) then
            call input_proc_bounds(procx_boundssg,procy_boundssg,procz_boundssg)
            procx_bounds=procx_boundssg; procy_bounds=procy_boundssg; procz_bounds=procz_boundssg
          endif
        else
          call input_proc_bounds(procx_bounds,procy_bounds,procz_bounds)
        endif

        close(lun_output)

      endif

      call mpibcast(procx_bounds,nprocx+1)
      call mpibcast(procy_bounds,nprocy+1)
      call mpibcast(procz_bounds,nprocz+1)
!
    endsubroutine rproc_bounds
!***********************************************************************
    subroutine wproc_bounds(file)
!
!  Export processor boundaries to file.
!
!  22-Feb-2012/PABourdin: adapted from io_dist
!  27-nov-2020/ccyang: make the file single
!
      character(len=*), intent(in) :: file
!
      integer :: ierr
!
!  Only one process is needed.
!
      if (.not. lroot) return
!
!  Write proc[xyz]_bounds.
!
      open (lun_output, FILE=file, FORM='unformatted', IOSTAT=ierr, status='replace')
      if (ierr /= 0) call fatal_error("wproc_bounds", "Cannot open "//trim(file)//" for writing")

      write (lun_output) procx_bounds
      write (lun_output) procy_bounds
      write (lun_output) procz_bounds
      close (lun_output)
!
    endsubroutine wproc_bounds
!***********************************************************************
    subroutine input_proc_bounds_double(procx_bounds,procy_bounds,procz_bounds)
!
!   Import processor boundaries from file.in double precision
!
!   23-oct-13/MR: derived from rproc_bounds
!
      real(KIND=rkind8), dimension(0:nprocx), intent(OUT):: procx_bounds
      real(KIND=rkind8), dimension(0:nprocy), intent(OUT):: procy_bounds
      real(KIND=rkind8), dimension(0:nprocz), intent(OUT):: procz_bounds
!
      read(lun_input) procx_bounds
      read(lun_input) procy_bounds
      read(lun_input) procz_bounds
!
    endsubroutine input_proc_bounds_double
!***********************************************************************
    subroutine input_proc_bounds_single(procx_bounds,procy_bounds,procz_bounds)
!
!   Import processor boundaries from file.in single precision
!
!   23-oct-13/MR: derived from rproc_bounds
!
      real(KIND=rkind4), dimension(0:nprocx), intent(OUT):: procx_bounds
      real(KIND=rkind4), dimension(0:nprocy), intent(OUT):: procy_bounds
      real(KIND=rkind4), dimension(0:nprocz), intent(OUT):: procz_bounds
!
      read(lun_input) procx_bounds
      read(lun_input) procy_bounds
      read(lun_input) procz_bounds
!
    endsubroutine input_proc_bounds_single
!***********************************************************************
    subroutine output_ode(file)
!
!  Write ODE snapshot file with time.
!
!  08-Sep-2023/MR: coded
!     
      use FarrayManager, only: farray_retrieve_metadata_ode
      
      character (len=*) :: file
!     
      character(len=7) :: file_
      character(LEN=30), dimension(n_odevars) :: names
      integer, dimension(n_odevars) :: lengs
      integer :: num,lun
!     
      if (.not. lroot) return
!     
      if (file(1:3)=='VAR') then
        file_='ODE'//file(4:)
      else
        file_='ode.dat'
      endif
      
      lun = lun_output+20
      open(lun,FILE=trim(directory_collect)//'/'//trim(file_),FORM='unformatted')
      if (n_odevars > 0) then
        num = farray_retrieve_metadata_ode(names,lengs)
        write(lun) n_odevars,num
        write(lun) names(n_odevars-num+1:)
        write(lun) lengs(n_odevars-num+1:)
        write(lun) f_ode
        write(lun) real(t)
      endif
      close(lun)
      
      if (ip<=10) print*,'written ODE snapshot ', trim (file_)
    
    endsubroutine output_ode
!***********************************************************************
    subroutine input_ode(file)
!
!  Read ODE snapshot file with time.
!
!  08-Sep-2023/MR: coded
!
      use FarrayManager, only: farray_retrieve_metadata_ode
      use File_Io, only: file_exists
      use General, only: rtoa
      use Mpicomm, only: mpibcast

      character (len=*), intent(IN) :: file

      character (len=7) :: file_
      integer :: n_in, num_in, num
      character(LEN=30), dimension(n_odevars) :: names
      character(LEN=30), dimension(:), allocatable :: names_in
      integer, dimension(n_odevars) :: lengs
      integer, dimension(:), allocatable :: lengs_in
      real :: t_in
        
      if (lroot) then
!
        if (file(1:3)=='VAR') then
          file_='ODE'//file(4:)
        else
          file_='ode.dat'
        endif
!
        if (file_exists(trim(directory_collect)//'/'//trim(file_))) then
          if (ip<=8) print*, 'read ODE snapshot', trim (file_)
          open(lun_input,FILE=trim(directory_collect)//'/'//trim(file_),FORM='unformatted')
          num = farray_retrieve_metadata_ode(names,lengs)
        
          read(lun_input) n_in, num_in
          if (n_in /= n_odevars) call fatal_error("input_ode", &
                                 "numbers of slots differ between "//trim(file)//" and 'f_ode' array")
          if (num_in /= num) call fatal_error("input_ode", &
                                  "numbers of variables differ between "//trim(file)//" and 'f_ode' array")
          if (n_in > 0) then
            allocate(lengs_in(num), names_in(num))
            read(lun_input) names_in
            if (any(names_in /= names(n_odevars-num+1:))) call fatal_error("input_ode","variable names differ")
            read(lun_input) lengs_in
            if (any(lengs_in /= lengs(n_odevars-num+1:))) call fatal_error("input_ode","variable lengths differ")
            read(lun_input) f_ode
          endif
        
          read(lun_input) t_in
          if (t_in/=real(t)) call warning("input_ode","times differ between "//trim(file_)//" and " &
              //trim(file)//": "//trim(rtoa(t_in))//"/="//trim(rtoa(real(t))))
          close(lun_input) 
        else
          call warning('input_ode','no ODE data available')
        endif
      endif

      call mpibcast(f_ode,n_odevars)
  
    endsubroutine input_ode
!***********************************************************************

