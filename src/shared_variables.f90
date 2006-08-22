! $Id: shared_variables.f90,v 1.8 2006-08-22 12:21:48 mee Exp $ 
!
!  This module is an interface to allow modules
!  to register pointers to their internal variables so that
!  other modules my then request them by name.
!
!  This uses a linked list of pointers and is neither efficient
!  nor easy to use.  THIS IS ON PURPOSE (a deterrent)!
!
!  Shared variable should always be avoided for portability 
!  and generality reasons.  This module allows for the possibility
!  when needs must but trys to protect agains screw ups that
!  can derive from shared quantities.
!
!  When used modules should call the get and put routines
!  at initialize_ time for optimal performance.
!
!  Variables to be added to the list must have the target property
!  And a pointer must be provided when getting a variable from
!  the list. 
!
!  Currently only scalar and 1D reals and integers may be used
!  2D could perhaps be added but 3D almost certainly should not
!  be shared this way.
!
module SharedVariables 
!
  use Messages
!
  implicit none
!
  private
!
  public :: initialize_shared_variables
  public :: put_shared_variable
  public :: get_shared_variable
!
  interface get_shared_variable
    module procedure get_variable_real0d
    module procedure get_variable_real1d
    module procedure get_variable_int0d
    module procedure get_variable_int1d
!    module procedure get_variable_char
  endinterface
!
  interface put_shared_variable
    module procedure put_variable_real0d
    module procedure put_variable_real1d
    module procedure put_variable_int0d
    module procedure put_variable_int1d
!    module procedure put_variable_char
  endinterface
!
! Used internally to keep track ot the type of data
! stored in the shared variables list.
!
  integer, parameter :: iSHVAR_TYPE_REAL0D=1
  integer, parameter :: iSHVAR_TYPE_REAL1D=2
  integer, parameter :: iSHVAR_TYPE_REAL2D=3
  integer, parameter :: iSHVAR_TYPE_INT0D=10
  integer, parameter :: iSHVAR_TYPE_INT1D=11
  integer, parameter :: iSHVAR_TYPE_INT2D=12
  integer, parameter :: iSHVAR_TYPE_CHAR0D=20
!
! Some possible error codes when getting a variable
! (if the user doesn't even ask for the error code
!  they get a fatal_error)
!
  integer, public, parameter :: iSHVAR_ERR_NOSUCHVAR=1
  integer, public, parameter :: iSHVAR_ERR_WRONGTYPE=2
  integer, public, parameter :: iSHVAR_ERR_DUPLICATE=3
  integer, public, parameter :: iSHVAR_ERR_NOTASSOCIATED=4
!
! Store pointers to variables in a general linked
! list structure.  Shame we can't have (void *) pointers.
!
  type shared_variable_list
!
! Shared variable metadata
!
    character (len=30) :: varname
    integer            :: vartype
!
! Possible data types
!
    real,                  pointer :: real0d
    real, dimension(:),    pointer :: real1d
    integer,               pointer :: int0d
    integer, dimension(:), pointer :: int1d
!    character (len=*),     pointer :: char0D
!
! Linked list link to next list element 
!
    type (shared_variable_list), pointer :: next
  endtype shared_variable_list
!
! The head of the list (initially empty)
!
  type (shared_variable_list), pointer :: thelist
 
  contains

!***********************************************************************
    subroutine initialize_shared_variables(lreloading)
!
      logical :: lreloading
!
      if (lreloading) then
        call free_list(thelist)
      else
        NULLIFY(thelist)
      endif
!
    endsubroutine initialize_shared_variables
!***********************************************************************
    subroutine get_variable_real0d(varname,variable,ierr) 
      character (len=*) :: varname
      real, pointer :: variable
      type (shared_variable_list), pointer :: item
      integer, optional :: ierr
!
      if (present(ierr)) ierr=0
!
      item=>thelist
      do while (associated(item))
        if (item%varname==varname) then
          if (item%vartype==iSHVAR_TYPE_REAL0D) then
            variable=>item%real0D
            if (.not.associated(item%real0D)) then
              if (present(ierr)) then
                ierr=iSHVAR_ERR_NOTASSOCIATED
                return
              endif
              print*,"Getting shared variable: ",varname
              call fatal_error("get_variable","Data pointer is not associated.")
            endif
            return
          else
            nullify(variable)
            if (present(ierr)) then
              ierr=iSHVAR_ERR_WRONGTYPE
              return
            endif
            print*,"Getting shared variable: ",varname
            call fatal_error("get_variable","Shared variable has the wrong type!")
          endif
        endif
        item=>item%next
      enddo
      nullify(variable)
      if (present(ierr)) then
        ierr=iSHVAR_ERR_NOSUCHVAR
        return
      endif
      print*,"Getting shared variable: ",varname
      call fatal_error("get_variable","Shared variable does not exist!")
!    
    endsubroutine get_variable_real0d
!***********************************************************************
    subroutine get_variable_real1d(varname,variable,ierr) 
      character (len=*) :: varname
      real, dimension(:), pointer :: variable
      type (shared_variable_list), pointer :: item
      integer, optional :: ierr
!
      if (present(ierr)) ierr=0
!
      item=>thelist
      do while (associated(item))
        if (item%varname==varname) then
          if (item%vartype==iSHVAR_TYPE_REAL1D) then
            variable=>item%real1D
            if (.not.associated(item%real1D)) then
              if (present(ierr)) then
                ierr=iSHVAR_ERR_NOTASSOCIATED
                return
              endif
              print*,"Getting shared variable: ",varname
              call fatal_error("get_variable","Data pointer is not associated.")
            endif
            return
          else
            nullify(variable)
            if (present(ierr)) then
              ierr=iSHVAR_ERR_WRONGTYPE
              return
            endif
            print*,"Getting shared variable: ",varname
            call fatal_error("get_variable","Shared variable has the wrong type!")
          endif
        endif
        item=>item%next
      enddo
      nullify(variable)
      if (present(ierr)) then
        ierr=iSHVAR_ERR_NOSUCHVAR
        return
      endif
      print*,"Getting shared variable: ",varname
      call fatal_error("get_variable","Shared variable does not exist!")
!    
    endsubroutine get_variable_real1d
!***********************************************************************
    subroutine get_variable_int0d(varname,variable,ierr) 
      character (len=*) :: varname
      integer, pointer :: variable
      type (shared_variable_list), pointer :: item
      integer, optional :: ierr
!
      if (present(ierr)) ierr=0
!
      item=>thelist
      do while (associated(item))
        if (item%varname==varname) then
          if (item%vartype==iSHVAR_TYPE_INT0D) then
            variable=>item%int0D
            if (.not.associated(item%int0D)) then
              if (present(ierr)) then
                ierr=iSHVAR_ERR_NOTASSOCIATED
                return
              endif
              print*,"Getting shared variable: ",varname
              call fatal_error("get_variable","Data pointer is not associated.")
            endif
            return
          else
            nullify(variable)
            if (present(ierr)) then
              ierr=iSHVAR_ERR_WRONGTYPE
              return
            endif
            print*,"Getting shared variable: ",varname
            call fatal_error("get_variable","Shared variable has the wrong type!")
          endif
        endif
        item=>item%next
      enddo
      nullify(variable)
      if (present(ierr)) then
        ierr=iSHVAR_ERR_NOSUCHVAR
        return
      endif
      print*,"Getting shared variable: ",varname
      call fatal_error("get_variable","Shared variable does not exist!")
!    
    endsubroutine get_variable_int0d
!***********************************************************************
    subroutine get_variable_int1d(varname,variable,ierr) 
      character (len=*) :: varname
      integer, dimension(:), pointer :: variable
      type (shared_variable_list), pointer :: item
      integer, optional :: ierr
!
      if (present(ierr)) ierr=0
!
      item=>thelist
      do while (associated(item))
        if (item%varname==varname) then
          if (item%vartype==iSHVAR_TYPE_INT1D) then
            variable=>item%int1D
            if (.not.associated(item%int1D)) then
              if (present(ierr)) then
                ierr=iSHVAR_ERR_NOTASSOCIATED
                return
              endif
              print*,"Getting shared variable: ",varname
              call fatal_error("get_variable","Data pointer is not associated.")
            endif
            return
          else
            nullify(variable)
            if (present(ierr)) then
              ierr=iSHVAR_ERR_WRONGTYPE
              return
            endif
            print*,"Getting shared variable: ",varname
            call fatal_error("get_variable","Shared variable has the wrong type!")
          endif
        endif
        item=>item%next
      enddo
      nullify(variable)
      if (present(ierr)) then
        ierr=iSHVAR_ERR_NOSUCHVAR
        return
      endif
      print*,"Getting shared variable: ",varname
      call fatal_error("get_variable","Shared variable does not exist!")
!    
    endsubroutine get_variable_int1d
!***********************************************************************
    subroutine put_variable_int0d(varname,variable,ierr) 
      character (len=*) :: varname
      integer, target :: variable
      type (shared_variable_list), pointer :: new
      integer, optional :: ierr
!
      if (present(ierr)) ierr=0
!
      if (variable_exists(varname)) then
        if (present(ierr)) then
          ierr=iSHVAR_ERR_DUPLICATE
          return
        endif
        print*,"Setting shared variable: ",varname
        call fatal_error("get_variable","Shared variable name already exists!")
      endif
!
      call new_item_atstart(thelist,new=new)
      new%varname=varname
      new%vartype=iSHVAR_TYPE_INT0D
      new%int0D=>variable
!    
    endsubroutine put_variable_int0d
!***********************************************************************
    subroutine put_variable_int1d(varname,variable,ierr) 
      character (len=*) :: varname
      integer, dimension(:), target :: variable
      type (shared_variable_list), pointer :: new
      integer, optional :: ierr
!
      if (present(ierr)) ierr=0
!
      if (variable_exists(varname)) then
        if (present(ierr)) then
          ierr=iSHVAR_ERR_DUPLICATE
          return
        endif
        print*,"Setting shared variable: ",varname
        call fatal_error("get_variable","Shared variable name already exists!")
      endif
!
      call new_item_atstart(thelist,new=new)
      new%varname=varname
      new%vartype=iSHVAR_TYPE_INT1D
      new%int1D=>variable
!    
    endsubroutine put_variable_int1d
!***********************************************************************
    subroutine put_variable_real0d(varname,variable,ierr) 
!
      character (len=*) :: varname
      real, target :: variable
      integer, optional :: ierr
!
      type (shared_variable_list), pointer :: new
!
      if (present(ierr)) ierr=0
!
      if (variable_exists(varname)) then
        if (present(ierr)) then
          ierr=iSHVAR_ERR_DUPLICATE
          return
        endif
        print*,"Setting shared variable: ",varname
        call fatal_error("get_variable","Shared variable name already exists!")
      endif
!
      call new_item_atstart(thelist,new=new)
      new%varname=varname
      new%vartype=iSHVAR_TYPE_REAL0D
      new%real0D=>variable
!    
    endsubroutine put_variable_real0d
!***********************************************************************
    subroutine put_variable_real1d(varname,variable,ierr) 
      character (len=*) :: varname
      real, dimension(:), target :: variable
      type (shared_variable_list), pointer :: new
      integer, optional :: ierr
!
      if (present(ierr)) ierr=0
!
      if (variable_exists(varname)) then
        if (present(ierr)) then
          ierr=iSHVAR_ERR_DUPLICATE
          return
        endif
        print*,"Setting shared variable: ",varname
        call fatal_error("get_variable","Shared variable name already exists!")
      endif
!
      call new_item_atstart(thelist,new=new)
      new%varname=varname
      new%vartype=iSHVAR_TYPE_REAL1D
      new%real1D=>variable
!    
    endsubroutine put_variable_real1d
!***********************************************************************
    function variable_exists(varname) 
      character (len=*) :: varname
      logical :: variable_exists
      type (shared_variable_list), pointer :: item
!
      item=>thelist
      do while (associated(item))
        if (item%varname==varname) then
          variable_exists=.true.
          return
        endif
        item=>item%next
      enddo
      variable_exists=.false.
      return
!    
    endfunction variable_exists
!***********************************************************************
    subroutine free_list(list) 
      type (shared_variable_list), pointer :: list
      type (shared_variable_list), pointer :: next
 
      do while (associated(list))
        next=>list%next
        deallocate(list)
        list=>next
      enddo
      nullify(list)
    endsubroutine free_list
!***********************************************************************
    subroutine new_item_atstart(list,new) 
      type (shared_variable_list), pointer :: list
      type (shared_variable_list), optional, pointer :: new
      type (shared_variable_list), pointer :: new_

      allocate(new_)
      new_%next=>list 
      list=>new_
      if (present(new)) new=>new_ 
    endsubroutine new_item_atstart
!***********************************************************************
endmodule SharedVariables
