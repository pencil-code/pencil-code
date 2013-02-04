! $Id$
!
!  This module is an interface to allow modules
!  to register pointers to their internal variables so that
!  other modules my then request them by name.
!
!  This uses a linked list of pointers and is neither efficient
!  nor easy to use.  THIS IS ON PURPOSE (a deterrent)!
!
!  The main point of using shared variables is that it avoids
!  unnecesary dependencies between modules!
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
!  Currently only scalar, 1D, and 2D reals and integers may be used
!  2D could perhaps be added but 3D almost certainly should not
!  be shared this way.
!
!  19-jul-06/tony: coded
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
  public :: sharedvars_error_string
  public :: sharedvars_clean_up
!
  interface get_shared_variable
    module procedure get_variable_real0d
    module procedure get_variable_real1d
    module procedure get_variable_real2d
    module procedure get_variable_int0d
    module procedure get_variable_int1d
    module procedure get_variable_logical0d
    module procedure get_variable_logical1d
  endinterface
!
  interface put_shared_variable
    module procedure put_variable_real0d
    module procedure put_variable_real1d
    module procedure put_variable_real2d
    module procedure put_variable_int0d
    module procedure put_variable_int1d
    module procedure put_variable_logical0d
    module procedure put_variable_logical1d
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
  integer, parameter :: iSHVAR_TYPE_LOG0D=20
  integer, parameter :: iSHVAR_TYPE_LOG1D=21
  integer, parameter :: iSHVAR_TYPE_LOG2D=22
  integer, parameter :: iSHVAR_TYPE_CHAR0D=30
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
    real, dimension(:,:),  pointer :: real2d
    integer,               pointer :: int0d
    integer, dimension(:), pointer :: int1d
    logical,               pointer :: log0d
    logical, dimension(:), pointer :: log1d
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
!
  contains
!
!***********************************************************************
    subroutine initialize_shared_variables(lreloading)
!
!  Comment me.
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
!
!  Comment me.
!
      character (len=*) :: varname
      real, pointer :: variable
      integer, optional :: ierr
      type (shared_variable_list), pointer :: item
!
      intent(in)  :: varname
      intent(out) :: ierr
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
              print*, 'Getting shared variable: ',varname
              call fatal_error('get_variable', 'Data pointer is not associated.')
            endif
            return
          else
            nullify(variable)
            if (present(ierr)) then
              ierr=iSHVAR_ERR_WRONGTYPE
              return
            endif
            print*, 'Getting shared variable: ',varname
            call fatal_error('get_variable', 'Shared variable has the wrong type!')
          endif
        endif
        item=>item%next
      enddo
!
      nullify(variable)
!
      if (present(ierr)) then
        ierr=iSHVAR_ERR_NOSUCHVAR
        return
      endif
!
      print*, 'Getting shared variable: ',varname
      call fatal_error('get_variable', 'Shared variable does not exist!')
!
    endsubroutine get_variable_real0d
!***********************************************************************
    subroutine get_variable_real1d(varname,variable,ierr)
!
!  Comment me.
!
      character (len=*) :: varname
      real, dimension(:), pointer :: variable
      integer, optional :: ierr
      type (shared_variable_list), pointer :: item
!
      intent(in)  :: varname
      intent(out) :: ierr
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
              print*, 'Getting shared variable: ',varname
              call fatal_error('get_variable', 'Data pointer is not associated.')
            endif
            return
          else
            nullify(variable)
            if (present(ierr)) then
              ierr=iSHVAR_ERR_WRONGTYPE
              return
            endif
            print*, 'Getting shared variable: ',varname
            call fatal_error('get_variable', 'Shared variable has the wrong type!')
          endif
        endif
        item=>item%next
      enddo
!
      nullify(variable)
!
      if (present(ierr)) then
        ierr=iSHVAR_ERR_NOSUCHVAR
        return
      endif
!
      print*, 'Getting shared variable: ',varname
      call fatal_error('get_variable', 'Shared variable does not exist!')
!
    endsubroutine get_variable_real1d
!***********************************************************************
    subroutine get_variable_real2d(varname,variable,ierr)
!
!  get 2-D array
!
      character (len=*) :: varname
      real, dimension(:,:), pointer :: variable
      integer, optional :: ierr
      type (shared_variable_list), pointer :: item
!
      intent(in)  :: varname
      intent(out) :: ierr
!
      if (present(ierr)) ierr=0
!
      item=>thelist
      do while (associated(item))
        if (item%varname==varname) then
          if (item%vartype==iSHVAR_TYPE_REAL2D) then
            variable=>item%real2D
            if (.not.associated(item%real2D)) then
              if (present(ierr)) then
                ierr=iSHVAR_ERR_NOTASSOCIATED
                return
              endif
              print*, 'Getting shared variable: ',varname
              call fatal_error('get_variable', 'Data pointer is not associated.')
            endif
            return
          else
            nullify(variable)
            if (present(ierr)) then
              ierr=iSHVAR_ERR_WRONGTYPE
              return
            endif
            print*, 'Getting shared variable: ',varname
            call fatal_error('get_variable', 'Shared variable has the wrong type!')
          endif
        endif
        item=>item%next
      enddo
!
      nullify(variable)
!
      if (present(ierr)) then
        ierr=iSHVAR_ERR_NOSUCHVAR
        return
      endif
!
      print*, 'Getting shared variable: ',varname
      call fatal_error('get_variable', 'Shared variable does not exist!')
!
    endsubroutine get_variable_real2d
!***********************************************************************
    subroutine get_variable_int0d(varname,variable,ierr)
!
!  Comment me.
!
      character (len=*) :: varname
      integer, pointer :: variable
      integer, optional :: ierr
      type (shared_variable_list), pointer :: item
!
      intent(in)  :: varname
      intent(out) :: ierr
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
              print*, 'Getting shared variable: ',varname
              call fatal_error('get_variable', 'Data pointer is not associated.')
            endif
            return
          else
            nullify(variable)
            if (present(ierr)) then
              ierr=iSHVAR_ERR_WRONGTYPE
              return
            endif
            print*, 'Getting shared variable: ',varname
            call fatal_error('get_variable', 'Shared variable has the wrong type!')
          endif
        endif
        item=>item%next
      enddo
!
      nullify(variable)
!
      if (present(ierr)) then
        ierr=iSHVAR_ERR_NOSUCHVAR
        return
      endif
!
      print*, 'Getting shared variable: ',varname
      call fatal_error('get_variable', 'Shared variable does not exist!')
!
    endsubroutine get_variable_int0d
!***********************************************************************
    subroutine get_variable_int1d(varname,variable,ierr)
!
!  Comment me.
!
      character (len=*) :: varname
      integer, dimension(:), pointer :: variable
      integer, optional :: ierr
      type (shared_variable_list), pointer :: item
!
      intent(in)  :: varname
      intent(out) :: ierr
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
              print*, 'Getting shared variable: ',varname
              call fatal_error('get_variable', 'Data pointer is not associated.')
            endif
            return
          else
            nullify(variable)
            if (present(ierr)) then
              ierr=iSHVAR_ERR_WRONGTYPE
              return
            endif
            print*, 'Getting shared variable: ',varname
            call fatal_error('get_variable', 'Shared variable has the wrong type!')
          endif
        endif
        item=>item%next
      enddo
!
      nullify(variable)
!
      if (present(ierr)) then
        ierr=iSHVAR_ERR_NOSUCHVAR
        return
      endif
!
      print*, 'Getting shared variable: ',varname
      call fatal_error('get_variable', 'Shared variable does not exist!')
!
    endsubroutine get_variable_int1d
!***********************************************************************
    subroutine get_variable_logical0d(varname,variable,ierr)
!
!  Comment me.
!
      character (len=*) :: varname
      logical, pointer :: variable
      integer, optional :: ierr
      type (shared_variable_list), pointer :: item
!
      intent(in)  :: varname
      intent(out) :: ierr
!
      if (present(ierr)) ierr=0
!
      item=>thelist
      do while (associated(item))
        if (item%varname==varname) then
          if (item%vartype==iSHVAR_TYPE_LOG0D) then
            variable=>item%log0D
            if (.not.associated(item%log0D)) then
              if (present(ierr)) then
                ierr=iSHVAR_ERR_NOTASSOCIATED
                return
              endif
              print*, 'Getting shared variable: ', varname
              call fatal_error('get_variable', 'Data pointer is not associated.')
            endif
            return
          else
            nullify(variable)
            if (present(ierr)) then
              ierr=iSHVAR_ERR_WRONGTYPE
              return
            endif
            print*, 'Getting shared variable: ', varname
            call fatal_error('get_variable', 'Shared variable has the wrong type!')
          endif
        endif
        item=>item%next
      enddo
!
      nullify(variable)
!
      if (present(ierr)) then
        ierr=iSHVAR_ERR_NOSUCHVAR
        return
      endif
!
      print*, 'Getting shared variable: ' ,varname
      call fatal_error('get_variable', 'Shared variable does not exist!')
!
    endsubroutine get_variable_logical0d
!***********************************************************************
    subroutine get_variable_logical1d(varname,variable,ierr)
!
!  Comment me.
!
      character (len=*) :: varname
      logical, dimension(:), pointer :: variable
      integer, optional :: ierr
      type (shared_variable_list), pointer :: item
!
      intent(in)  :: varname
      intent(out) :: ierr
!
      if (present(ierr)) ierr=0
!
      item=>thelist
      do while (associated(item))
        if (item%varname==varname) then
          if (item%vartype==iSHVAR_TYPE_LOG1D) then
            variable=>item%log1D
            if (.not.associated(item%log1D)) then
              if (present(ierr)) then
                ierr=iSHVAR_ERR_NOTASSOCIATED
                return
              endif
              print*, 'Getting shared variable: ', varname
              call fatal_error('get_variable', 'Data pointer is not associated.')
            endif
            return
          else
            nullify(variable)
            if (present(ierr)) then
              ierr=iSHVAR_ERR_WRONGTYPE
              return
            endif
            print*, 'Getting shared variable: ', varname
            call fatal_error('get_variable', 'Shared variable has the wrong type!')
          endif
        endif
        item=>item%next
      enddo
!
      nullify(variable)
!
      if (present(ierr)) then
        ierr=iSHVAR_ERR_NOSUCHVAR
        return
      endif
!
      print*, 'Getting shared variable: ' ,varname
      call fatal_error('get_variable', 'Shared variable does not exist!')
!
    endsubroutine get_variable_logical1d
!***********************************************************************
    subroutine put_variable_int0d(varname,variable,ierr)
!
!  Comment me.
!
      character (len=*) :: varname
      integer, target :: variable
      integer, optional :: ierr
      type (shared_variable_list), pointer :: new
!
      intent(in)  :: varname
      intent(out) :: ierr
!
      if (present(ierr)) ierr=0
!
      new=>find_variable(varname)
      if (associated(new)) then
        if (associated(new%int0d,target=variable)) return
        if (present(ierr)) then
          ierr=iSHVAR_ERR_DUPLICATE
          return
        endif
        print*, 'Setting shared variable: ',varname
        call fatal_error('get_variable', 'Shared variable name already exists!')
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
!
!  Comment me.
!
      character (len=*) :: varname
      integer, dimension(:), target :: variable
      integer, optional :: ierr
      type (shared_variable_list), pointer :: new
!
      intent(in)  :: varname
      intent(out) :: ierr
!
      if (present(ierr)) ierr=0
!
      new=>find_variable(varname)
      if (associated(new)) then
        if (associated(new%int1d,target=variable)) return
        if (present(ierr)) then
          ierr=iSHVAR_ERR_DUPLICATE
          return
        endif
        print*, 'Setting shared variable: ',varname
        call fatal_error('get_variable', 'Shared variable name already exists!')
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
!  Comment me.
!
      character (len=*) :: varname
      real, target :: variable
      integer, optional :: ierr
      type (shared_variable_list), pointer :: new
!
      intent(in)  :: varname
      intent(out) :: ierr
!
      if (present(ierr)) ierr=0
!
      new=>find_variable(varname)
      if (associated(new)) then
        if (associated(new%real0D,target=variable)) return
        if (present(ierr)) then
          ierr=iSHVAR_ERR_DUPLICATE
          return
        endif
        print*, 'Setting shared variable: ',varname
        call fatal_error('get_variable', 'Shared variable name already exists!')
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
!
!  put 1-D array into shared variable
!
      character (len=*) :: varname
      real, dimension(:), target :: variable
      integer, optional :: ierr
      type (shared_variable_list), pointer :: new
!
      intent(in)  :: varname
      intent(out) :: ierr
!
      if (present(ierr)) ierr=0
!
      new=>find_variable(varname)
      if (associated(new)) then
        if (associated(new%real1D,target=variable)) return
        if (present(ierr)) then
          ierr=iSHVAR_ERR_DUPLICATE
          return
        endif
        print*, 'Setting shared variable: ',varname
        call fatal_error('get_variable', 'Shared variable name already exists!')
      endif
!
      call new_item_atstart(thelist,new=new)
      new%varname=varname
      new%vartype=iSHVAR_TYPE_REAL1D
      new%real1D=>variable
!
    endsubroutine put_variable_real1d
!***********************************************************************
    subroutine put_variable_real2d(varname,variable,ierr)
!
!  put 2-D array into shared variable
!
      character (len=*) :: varname
      real, dimension(:,:), target :: variable
      integer, optional :: ierr
      type (shared_variable_list), pointer :: new
!
      intent(in)  :: varname
      intent(out) :: ierr
!
      if (present(ierr)) ierr=0
!
      new=>find_variable(varname)
      if (associated(new)) then
        if (associated(new%real2D,target=variable)) return
        if (present(ierr)) then
          ierr=iSHVAR_ERR_DUPLICATE
          return
        endif
        print*, 'Setting shared variable: ',varname
        call fatal_error('get_variable', 'Shared variable name already exists!')
      endif
!
      call new_item_atstart(thelist,new=new)
      new%varname=varname
      new%vartype=iSHVAR_TYPE_REAL2D
      new%real2D=>variable
!
    endsubroutine put_variable_real2d
!***********************************************************************
    subroutine put_variable_logical0d(varname,variable,ierr)
!
!  Comment me.
!
      character (len=*) :: varname
      logical, target :: variable
      integer, optional :: ierr
      type (shared_variable_list), pointer :: new
!
      intent(in)  :: varname
      intent(out) :: ierr
!
      if (present(ierr)) ierr=0
!
      new=>find_variable(varname)
      if (associated(new)) then
        if (associated(new%log0D,target=variable)) return
        if (present(ierr)) then
          ierr=iSHVAR_ERR_DUPLICATE
          return
        endif
        print*, 'Setting shared variable: ',varname
        call fatal_error('get_variable', 'Shared variable name already exists!')
      endif
!
      call new_item_atstart(thelist,new=new)
      new%varname=varname
      new%vartype=iSHVAR_TYPE_LOG0D
      new%log0D=>variable
!
    endsubroutine put_variable_logical0d
!***********************************************************************
    subroutine put_variable_logical1d(varname,variable,ierr)
!
!  Comment me.
!
      character (len=*) :: varname
      logical, dimension(:), target :: variable
      integer, optional :: ierr
      type (shared_variable_list), pointer :: new
!
      intent(in)  :: varname
      intent(out) :: ierr
!
      if (present(ierr)) ierr=0
!
      new=>find_variable(varname)
      if (associated(new)) then
        if (associated(new%log1D,target=variable)) return
        if (present(ierr)) then
          ierr=iSHVAR_ERR_DUPLICATE
          return
        endif
        print*, 'Setting shared variable: ',varname
        call fatal_error('get_variable', 'Shared variable name already exists!')
      endif
!
      call new_item_atstart(thelist,new=new)
      new%varname=varname
      new%vartype=iSHVAR_TYPE_LOG1D
      new%log1D=>variable
!
    endsubroutine put_variable_logical1d
!***********************************************************************
    function find_variable(varname)
!
!  Comment me.
!
      character (len=*) :: varname
      type (shared_variable_list), pointer :: find_variable
!
      intent(in)  :: varname
!
      find_variable=>thelist
      do while (associated(find_variable))
        if (find_variable%varname==varname) then
          return
        endif
        find_variable=>find_variable%next
      enddo
!
      nullify(find_variable)
!
      return
!
    endfunction find_variable
!***********************************************************************
    subroutine free_list(list)
!
!  Comment me.
!
      type (shared_variable_list), pointer :: list
      type (shared_variable_list), pointer :: next
!
      do while (associated(list))
        next=>list%next
        deallocate(list)
        list=>next
      enddo
!
      nullify(list)
!
    endsubroutine free_list
!***********************************************************************
    subroutine new_item_atstart(list,new)
!
!  Comment me.
!
      type (shared_variable_list), pointer :: list
      type (shared_variable_list), optional, pointer :: new
!
      type (shared_variable_list), pointer :: new_
!
      allocate(new_)
      new_%next=>list
      list=>new_
      if (present(new)) new=>new_
!
    endsubroutine new_item_atstart
!***********************************************************************
    function sharedvars_error_string(ierr) result(string)
!
!  Tell what went wrong
!
!  29-aug-06/wolf: adapted
!
      use Cparam, only: labellen
!
      character(len=labellen) :: string
      integer                 :: ierr
!
      intent(in)  :: ierr
!
      select case (ierr)
      case (iSHVAR_ERR_NOSUCHVAR);     string='SHVAR_ERR_NOSUCHVAR'
      case (iSHVAR_ERR_WRONGTYPE);     string='SHVAR_ERR_WRONGTYPE'
      case (iSHVAR_ERR_DUPLICATE);     string='SHVAR_ERR_DUPLICATE'
      case (iSHVAR_ERR_NOTASSOCIATED); string='SHVAR_ERR_NOTASSOCIATED'
      case default;                   string='Undocumented ierr!'
      endselect
!
    endfunction sharedvars_error_string
!***********************************************************************
    subroutine sharedvars_clean_up()
!
!  Free any allocated memory, so G95 does not complain after the run
!
!  18-may-2007/wolf: coded
!
      call free_list(thelist)
!
    endsubroutine sharedvars_clean_up
!***********************************************************************
endmodule SharedVariables
