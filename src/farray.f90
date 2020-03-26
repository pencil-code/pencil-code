! $Id$
!
!  This module allocates and manages indices in the f-array
!  in a controlled way.  This includes handling different
!  types of variable that may be stored in the f-array like
!  PDE variables, auxiliaries, communicated auxiliaries
!  and global variables.
!
!  The code follows the similar principles to the SharedVariables
!  module but only handles integer indices.
!
module FArrayManager
!
  use Cparam, only: mvar,maux,mglobal,maux_com,mscratch
  use Cdata, only: nvar,naux,naux_com,datadir,lroot,lwrite_aux,lreloading
  use HDF5_IO
  use Messages
!
  implicit none
!
  private
!
  public :: farray_register_variable
  public :: farray_register_pde
  public :: farray_register_auxiliary
  public :: farray_register_global
  public :: farray_use_variable
  public :: farray_index_append
  public :: farray_index_reset
  public :: farray_use_pde
  public :: farray_use_auxiliary
  public :: farray_use_global
  public :: farray_size_by_name
  public :: farray_type_by_name
!
  public :: farray_acquire_scratch_area
  public :: farray_release_scratch_area
!
  public :: farray_clean_up
!
! Types of variable that may be stored in the f-array
!
  integer, public, parameter :: iFARRAY_TYPE_NOSUCHTYPE=0
  integer, public, parameter :: iFARRAY_TYPE_PDE=1
  integer, public, parameter :: iFARRAY_TYPE_COMM_AUXILIARY=2
  integer, public, parameter :: iFARRAY_TYPE_AUXILIARY=3
  integer, public, parameter :: iFARRAY_TYPE_GLOBAL=4
  integer, public, parameter :: iFARRAY_TYPE_SCRATCH=5
!
! Some possible error codes when getting a variable
! (if the user doesn't even ask for the error code
!  they get a fatal_error)
!
  integer, public, parameter :: iFARRAY_ERR_NOSUCHVAR=1
  integer, public, parameter :: iFARRAY_ERR_NOSUCHCOMPONENT=2
  integer, public, parameter :: iFARRAY_ERR_WRONGTYPE=3
  integer, public, parameter :: iFARRAY_ERR_WRONGSIZE=4
  integer, public, parameter :: iFARRAY_ERR_DUPLICATE=5
  integer, public, parameter :: iFARRAY_ERR_INDEXMISMATCH=6
  integer, public, parameter :: iFARRAY_ERR_OUTOFSPACE=7
!
  type pp
    integer, pointer :: p
  endtype pp
!
! Store pointers to variables in a general linked
! list structure.  Shame we can't have (void *) pointers.
!
  type farray_contents_list
!
! f-array variable metadata
!
    character (len=30) :: varname
    integer            :: ncomponents
    integer            :: vartype
    type(pp), dimension(:), pointer :: ivar
!
! Linked list link to next list element
!
    type (farray_contents_list), pointer :: next
    type (farray_contents_list), pointer :: previous
  endtype farray_contents_list
!
! The head of the list (initially empty)
!
  type (farray_contents_list), pointer :: thelist
!
! Keep track of which spaces are currently in use.
!
  logical, dimension(mscratch+1) :: scratch_used
!
! Type counters (need to move nvar, naux and naux_com here too)
!
  integer :: nscratch=0
  integer :: nglobal=0

  contains

!***********************************************************************
    subroutine farray_register_pde(varname,ivar,vector,array,ierr)
!
!  Register a PDE variable in the f array.
!
      character (len=*), intent(in) :: varname
      integer, intent(out)  :: ivar
      integer, optional, intent(in) :: vector, array
      integer, optional, intent(out) :: ierr
!
      integer, parameter :: vartype = iFARRAY_TYPE_PDE
!
      call farray_register_variable(varname,ivar,vartype,vector=vector,array=array,ierr=ierr)
!
    endsubroutine farray_register_pde
!***********************************************************************
    subroutine farray_register_global(varname,ivar,vector,array,ierr)
!
!  Register a global variable in the f array.
!
      character (len=*), intent(in) :: varname
      integer, intent(out)  :: ivar
      integer, optional, intent(in) :: vector, array
      integer, optional, intent(out) :: ierr
!
      integer, parameter :: vartype = iFARRAY_TYPE_GLOBAL
!
      call farray_register_variable(varname,ivar,vartype,vector=vector,array=array,ierr=ierr)
!
    endsubroutine farray_register_global
!***********************************************************************
    subroutine farray_register_auxiliary(varname,ivar,communicated,vector,array,aux,ierr)
!
!  Register an auxiliary variable in the f array.
!
      use General, only: loptest
!
      character (len=*), intent(in) :: varname
      integer, intent(out)  :: ivar
      logical, optional, intent(in) :: communicated, aux
      integer, optional, intent(in) :: vector, array
      integer, optional, intent(out) :: ierr
!
      integer :: vartype
!
      if (loptest(communicated)) then
        vartype = iFARRAY_TYPE_COMM_AUXILIARY
      else
        vartype = iFARRAY_TYPE_AUXILIARY
      endif
!
      call farray_register_variable(varname,ivar,vartype,vector=vector,array=array,aux=aux,ierr=ierr)
!
    endsubroutine farray_register_auxiliary
!***********************************************************************
    subroutine farray_register_variable(varname,ivar,vartype,vector,array,aux,ierr)
!
! 12-may-12/MR: avoid writing of auxiliary variable index into index.pro if
!               variable not written into var.dat
! 16-jan-17/MR: avoid fatal error due to try of re-registering of an already existing variable at reloading.
!
      use General, only: ioptest

      character (len=*), intent(in) :: varname
      integer, target, intent(out)  :: ivar
      integer, intent(in)  :: vartype
      integer, optional, intent(in) :: vector, array
      logical, optional, intent(in) :: aux
      integer, optional, intent(out) :: ierr
!
      type (farray_contents_list), pointer :: item, new
      integer :: ncomponents, narray, nvars
!
      if (vartype==iFARRAY_TYPE_SCRATCH) then
        call fatal_error("farray_register_variable", &
          "Registering "//trim(varname)//" as scratch variable fails. Use the"// &
          "farray_acquire_scratch_area routine.")
      endif
!
      if (present(ierr)) ierr=0
      ncomponents=ioptest(vector,1)
      narray=ioptest(array,1)
      nvars=ncomponents*narray
!
      item => find_by_name(varname)
!
! Already existing variable.
!
      if (associated(item)) then
        if (item%ncomponents/=ncomponents) then
          if (present(ierr)) then
            ierr=iFARRAY_ERR_WRONGSIZE
            ivar=0
            return
          endif
          call fatal_error("farray_register_variable", &
            "Registering "//trim(varname)//" fails: Name already exists but with a different "//&
            "number of components!")
        endif
        if (item%vartype/=vartype) then
          if (present(ierr)) then
            ierr=iFARRAY_ERR_WRONGTYPE
            ivar=0
            return
          endif
          call fatal_error("farray_register_variable", &
                           "Registering "//trim(varname)//" fails: Name already exists but with a different variable type!")
        endif
        if (.not.lreloading) then
          if (present(ierr)) then
            ierr=iFARRAY_ERR_DUPLICATE
            ivar=0
            return
          endif
          call fatal_error("farray_register_variable","Registering "//trim(varname)//" fails: Name already exists!")
        endif
      else
!
! New variable.
!
        select case (vartype)
          case (iFARRAY_TYPE_PDE)
            if (nvar+nvars>mvar) then
              if (present(ierr)) then
                ierr=iFARRAY_ERR_OUTOFSPACE
                ivar=0
                return
              endif
              call fatal_error("farray_register_variable", &
            "Registering "//trim(varname)//" fails: Insufficient mvar variables allocated.  This probably means "//&
            "the MVAR CONTRIBUTION header is incorrect in one of the physics "// &
            "modules or in cparam.local.")
            endif
          case (iFARRAY_TYPE_AUXILIARY)
            if (naux+nvars>maux) then
              if (present(ierr)) then
                ierr=iFARRAY_ERR_OUTOFSPACE
                ivar=0
                return
              endif
              call fatal_error("farray_register_variable", &
            "Registering "//trim(varname)//" fails: Insufficient maux variables allocated.  This either means "//&
            "the MAUX CONTRIBUTION header is incorrect in one of the physics "// &
            "modules. Or that there you are using some code that can, depending "// &
            "on runtime parameters, require extra auxiliary variables.  For the "// &
            "latter try adding an MAUX CONTRIBUTION header to cparam.local.")
            endif
          case (iFARRAY_TYPE_COMM_AUXILIARY)
            if (naux_com+nvars>maux_com) then
              if (present(ierr)) then
                ierr=iFARRAY_ERR_OUTOFSPACE
                ivar=0
                return
              endif
              call fatal_error("farray_register_variable", &
            "Registering "//trim(varname)//" fails: Insufficient maux_com variables allocated.  This either means "//&
            "the COMMUNICATED AUXILIARIES header is incorrect in one of the physics "// &
            "modules. Or that there you are using some code that can, depending "// &
            "on runtime parameters, require extra auxiliary variables.  For the "// &
            "latter try adding an MAUX CONTRIBUTION and COMMUNICATED AUXILIARIES "// &
            "headers to cparam.local.")
            endif
          case (iFARRAY_TYPE_GLOBAL)
            if (nglobal+nvars>mglobal) then
              if (present(ierr)) then
                ierr=iFARRAY_ERR_OUTOFSPACE
                ivar=0
                return
              endif
              call fatal_error("farray_register_variable", &
            "Registering "//trim(varname)//" fails: Insufficient mglobal variables allocated.  This either means "//&
            "the MGLOBAL CONTRIBUTION header is incorrect in one of the physics "// &
            "modules. Or that there you are using some code that can, depending "// &
            "on runtime parameters, require extra global variables.  For the "// &
            "latter try adding an MGLOBAL CONTRIBUTION header to cparam.local.")
            endif
          case default
            call fatal_error("farray_register_variable", &
            "Registering "//trim(varname)//" fails: Invalid vartype set")
        endselect
!
        call new_item_atstart(thelist,new=new)
        new%varname     = varname
        new%vartype     = vartype
        new%ncomponents = ncomponents
        allocate(new%ivar(nvars))
        new%ivar(1)%p => ivar
!
        select case (vartype)
          case (iFARRAY_TYPE_PDE)
            ivar=nvar+1
            nvar=nvar+nvars
          case (iFARRAY_TYPE_COMM_AUXILIARY)
            ivar=mvar+naux_com+1
            naux=naux+nvars
            naux_com=naux_com+nvars
          case (iFARRAY_TYPE_AUXILIARY)
            ivar=mvar+maux_com+(naux-naux_com)+1
            naux=naux+nvars
          case (iFARRAY_TYPE_GLOBAL)
            ivar=mvar+maux+nglobal+1
            nglobal=nglobal+nvars
        endselect
!
        call save_analysis_info(new)
!
!  write varname and index into index.pro file (for idl)
!  except for auxiliary variables which are not written into var.dat
!
        if ( .not.lwrite_aux .and. (vartype==iFARRAY_TYPE_COMM_AUXILIARY .or. &
                                    vartype==iFARRAY_TYPE_AUXILIARY )) return
!
        call farray_index_append('i'//varname,ivar,vector=vector,array=array,aux=aux)
!
      endif
!
    endsubroutine farray_register_variable
!***********************************************************************
    subroutine farray_index_append(varname,ivar,vector,array,aux)
!
! 14-Oct-2018/PAB: coded
!
      use General, only: itoa
!
      character (len=*), intent(in) :: varname
      integer, intent(in) :: ivar
      integer, optional, intent(in) :: vector, array
      logical, optional, intent(in) :: aux
!
      character (len=len(varname)) :: component
      integer :: pos, l
      logical :: use_aux
!
      use_aux=.false.
      if (present(aux)) use_aux=aux
!
      call index_append(trim(varname),ivar,vector=vector,array=array)
      if (.not. present (array) .and. present (vector)) then
        ! expand vectors: iuu => (iux,iuy,iuz), iaa => (iax,iay,iaz), etc.
        component = trim(varname)
        l = len(trim(component))
        if (vector == 3 .and. use_aux.eqv..false.) then
          if (l == 3) then
            ! double endings: iuu, iaa, etc.
            if (component(2:2) == component(3:3)) l = 2
          endif
          call index_append(trim(component(1:l))//'x',ivar)
          call index_append(trim(component(1:l))//'y',ivar+1)
          call index_append(trim(component(1:l))//'z',ivar+2)
        elseif (vector >= 2) then
          do pos=1, vector
            call index_append(trim(component(1:l))//trim(itoa(pos)),ivar+pos-1)
          enddo
        endif
      endif
!
    endsubroutine farray_index_append
!***********************************************************************
    subroutine farray_index_reset()
!
! 14-oct-18/PAB: coded
!
      call index_reset()
!
    endsubroutine farray_index_reset
!***********************************************************************
    subroutine farray_acquire_scratch_area(varname,ivar,vector,ierr)
!
      character (len=*) :: varname
      integer           :: ivar
      integer           :: i
      type (farray_contents_list), pointer :: item, new
      integer           :: ncomponents
      integer, optional :: ierr
      integer, optional :: vector
      integer           :: memstat, consecutive_free, iscratch
!
      intent(in)  :: varname,vector
      intent(out) :: ivar,ierr
!
      ncomponents=1
      if (present(ierr)) ierr=0
      if (present(vector)) ncomponents=vector
!
      item => find_by_name(varname,only_scratch=.true.)
      if (associated(item)) then
        if (item%ncomponents/=ncomponents) then
          if (present(ierr)) then
            ierr=iFARRAY_ERR_WRONGSIZE
            ivar=0
            return
          endif
          print*,"Acquiring f-array scratch area: ",varname
          call fatal_error("farray_acquire_scratch_area", &
            "f array variable name already exists and with a different "//&
            "number of components!")
        endif
        if (present(ierr)) then
          ierr=iFARRAY_ERR_DUPLICATE
          ivar=0
          return
        endif
        print*,"Acquiring f-array scratch area: ",varname
        call fatal_error("acquire_scratch_area","f-array scratch area name already exists!")
      endif
!
      if (nscratch+ncomponents>mscratch) then
        if (present(ierr)) then
          ierr=iFARRAY_ERR_OUTOFSPACE
          ivar=0
          return
        endif
        print*,"Acquiring f-array scratch area: ",varname
        call fatal_error("farray_acquire_scratch_area", &
      "There are insufficient mscratch slots allocated. This means that"// &
      "the total number of 3D scratch arrays used at any one time is more"// &
      "than the number currently allocated.  Try adding "// &
      "an MSCRATCH CONTRIBUTION header to cparam.local or increasing the value.")
      endif
!
      consecutive_free=0
      iscratch=0
!
      if (iscratch+ncomponents>mscratch) then
        if (present(ierr)) then
          ierr=iFARRAY_ERR_OUTOFSPACE
          ivar=0
          return
        endif
        print*,"Acquiring f-array scratch area: ",varname
        call fatal_error("farray_acquire_scratch_area", &
        "There are insufficient consecutive mscratch slots available. "// &
        "The scratch area has become fragmented!")
      endif
!
      do i=1,mscratch
        if (consecutive_free==0) iscratch=i
        if (scratch_used(i)) then
          consecutive_free=0
        else
          consecutive_free=consecutive_free+1
          if (consecutive_free>=ncomponents) exit
        endif
      enddo
!
      call new_item_atstart(thelist,new=new)
      new%varname=varname
      new%vartype=iFARRAY_TYPE_SCRATCH
      new%ncomponents=ncomponents
!
      allocate(new%ivar(ncomponents),stat=memstat)
      allocate(new%ivar(1)%p)
!
      ivar=mvar+maux+mglobal+iscratch
      new%ivar(1)%p=ivar
      if (ncomponents==1) then
        scratch_used(iscratch)=.true.
      else
        scratch_used(iscratch:iscratch+ncomponents-1)=.false.
      endif
!
      nscratch=nscratch+ncomponents
!
    endsubroutine farray_acquire_scratch_area
!***********************************************************************
    subroutine farray_release_scratch_area(varname,ivar,vector,ierr)
!
      character (len=*) :: varname
      integer           :: ivar
      type (farray_contents_list), pointer :: item
      integer :: ncomponents
      integer, optional :: ierr
      integer, optional :: vector
      integer :: iscratch
!
      intent(in)  :: varname,ivar,vector
      intent(out) :: ierr
!
      ncomponents=1
      if (present(ierr)) ierr=0
      if (present(vector)) ncomponents=vector
!
      item => find_by_name(varname,only_scratch=.true.)
      if (.not.associated(item)) then
        if (present(ierr)) then
          ierr=iFARRAY_ERR_NOSUCHVAR
          return
        endif
        print*,"Releasing f-array scratch area: ",varname
        call fatal_error("farray_release_scratch_area", &
                       "scratch area not found")
      endif
!
      if (item%ivar(1)%p/=ivar) then
        if (present(ierr)) then
          ierr=iFARRAY_ERR_INDEXMISMATCH
          return
        endif
        print*,"Releasing f-array scratch area: ",varname
        call fatal_error("farray_release_scratch_area", &
            "scratch area found but not"// &
            "with the index specified")
      endif
!
      deallocate(item%ivar(1)%p)
      deallocate(item%ivar)
      iscratch=ivar-mvar-maux-mglobal

      if (ncomponents==1) then
        scratch_used(iscratch)=.true.
      else
        scratch_used(iscratch:iscratch+ncomponents-1)=.false.
      endif
      nscratch=nscratch-item%ncomponents
!
      call delete_item_from_list(item)
!
    endsubroutine farray_release_scratch_area
!***********************************************************************
    subroutine save_analysis_info(item)
!
! Process an item/slot in the farray and set/write any metadata required
! by internal or external analysis tools eg. IDL, OpenDX...
!
      use Cdata, only: varname
      use General, only: itoa
!
      type (farray_contents_list), pointer :: item
      integer :: i
!
!  Put variable name in array for
!  use by analysis tool output
!
      if (item%ncomponents>1) then
        do i=0,item%ncomponents-1
          varname(item%ivar(1)%p+i) = trim(item%varname)//trim(itoa(i+1))
        enddo
      else
        varname(item%ivar(1)%p) = trim(item%varname)
      endif
!
    endsubroutine save_analysis_info
!***********************************************************************
    subroutine farray_use_pde(varname,ivar,vector,ierr)
!
      character (len=*) :: varname
      integer, pointer  :: ivar
      integer, parameter :: vartype = iFARRAY_TYPE_PDE
      integer, optional :: ierr
      integer, optional :: vector
!
      intent(in)  :: varname,vector
      intent(out) :: ierr
!
      if (present(ierr).and.present(vector)) then
        call farray_use_variable(varname,ivar,vartype,vector=vector,ierr=ierr)
      elseif (present(ierr)) then
        call farray_use_variable(varname,ivar,vartype,ierr=ierr)
      elseif (present(vector)) then
        call farray_use_variable(varname,ivar,vartype,vector=vector)
      else
        call farray_use_variable(varname,ivar,vartype)
      endif
!
    endsubroutine farray_use_pde
!***********************************************************************
    subroutine farray_use_global(varname,ivar,vector,ierr)
!
      character (len=*) :: varname
      integer, pointer  :: ivar
      integer, parameter :: vartype = iFARRAY_TYPE_GLOBAL
      integer, optional :: ierr
      integer, optional :: vector
!
      intent(in)  :: varname,vector
      intent(out) :: ierr
!
      if (present(ierr).and.present(vector)) then
        call farray_use_variable(varname,ivar,vartype,vector=vector,ierr=ierr)
      elseif (present(ierr)) then
        call farray_use_variable(varname,ivar,vartype,ierr=ierr)
      elseif (present(vector)) then
        call farray_use_variable(varname,ivar,vartype,vector=vector)
      else
        call farray_use_variable(varname,ivar,vartype)
      endif
!
    endsubroutine farray_use_global
!***********************************************************************
    subroutine farray_use_auxiliary(varname,ivar,communicated,vector,ierr)
!
      character (len=*) :: varname
      integer, pointer  :: ivar
      integer :: vartype
      integer, optional :: ierr
      integer, optional :: vector
      logical, optional :: communicated
!
      intent(in)  :: varname,communicated,vector
      intent(out) :: ierr
!
      if (present(communicated)) then
        if (communicated) then
          vartype=iFARRAY_TYPE_COMM_AUXILIARY
        else
          vartype=iFARRAY_TYPE_AUXILIARY
        endif
      endif
!
      if (present(ierr).and.present(vector)) then
        call farray_use_variable(varname,ivar,vartype,vector=vector,ierr=ierr)
      elseif (present(ierr)) then
        call farray_use_variable(varname,ivar,vartype,ierr=ierr)
      elseif (present(vector)) then
        call farray_use_variable(varname,ivar,vartype,vector=vector)
      else
        call farray_use_variable(varname,ivar,vartype)
      endif
!
    endsubroutine farray_use_auxiliary
!***********************************************************************
    subroutine farray_use_variable(varname,ivar,vartype,component,vector,ierr)
!
      character (len=*) :: varname
      integer, pointer  :: ivar
      type (farray_contents_list), pointer :: item
      integer, optional :: ierr
      integer, optional :: vector
      integer, optional :: vartype
      integer, optional :: component
!
      intent(in)  :: varname,vartype,component,vector
      intent(out) :: ierr
!
      if (present(ierr)) ierr=0
!
      item => find_by_name(varname)
      if (associated(item)) then
        if (present(vector)) then
          if (item%ncomponents/=vector) then
            if (present(ierr)) then
              ierr=iFARRAY_ERR_WRONGSIZE
              nullify(ivar)
              return
            endif
            print*,"Using f-array variable: ",varname
            call fatal_error("farray_use_variable","f array variable name already exists but with the wrong number of components!")
          endif
        endif
        if (present(vartype)) then
          if (item%vartype/=vartype) then
            if (present(ierr)) then
              ierr=iFARRAY_ERR_WRONGTYPE
              nullify(ivar)
              return
            endif
            print*,"Using f-array variable: ",varname
            call fatal_error("farray_use_variable","f array variable name already exists but with the wrong variable type!")
          endif
        endif
!
        if (present(component)) then
          if ((component<=0).or.(component>item%ncomponents)) then
            if (present(ierr)) then
              ierr=iFARRAY_ERR_NOSUCHCOMPONENT
              nullify(ivar)
              return
            endif
            print*,"Using f-array variable: ",varname
            call fatal_error("farray_use_variable","f-array variable not found!")
          endif
          ivar=>item%ivar(component)%p
        else
          ivar=>item%ivar(1)%p
          return
        endif
      endif
!
      if (present(ierr)) then
        ierr=iFARRAY_ERR_NOSUCHVAR
        nullify(ivar)
        return
      endif
      print*,"Using f-array variable: ",varname
      call fatal_error("farray_use_variable","f-array variable not found!")
!
    endsubroutine farray_use_variable
!***********************************************************************
    function variable_exists(varname,include_scratch)
!
      character (len=*) :: varname
      logical, optional :: include_scratch
      logical :: variable_exists
      type (farray_contents_list), pointer :: item
!
      intent(in) :: varname,include_scratch
!
      item=>thelist
      do while (associated(item))
        if (item%vartype==iFARRAY_TYPE_SCRATCH) then
          if (present(include_scratch)) then
            if (include_scratch) then
              if (item%varname==varname) then
                variable_exists=.true.
                return
              endif
            endif
          endif
        else
          if (item%varname==varname) then
            variable_exists=.true.
            return
          endif
        endif
        item=>item%next
      enddo
      variable_exists=.false.
      return
!
    endfunction variable_exists
!***********************************************************************
    function find_by_name(varname, include_scratch, only_scratch)
!
      character (len=*) :: varname
      logical, optional :: include_scratch
      logical, optional :: only_scratch
      type (farray_contents_list), pointer :: find_by_name
!
      intent(in) :: varname,include_scratch,only_scratch
!
      find_by_name=>thelist
      do while (associated(find_by_name))
!
! If we only want to search the scratch areas
!
        if (present(only_scratch)) then
          if (only_scratch) then
            if (find_by_name%vartype/=iFARRAY_TYPE_SCRATCH) then
              find_by_name=>find_by_name%next
              cycle
            endif
            if (find_by_name%varname==varname) then
              return
            endif
          endif
        endif
!
! Otherwise check everything
!
        if (find_by_name%varname==varname) then
          if (find_by_name%vartype==iFARRAY_TYPE_SCRATCH) then
            if (present(include_scratch)) then
              if (include_scratch) then
                return
              endif
            endif
          else
            return
          endif
        endif
        find_by_name=>find_by_name%next
      enddo
      NULLIFY(find_by_name)
      return
!
    endfunction find_by_name
!***********************************************************************
    function farray_size_by_name(varname)
!
      character (len=*) :: varname
      integer :: farray_size_by_name
      type (farray_contents_list), pointer :: item
!
      intent(in) :: varname
!
      item=>find_by_name(varname)
      if (associated(item)) then
        farray_size_by_name=item%ncomponents
        return
      endif
      farray_size_by_name=0
      return
!
    endfunction farray_size_by_name
!***********************************************************************
    function farray_type_by_name(varname)
!
      character (len=*) :: varname
      integer :: farray_type_by_name
      type (farray_contents_list), pointer :: item
!
      intent(in) :: varname
!
      item=>find_by_name(varname)
      if (associated(item)) then
        farray_type_by_name=item%vartype
        return
      endif
      farray_type_by_name=iFARRAY_TYPE_NOSUCHTYPE
      return
!
    endfunction farray_type_by_name
!!***********************************************************************
!    function farray_index_by_name(varname,component,ierr)
!      character (len=*) :: varname
!      integer, pointer :: farray_find_by_name
!      integer, optional :: component
!      integer, optional :: ierr
!      type (farray_contents_list), pointer :: item
!!
!      intent(in) :: varname,component,ierr
!!
!      item=>thelist
!      do while (associated(item))
!        if (item%varname==varname) then
!          if (present(component)) then
!            if ((component<=0).or.(component>item%ncomponents) then
!            endif
!            farray_find_by_name=>item%ivar(component)
!          else
!            farray_find_by_name=>item%ivar(1)
!          return
!        endif
!        item=>item%next
!      enddo
!      NULLIFY(farray_find_by_name)
!      return
!!
!    endfunction farray_index_by_name
!***********************************************************************
    subroutine free_list(list)
!
!  Free any memory allocated by list
!
!  18-may-2007/wolf: documented and added list%ivar deallocation
!
      type (farray_contents_list), pointer :: list
      type (farray_contents_list), pointer :: next
!
      do while (associated(list))
        next => list%next
        deallocate(list%ivar)
        deallocate(list)
        list => next
      enddo
      nullify(list)
!
    endsubroutine free_list
!***********************************************************************
    subroutine delete_item_from_list(item)
!
      type (farray_contents_list), pointer :: item
!
      item%next%previous => item%previous
      item%previous%next => item%next
      deallocate(item)
      nullify(item)
!
    endsubroutine delete_item_from_list
!***********************************************************************
    subroutine new_item_atstart(list,new)
!
!  Insert new item at beginning of list [like Perl's unshift(), or Lisp's
!  (cons 'item list)]
!
!  18-may-2007/wolf: documented
!
      type (farray_contents_list), pointer :: list
      type (farray_contents_list), optional, pointer :: new
      type (farray_contents_list), pointer :: new_
!
      allocate(new_)
      new_%next => list
      nullify(new_%previous)
      list => new_
      if (present(new)) new => new_
!
    endsubroutine new_item_atstart
!***********************************************************************
    subroutine farray_clean_up()
!
!  Free any memory allocated for farrays, so G95 doesn't warn about still
!  allocated memory.
!
!  18-may-2007/wolf: coded
!
      call free_list(thelist)
!
    endsubroutine farray_clean_up
!***********************************************************************
endmodule FArrayManager
