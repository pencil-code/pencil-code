! $Id: farray.f90,v 1.2 2006-08-22 12:08:32 mee Exp $ 
!
!  This module allocates and manages indices in the f-array
!  in a controlled way.  THis includes handling different 
!  types of variable that may be stored in the f-array like
!  PDE variables, auxilliaries, communicated auxilliaries
!  and global variables.  
!
!  The code follows the similar principles to the SharedVariables
!  module but only handles integer indices.
!
module FArrayManager 
!
  use Cparam, only: mvar,maux,mglobal,maux_com
  use Cdata, only: nvar,naux,naux_com
  use Messages
!
  implicit none
!
  private
!
  public :: farray_register_variable
  public :: farray_register_pde
  public :: farray_register_auxilliary
  public :: farray_register_global
  public :: farray_use_variable
  public :: farray_use_pde
  public :: farray_use_auxilliary
  public :: farray_use_global
  public :: farray_size_by_name
  public :: farray_type_by_name
!
! Types of variable that may be stored in the f-array
!
  integer, public, parameter :: iFARRAY_TYPE_NOSUCHTYPE=0
  integer, public, parameter :: iFARRAY_TYPE_PDE=1
  integer, public, parameter :: iFARRAY_TYPE_COMM_AUXILLIARY=2
  integer, public, parameter :: iFARRAY_TYPE_AUXILLIARY=3
  integer, public, parameter :: iFARRAY_TYPE_GLOBAL=4
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

  type pp
    integer, pointer :: p
  end type pp

!
! Store pointers to variables in a general linked
! list structure.  Shame we can't have (void *) pointers.
!
  type farray_contents_list
!
! Shared variable metadata
!
    character (len=30) :: varname
    integer            :: ncomponents
    integer            :: vartype
    type(pp), dimension(:), pointer :: ivar
!
! Linked list link to next list element 
!
    type (farray_contents_list), pointer :: next
  endtype farray_contents_list
!
! The head of the list (initially empty)
!
  type (farray_contents_list), pointer :: thelist
 
  contains

!***********************************************************************
    subroutine farray_register_pde(varname,ivar,vector,ierr) 
      character (len=*) :: varname
      integer, target  :: ivar
      type (farray_contents_list), pointer :: item
      integer :: ncomponents
      integer, parameter :: vartype = iFARRAY_TYPE_PDE
      integer, optional :: ierr
      integer, optional :: vector

      if (present(ierr).and.present(vector)) then
        call farray_register_variable(varname,ivar,vartype,vector=vector,ierr=ierr)
      elseif (present(ierr)) then
        call farray_register_variable(varname,ivar,vartype,ierr=ierr)
      elseif (present(vector)) then
        call farray_register_variable(varname,ivar,vartype,vector=vector)
      else
        call farray_register_variable(varname,ivar,vartype)
      endif
 
    endsubroutine farray_register_pde
!***********************************************************************
    subroutine farray_register_global(varname,ivar,vector,ierr) 
      character (len=*) :: varname
      integer, target  :: ivar
      type (farray_contents_list), pointer :: item
      integer :: ncomponents
      integer, parameter :: vartype = iFARRAY_TYPE_GLOBAL
      integer, optional :: ierr
      integer, optional :: vector

      if (present(ierr).and.present(vector)) then
        call farray_register_variable(varname,ivar,vartype,vector=vector,ierr=ierr)
      elseif (present(ierr)) then
        call farray_register_variable(varname,ivar,vartype,ierr=ierr)
      elseif (present(vector)) then
        call farray_register_variable(varname,ivar,vartype,vector=vector)
      else
        call farray_register_variable(varname,ivar,vartype)
      endif
 
    endsubroutine farray_register_global
!***********************************************************************
    subroutine farray_register_auxilliary(varname,ivar,communicated,vector,ierr) 
      character (len=*) :: varname
      integer, target  :: ivar
      type (farray_contents_list), pointer :: item
      integer :: ncomponents
      integer :: vartype
      integer, optional :: ierr
      integer, optional :: vector
      logical, optional :: communicated

      if (present(communicated)) then
        if (communicated) then 
          vartype=iFARRAY_TYPE_COMM_AUXILLIARY
        else 
          vartype=iFARRAY_TYPE_AUXILLIARY
        endif
      endif

      if (present(ierr).and.present(vector)) then
        call farray_register_variable(varname,ivar,vartype,vector=vector,ierr=ierr)
      elseif (present(ierr)) then
        call farray_register_variable(varname,ivar,vartype,ierr=ierr)
      elseif (present(vector)) then
        call farray_register_variable(varname,ivar,vartype,vector=vector)
      else
        call farray_register_variable(varname,ivar,vartype)
      endif
 
    endsubroutine farray_register_auxilliary
!***********************************************************************
    subroutine farray_register_variable(varname,ivar,vartype,vector,ierr) 
      character (len=*) :: varname
      integer, target   :: ivar
      integer           :: vartype, i
      type (farray_contents_list), pointer :: item, new
      integer :: ncomponents
      integer, optional :: ierr
      integer, optional :: vector
      integer :: memstat
!
      
      ncomponents=1
      if (present(ierr)) ierr=0
      if (present(vector)) ncomponents=vector
!
      item => find_by_name(varname)
      if (associated(item)) then
        if (item%ncomponents/=ncomponents) then
          if (present(ierr)) then
            ierr=iFARRAY_ERR_WRONGSIZE
            ivar=0
            return
          endif
          print*,"Setting f-array variable: ",varname
          call fatal_error("farray_register_variable","f array variable name already exists but with a different number of components!")
        endif
        if (item%vartype/=vartype) then
          if (present(ierr)) then
            ierr=iFARRAY_ERR_WRONGTYPE
            ivar=0
            return
          endif
          print*,"Setting f-array variable: ",varname
          call fatal_error("farray_register_variable","f array variable name already exists but with a variable type!")
        endif
        if (present(ierr)) then
          ierr=iFARRAY_ERR_DUPLICATE
          ivar=0
          return
        endif
        print*,"Setting shared variable: ",varname
        call fatal_error("farray_register_variable","Shared variable name already exists!")
      endif
!
      call new_item_atstart(thelist,new=new)
      new%varname=varname
      new%vartype=vartype
      new%ncomponents=ncomponents
      allocate(new%ivar(ncomponents),stat=memstat)
      new%ivar(1)%p=>ivar

      select case (vartype)
        case (iFARRAY_TYPE_PDE)
          ivar=nvar+1
          nvar=nvar+ncomponents
        case (iFARRAY_TYPE_COMM_AUXILLIARY)
          ivar=mvar+naux_com+1
          naux=naux+ncomponents
          naux_com=naux_com+ncomponents
        case (iFARRAY_TYPE_AUXILLIARY)
          ivar=mvar+naux+1
          naux=naux+ncomponents
!        case (iFARRAY_TYPE_GLOBAL)
!          ivar=mvar+maux+nglobal
!          nglobal=nglobal+ncomponents
      endselect
!
! Keep a list of component indices
!
 !     do i=2,ncomponents
 !       new%ivar(i)%p=new%ivar(i-1)%p+1
 !     enddo
!
      call save_analysis_info(new)
!
!      ivar=>new%ivar(1)
!    
    endsubroutine farray_register_variable
!***********************************************************************
    subroutine save_analysis_info(item) 
!
      use Cdata, only: varname
      use General, only: chn
!
      type (farray_contents_list), pointer :: item
      integer :: i
      character (len=5) :: istr

      do i=1,item%ncomponents
        call chn(i,istr)
!
!  Put variable name in array
!
        if (item%ncomponents>1) then      
          varname(item%ivar(1)%p+i-1) = item%varname//trim(istr)
        else
          varname(item%ivar(1)%p) = item%varname
        endif
      enddo
    endsubroutine save_analysis_info
!***********************************************************************
    subroutine farray_use_pde(varname,ivar,vector,ierr) 
      character (len=*) :: varname
      integer, pointer  :: ivar
      type (farray_contents_list), pointer :: item
      integer :: ncomponents
      integer, parameter :: vartype = iFARRAY_TYPE_PDE
      integer, optional :: ierr
      integer, optional :: vector

      if (present(ierr).and.present(vector)) then
        call farray_use_variable(varname,ivar,vartype,vector=vector,ierr=ierr)
      elseif (present(ierr)) then
        call farray_use_variable(varname,ivar,vartype,ierr=ierr)
      elseif (present(vector)) then
        call farray_use_variable(varname,ivar,vartype,vector=vector)
      else
        call farray_use_variable(varname,ivar,vartype)
      endif
 
    endsubroutine farray_use_pde
!***********************************************************************
    subroutine farray_use_global(varname,ivar,vector,ierr) 
      character (len=*) :: varname
      integer, pointer  :: ivar
      type (farray_contents_list), pointer :: item
      integer :: ncomponents
      integer, parameter :: vartype = iFARRAY_TYPE_GLOBAL
      integer, optional :: ierr
      integer, optional :: vector

      if (present(ierr).and.present(vector)) then
        call farray_use_variable(varname,ivar,vartype,vector=vector,ierr=ierr)
      elseif (present(ierr)) then
        call farray_use_variable(varname,ivar,vartype,ierr=ierr)
      elseif (present(vector)) then
        call farray_use_variable(varname,ivar,vartype,vector=vector)
      else
        call farray_use_variable(varname,ivar,vartype)
      endif
 
    endsubroutine farray_use_global
!***********************************************************************
    subroutine farray_use_auxilliary(varname,ivar,communicated,vector,ierr) 
      character (len=*) :: varname
      integer, pointer  :: ivar
      type (farray_contents_list), pointer :: item
      integer :: ncomponents
      integer :: vartype
      integer, optional :: ierr
      integer, optional :: vector
      logical, optional :: communicated

      if (present(communicated)) then
        if (communicated) then
          vartype=iFARRAY_TYPE_COMM_AUXILLIARY
        else
          vartype=iFARRAY_TYPE_AUXILLIARY
        endif
      endif

      if (present(ierr).and.present(vector)) then
        call farray_use_variable(varname,ivar,vartype,vector=vector,ierr=ierr)
      elseif (present(ierr)) then
        call farray_use_variable(varname,ivar,vartype,ierr=ierr)
      elseif (present(vector)) then
        call farray_use_variable(varname,ivar,vartype,vector=vector)
      else
        call farray_use_variable(varname,ivar,vartype)
      endif
 
    endsubroutine farray_use_auxilliary
!***********************************************************************
    subroutine farray_use_variable(varname,ivar,vartype,component,vector,ierr) 
      character (len=*) :: varname
      integer, pointer  :: ivar
      integer           :: icomp
      type (farray_contents_list), pointer :: item
      integer :: ncomponents
      integer, optional :: ierr
      integer, optional :: vector
      integer, optional :: vartype
      integer, optional :: component
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

        if (present(component)) then
          if ((component<=0).or.(component>item%ncomponents)) then
            if (present(ierr)) then
              ierr=iFARRAY_ERR_NOSUCHCOMPONENT
              nullify(ivar)
              return
            endif
            print*,"Using f-array variable: ",varname
            call fatal_error("farray_use_variable","F-array variable not found!")
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
      print*,"Using shared variable: ",varname
      call fatal_error("farray_use_variable","F-array variable not found!")
!    
    endsubroutine farray_use_variable
!***********************************************************************
    function variable_exists(varname) 
      character (len=*) :: varname
      logical :: variable_exists
      type (farray_contents_list), pointer :: item
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
    function find_by_name(varname) 
      character (len=*) :: varname
      type (farray_contents_list), pointer :: find_by_name
!
      find_by_name=>thelist
      do while (associated(find_by_name))
        if (find_by_name%varname==varname) then
          return
        endif
        find_by_name=>find_by_name%next
      enddo
      NULLIFY(find_by_name)
      return
!    
    endfunction find_by_name
!***********************************************************************
    function farray_size_by_name(varname) 
      character (len=*) :: varname
      integer :: farray_size_by_name
      type (farray_contents_list), pointer :: item
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
      character (len=*) :: varname
      integer, pointer :: farray_type_by_name
      type (farray_contents_list), pointer :: item
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
      type (farray_contents_list), pointer :: list
      type (farray_contents_list), pointer :: next
 
      do while (associated(list))
        next=>list%next
        deallocate(list)
        list=>next
      enddo
      nullify(list)
    endsubroutine free_list
!***********************************************************************
    subroutine new_item_atstart(list,new) 
      type (farray_contents_list), pointer :: list
      type (farray_contents_list), optional, pointer :: new
      type (farray_contents_list), pointer :: new_

      allocate(new_)
      new_%next=>list 
      list=>new_
      if (present(new)) new=>new_ 
    endsubroutine new_item_atstart
!***********************************************************************
endmodule FArrayManager
