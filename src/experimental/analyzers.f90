
module Analyzers
  implicit none

  abstract interface
    function AnalyzerTemplate(dataset,xdim1,xdim2,tlen,resultlen) result(analysis)
      implicit none
      integer,intent(in) :: xdim1,xdim2,tlen,resultlen
      real, dimension(xdim1,xdim2,tlen), intent(in) :: dataset
      real, dimension(xdim1,xdim2,resultlen) :: analysis
    end function AnalyzerTemplate
  end interface

  integer, parameter :: npossibleanalyzers = 2
  character(len=30), dimension(npossibleanalyzers), parameter :: possibleanalyzers = &
                                       [ character(len=30) :: 'sumt', 'sumt2' ]
  
  contains

  subroutine getAnalyzer(analyzername,analyzer,resultlen)
    implicit none
    character(len=30), intent(in) :: analyzername
    procedure(AnalyzerTemplate), pointer , intent(inout):: analyzer
    integer, intent(inout) :: resultlen
    if (trim(analyzername) == 'sumt') then
      resultlen = 1
      analyzer => sumt
    else if (trim(analyzername) == 'sumt2') then
      resultlen = 2
      analyzer => sumt2
    else
      resultlen = 1
      analyzer => identity
    end if
  end subroutine getAnalyzer

  function sumt(dataset, xdim1, xdim2, tlen, resultlen) result(analysis)
    implicit none
    integer,intent(in) :: xdim1,xdim2,tlen,resultlen
    real, dimension(xdim1,xdim2,tlen), intent(in) :: dataset
    real, dimension(xdim1,xdim2,resultlen) :: analysis
    !write(*,*) 'analyzing...sumt', xdim1, xdim2, tlen, resultlen
    analysis(1:xdim1,1:xdim2,1) = sum(dataset,dim=3)
  end function

  function sumt2(dataset, xdim1, xdim2, tlen, resultlen) result(analysis)
    implicit none
    integer,intent(in) :: xdim1,xdim2,tlen,resultlen
    real, dimension(xdim1,xdim2,tlen), intent(in) :: dataset
    real, dimension(xdim1,xdim2,resultlen) :: analysis
    write(*,*) 'analyzing...sumt2', xdim1, xdim2, tlen, resultlen
    analysis(1:xdim1,1:xdim2,1) = sum(dataset,dim=3)
    analysis(1:xdim1,1:xdim2,2) = sum(dataset,dim=3)
  end function

  function identity(dataset, xdim1, xdim2, tlen, resultlen) result(analysis)
    implicit none
    integer,intent(in) :: xdim1,xdim2,tlen,resultlen
    real, dimension(xdim1,xdim2,tlen), intent(in) :: dataset
    real, dimension(xdim1,xdim2,resultlen) :: analysis
    write(*,*) 'analyzing...identity', xdim1, xdim2, tlen, resultlen
    analysis(:,:,1:resultlen) = dataset(:,:,1:resultlen)
  end function

end module Analyzers
