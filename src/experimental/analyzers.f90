
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

  integer, parameter :: npossibleanalyzers = 3
  character(len=30), dimension(npossibleanalyzers), parameter :: possibleanalyzers = &
                                       [ character(len=30) :: 'sum', 'average', 'identity']
  
  contains

  subroutine getAnalyzer(analyzername,analyzer,resultlen)
    implicit none
    character(len=30), intent(in) :: analyzername
    procedure(AnalyzerTemplate), pointer , intent(inout):: analyzer
    integer, intent(inout) :: resultlen
    if (trim(analyzername) == 'sum') then
      resultlen = 1
      analyzer => analyzer_sum
    else if (trim(analyzername) == 'average') then
      resultlen = 1
      analyzer => analyzer_average
    else if (trim(analyzername) == 'identity') then
      resultlen = 1
      analyzer => analyzer_identity
    else
      resultlen = -1
      nullify(analyzer)
    end if
  end subroutine getAnalyzer

  function analyzer_sum(dataset, xdim1, xdim2, tlen, resultlen) result(analysis)
    implicit none
    integer,intent(in) :: xdim1,xdim2,tlen,resultlen
    real, dimension(xdim1,xdim2,tlen), intent(in) :: dataset
    real, dimension(xdim1,xdim2,resultlen) :: analysis
    analysis(1:xdim1,1:xdim2,1) = sum(dataset,dim=3)
  end function

  function analyzer_average(dataset, xdim1, xdim2, tlen, resultlen) result(analysis)
    implicit none
    integer,intent(in) :: xdim1,xdim2,tlen,resultlen
    real, dimension(xdim1,xdim2,tlen), intent(in) :: dataset
    real, dimension(xdim1,xdim2,resultlen) :: analysis
    analysis(1:xdim1,1:xdim2,1) = sum(dataset,dim=3)/tlen
  end function

  function analyzer_identity(dataset, xdim1, xdim2, tlen, resultlen) result(analysis)
    implicit none
    integer,intent(in) :: xdim1,xdim2,tlen,resultlen
    real, dimension(xdim1,xdim2,tlen), intent(in) :: dataset
    real, dimension(xdim1,xdim2,resultlen) :: analysis
    write(*,*) 'analyzing...identity', xdim1, xdim2, tlen, resultlen
    analysis(:,:,1:resultlen) = dataset(:,:,1:resultlen)
  end function

end module Analyzers
