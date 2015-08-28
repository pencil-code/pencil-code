
module Analyzers
  use EMD, only: analyzer_emd
  implicit none

  abstract interface
    function AnalyzerTemplate(dataset,xdim1,xdim2,tlen,resultlen1,resultlen2) result(analysis)
      implicit none
      integer,intent(in) :: xdim1,xdim2,tlen,resultlen1,resultlen2
      real(kind=8), dimension(xdim1,xdim2,tlen), intent(in) :: dataset
      real(kind=8), dimension(xdim1,xdim2,resultlen1,resultlen2) :: analysis
    end function AnalyzerTemplate
  end interface

  integer, parameter :: npossibleanalyzers = 4
  character(len=30), dimension(npossibleanalyzers), parameter :: possibleanalyzers = &
                                       [ character(len=30) :: 'sum', 'average', 'identity', 'emd']
  
  contains

  subroutine getAnalyzer(analyzername,analyzer,resultlen1,resultlen2)
    implicit none
    character(len=30), intent(in) :: analyzername
    procedure(AnalyzerTemplate), pointer , intent(inout):: analyzer
    integer, intent(inout) :: resultlen1,resultlen2
    if (trim(analyzername) == 'sum') then
      resultlen1  = 1
      resultlen2  = 1
      analyzer => analyzer_sum
    else if (trim(analyzername) == 'average') then
      resultlen1  = 1
      resultlen2  = 1
      analyzer => analyzer_average
    else if (trim(analyzername) == 'identity') then
      resultlen1  = resultlen1
      resultlen2  = 1
      analyzer => analyzer_identity
    else if (trim(analyzername) == 'emd') then
      resultlen1  = resultlen1
      resultlen2  = 10
      analyzer => analyzer_emd
    else
      resultlen1  = -1
      resultlen2  = -1
      nullify(analyzer)
    end if
  end subroutine getAnalyzer

  function analyzer_sum(dataset, xdim1, xdim2, tlen, resultlen1, resultlen2) result(analysis)
    implicit none
    integer,intent(in) :: xdim1,xdim2,tlen,resultlen1,resultlen2
    real(kind=8), dimension(xdim1,xdim2,tlen), intent(in) :: dataset
    real(kind=8), dimension(xdim1,xdim2,resultlen1,resultlen2) :: analysis
    analysis(1:xdim1,1:xdim2,1,1) = sum(dataset,dim=3)
  end function

  function analyzer_average(dataset, xdim1, xdim2, tlen, resultlen1, resultlen2) result(analysis)
    implicit none
    integer,intent(in) :: xdim1,xdim2,tlen,resultlen1,resultlen2
    real(kind=8), dimension(xdim1,xdim2,tlen), intent(in) :: dataset
    real(kind=8), dimension(xdim1,xdim2,resultlen1,resultlen2) :: analysis
    analysis(1:xdim1,1:xdim2,1,1) = sum(dataset,dim=3)/tlen
  end function

  function analyzer_identity(dataset, xdim1, xdim2, tlen, resultlen1, resultlen2) result(analysis)
    implicit none
    integer,intent(in) :: xdim1,xdim2,tlen,resultlen1,resultlen2
    real(kind=8), dimension(xdim1,xdim2,tlen), intent(in) :: dataset
    real(kind=8), dimension(xdim1,xdim2,resultlen1,resultlen2) :: analysis
    write(*,*) 'analyzing...identity', xdim1, xdim2, tlen, resultlen1, resultlen2
    analysis(:,:,1:tlen,1) = dataset(:,:,1:tlen)
  end function

end module Analyzers
