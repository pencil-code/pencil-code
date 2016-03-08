
module Analyzers
  !use EMD, only: analyzer_emd
  implicit none

  abstract interface
    function AnalyzerTemplate(dataset,xdim1,xdim2,xdim3,tlen,resultlen) result(analysis)
      implicit none
      integer,intent(in) :: xdim1,xdim2,xdim3,tlen,resultlen
      real, dimension(xdim1,xdim2,xdim3,tlen), intent(in) :: dataset
      real, dimension(xdim1,xdim2,xdim3,resultlen) :: analysis
    end function AnalyzerTemplate
  end interface

  integer, parameter :: npossibleanalyzers = 2
  character(len=30), dimension(npossibleanalyzers), parameter :: possibleanalyzers = &
                                       [ character(len=30) :: 'mean', 'data']
!                                       [ character(len=30) :: 'sum', 'average', 'identity', 'emd']
  
  contains

  subroutine getAnalyzer(analyzername,analyzer,resultlen)
    implicit none
    character(len=30), intent(in) :: analyzername
    procedure(AnalyzerTemplate), pointer , intent(inout):: analyzer
    integer, intent(inout) :: resultlen
    if (trim(analyzername) == 'mean') then
      resultlen  = 1
      analyzer => analyzer_mean
    else if (trim(analyzername) == 'data') then
      resultlen  = resultlen
      analyzer => analyzer_data
    !else if (trim(analyzername) == 'emd') then
    !  resultlen  = resultlen
    !  analyzer => analyzer_emd
    else
      resultlen  = -1
      nullify(analyzer)
    end if
  end subroutine getAnalyzer

  function analyzer_mean(dataset, xdim1, xdim2, xdim3, tlen, resultlen) result(analysis)
    implicit none
    integer,intent(in) :: xdim1,xdim2,xdim3,tlen,resultlen
    real, dimension(xdim1,xdim2,xdim3,tlen), intent(in) :: dataset
    real, dimension(xdim1,xdim2,xdim3,resultlen) :: analysis
    analysis(1:xdim1,1:xdim2,1:xdim3,1) = sum(dataset,dim=4)/tlen
  end function

  function analyzer_data(dataset, xdim1, xdim2, xdim3, tlen, resultlen) result(analysis)
    implicit none
    integer,intent(in) :: xdim1,xdim2,xdim3,tlen,resultlen
    real, dimension(xdim1,xdim2,xdim3,tlen), intent(in) :: dataset
    real, dimension(xdim1,xdim2,xdim3,resultlen) :: analysis
    analysis(1:xdim1,1:xdim2,1:xdim3,1:resultlen) = dataset(1:xdim1,1:xdim2,1:xdim3,1:tlen)
  end function

end module Analyzers
