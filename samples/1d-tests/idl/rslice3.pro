; $Id$
;
;  show two different slices
;
print,'for 1-d fields in the x-direction'
print,'type ".r start" first"
;
file_slice1='data/proc0/slice_ux.xz'
file_slice2='data/proc0/slice_lnrho.xz'
file_slice3='data/proc0/slice_ss.xz'
;
t=zero
xz_slice1=one*fltarr(nx,nz)
xz_slice2=one*fltarr(nx,nz)
xz_slice3=one*fltarr(nx,nz)
;
close,1 & openr,1,file_slice1,/f77
close,2 & openr,2,file_slice2,/f77
close,3 & openr,3,file_slice3,/f77
print,'files opened ok'
;
;  set default waiting time w
;  real all snapshots and plot
;
n=0
ps=0
default,w,.2
!p.charsize=2
!p.multi=[0,1,3]
!x.margin=[15,2]
xxx=x[l1:l2]
;
Ma=0.01 & yru=.011 & yrlr=yru & yrs=.00005
Ma=0.05 & yru=.051 & yrlr=yru & yrs=.00005
Ma=0.1 & yru=.1 & yrlr=yru & yrs=.005
Ma=1.0 & yru=1. & yrlr=yru & yrs=1.
;
while not eof(1) do begin
  readu,1,xz_slice1,t
  readu,2,xz_slice2,t
  readu,3,xz_slice3,t
  print,t
  ;plot_io,exp(xz_slice1(*,n));,yr=[-1,1]*1e+1
  plot,xxx,xz_slice1(*,n),ytit='u',yst=3,yr=[-.1,1]*yru,ps=ps
  plot,xxx,xz_slice2(*,n),ytit='lnrho',yst=3,yr=[-.1,2.4]*yrlr,ps=ps
  plot,xxx,xz_slice3(*,n),ytit='s',yst=3,yr=[-.7,.4]*yrs,ps=ps
  wait,w
end
;
!p.multi=0
close,1
close,2
close,3
END
