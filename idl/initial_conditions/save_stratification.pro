pro,save_stratification,z,lnrho,ss,zresolution=zresolution,filename=filename

default,filename,'stratification.dat'

if keyword_set(zresolution) then begin
  zout=z
  lnrhoout=lnrho
  ssout=ss
  zresolution=(size(z))[1]
end else begin
  zout=rebin(zout,zresolution)
  lnrhoout=rebin(lnrhoout,zresolution)
  ssout=rebin(ssout,zresolution)
end

get_lun,lun
openw,lun,filename
for i=1,zresolution do printf,lun,zout[i],lnrhoout[i],ssout[i]
close,lun
free_lun,lun
