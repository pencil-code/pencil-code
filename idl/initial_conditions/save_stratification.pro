;$Id$
pro save_stratification,z,lnrho,ss,zresolution=zresolution,filename=filename
;
;  seems like a rather special routine to me
;
default,filename,'stratification.dat'

if keyword_set(zresolution) then begin
  zout=congrid(reform(z,n_elements(z)),zresolution,1,1)
  lnrhoout=congrid(reform(lnrho,n_elements(lnrho)),zresolution,1,1)
  ssout=congrid(ss,zresolution,1,1)
end else begin
  zout=z
  lnrhoout=lnrho
  ssout=ss
  zresolution=(size(z))[1]
end

get_lun,lun
openw,lun,filename
for i=0,zresolution-1 do printf,lun,zout[i],lnrhoout[i],ssout[i]
close,lun
free_lun,lun

end
