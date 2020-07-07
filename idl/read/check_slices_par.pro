function check_slices_par, field, readdir, switches
;
; Checks whether <field> is in video.in. 
; Reads slice switches from <readdir>/slice_position.dat into structure <switches>.
; Returns 1 on success, else 0.
;
;  1-jun-18/MR: copy+pasted from 
; 17-mar-20/MR: added xyz line
;  7-jul-20/PAB: added regex hook, added check for existence of orig field
;
fields=rstringlist('video.in')
if (fields[0] eq '') then begin
  print, 'No entries in video.in.'
  return, 0
endif
;
if (is_defined(field)) then begin

  if (field ne '') then begin

    pos=stregex(field,'[1-9][0-9]*[xyz]?$')
    if (pos lt 0) then pos=stregex(field,'[xyz]$')

    field_base = field
    if (pos ge 0) then field_base = strtrim(strmid(field,0,pos),2)
    if (strlen(field_base) lt 2) then field_base=field
;
;AB: Matthias, please check; without the "and ..." it would not work for us.
;
    if (((where(field_base eq fields))[0] eq -1) and ((where(field eq fields))[0] eq -1)) then begin
      print, 'Field "'+strtrim(field_base,2)+'" not in video.in!!!'
      return, 0
    endif

  endif else $
    field=fields[0]

endif else $
  field=fields[0]
;
; Read slice switches
;
cdum=''
switches=create_struct('xyread',0,'xzread',0,'yzread',0,'xy2read',0,'xy3read',0,'xy4read',0,'xz2read',0)
on_ioerror, noposition
openr, lun_1, readdir+'/slice_position.dat', /get_lun
on_ioerror, NULL
readf, lun_1, cdum & switches.xyread  = strmid(strtrim(cdum,2),0,1) eq 'T' ? 1 : 0
readf, lun_1, cdum & switches.xy2read = strmid(strtrim(cdum,2),0,1) eq 'T' ? 1 : 0
readf, lun_1, cdum & switches.xy3read = strmid(strtrim(cdum,2),0,1) eq 'T' ? 1 : 0
readf, lun_1, cdum & switches.xy4read = strmid(strtrim(cdum,2),0,1) eq 'T' ? 1 : 0
readf, lun_1, cdum & switches.xzread  = strmid(strtrim(cdum,2),0,1) eq 'T' ? 1 : 0
readf, lun_1, cdum & switches.xz2read = strmid(strtrim(cdum,2),0,1) eq 'T' ? 1 : 0
readf, lun_1, cdum & switches.yzread  = strmid(strtrim(cdum,2),0,1) eq 'T' ? 1 : 0
close, 1
free_lun, lun_1
return, 1
;
noposition:
print, 'No slice_position.dat found in "'+readdir+'"!!!'
print, 'For backwards compatibility it is assumed that XY, XY2, XZ and YZ slices do exist.'

switches.xyread = 1
switches.xy2read = 1
switches.xy3read = 0
switches.xy4read = 0
switches.xzread = 1
switches.xz2read = 0
switches.yzread = 1

return, 1
;
end
