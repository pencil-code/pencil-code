function check_slices_par, field, stride, readdir, switches
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
if (strtrim(fields[0],2) eq '') then begin
  message, 'No entries in "video.in"!', /warn
  return, 0
endif
;
default, field, fields[0]
;
if ((is_str (field) le 0) or (strtrim(field,2) eq '')) then begin
  message, '"field" must be a non-empty string!', /warn
  return, 0
endif
;
field_base = field
pos = stregex(field,'[1-9][0-9]*[xyz]?$')
if (pos lt 0) then pos = stregex(field,'[xyz]$')
if (pos le 1) then pos = 0
if (pos ge 0) then field_base = strtrim(strmid(field,0,pos),2)
;
;AB: Matthias, please check; without the "or ..." it would not work for us.
;
if (not (any (field_base eq fields) or any (field eq fields))) then begin
  print, 'Field "'+strtrim(field_base,2)+'" is not listed in "video.in"!'
  return, 0
endif
;
; Read slice switches
;
cdum=''
switches=create_struct('xyread',0,'xzread',0,'yzread',0,'xy2read',0,'xy3read',0,'xy4read',0,'xz2read',0,'rread',0)

for iread=0,1 do begin

  on_ioerror, noposition
  openr, lun_1, readdir+'/slice_position.dat', /get_lun
  on_ioerror, cont
  readf, lun_1, cdum & switches.xyread  = strmid(strtrim(cdum,2),0,1) eq 'T' ? 1 : 0
  readf, lun_1, cdum & switches.xy2read = strmid(strtrim(cdum,2),0,1) eq 'T' ? 1 : 0
  readf, lun_1, cdum & switches.xy3read = strmid(strtrim(cdum,2),0,1) eq 'T' ? 1 : 0
  readf, lun_1, cdum & switches.xy4read = strmid(strtrim(cdum,2),0,1) eq 'T' ? 1 : 0
  readf, lun_1, cdum & switches.xzread  = strmid(strtrim(cdum,2),0,1) eq 'T' ? 1 : 0
  readf, lun_1, cdum & switches.xz2read = strmid(strtrim(cdum,2),0,1) eq 'T' ? 1 : 0
  readf, lun_1, cdum & switches.yzread  = strmid(strtrim(cdum,2),0,1) eq 'T' ? 1 : 0
  readf, lun_1, cdum & switches.rread  = strmid(strtrim(cdum,2),0,1) eq 'T' ? 1 : 0
cont:
  close, 1
  free_lun, lun_1
  on_ioerror,  NULL
  return, 1

noposition:
  if (iread eq 0) then begin
    close, 1
    print, 'No slice_position.dat found in "'+readdir+'"!!!'
    if (n_params() eq 1) then stride=1
    spawn, 'echo "- calling read_videofiles."; read_videofiles '+strtrim(field,2)+' '+strtrim(stride,2)+' >& /dev/null'
  endif
endfor

print, 'For backwards compatibility it is assumed that XY, XY2, XZ and YZ slices do exist.'

switches.xyread = 1
switches.xy2read = 1
switches.xy3read = 0
switches.xy4read = 0
switches.xzread = 1
switches.xz2read = 0
switches.yzread = 1
switches.rread = 0

return, 1
;
end
