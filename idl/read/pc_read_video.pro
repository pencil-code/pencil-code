;;
;; $Id$
;;
;; NAME:
;;      PC_READ_VIDEO
;;
;; PURPOSE:
;;     Read video slices and put in structure
;;
;; PARAMETERS:
;;     Field needs to be the name as it appears in the name of the data file
;;     to be read: i.e. for scalar quantities file name is slice_<field>.*,
;;     but for vector quantities slice_<variable>[1-3].* where variable is
;;     uu, aa, bb etc. Then field has to be <variable>[1-3].
;;
;; MODIFICATION HISTORY:
;;     Written by: Anders Johansen (johansen@mpia.de) on 28.06.2007
;;
pro pc_read_video, field=field, object=object, nt=nt, njump=njump, $
    dim=dim, datadir=datadir, proc=proc, swap_endian=swap_endian, $
    xy2read=xy2read, xyread=xyread, xzread=xzread, yzread=yzread, $
    xz2read=xz2read, print=print, mask=mask, fail=fail
COMPILE_OPT IDL2,HIDDEN
common pc_precision, zero, one, precision, data_type, data_bytes, type_idl
;
; Default values.
;
datadir = pc_get_datadir(datadir)
default, proc, -1
default, nt, 100
default, njump, 0
default, swap_endian, 0
default, xyread, 0
default, xy2read, 0
default, xy3read, 0
default, xy4read, 0
default, xzread, 0
default, xz2read, 0
default, yzread, 0
default, print, 1

fields=rstringlist('video.in')
if fields[0] eq '' then return


if is_defined(field) then begin

  if field ne '' then begin

    pos=stregex(field,'[1-9][0-9]*')
    field_base = pos ge 0 ? strtrim(strmid(field,0,pos)) : field

    if (where(field_base eq fields))[0] eq -1 then begin
      print, 'Field "'+strtrim(field_base,2)+'" not in video.in!!!'
      return
    endif
 
  endif else $
    field=fields[0]

endif else $
  field=fields[0]
;
; Read either from global data directory or from single processor directory.
;
if (proc eq -1) then begin
  readdir=datadir
endif else begin
  readdir=datadir+'/proc'+strtrim(proc,2)
  if not file_test(readdir,/directory) then begin
    print, '"'+readdir+'" does not exist or is not a directory!!!'
    return
  endif
endelse
;
; Read slice switches
;
cdum=''
on_ioerror, noposition
openr, lun_1, readdir+'/slice_position.dat', /get_lun
on_ioerror, NULL
readf, lun_1, cdum & xyread  = strmid(strtrim(cdum,2),0,1) eq 'T' ? 1 : 0
readf, lun_1, cdum & xy2read = strmid(strtrim(cdum,2),0,1) eq 'T' ? 1 : 0
readf, lun_1, cdum & xy3read = strmid(strtrim(cdum,2),0,1) eq 'T' ? 1 : 0
readf, lun_1, cdum & xy4read = strmid(strtrim(cdum,2),0,1) eq 'T' ? 1 : 0
readf, lun_1, cdum & xzread  = strmid(strtrim(cdum,2),0,1) eq 'T' ? 1 : 0
readf, lun_1, cdum & xz2read = strmid(strtrim(cdum,2),0,1) eq 'T' ? 1 : 0
readf, lun_1, cdum & yzread  = strmid(strtrim(cdum,2),0,1) eq 'T' ? 1 : 0
close, 1
free_lun, lun_1
;
; Combine with mask
;
if arg_present(mask) then begin
  xyread = xyread and mask[0]
  xy2read = xy2read and mask[1]
  xy3read = xy3read and mask[2]
  xy4read = xy4read and mask[3]
  xzread = xzread and mask[4]
  xzread2 = xz2read and mask[5]
  xyread = xyread and mask[6]
endif

;
; Read dimensions and set precision.
;
pc_read_dim, obj=dim, datadir=readdir, /quiet
;
; Define filenames of slices, declare arrays for reading and storing, open slice files.
;
on_ioerror, data_file_missing
if (xyread) then begin
  tag='xy'
  file_slice1=readdir+'/slice_'+field+'.xy'
  xy_tmp =fltarr(dim.nx,dim.ny)*one
  xy =fltarr(dim.nx,dim.ny,nt)*one
  openr, lun_1, file_slice1, /f77, /get_lun, swap_endian=swap_endian
endif
if (xy2read) then begin
  tag='xy2'
  file_slice2=readdir+'/slice_'+field+'.xy2'
  xy2_tmp=fltarr(dim.nx,dim.ny)*one
  xy2=fltarr(dim.nx,dim.ny,nt)*one
  openr, lun_2, file_slice2, /f77, /get_lun, swap_endian=swap_endian
endif
if (xy3read) then begin
  tag='xy3'
  file_slice3=readdir+'/slice_'+field+'.xy3'
  xy3_tmp=fltarr(dim.nx,dim.ny)*one
  xy3=fltarr(dim.nx,dim.ny,nt)*one
  openr, lun_3, file_slice3, /f77, /get_lun, swap_endian=swap_endian
endif
if (xy4read) then begin
  tag='xy4'
  file_slice4=readdir+'/slice_'+field+'.xy4'
  xy4_tmp=fltarr(dim.nx,dim.ny)*one
  xy4=fltarr(dim.nx,dim.ny,nt)*one
  openr, lun_4, file_slice4, /f77, /get_lun, swap_endian=swap_endian
endif
if (xzread) then begin
  tag='xz'
  file_slice5=readdir+'/slice_'+field+'.xz'
  xz_tmp =fltarr(dim.nx,dim.nz)*one
  xz =fltarr(dim.nx,dim.nz,nt)*one
  openr, lun_5, file_slice5, /f77, /get_lun, swap_endian=swap_endian
endif
if (xz2read) then begin
  tag='xz2'
  file_slice6=readdir+'/slice_'+field+'.xz2'
  xz2_tmp=fltarr(dim.nx,dim.nz)*one
  xz2=fltarr(dim.nx,dim.nz,nt)*one
  openr, lun_6, file_slice6, /f77, /get_lun, swap_endian=swap_endian
endif
if (yzread) then begin
  tag='yz'
  file_slice7=readdir+'/slice_'+field+'.yz'
  yz_tmp =fltarr(dim.ny,dim.nz)*one
  yz =fltarr(dim.ny,dim.nz,nt)*one
  openr, lun_7, file_slice7, /f77, /get_lun, swap_endian=swap_endian
endif
on_ioerror, NULL
;
t_tmp=one & t=fltarr(nt)*one
;
; Possible to skip njump slices.
;
if (njump gt 0) then begin
  dummy=zero
  for ijump=0,njump-1 do begin
    if (xyread)  then begin 
      tag='xy'
      readu, lun_1, dummy1
      if (eof(lun_1)) then goto, jumpfail
    endif
    if (xy2read) then begin 
      tag='xy2'
      readu, lun_2, dummy1
      if (eof(lun_2)) then goto, jumpfail
    endif
    if (xy3read) then begin 
      tag='xy3'
      readu, lun_3, dummy1
      if (eof(lun_3)) then goto, jumpfail
    endif
    if (xy4read) then begin 
      tag='xy4'
      readu, lun_4, dummy1
      if (eof(lun_4)) then goto, jumpfail
    endif
    if (xzread)  then begin 
      tag='xz'
      readu, lun_5, dummy1
      if (eof(lun_5)) then goto, jumpfail
    endif
    if (xz2read) then begin 
      tag='xz2'
      readu, lun_6, dummy1
      if (eof(lun_6)) then goto, jumpfail
    endif
    if (yzread)  then begin 
      tag='yz'
      readu, lun_7, dummy1
      if (eof(lun_7)) then goto, jumpfail
    endif
  endfor
endif
;
; Read slices at nt times.
;
for it=0,nt-1 do begin
;
; Stop if end of file reached.
;
  if (xyread)  then begin 
    if (eof(lun_1)) then break
    readu, lun_1, xy_tmp, t_tmp, slice_zpos
    xy [*,*,it]=xy_tmp
  endif
  if (xy2read) then begin
    if (eof(lun_2)) then break
    readu, lun_2, xy2_tmp, t_tmp, slice_z2pos 
    xy2[*,*,it]=xy2_tmp
  endif
  if (xy3read) then begin
    if (eof(lun_3)) then break
    readu, lun_3, xy3_tmp, t_tmp, slice_z3pos 
    xy3[*,*,it]=xy3_tmp
  endif
  if (xy4read) then begin
    if (eof(lun_4)) then break
    readu, lun_4, xy4_tmp, t_tmp, slice_z4pos
    xy4[*,*,it]=xy4_tmp
  endif
  if (xzread)  then begin
    if (eof(lun_5)) then break
    readu, lun_5, xz_tmp, t_tmp, slice_ypos 
    xz [*,*,it]=xz_tmp
  endif
  if (xz2read) then begin
    if (eof(lun_6)) then break
    readu, lun_6, xz2_tmp, t_tmp, slice_y2pos 
    xz2[*,*,it]=xz2_tmp
  endif
  if (yzread)  then begin
    if (eof(lun_7)) then break
    readu, lun_7, yz_tmp, t_tmp, slice_xpos 
    yz [*,*,it]=yz_tmp
  endif
  t[it]=t_tmp
endfor
;
; Close files.
;
if (xyread)  then begin
  close, lun_1
  free_lun, lun_1
endif
if (xy2read) then begin
  close, lun_2
  free_lun, lun_2
endif
if (xy3read) then begin
  close, lun_3
  free_lun, lun_3
endif
if (xy4read) then begin
  close, lun_4
  free_lun, lun_4
endif
if (xzread)  then begin
  close, lun_5
  free_lun, lun_5
endif
if (xz2read) then begin
  close, lun_6
  free_lun, lun_6
endif
if (yzread)  then begin
  close, lun_7
  free_lun, lun_7
endif
;
; Truncate the data if less than nt slices were read.
;
if it lt nt then begin
  it -= 1
  if xyread  then xy = xy[*,*,0:it]
  if xy2read then xy2 = xy2[*,*,0:it]
  if xy3read then xy3 = xy3[*,*,0:it]
  if xy4read then xy4 = xy4[*,*,0:it]
  if xzread  then xz = xz[*,*,0:it]
  if xz2read then xz2 = xz2[*,*,0:it]
  if yzread  then yz = yz[*,*,0:it]
  t = t[0:it]
endif
;
; Build structure of all the variables.
;
object = create_struct(name=objectname,'t',t)
if xyread  then object=create_struct(object,'xy',xy)
if xy2read then object=create_struct(object,'xy2',xy2)
if xy3read then object=create_struct(object,'xy3',xy3)
if xy4read then object=create_struct(object,'xy4',xy4)
if xzread  then object=create_struct(object,'xz',xz)
if xz2read then object=create_struct(object,'xz2',xz2)
if yzread  then object=create_struct(object,'yz',yz)
;
; If requested print a summary.
;
if keyword_set(print) then begin
  print, 'field="', field, '", datadir="', readdir, '"'
  print, 'min(t)  , max(t)   = ', min(t), ',', max(t)
  if xyread  then print, 'min(xy) , max(xy)  = ', min(xy), ', ', max(xy)
  if xy2read then print, 'min(xy2), max(xy2) = ', min(xy2), ', ', max(xy2)
  if xy3read then print, 'min(xy3), max(xy3) = ', min(xy3), ', ', max(xy3)
  if xy4read then print, 'min(xy4), max(xy4) = ', min(xy4), ', ', max(xy4)
  if xzread  then print, 'min(xz) , max(xz)  = ', min(xz), ', ', max(xz)
  if xz2read then print, 'min(xz2), max(xz2) = ', min(xz2), ', ', max(xz2)
  if yzread  then print, 'min(yz) , max(yz)  = ', min(yz), ', ', max(yz)
endif

mask = [xyread,xy2read,xy3read,xy4read,xzread,xz2read,yzread]
if keyword_set(fail) then fail=0
return
;
jumpfail:
print, 'Not enough data for requested jump over ',strtrim(string(njump),2),' first times'
print, 'in slice_"'+field+'.'+tag+'"!!!' 
if keyword_set(fail) then fail=1
return
;
noposition:
print, 'No slice_position.dat found in "'+readdir+'"!!!'
if keyword_set(fail) then fail=1
return
;
data_file_missing:
print, 'No data file slice_"'+field+'.'+tag+'" found in "'+readdir+'"!!!'
if keyword_set(fail) then fail=1
;
end
