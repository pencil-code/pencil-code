;;
;; $Id$
;;
;; NAME:
;;      PC_READ_VIDEO
;;
;; PURPOSE:
;;     Read video slices and put in structure
;;
;; MODIFICATION HISTORY:
;;     Written by: Anders Johansen (johansen@mpia.de) on 28.06.2007
;;
pro pc_read_video, field=field, object=object, nt=nt, njump=njump, $
    dim=dim, datadir=datadir, proc=proc, swap_endian=swap_endian, $
    xy2read=xy2read, xyread=xyread, xzread=xzread, yzread=yzread, $
    xz2read=xz2read, print=print, quiet=quiet, help=help
COMPILE_OPT IDL2,HIDDEN
common pc_precision, zero, one, precision, data_type, data_bytes, type_idl
;
; Default values.
;
default, field, 'lnrho'
datadir = pc_get_datadir(datadir)
default, proc, -1
default, nt, 100
default, njump, 0
default, swap_endian, 0
default, xy2read, 1
default, xyread, 1
default, xzread, 1
default, yzread, 1
default, xz2read, 0
default, print, 1
default, quiet, 0
default, help, 0
;
; Read either from global data directory or from single processor directory.
;
if (proc eq -1) then begin
  readdir=datadir
endif else begin
  readdir=datadir+'/proc'+strtrim(proc,2)
endelse
;
; Read dimensions and set precision.
;
pc_read_dim, obj=dim, datadir=readdir, /quiet
;
; Define filenames of slices.
;
if (xy2read) then file_slice1=readdir+'/slice_'+field+'.xy2'
if (xyread) then  file_slice2=readdir+'/slice_'+field+'.xy'
if (xzread) then  file_slice3=readdir+'/slice_'+field+'.xz'
if (yzread) then  file_slice4=readdir+'/slice_'+field+'.yz'
if (xz2read) then file_slice5=readdir+'/slice_'+field+'.xz2'
;
; Read data with proper precision.
;
if (xy2read) then xy2_tmp=fltarr(dim.nx,dim.ny)*one
if (xyread) then xy_tmp =fltarr(dim.nx,dim.ny)*one
if (xzread) then xz_tmp =fltarr(dim.nx,dim.nz)*one
if (yzread) then yz_tmp =fltarr(dim.ny,dim.nz)*one
if (xz2read) then xz2_tmp=fltarr(dim.nx,dim.nz)*one
t_tmp  =one
;
; Put data in arrays.
;
xy2=0
xy=0
xz=0
yz=0
xz2=0

if (xy2read) then xy2=fltarr(dim.nx,dim.ny,nt)*one
if (xyread) then xy =fltarr(dim.nx,dim.ny,nt)*one
if (xzread) then xz =fltarr(dim.nx,dim.nz,nt)*one
if (xzread) then yz =fltarr(dim.ny,dim.nz,nt)*one
if (xz2read) then xz2=fltarr(dim.nx,dim.nz,nt)*one
t  =fltarr(nt)*one
;
; Open slice files.
;
if (xy2read) then begin
  openr, lun_1, file_slice1, /f77, /get_lun, swap_endian=swap_endian
endif
if (xyread) then begin
  openr, lun_2, file_slice2, /f77, /get_lun, swap_endian=swap_endian
endif
if (xzread) then begin
  openr, lun_3, file_slice3, /f77, /get_lun, swap_endian=swap_endian
endif
if (yzread) then begin
  openr, lun_4, file_slice4, /f77, /get_lun, swap_endian=swap_endian
endif
if (xz2read) then begin
  openr, lun_5, file_slice5, /f77, /get_lun, swap_endian=swap_endian
endif
;
; Read slices at nt times.
;
for it=0,nt-1 do begin
;
; Possible to skip njump slices.
;
  if (njump gt 0) then begin
    dummy=zero
    for ijump=0,njump-1 do begin
      if (xy2read) then readu, lun_1, dummy1
      if (xyread)  then readu, lun_2, dummy1
      if (xzread)  then readu, lun_3, dummy1
      if (yzread)  then readu, lun_4, dummy1
      if (xz2read) then readu, lun_5, dummy1
      if (xy2read) then begin
        if (eof(lun_1)) then break
      endif
      if (xyread) then begin
        if (eof(lun_2)) then break
      endif
      if (xzread) then begin
        if (eof(lun_3)) then break
      endif
      if (yzread) then begin
        if (eof(lun_4)) then break
      endif
      if (xz2read) then begin
        if (eof(lun_5)) then break
      endif
    endfor
  endif
;
; Stop if end of file reached.
;
  if (xy2read) then begin
    if (eof(lun_1)) then break
  endif
  if (xyread) then begin
    if (eof(lun_2)) then break
  endif
  if (xzread) then begin
    if (eof(lun_3)) then break
  endif
  if (yzread) then begin
    if (eof(lun_4)) then break
  endif
  if (xz2read) then begin
    if (eof(lun_5)) then break
  endif
;
  if (xy2read) then readu, lun_1, xy2_tmp, t_tmp, slice_z2pos
  if (xyread)  then readu, lun_2, xy_tmp, t_tmp, slice_zpos
  if (xzread)  then readu, lun_3, xz_tmp, t_tmp, slice_ypos
  if (yzread)  then readu, lun_4, yz_tmp, t_tmp, slice_xpos
  if (xz2read) then readu, lun_5, xz2_tmp, t_tmp, slice_y2pos
  if (xy2read) then xy2[*,*,it]=xy2_tmp
  if (xyread) then xy [*,*,it]=xy_tmp
  if (xzread) then xz [*,*,it]=xz_tmp
  if (yzread) then yz [*,*,it]=yz_tmp
  if (xz2read) then xz2[*,*,it]=xz2_tmp
  t[it]=t_tmp
endfor
;
; Close files.
;
if (xy2read) then begin
  close, lun_1
  free_lun, lun_1
endif
if (xyread)  then begin
  close, lun_2
  free_lun, lun_2
endif
if (xzread)  then begin
  close, lun_3
  free_lun, lun_3
endif
if (yzread)  then begin
  close, lun_4
  free_lun, lun_4
endif
if (xz2read) then begin
  close, lun_5
  free_lun, lun_5
endif
;
; Truncate the data if less than nt slices were read.
;
if it lt nt then begin
  it -= 1
if (xyread) then  xy = xy[*,*,0:it]
if (xzread) then  xz = xz[*,*,0:it]
if (yzread) then  yz = yz[*,*,0:it]
if (xy2read) then  xy2 = xy2[*,*,0:it]
if (xz2read) then  xz2 = xz2[*,*,0:it]
  t = t[0:it]
endif
;
; Build structure of all the variables.
;
object = create_struct(name=objectname,['t','xy','xy2','xz','yz','xz2'], $
    t, xy, xy2, xz, yz,xz2)
;
; If requested print a summary.
;
if (keyword_set(print)) then begin
  print, 'field=''', field, ''', datadir=''', readdir, ''''
  print, 'min(t)  , max(t)   = ', min(t), ',', max(t)
  if (xyread)  then print, 'min(xy) , max(xy)  = ', min(xy), ', ', max(xy)
  if (xy2read) then print, 'min(xy2), max(xy2) = ', min(xy2), ', ', max(xy2)
  if (xzread)  then print, 'min(xz) , max(xz)  = ', min(xz), ', ', max(xz)
  if (yzread)  then print, 'min(yz) , max(yz)  = ', min(yz), ', ', max(yz)
  if (xz2read) then print, 'min(xz2), max(xz2) = ', min(xz2), ', ', max(xz2)
endif
;
end
