;
; $Id$
;
; NAME:
;      PC_READ_VIDEO
;+
; PURPOSE:
;     Read video slices and put in structure
;
; PARAMETERS:
;     Field needs to be the name as it appears in the name of the data file
;     to be read: i.e. for scalar quantities file name is slice_<field>.*,
;     but for vector quantities slice_<variable>[1-3].* where variable is
;     uu, aa, bb etc. Then field has to be <variable>[1-3].
;-
; MODIFICATION HISTORY:
;     Written by: Anders Johansen (johansen@mpia.de) on 28.06.2007
;
pro pc_read_video, field=field, object=object, nt=nt, njump=njump, stride=stride, $
    dim=dim, datadir=datadir, proc=proc, swap_endian=swap_endian, help=help, $
    xy2read=xy2read, xyread=xyread, xzread=xzread, yzread=yzread, xz2read=xz2read, $
    rread=rread, print=print, mask=mask, fail=fail, old_format=old_form, single=single
COMPILE_OPT IDL2,HIDDEN
common pc_precision, zero, one, precision, data_type, data_bytes, type_idl
;
  if (keyword_set(help)) then begin
    doc_library, 'pc_read_video'
    return
  endif
;
; Default values.
;
planes = [ 'xy', 'xz', 'yz', 'xy2', 'xy3', 'xy4', 'xz2', 'r' ]
datadir = pc_get_datadir(datadir)
default, proc, -1
default, nt, 100
default, njump, 0
default, stride, 1
default, swap_endian, 0
default, xyread, 1
default, xy2read, 1
default, xy3read, 1
default, xy4read, 1
default, xzread, 1
default, xz2read, 1
default, yzread, 1
default, rread, 1
default, print, 1
default, single, 0
;
; Read dimensions and set precision.
;
if (not is_defined(dim)) then pc_read_dim, obj=dim, datadir=readdir, /quiet
pc_set_precision, datadir=datadir, dim=dim, /quiet
;
; Load HDF5 data, if available
;
  if (not keyword_set (old_format) and file_test (datadir+'/slices', /directory)) then begin
    field=field[0]
    num_planes = n_elements (planes)
    default, mask, replicate (1B, num_planes)
    mask = mask and [ xyread, xzread, yzread, xy2read, xy3read, xy4read, xz2read ]
    mask = mask and file_test (datadir+'/slices/'+field+'_'+planes+'.h5')
    object = { field:field }
    for i = 0, num_planes-1 do begin
      if (mask[i]) then object = create_struct (object, planes[i], pc_read_slice (field, planes[i], datadir=datadir, time=t, coord=coord, pos=pos, single=single))
    endfor
    default, t, single ? 0. : zero
    default, coord, single ? !Values.F_NaN : !Values.D_NaN*zero
    default, pos, -1
    object = create_struct (object, 't', t, 'coordinate', coord, 'position', pos, 'num_planes', num_planes, 'num_snapshots', n_elements (t))
    if (keyword_set (print)) then begin
      ; print summary
      print, 'field="', field, '", datadir="', datadir, '"'
      print, 'min(t)  , max(t)   = ', min(t), ',', max(t)
      for i = 0, num_planes-1 do begin
        if (mask[i]) then print, 'min('+planes[i]+') , max('+planes[i]+')  = ', minmax(object.(i+1))
      endfor
    endif
    return
  endif
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
if not check_slices_par(field, stride, readdir, s) then goto, slices_par_missing
;
; Combine with switch preselection
;
s.xyread = s.xyread and xyread 
s.xy2read = s.xy2read and xy2read 
s.xy3read = s.xy3read and xy3read 
s.xy4read = s.xy4read and xy4read 
s.xzread = s.xzread and xzread 
s.xz2read = s.xz2read and xz2read 
s.yzread = s.yzread and yzread 
if (is_defined(s.rread)) then s.rread = s.rread and rread 
;
; Combine with mask
;
if arg_present(mask) then begin
  if (size(mask) ge 1) then s.xyread = s.xyread and mask[0]
  if (size(mask) ge 2) then s.xzread = s.xzread and mask[1]
  if (size(mask) ge 3) then s.yzread = s.yzread and mask[2]
  if (size(mask) ge 4) then s.xy2read = s.xy2read and mask[3]
  if (size(mask) ge 5) then s.xy3read = s.xy3read and mask[4]
  if (size(mask) ge 6) then s.xy4read = s.xy4read and mask[5]
  if (size(mask) ge 7) then s.xz2read = s.xz2read and mask[6]
  if (is_defined(s.rread) and size(mask) ge 8) then s.rread = s.rread and mask[7]
endif
;
; Define filenames of slices, declare arrays for reading and storing, open slice files,
; call read_videofiles if necessary.
;
present=0
if (s.xyread) then begin
  tag='xy'
  file_slice=readdir+'/slice_'+field+'.xy'
  xy_tmp =make_array(dim.nx,dim.ny, type=type_idl)
  xy = make_array(dim.nx,dim.ny,nt, type= single ? 4 : type_idl)
  if (proc eq -1) and (not file_test(file_slice)) then $
    spawn, 'echo Data missing - calling read_videofiles.; read_videofiles '+strtrim(field,2)+' '+strtrim(stride,2)+' >& /dev/null' 
  if (file_test(file_slice)) then begin 
    openr, lun_1, file_slice, /f77, /get_lun, swap_endian=swap_endian
    present=1
  endif
endif
if (s.xy2read) then begin
  tag='xy2'
  file_slice=readdir+'/slice_'+field+'.xy2'
  xy2_tmp=make_array(dim.nx,dim.ny, type=type_idl)
  xy2= make_array(dim.nx,dim.ny,nt, type= single ? 4 : type_idl)
  if (proc eq -1) and (not file_test(file_slice)) then $
    spawn, 'echo Data missing - calling read_videofiles.; read_videofiles '+strtrim(field,2)+' '+strtrim(stride,2)+' >& /dev/null'
  if (file_test(file_slice)) then begin 
    openr, lun_2, file_slice, /f77, /get_lun, swap_endian=swap_endian
    present=1
  endif
endif
if (s.xy3read) then begin
  tag='xy3'
  file_slice=readdir+'/slice_'+field+'.xy3'
  xy3_tmp=make_array(dim.nx,dim.ny, type=type_idl)
  xy3= make_array(dim.nx,dim.ny,nt, type= single ? 4 : type_idl)
  if (proc eq -1) and (not file_test(file_slice)) then $
    spawn, 'echo Data missing - calling read_videofiles.; read_videofiles '+strtrim(field,2)+' '+strtrim(stride,2)+' >& /dev/null'
  if (file_test(file_slice)) then begin 
    openr, lun_3, file_slice, /f77, /get_lun, swap_endian=swap_endian
    present=1
  endif
endif
if (s.xy4read) then begin
  tag='xy4'
  file_slice=readdir+'/slice_'+field+'.xy4'
  xy4_tmp=make_array(dim.nx,dim.ny, type=type_idl)
  xy4= make_array(dim.nx,dim.ny,nt,type= single ? 4 : type_idl)
  if (proc eq -1) and (not file_test(file_slice)) then $
    spawn, 'echo Data missing - calling read_videofiles.; read_videofiles '+strtrim(field,2)+' '+strtrim(stride,2)+' >& /dev/null'
  if (file_test(file_slice)) then begin 
    openr, lun_4, file_slice, /f77, /get_lun, swap_endian=swap_endian
    present=1
  endif
endif
if (s.xzread) then begin
  tag='xz'
  file_slice=readdir+'/slice_'+field+'.xz'
  xz_tmp =make_array(dim.nx,dim.nz, type=type_idl)
  xz = make_array(dim.nx,dim.nz,nt, type= single ? 4 : type_idl)
  if (proc eq -1) and (not file_test(file_slice)) then $
    spawn, 'echo Data missing - calling read_videofiles.; read_videofiles '+strtrim(field,2)+' '+strtrim(stride,2)+' >& /dev/null'
  if (file_test(file_slice)) then begin 
    openr, lun_5, file_slice, /f77, /get_lun, swap_endian=swap_endian
    present=1
  endif
endif
if (s.xz2read) then begin
  tag='xz2'
  file_slice=readdir+'/slice_'+field+'.xz2'
  xz2_tmp=make_array(dim.nx,dim.nz, type=type_idl)
  xz2=make_array(dim.nx,dim.nz,nt, type= single ? 4 : type_idl)
  if (proc eq -1) and (not file_test(file_slice)) then $
    spawn, 'echo Data missing - calling read_videofiles.; read_videofiles '+strtrim(field,2)+' '+strtrim(stride,2)+' >& /dev/null'
  if (file_test(file_slice)) then begin 
    openr, lun_6, file_slice, /f77, /get_lun, swap_endian=swap_endian
    present=1
  endif
endif
if (s.yzread) then begin
  tag='yz'
  file_slice=readdir+'/slice_'+field+'.yz'
  yz_tmp =make_array(dim.ny,dim.nz, type=type_idl)
  yz =make_array(dim.ny,dim.nz,nt, type= single ? 4 : type_idl)
  if (proc eq -1) and (not file_test(file_slice)) then $
    spawn, 'echo Data missing - calling read_videofiles.; read_videofiles '+strtrim(field,2)+' '+strtrim(stride,2)+' >& /dev/null'
  if (file_test(file_slice)) then begin 
    openr, lun_7, file_slice, /f77, /get_lun, swap_endian=swap_endian
    present=1
  endif
endif
if (s.rread) then begin
  tag='r'
  file_slice=readdir+'/slice_'+field+'.r'
  pc_read_param, obj=param, /run_, datadir=readdir, /quiet
  r_tmp =make_array(param.nth_rslice,param.nph_rslice, type=type_idl)
  r =make_array(param.nth_rslice,param.nph_rslice,nt, type= single ? 4 : type_idl)
  if (proc eq -1) and (not file_test(file_slice)) then $
    spawn, 'echo Data missing - calling read_videofiles.; read_videofiles '+strtrim(field,2)+' '+strtrim(stride,2)+' >& /dev/null'
  if (file_test(file_slice)) then begin
    openr, lun_8, file_slice, /f77, /get_lun, swap_endian=swap_endian
    present=1
  endif
endif
if not present then goto, data_file_missing
;
t_tmp=one & t=make_array(nt, type= single ? 4 : type_idl)
;
; Possible to skip njump slices.
;
if (njump gt 0) then begin
  dummy1=zero
  for ijump=0,njump-1 do begin
    if is_defined(lun_1) then begin 
      tag='xy'
      readu, lun_1, dummy1
      if (eof(lun_1)) then goto, jumpfail
    endif
    if is_defined(lun_2) then begin 
      tag='xy2'
      readu, lun_2, dummy1
      if (eof(lun_2)) then goto, jumpfail
    endif
    if is_defined(lun_3) then begin 
      tag='xy3'
      readu, lun_3, dummy1
      if (eof(lun_3)) then goto, jumpfail
    endif
    if is_defined(lun_4) then begin 
      tag='xy4'
      readu, lun_4, dummy1
      if (eof(lun_4)) then goto, jumpfail
    endif
    if is_defined(lun_5) then begin 
      tag='xz'
      readu, lun_5, dummy1
      if (eof(lun_5)) then goto, jumpfail
    endif
    if is_defined(lun_6) then begin 
      tag='xz2'
      readu, lun_6, dummy1
      if (eof(lun_6)) then goto, jumpfail
    endif
    if is_defined(lun_7)  then begin 
      tag='yz'
      readu, lun_7, dummy1
      if (eof(lun_7)) then goto, jumpfail
    endif
    if is_defined(lun_8)  then begin
      tag='r'
      readu, lun_8, dummy1
      if (eof(lun_8)) then goto, jumpfail
    endif
  endfor
endif
;
; Read slices at nt times.
; Stride is not (again) applied if data are already collected from all processors.
;
if proc eq -1 then stride=1

for it=0,nt-1,stride do begin
;
; Stop if end of file reached.
;
  if is_defined(lun_1) then begin
    if (eof(lun_1)) then break
    readu, lun_1, xy_tmp, t_tmp, slice_zpos
    xy [*,*,it]=xy_tmp
  endif
  if is_defined(lun_2) then begin
    if (eof(lun_2)) then break
    readu, lun_2, xy2_tmp, t_tmp, slice_z2pos 
    xy2[*,*,it]=xy2_tmp
  endif
  if is_defined(lun_3) then begin
    if (eof(lun_3)) then break
    readu, lun_3, xy3_tmp, t_tmp, slice_z3pos 
    xy3[*,*,it]=xy3_tmp
  endif
  if is_defined(lun_4) then begin
    if (eof(lun_4)) then break
    readu, lun_4, xy4_tmp, t_tmp, slice_z4pos
    xy4[*,*,it]=xy4_tmp
  endif
  if is_defined(lun_5) then begin
    if (eof(lun_5)) then break
    readu, lun_5, xz_tmp, t_tmp, slice_ypos 
    xz [*,*,it]=xz_tmp
  endif
  if is_defined(lun_6) then begin
    if (eof(lun_6)) then break
    readu, lun_6, xz2_tmp, t_tmp, slice_y2pos 
    xz2[*,*,it]=xz2_tmp
  endif
  if is_defined(lun_7) then begin
    if (eof(lun_7)) then break
    readu, lun_7, yz_tmp, t_tmp, slice_xpos 
    yz [*,*,it]=yz_tmp
  endif
  if is_defined(lun_8) then begin
    if (eof(lun_8)) then break
    readu, lun_8, r_tmp, t_tmp, slice_rpos
    r [*,*,it]=r_tmp
  endif
  t[it]=t_tmp
endfor
;
; Close files.
;
if is_defined(lun_1) then free_lun, lun_1
if is_defined(lun_2) then free_lun, lun_2
if is_defined(lun_3) then free_lun, lun_3
if is_defined(lun_4) then free_lun, lun_4
if is_defined(lun_5) then free_lun, lun_5
if is_defined(lun_6) then free_lun, lun_6
if is_defined(lun_7) then free_lun, lun_7
if is_defined(lun_8) then free_lun, lun_8
;
; Truncate the data if less than nt slices were read.
;
if it lt nt then begin
  it -= 1
  if s.xyread  then xy = xy[*,*,0:it]
  if s.xy2read then xy2 = xy2[*,*,0:it]
  if s.xy3read then xy3 = xy3[*,*,0:it]
  if s.xy4read then xy4 = xy4[*,*,0:it]
  if s.xzread  then xz = xz[*,*,0:it]
  if s.xz2read then xz2 = xz2[*,*,0:it]
  if s.yzread  then yz = yz[*,*,0:it]
  if s.rread  then r = r[*,*,0:it]
  t = t[0:it]
endif
;
; Build structure of all the variables.
;
object = create_struct(name=objectname,'t',t)
if s.xyread  then object=create_struct(object,'xy',xy)
if s.xzread  then object=create_struct(object,'xz',xz)
if s.yzread  then object=create_struct(object,'yz',yz)
if s.xy2read then object=create_struct(object,'xy2',xy2)
if s.xy3read then object=create_struct(object,'xy3',xy3)
if s.xy4read then object=create_struct(object,'xy4',xy4)
if s.xz2read then object=create_struct(object,'xz2',xz2)
if s.rread then object=create_struct(object,'r',r)
;
; If requested print a summary.
;
if keyword_set(print) then begin
  print, 'field="', field, '", datadir="', readdir, '"'
  print, 'min(t)  , max(t)   = ', min(t), ',', max(t)
  if s.xyread  then print, 'min(xy) , max(xy)  = ', min(xy), ', ', max(xy)
  if s.xzread  then print, 'min(xz) , max(xz)  = ', min(xz), ', ', max(xz)
  if s.yzread  then print, 'min(yz) , max(yz)  = ', min(yz), ', ', max(yz)
  if s.xy2read then print, 'min(xy2), max(xy2) = ', min(xy2), ', ', max(xy2)
  if s.xy3read then print, 'min(xy3), max(xy3) = ', min(xy3), ', ', max(xy3)
  if s.xy4read then print, 'min(xy4), max(xy4) = ', min(xy4), ', ', max(xy4)
  if s.xz2read then print, 'min(xz2), max(xz2) = ', min(xz2), ', ', max(xz2)
  if s.rread   then print, 'min(r),   max(r)   = ', min(r), ', ', max(r)
endif

mask = [s.xyread,s.xzread,s.yzread,s.xy2read,s.xy3read,s.xy4read,s.xz2read,s.rread]
if arg_present(fail) then fail=0
return
;
jumpfail:
print, 'Not enough data for requested jump over ',strtrim(string(njump),2),' first times'
print, 'in slice_"'+field+'.'+tag+'"!!!' 
if arg_present(fail) then fail=1
return
;
data_file_missing:
print, 'No data file "slice_'+field+'.'+tag+'" found in "'+readdir+'"!!!'
if arg_present(fail) then fail=1
return
;
slices_par_missing:
print, 'No slices parameters found!!!'
if arg_present(fail) then fail=1

end
