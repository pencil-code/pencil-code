;
; $Id$
;
;   Read a slice from a HDF5 file
;
;  Author: Philippe Bourdin
;  $Date: 2019-04-30 17:16:52 $
;  $Revision: 1.0 $
;
;  04-Apr-2019/PABourdin: coded for Zafira
;  09-Oct-2019/MR: removed premature reforms
;
function pc_read_slice, quantity, plane, time=time, coordinate=coord, position=pos, datadir=datadir, first=first, last=last, skip=skip, single=single

  ; load one HDF5 slice, defined by quantity and plane, from directory datadir/'slices'

  datadir = pc_get_datadir (datadir)
  default, first, 1
  default, last, pc_read ('last', filename=quantity+'_'+plane+'.h5', datadir=datadir+'/slices')
  default, skip, 0
  default, single, 0

  if first gt last then first=last
  if skip lt 0 then skip=0

  step = skip + 1
  num = 1 + (last - first) / step

  if (keyword_set(single)) then begin
    data = fltarr ([h5_get_size ('1/data'), num]) 
    if arg_present(time) then time = fltarr (num)
  endif else begin
    data = dblarr ([h5_get_size ('1/data'), num])
    if arg_present(time) then time = dblarr (num)
  endelse

  if arg_present(coord) then coord = time 
  if arg_present(pos) then pos = lonarr (num)

  ; iterate over slices
  cut_upto=-1
  for slice = first, last, step do begin

    catch, error_status
    ;This statement begins the error handler:
    if error_status ne 0 then begin
      print, 'Error: slice '+strtrim (slice, 2)+' not found!!! Reading terminated.'
      cut_upto=index-1
      break
    endif

    index = (slice - first) / step
    group = strtrim (slice, 2)
    data[*,*,index] = pc_read (group+'/data',single=single)
    if arg_present(time) then time[index] = pc_read (group+'/time')
    if arg_present(coord) then coord[index] = pc_read (group+'/coordinate')
    if arg_present(pos) then pos[index] = pc_read (group+'/position')
  end
  catch, /cancel
  h5_close_file

  if cut_upto ge 0 then begin
    data=data[*,*,0:cut_upto]
    if arg_present(time) then time=time[0:cut_upto]
    if arg_present(coord) then coord=coord[0:cut_upto]
    if arg_present(pos) then pos=pos[0:cut_upto]
  endif

  return, reform(data,/overwrite)
end

