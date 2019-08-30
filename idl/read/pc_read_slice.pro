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
;
;
function pc_read_slice, quantity, plane, time=time, coordinate=coord, position=pos, datadir=datadir, first=first, last=last, skip=skip

  if (size (datadir, /type) eq 0) then datadir = pc_get_datadir (datadir)

  ; load one HDF5 slice, if datadir ends with '/slices'
  last_slice = pc_read ('last', filename=quantity+'_'+plane+'.h5', datadir=datadir+'/slices')

  default, first, 1
  default, last, last_slice
  default, skip, 0
  step = skip + 1
  num = 1 + (last - first) / step

  data = reform (dblarr ([ h5_get_size ('1/data'), num ]))
  time = reform (dblarr (num))
  coord = reform (dblarr (num))
  pos = reform (lonarr (num))

  ; iterate over slices
  for slice = first, last, step do begin
    index = (slice - first) / step
    group = strtrim (slice, 2)
    data[*,*,index] = pc_read (group+'/data')
    time[index] = pc_read (group+'/time')
    coord[index] = pc_read (group+'/coordinate')
    pos[index] = pc_read (group+'/position')
  end
  h5_close_file

  return, data
end

