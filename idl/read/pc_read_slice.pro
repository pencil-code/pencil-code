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
function pc_read_slice, quantity, plane, time=time, coordinate=coord, position=pos, datadir=datadir

  if (size (datadir, /type) eq 0) then datadir = pc_get_datadir (datadir)

  ; load one HDF5 slice, if datadir ends with '/slices'
  last = pc_read ('last', filename=quantity+'_'+plane+'.h5', datadir=datadir+'/slices')
  data = reform (dblarr ([ h5_get_size ('1/data'), last ]))
  time = reform (dblarr (last))
  coord = reform (dblarr (last))
  pos = reform (lonarr (last))

  ; iterate over slices
  for slice = 1, last do begin
    group = strtrim (slice, 2)
    data[*,*,slice-1] = pc_read (group+'/data')
    time[slice-1] = pc_read (group+'/time')
    coord[slice-1] = pc_read (group+'/coordinate')
    pos[slice-1] = pc_read (group+'/position')
  end
  h5_close_file

  return, data
end

