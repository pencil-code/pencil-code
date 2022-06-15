function read_profile, filename, z, datadir=datadir

if (size (z, /type) eq 0) then message, 'ERROR: "z" is undefined.'
if (size (datadir, /type) ne 7) then path = 'data/'
if (size (datadir, /type) eq 7) then if (datadir eq '') then path = 'data/'
if (size (path, /type) ne 7) then path = datadir
path += '/../driver/'

len_double = 8

; determine the number of data points in the profile
n_data = ((file_info(path+filename)).size - 2*2*4) / (len_double * 2)

data   = dblarr(n_data)
data_z = dblarr(n_data)


; read profile
OPENU, lun, path+filename, /F77, /GET_LUN
READU, lun, data
READU, lun, data_z
CLOSE, lun
FREE_LUN, lun

; interpolate logarthmic data to Pencil grid profile
profile_out = interpolate_profile(data, data_z, n_data, z)

return,profile_out 



end
