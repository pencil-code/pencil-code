function hdf5_contains, label, group=group, name=name

	common hdf5_file_info, file_id, file_name

	if (size (file_id, /type) eq 0) then file_id = !Values.D_NaN

	if (finite (file_id, /NaN)) then begin
		print, "ERROR: no HDF5 file is open!"
		stop
		return, !Values.D_NaN
	end

	pos = strpos (label, '/', /reverse_search)
	if (pos lt 0) then begin
		group = '/'
		name = label
	end else begin
		group = strmid (label, 0, pos)
		name = strmid (label, pos+1)
		exists = hdf5_contains (group)
		if (not exists) then return, exists
	end

	list = hdf5_content (group)
	found = total (list eq name)

	return, (found gt 0)
end
