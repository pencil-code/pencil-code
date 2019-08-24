function h5_contains, label, group=group, name=name

	common h5_file_info, file_id, file_name, group_name, group_content

	if (size (file_id, /type) eq 0) then file_id = !Values.D_NaN

	if (finite (file_id, /NaN)) then begin
		print, "ERROR: no HDF5 file is open!"
		stop
		return, !Values.D_NaN
	end

	if (strmid (label, 0, 1) eq '/') then label = strmid (label, 1)

	if (label eq '/') then return, (label eq '')
	pos = strpos (label, '/', /reverse_search)
	if (pos lt 0) then begin
		group = '/'
		name = label
	end else begin
		group = strmid (label, 0, pos)
		name = strmid (label, pos+1)
help, label, group, name
		exists = h5_contains (group)
		if (not exists) then return, exists
	end

	found = total (h5_content (group) eq name)

	return, (found gt 0)
end
