function h5_get_type, label, exists=exists

	common h5_file_info, file_id, file_name, group_name, group_content

	if (size (file_id, /type) eq 0) then file_id = !Values.D_NaN

	if (finite (file_id, /NaN)) then begin
		print, "ERROR: no HDF5 file is open!"
		stop
		return, !Values.D_NaN
	end

	if (not keyword_set (exists)) then begin
		if (not h5_contains (label)) then return, 0
		exists = 1
	end

	num = h5_get_size (label)
	dims = n_elements (num)
	if ((dims eq 1) and (num eq 0)) then begin
		idl_type = size (h5_read (label), /type)
	end else begin
		idl_type = size (h5_read (label, start=replicate (0, dims), count=replicate (1, dims)), /type)
	end

	return, idl_type
end
