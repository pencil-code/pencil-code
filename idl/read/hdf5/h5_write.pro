pro h5_write, label, data, transpose=transpose, filename=filename, truncate=truncate, close=close

	common h5_file_info, file_id, file_name, group_name, group_content

	if (size (filename, /type) eq 7) then h5_open_file, filename, /write, truncate=truncate

	exists = h5_contains (label, group=group, name=name)

	if (exists) then begin
		print, "WARNING: can not overwrite existing '"+label+"' in '"+file_name+"'!"
		return
	end

	target_size = size (data, /dimensions)
	if (keyword_set (transpose) and (size (data, /n_dimensions) ge 2)) then target_size = reverse (target_size)

	datatype_id = h5t_idl_create (data)
	if (size (data, /n_dimensions) eq 0) then begin
		dataspace_id = h5s_create_scalar ()
	end else begin
		dataspace_id = h5s_create_simple (target_size)
	end
	group_id = h5g_open (file_id, group)
	dataset_id = h5d_create (group_id, name, datatype_id, dataspace_id)
	if (keyword_set (transpose) and (size (data, /n_dimensions) ge 2)) then begin
		h5d_write, dataset_id, reform (transpose (data), target_size)
	end else begin
		h5d_write, dataset_id, data
	end
	h5d_close, dataset_id
	h5g_close, group_id
	h5s_close, dataspace_id
	h5t_close, datatype_id

	if (keyword_set (close)) then h5_close_file
end

