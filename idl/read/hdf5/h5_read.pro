function h5_read, label, start=start, count=count, transpose=transpose, filename=filename, close=close

	common h5_file_info, file_id, file_name, group_name, group_content

	if (size (file_id, /type) eq 0) then file_id = !Values.D_NaN

	if (size (filename, /type) eq 7) then h5_open_file, filename

	if (finite (file_id, /NaN)) then begin
		print, "ERROR: no HDF5 file is open!"
		stop
		return, !Values.D_NaN
	end

	dataset = H5D_OPEN (file_id, label)

	num = n_elements (start) < n_elements (count)
	if (num ge 1) then begin
		start = start[0:num-1]
		count = count[0:num-1]
		dataspace = H5D_GET_SPACE (dataset)
		H5S_SELECT_HYPERSLAB, dataspace, start, count, /reset
		mem_space = H5S_CREATE_SIMPLE (count)
		data = H5D_READ (dataset, file_space=dataspace, memory_space=mem_space)
		H5S_CLOSE, dataspace
		H5S_CLOSE, mem_space
	end else begin
		data = H5D_READ (dataset)
	end
	H5D_CLOSE, dataset

	if (keyword_set (transpose) and (size (data, /n_dimensions) ge 2)) then begin
		target_size = reverse (size (data, /dimensions))
		data = reform (transpose (data), target_size)
	end

	if ((size (data, /n_dimensions) eq 1) and (n_elements (data) eq 1)) then begin
		data = data[0]
	end

	if (keyword_set (close)) then h5_close_file

	return, data
end
