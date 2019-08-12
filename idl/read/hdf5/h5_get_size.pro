function h5_get_size, label, transpose=transpose

	common h5_file_info, file_id, file_name, group_name, group_content

	if (size (file_id, /type) eq 0) then file_id = !Values.D_NaN

	if (finite (file_id, /NaN)) then begin
		print, "ERROR: no HDF5 file is open!"
		stop
		return, !Values.D_NaN
	end

	dataset = H5D_OPEN (file_id, label)

	dataspace = H5D_GET_SPACE (dataset)
	H5S_SELECT_ALL, dataspace
	dims = H5S_GET_SIMPLE_EXTENT_DIMS (dataspace)

	if (keyword_set (transpose) and (n_elements (dims) ge 2)) then dims = reverse (dims)

	H5S_CLOSE, dataspace
	H5D_CLOSE, dataset

	return, dims
end
