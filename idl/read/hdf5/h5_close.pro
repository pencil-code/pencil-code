pro hdf5_close

	common hdf5_file_info, file_id, file_name

	if (size (file_id, /type) eq 0) then file_id = !Values.D_NaN

	if (finite (file_id, /NaN)) then begin
		print, "WARNING: no opened HDF5 file!"
		return
	end

	H5F_CLOSE, file_id
	file_id = !Values.D_NaN
	file_name = ''
end
