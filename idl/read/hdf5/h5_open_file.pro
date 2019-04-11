pro h5_open_file, file, write=write, truncate=truncate

	common h5_file_info, file_id, file_name

	if (size (file_id, /type) eq 0) then file_id = !Values.D_NaN

	if (not finite (file_id, /NaN)) then begin
		H5F_CLOSE, file_id
		file_id = !Values.D_NaN
	end

	if (file_test (file, /regular) and keyword_set (truncate)) then begin
		file_delete, file
	end

	if (not file_test (file, /regular) or keyword_set (truncate)) then begin
		if (keyword_set (write) or keyword_set (truncate)) then begin
			file_id = H5F_CREATE (file)
		end else begin
			print, "ERROR: file '"+file+"' does not exist!"
			stop
			return
		end
	end else begin
		if (not H5F_IS_HDF5 (file)) then begin
			print, "ERROR: '"+file+"' is not a HDF5 file!"
			stop
			return
		end

		file_name = file
		file_id = H5F_OPEN (file, write=write)
	end
end
