pro h5_open_file, filename, write=write, truncate=truncate

	common h5_file_info, file_id, file_name, group_name, group_content
	common pc_read_common, file

	if (size (file_id, /type) eq 0) then file_id = !Values.D_NaN

	if (not finite (file_id, /NaN)) then begin
		H5F_CLOSE, file_id
		file_id = !Values.D_NaN
	end

	if (file_test (filename, /regular) and keyword_set (truncate)) then begin
		file_delete, filename
	end

	if (not file_test (filename, /regular) or keyword_set (truncate)) then begin
		if (keyword_set (write) or keyword_set (truncate)) then begin
			file_id = H5F_CREATE (filename)
		end else begin
			print, "ERROR: file '"+filename+"' does not exist!"
			stop
			return
		end
	end else begin
		if (not H5F_IS_HDF5 (filename)) then begin
			print, "ERROR: '"+filename+"' is not a HDF5 file!"
			stop
			return
		end
		file_id = H5F_OPEN (filename, write=write)
	end

	file_name = filename
	file = filename
end
