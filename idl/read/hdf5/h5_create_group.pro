pro hdf5_create_group, group

	common hdf5_file_info, file_id, file_name

	exists = hdf5_contains (group, group=parent, name=child)
	if (exists) then return

	if (parent ne '/') then begin
		if (not hdf5_contains (parent)) then hdf5_create_group, parent
	end

	result = h5g_create (file_id, group)
end

