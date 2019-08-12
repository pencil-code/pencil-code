pro h5_create_group, group

	common h5_file_info, file_id, file_name, group_name, group_content

	exists = h5_contains (group, group=parent, name=child)
	if (exists) then return

	if (parent ne '/') then begin
		if (not h5_contains (parent)) then h5_create_group, parent
	end

	result = h5g_create (file_id, group)
end

