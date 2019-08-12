function h5_content, group, number=number, maximum=maximum

	common h5_file_info, file_id, file_name, group_name, group_content

	if (size (file_id, /type) eq 0) then file_id = !Values.D_NaN
	if (size (group_name, /type) ne 7) then group_name = ''

	if (finite (file_id, /NaN)) then begin
		print, "ERROR: no HDF5 file is open!"
		stop
		return, !Values.D_NaN
	end

	if (group eq group_name) then return, group_content

	object = h5g_get_objinfo (file_id, group)
	if (strupcase (object.type) ne 'GROUP') then return, [ '' ]

	num = h5g_get_nmembers (file_id, group)
	number = num
	if (keyword_set (maximum)) then num = num < (maximum > 1)
	list = strarr (num)
	for pos = 0L, num-1L do begin
		list[pos] = h5g_get_member_name (file_id, group, pos)
	end

	group_name = group
	group_content = list

	return, list
end
