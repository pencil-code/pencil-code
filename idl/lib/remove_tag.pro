function remove_tag, struct, tag

  if not is_struct(struct) then return, 0
  if not is_str(tag) then return, struct

  tags = tag_names(struct)
  inds = where(tags ne strupcase(tag),count)

  if count eq 0 then return, struct

  new_struct = create_struct(tags[inds[0]], struct.(inds[0]))

  ; Append the rest of the kept tags
  for i=1,count-1 do $
    new_struct = create_struct(new_struct, tags[inds[i]], struct.(inds[i]))

  undefine, struct
  return, new_struct
end
