; $Id$
;+
; NAME: REMOVE_TAG
;       
; PURPOSE:
;       Removes a tag from a structure by creating a new structure and (by default) removing the old one.
;       If old_struct is not a structure, returns 0.
;
; CATEGORY:
;       General helpers.
;
; CALLING SEQUENCE:
;       new_struct = REMOVE_TAG, old_struct, tag [, /KEEP]
;
; KEYWORDS:
;       KEEP: keep old structure
;
; MODIFICATION HISTORY:
;       Written by:    M. Rheinhardt, October 3rd, 2026
;-
function remove_tag, struct, tag, keep=keep

  if not is_struct(struct) then return, 0
  if not is_str(tag) then return, struct

  tags = tag_names(struct)
  inds = where(tags ne strupcase(tag),count)

  if count eq 0 then return, struct

  new_struct = create_struct(tags[inds[0]], struct.(inds[0]))

  ; Append the rest of the kept tags
  for i=1,count-1 do $
    new_struct = create_struct(new_struct, tags[inds[i]], struct.(inds[i]))

  if not keyword_set(keep) then undefine, struct
  return, new_struct
end
