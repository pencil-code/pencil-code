   function tag_exists, struct, field
;
;  Tests whether structure struct has an element named field.
;
     return, (where(tag_names(struct) eq field))(0) ne -1

   end
