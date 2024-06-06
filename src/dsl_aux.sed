# sed script for generation of field declarations for DSL code
# so far only PDE variables are considered (no auxiliaries)
# multiple register calls for the same field result in multiple declarations -> 
# run-time switches like ldensity_nolog would need to be taken into account in addition

/^ *call  *farray_register_auxiliary/ !d
# only array par
/call  *farray_register_auxiliary/ s/^[^!]*farray_register_auxiliary *( *[']\([0-9a-zA-Z_][0-9a-zA-Z_]*\)['] *,[^=]*array *= *\([a-zA-Z_0-9][a-zA-Z_0-9]*\)[^=]*$/auxiliary Field \U\1[\L\2] \\\/\\\/ \1\\n/
#
# only vector par
/call  *farray_register_auxiliary/ s/^[^!]*farray_register_auxiliary *( *[']\([0-9a-zA-Z_][0-9a-zA-Z_]*\)['] *,[^=]*vector *= *\([0-9][0-9]*\)[^=]*$/auxiliary Field \U\1X,\U\1Y,\U\1Z \\\/\\\/ \L\1\\n#define \U\1\L F\Lield\2(\U\1X,\U\1Y,\U\1Z)\\n/
#
# neither vector nor array
/call  *farray_register_auxiliary/ s/^[^!]*farray_register_auxiliary *( *[']\([0-9a-zA-Z_][0-9a-zA-Z_]*\)['] *,[^=]*$/auxiliary Field \U\1\L \\\/\\\/ \L\1\n/
#
# neither vector nor array, communicated
/call  *farray_register_auxiliary/ s/^[^!]*farray_register_auxiliary *( *[']\([0-9a-zA-Z_][0-9a-zA-Z_]*\)['] *,.*communicated *=[^=]*$/communicated auxiliary Field \U\1\L \\\/\\\/ \L\1\n/
#
# first vector then array par
/call  *farray_register_auxiliary/ s/^[^!]*farray_register_auxiliary *( *[']\([0-9a-zA-Z_][0-9a-zA-Z_]*\)['] *,[^=]*vector *= *\([0-9][0-9]*\) *,[^=]*array *= *\([a-zA-Z_0-9][a-zA-Z_0-9]*\).*$/auxiliary Field \U\1X[\L\3],\U\1Y[\L\3],\U\1Z[\L\3] \\\/\\\/ \L\1\\n\\\/\\\/#define \U\1\L F\Lield\2(\U\1X,\U\1Y,\U\1Z)\\n/
#
# first array then vector par
/call  *farray_register_auxiliary/ s/^[^!]*farray_register_auxiliary *( *[']\([0-9a-zA-Z_][0-9a-zA-Z_]*\)['] *,[^=]*array *= *\([a-zA-Z_0-9][a-zA-Z_0-9]*\) *,[^=]*vector *= *\([0-9][0-9]*\).*$/auxiliary Field \U\1X[\L\2],\U\1Y[\L\2],\U\1Z[\L\2] \\\/\\\/ \L\1\\n\\\/\\\/#define \U\1\L F\Lield\3(\U\1X,\U\1Y,\U\1Z)\\n/
#
q
