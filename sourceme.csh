#
#  Tries to set the PENCIL_HOME environment variable if it
#  doesn't exist yet, and then adds stuff to your PATH and IDL_PATH.
#  All the work is done by sourceme_quiet.csh.

set _sourceme_verbose

source sourceme_quiet.csh

unset _sourceme_verbose
