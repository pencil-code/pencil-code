#
#  This should be sourced before you start using the code.
#  You may or may not want to put the line
#  . $HOME/f90/pencil_modular/sourceme.sh
#  into your .bashrc
#
CDPATH="../:../../:../../../:$HOME/"
PATH=".:../bin:../../bin:../../../bin:../../../../bin:${HOME}/bin:${PATH}"
IDL_PATH="+../idl:+../../idl:+../../../idl:tmp:+$HOME/idl:<IDL_DEFAULT>"

export CDPATH PATH IDL_PATH

# alias local='cp -p $1 tmp.$$; \rm $1; mv tmp.$$ $1; chmod u+w $1'
local() {
  cp -p $1 tmp.$$; \rm $1; mv tmp.$$ $1; chmod u+w $1;
}
