#
#  This file tries to set the PENCIL_HOME environment variable if it
#  doesn't exist yet, and then adds stuff to your PATH and IDL_PATH.
#
#  If _sourceme_quiet is set, no output is printed, which enables you to
#  put the lines
#    export PENCIL_HOME=[...]
#    _sourceme_quiet=1; . $PENCIL_HOME/sourceme.sh; unset _sourceme_quiet
#  into your .bashrc
#

#  You may or may not want to put the lines
#    setenv PENCIL_HOME [...]
#    . $HOME/f90/pencil_modular/sourceme.sh
#  into your .bashrc

if [ -z $PENCIL_HOME ]; then
  unset _sourceme		# tabula rasa without PENCIL_HOME
  #
  # Try to identify position of the code's home directory:
  #
  for _dir in   . .. ../.. ../../.. ../../../.. \
                pencil pencil-code \
		f90/pencil f90/pencil-code \
		pencil_modular f90/pencil_modular ; do
    if ( [ -e $_dir/sourceme.csh ] && \
         [ -d $_dir/bin ]          && \
	 [ -d $_dir/doc ]          && \
	 [ -d $_dir/src ]          && \
	 [ -d $_dir/samples ]         \
       ); then
      PENCIL_HOME=`cd $_dir; echo $PWD`
      export PENCIL_HOME
      break
    fi
  done
  unset _dir

  if [ -z $PENCIL_HOME ]; then # no success
    echo "sourceme.sh: Cannot locate home directory of pencil code."
    echo "  Try sourcing me from the home directory itself, or set PENCIL_HOME"
    exit 0
  fi
fi

if [ -z $_sourceme_quiet ]; then echo "PENCIL_HOME = <$PENCIL_HOME>"; fi

if [ -z $_sourceme ]; then	# called for the first time?
  # CDPATH="./:../:../../:../../../:$HOME"
  if [ -d $PENCIL_HOME/bin ]; then
    #  Set shell path
    if [ -z $_sourceme_quiet ]; then echo "Adding $PENCIL_HOME/bin to PATH"; fi
    PATH=${PATH}:$PENCIL_HOME/bin
    
    IDL_PATH="./idl:../idl:+${PENCIL_HOME}/idl:./data:./tmp:${IDL_PATH=<IDL_DEFAULT>}"

    DXMACROS="${PENCIL_HOME}/dx/macros:${DXMACROS}"
    _sourceme="set"

    # export CDPATH PATH IDL_PATH
    export PATH DXMACROS IDL_PATH _sourceme
    
  else
    echo "Not adding $PENCIL_HOME/bin to PATH: not a directory"
  fi
fi
