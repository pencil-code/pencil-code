#
#  This file tries to set the PENCIL_HOME environment variable if it
#  doesn't exist yet, and then adds stuff to your PATH and IDL_PATH.
#
#  You may or may not want to put the lines
#    setenv PENCIL_HOME [...]
#    . $HOME/f90/pencil_modular/sourceme.sh
#  into your .bashrc
#

if [ -z $PENCIL_HOME ]; then
  #
  # Try to identify position of the code's home directory:
  #
  for _dir in   . .. ../.. ../../.. ../../../.. pencil pencil_modular \
                f90/pencil f90/pencil_modular ; do
    if ( [ -e $_dir/sourceme.csh ] && \
         [ -d $_dir/bin ]          && \
	 [ -d $_dir/doc ]          && \
	 [ -d $_dir/src ]          && \
         [ -d $_dir/runs ]            \
       ); then
      PENCIL_HOME=`cd $_dir; echo $PWD`
      break
    fi
  done
  unset _dir

  if [ -z $PENCIL_HOME ]; then # no success
    echo "Cannot locate home directory of pencil code."
    echo "Try sourcing me from the home directory itself, or set PENCIL_HOME"
    exit 0
  fi
fi

echo "PENCIL_HOME = <$PENCIL_HOME>"

if [ -z $_sourceme ]; then	# called for the first time?
  # CDPATH="./:../:../../:../../../:$HOME"
  if [ -d $PENCIL_HOME/bin ]; then
    #  Set shell path
    echo "Adding $PENCIL_HOME/bin to path"
    PATH=${PATH}:$PENCIL_HOME/bin
    IDL_PATH="./idl:../idl:+${PENCIL_HOME}:./tmp:${IDL_PATH=<IDL_DEFAULT>}"
    _sourceme="set"

    # export CDPATH PATH IDL_PATH
    export PATH IDL_PATH _sourceme
    
    # Now a separate shell script:
    # local() {
    #   cp -p $1 tmp.$$; \rm $1; mv tmp.$$ $1; chmod u+w $1;
    # }
  else
    echo "Not adding $PENCIL_HOME/bin to path: not a directory"
  fi
fi
