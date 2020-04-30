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
      unset cd   # some people are crazy enough to overload cd
      PENCIL_HOME=`cd $_dir; echo $PWD`
      export PENCIL_HOME
      break
    fi
  done
  unset _dir

  if [ -z $PENCIL_HOME ]; then # no success
    echo "sourceme.sh: Cannot locate home directory of pencil code."
    echo "  Try sourcing me from the home directory itself, or set PENCIL_HOME"
  fi
fi

if [ -z $_sourceme_quiet ]; then echo "PENCIL_HOME = <$PENCIL_HOME>"; fi

if [ -z $_sourceme ]; then	# called for the first time?
  # CDPATH="./:../:../../:../../../:$HOME"
  if ([ -n $PENCIL_HOME ] && [ -d $PENCIL_HOME/bin ]); then

    #  Set shell path
    if [ -z $_sourceme_quiet ]; then echo "Adding $PENCIL_HOME/{bin,utils{,/axel},remesh/bin} to PATH"; fi
    # PATH=${PATH}:$PENCIL_HOME/bin:$PENCIL_HOME/utils:$PENCIL_HOME/utils/axel:$PENCIL_HOME/remesh/bin
    PATH=${PATH}:$PENCIL_HOME/bin:$PENCIL_HOME/utils:$PENCIL_HOME/utils/axel:$PENCIL_HOME/utils/xiangyu:$PENCIL_HOME/remesh/bin:$PENCIL_HOME/src/scripts

    if ([ -d $PENCIL_HOME/src/astaroth/submodule/scripts ]); then
      export AC_HOME=$PENCIL_HOME/src/astaroth/submodule
      export PATH=${PATH}:$AC_HOME/scripts/
    fi

    #  Set path for DX macros
    DXMACROS="${PENCIL_HOME}/dx/macros${DXMACROS:+:$DXMACROS}"

    #  Set IDL path
    IDL_PATH="./idl:../idl:+${PENCIL_HOME}/idl:./data:./tmp:${IDL_PATH=<IDL_DEFAULT>}"

    # Set HDF5 path
    if [ -z $HDF5_HOME ]; then HDF5_HOME=`which h5fc 2>/dev/null | sed -e 's/\/bin\/h5fc$//'`; fi

#    #  Set Perl module path [no longer needed]
#    PERL5LIB="${PENCIL_HOME}/lib/perl${PERL5LIB:+:$PERL5LIB}"
#    export PERL5LIB
    #   Set PYTHONPATH
    if [ -z $PYTHONPATH ]; then
       PYTHONPATH="$PENCIL_HOME/python:$PWD/python"
    else
       PYTHONPATH="$PYTHONPATH:$PENCIL_HOME/python:$PWD/python"
    fi
    #  Set library path for linker
    if [ -z $LD_LIBRARY_PATH ]; then
      LD_LIBRARY_PATH="./src:./src/astaroth:./src/astaroth/submodule/build/src/core"
    else
      LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:./src:./src/astaroth:./src/astaroth/submodule/build/src/core"
    fi

    # Remember that sourceme has been successfully run
    _sourceme="set"

    # export CDPATH PATH IDL_PATH
    export PATH DXMACROS IDL_PATH PYTHONPATH _sourceme

  else
    if [ -n $PENCIL_HOME ]; then
      echo "Not adding $PENCIL_HOME/bin to PATH: not a directory"
    fi
  fi
fi
