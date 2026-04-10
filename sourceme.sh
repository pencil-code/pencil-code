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
  # "$(dirname "${BASH_SOURCE[0]}")" should always work in bash. The other
  # directories below are needed in case the user is using some other shell
  # (zsh, ash, dash, ...).
  #
  for _dir in "$(dirname "${BASH_SOURCE[0]}")" \
      . .. ../.. ../../.. ../../../.. \
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
    # remove first all paths, which contain "pencil-code" from PATH, then add (new) PC paths
    PATH=`echo $PATH | sed -e's/[^:]*pencil-code[^:]*://g' -e's/[^:]*pencil-code[^:]* *$//' -e's/:$//'`:$PENCIL_HOME/bin:$PENCIL_HOME/utils:$PENCIL_HOME/utils/axel:$PENCIL_HOME/utils/xiangyu:$PENCIL_HOME/remesh/bin:$PENCIL_HOME/src/scripts:./src/pre_and_post_processing

    export AC_HOME=$PENCIL_HOME/src/astaroth/submodule
    if ([ -d $PENCIL_HOME/src/astaroth/submodule/scripts ]); then
      PATH=${PATH}:$AC_HOME/scripts/
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
      LD_LIBRARY_PATH="./src:./src/astaroth:./src/astaroth/submodule/build/src/core:./src/astaroth/submodule/build/src/core/kernels:./src/astaroth/submodule/build/src/utils"
    else
      LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:./src:./src/astaroth:./src/astaroth/submodule/build/src/core:./src/astaroth/submodule/build/src/core/kernels:./src/astaroth/submodule/build/src/utils"
    fi

    # Remember that sourceme has been successfully run
    _sourceme="set"

    # export CDPATH PATH IDL_PATH
    export PATH DXMACROS IDL_PATH PYTHONPATH LD_LIBRARY_PATH _sourceme

  else
    if [ -n $PENCIL_HOME ]; then
      echo "Not adding $PENCIL_HOME/bin to PATH: not a directory"
    fi
  fi
fi

if [ -d .git ]; then
# 2025-Nov-11/Kishore: Matthias, I think it would be cleaner to check the config
# value by running
# `git config get pull.rebase` and checking that it does not return "false"
# I think the above has the advantage of also checking the value inherited from
# the global config, if any.
	if [[ `grep '^\srebase *= *false' .git/config` != "" ]]; then
	echo !!!WARNING - you have \"rebase = false\" settings in your .git/config!!!
	echo !!!Pull strategy should always be \"--rebase\" on all branches!!!
    fi
#
# Enforce basic pull policy to "rebase".
#
# Added -C flag to change the directory of the git command to $PENCIL_HOME
    git -C $PENCIL_HOME config pull.rebase true
#
#
# Enforce that all committers have set their email account and username in git before committing
# 2026-Apr-10/Kishore: commented the following as it overwrites the user's pre-existing hooks. In general, it is not safe (and very impolite) to automatically add hooks to a user's local git repo. @Touko, as far as I know, git already enforces that the name and email are set (I also tested this on a new user account on my machine). What is the problem that you are trying to solve?
# 2026-Apr-10/TP: @Kishore git does not enforce it since I was making commits without the user info set and was asked explicitly by Philippe to add my correct info, and the same held true for Fred and Matthias at some point. Could be that this is related to committing from HPC computers.
# Could we consider appending this snippet to the end of the pre-commit hook instead of overwriting?
# This could be of course done from the server side but would not know how to do it because of the mirroring shenanigans.
# #
#     touch $PENCIL_HOME/.git/hooks/pre-commit
#     echo "#!/bin/sh
# name=\$(git config user.name)
# email=\$(git config user.email)
# if [ -z \"\$name\" ] || [ -z \"\$email\" ]; then
#   echo \"Error: Git user.name and user.email must be set before committing!!\"
#   exit 1
# fi" > $PENCIL_HOME/.git/hooks/pre-commit
# 
#     chmod +x $PENCIL_HOME/.git/hooks/pre-commit
# 
fi
