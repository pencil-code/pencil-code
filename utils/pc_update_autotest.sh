#!/bin/sh
#
#  $Id$
#
#  If $PC_SET_VALIDATED is set, then update validated.dat
#
if [ $PC_SET_VALIDATED ]; then
  cd $PENCIL_HOME/misc/validation;
  svn up;
  echo "`date` (previous validated revision: `cat validated.dat`)" \
    >> `uname -n`.autotest;
  svn ci -m "updated by auto-test on `uname -n` by $USER";
fi

