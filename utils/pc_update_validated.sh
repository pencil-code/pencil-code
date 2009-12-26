#!/bin/sh
#
#  $Id$
#
#  If $PC_SET_VALIDATED is set, then update validated.dat
#
if [ $PC_SET_VALIDATED ]; then
  cd $PENCIL_HOME/misc/validation;
  awk '{print $4}' `uname -n`.auto-test |head -1 >validated.dat;
  svn ci -m "updated by auto-test on `uname -n` by $USER";
fi
