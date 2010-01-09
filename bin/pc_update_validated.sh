#!/bin/sh
#
#  $Id$
#
#  If $PC_SET_VALIDATED is set, then update validated.dat
#  Otherwise, don't do anything.
#  This script is used by auto-test and pc_auto-test.
#
if [ $PC_SET_VALIDATED ]; then
  cd $PENCIL_HOME/misc/validation;
  awk '{print $2}' /tmp/pc_current_revision.dat >validated.dat;
  svn ci -m "automatic validation completed: auto-test on `uname -n` by $USER";
fi
