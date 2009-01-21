## $Id$
##
## This file can be used by copy-snapshots to list the data directory on all
## the nodes.
##
## Usage: Copy file to data/copy-snapshots-eval-me.csh during a run.
## (the file must be executable; chmod 755 data/copy-snapshots-eval-me.csh)
##
## Only works with copy-snapshots_exp (set lcopysnapshots_exp=T in start.in)
##
set nodelist = (`echo $NODELIST | sed 's/:/ /g'`)

foreach node ($nodelist)
  echo "----------------- node: $node ----------------"
  $SSH $node ls -l $SCRATCH_DIR $SCRATCH_DIR/proc*
end
