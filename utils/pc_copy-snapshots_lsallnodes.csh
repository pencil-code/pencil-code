##
## This file can be used by copy-snapshots to list the data directory on all
## the nodes.
##
## Usage: Copy file to data/source-me.csh during a run.
##
set i=0
foreach node ($nodelist)
  echo "----------------- node: $node ----------------"
  $SSH $node ls -l $SCRATCH_DIR $SCRATCH_DIR/proc*
  set i=`expr $i + 1`
end
