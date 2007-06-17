#!/bin/csh
#
#alias gb 'cd $gt ; set gg=$gb ; set gb=$gt ; set gt=$gg ; echo $gt "->" $gb'
#alias gt 'set gt=$cwd; cd \!^ ; set gb=$cwd ; echo $gt "->" $gb'
#
du -sm *
mkdatadir_scratch
gt data1
pc_mkproctree 8
gb
move-VAR-to-data1
du -sm *
