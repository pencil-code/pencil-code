#!/bin/csh

#
# Name: pc_llsubmit
# Author: Anders Johansen
# Date: 20-aug-2004
# Description:
#   Shell script to easily submit jobs to the LoadLeveler queuing system.
#   Set non-default variables in small script jobscript.csh and include
#     source pc_llsubmit $argv
#   in the file. Then do e.g.
#     jobscript.csh start_run.csh
#   from the command line.
#
# Based on sleipsub (by theine)
#

#
# defaults
#
if (! $?job_name) set job_name = $argv
if (! $?job_type) set job_type = parallel
if (! $?total_tasks) set total_tasks = `perl -ne '$_ =~ /^\s*integer\b[^\\!]*ncpus\s*=\s*([0-9]*)/i && print $1' src/cparam.local`
if (! $?class) set class = qexp
if (! $?environment) set environment = "PENCIL_HOME=$PENCIL_HOME"
if (! $?input) set input = /dev/null
if (! $?output) set output = \$\(job_name\).\$\(jobid\).out
if (! $?error) set error = $output
if (! $?notification) set notification = never
#
# write temporary jobscript
#
echo "#\!/bin/sh" > pc_llsubmit.tmp
echo "# @ job_name = $job_name" >> pc_llsubmit.tmp
echo "# @ job_type = $job_type" >> pc_llsubmit.tmp
echo "# @ total_tasks = $total_tasks" >> pc_llsubmit.tmp
echo "# @ class = $class" >> pc_llsubmit.tmp
echo "# @ environment = $environment" >> pc_llsubmit.tmp
echo "# @ executable = $argv" >> pc_llsubmit.tmp
echo "# @ input = $input" >> pc_llsubmit.tmp
echo "# @ output = $output" >> pc_llsubmit.tmp
echo "# @ error =  $error" >> pc_llsubmit.tmp
echo "# @ notification = $notification" >> pc_llsubmit.tmp
echo "# @ resources = ConsumableCpus(1)" >> pc_llsubmit.tmp
echo "# @ queue" >> pc_llsubmit.tmp
#
# submit temporary jobscript and remove it afterwards
#
llsubmit pc_llsubmit.tmp && rm pc_llsubmit.tmp

# End of file pc_llsubmit
