#!/bin/false

get_temp_filename ()
{
  TMP_DIR=/tmp
  local newfile=${TMP_DIR}/pc_${USER}_`date +$H%M%s%d%m%Y`
#  if [ ! "$TMP_DIR" ]; then
#    if [ -d /tmp ]; then
#      $TMP_DIR=/tmp
#    else if [ -d /usr/tmp ]; then
#      $TMP_DIR=/usr/tmp
#    else if [ -d "$HOME" ]; then
#     [ ! -e "$HOME/.tmp"] && mkdir $HOME/.tmp && $TMP_DIR=$HOME/.tmp
#    fi
#    if [ ! "$tmp_dir" ]; then
#      echo "CANNOT FIND DIRECTORY FOR TEMPORARY FILES: need to set TMP_DIR in your pencil-code config."
#      exit 1
#    fi
#  fi


#  while [ -e "$newfile" ]; do 
    newfile="${TMP_DIR}/pc_${USER}_`date +$H%M%s%d%m%Y`"
#  done
  echo $newfile
}

usage_from_header ()
{
  # Knicked from one of Wolfgang's perl scripts. 
  # Paragraph _must_ contain `Description:' or `Usage:'
  # Drop `Author:', etc. (anything before `Description:' or `Usage:')
  # Don't print comment sign:
  cat $1 | perl -e 'local $/ = ""; while (<STDIN>) { next unless /^\s*\#\s*(Description|Usage):/m; s/.*?\n(\s*\#\s*(Description|Usage):\s*\n.*)/$1/s; s/^\s*# ?//mg; last; } print ($_ or "<No usage information found>\n");'
}

split_path ()
{
  local file=$1
}

cvs_get_status ()
{
  local cvsfile=$1
  local temp_file=`get_temp_filename`
  local restorecwd=`pwd`

  if [ -e "$cvsfile" ]; then
    if [ -h $cvsfile ]; then
      local oldfile=$cvsfile
      local linkedfile=`readlink $oldfile`
      local linkedpath=`echo $linkedfile | sed -e 's@\(.*\)\/\([^\/]*\)$@\1@'`
      local cvsfile=`echo $linkedfile | sed -e 's@\(.*\)\/\([^\/]*\)$@\2@'`
      cd $linkedpath
    fi
    cvs status $cvsfile > $temp_file
    cvs_status=`grep "Status:" $temp_file | sed -e 's/.*Status: \(.*\)$/\1/'`
    cvs_repository_revision=`grep "Repository revision:" $temp_file | sed -e 's/.*Repository Revision: \(.*\)$/\1/'`
    rm -f $temp_file
  else
    cvs_status="Missing File"
    [ $cvs_repository_version ] unset cvs_repository_version 
  fi
  cd $restorecwd
}
