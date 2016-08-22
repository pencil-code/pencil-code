#/bin/sh

# This script emulates "git pull --rebase --autostash" for older versions of
# Git. Essentially the entire code is taken from
#
#     https://github.com/git/git/blob/master/git-rebase.sh
#     https://github.com/git/git/blob/master/git-sh-setup.sh
#
# Author: Tobias Heinemann <tobias.heinemann@gmail.com>
# Date: Fri Aug 12 2016

die () {
    die_with_status 1 "$@"
}

die_with_status () {
    status=$1
    shift
    printf >&2 '%s\n' "$*"
    exit "$status"
}

git_dir_init () {
	GIT_DIR=$(git rev-parse --git-dir) || exit
	if [ -z "$SUBDIRECTORY_OK" ]
	then
		test -z "$(git rev-parse --show-cdup)" || {
			exit=$?
			echo >&2 "You need to run this command from the toplevel of the working tree."
			exit $exit
		}
	fi
	test -n "$GIT_DIR" && GIT_DIR=$(cd "$GIT_DIR" && pwd) || {
		echo >&2 "Unable to determine absolute path of git directory"
		exit 1
	}
	: "${GIT_OBJECT_DIRECTORY="$(git rev-parse --git-path objects)"}"
}

require_clean_work_tree () {
    git rev-parse --verify HEAD >/dev/null || exit 1
    git update-index -q --ignore-submodules --refresh
    err=0

    if ! git diff-files --quiet --ignore-submodules
    then
	echo >&2 "Cannot $1: You have unstaged changes."
	err=1
    fi

    if ! git diff-index --cached --quiet --ignore-submodules HEAD --
    then
	if [ $err = 0 ]
	then
	    echo >&2 "Cannot $1: Your index contains uncommitted changes."
	else
	    echo >&2 "Additionally, your index contains uncommitted changes."
	fi
	err=1
    fi

    if [ $err = 1 ]
    then
	test -n "$2" && echo >&2 "$2"
	exit 1
    fi
}

# This sets GIT_DIR
git_dir_init

# Temporary directory for autostashing
state_dir="$GIT_DIR"/pc-autostash

# Read rebase.autostash option from Git config -- default to true
autostash=$(git config --bool rebase.autostash || echo true)

# Stash uncommited changes
if test "$autostash" = true && ! (require_clean_work_tree) 2>/dev/null
then
    stash_sha1=$(git stash create "autostash") ||
    die "Cannot autostash"
    mkdir -p "$state_dir" &&
    echo $stash_sha1 >"$state_dir/autostash" &&
    stash_abbrev=$(git rev-parse --short $stash_sha1) &&
    echo "Created autostash: $stash_abbrev" &&
    git reset --hard
fi

# Do the actual "git pull --rebase"
git fetch
git rebase

# Apply previously saved stash
if test -f "$state_dir/autostash"
then
    stash_sha1=$(cat "$state_dir/autostash")
    if git stash apply $stash_sha1 2>&1 >/dev/null
    then
	echo "Applied autostash."
    else
	git stash store -m "autostash" -q $stash_sha1 ||
	die "Cannot store \$stash_sha1"
	echo Applying autostash resulted in conflicts.
        echo Your changes are safe in the stash.
        echo You can run "git stash pop" or "git stash drop" at any time.
    fi
fi

# Remove temporary directory
rm -rf "$state_dir"
