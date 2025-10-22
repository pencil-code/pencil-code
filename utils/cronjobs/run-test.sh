#!/bin/bash

# Sample /etc/crontab entries:
### default extended test:
# 50 2		* * *   pencil	/home/pencil/run-test.sh
### same as above, but different output filenames:
# 50 2		* * *   pencil	/home/pencil/run-test.sh 3 "extended"
### normal test:
# 10 */2	* * *   pencil	/home/pencil/run-test.sh 2 "normal" "normal-previous"
### basic test, minutely from 5 to 23:
# 0 5-23	* * *   pencil	/home/pencil/run-test.sh 0,1 "basic"

# settings
TEST="auto-test"
LAST="previous"
LOGFILE="index.html"
HOME_DIR="/home/pencil"
TARGET_DIR="/var/www/pencil-code/tests"
UPDATE_DIR="/home/pencil/UPDATE"
SVN_CHECKOUT="https://pencil-code.org/svn/trunk/"

# set maximum test level, if given explicitly
if [ -n "$1" ] ; then
	LEVEL="$1"
fi
# set output filename, if given explicitly
if [ -n "$2" ] ; then
	TEST="$2"
fi
# set backup filename, if given explicitly
if [ -n "$3" ] ; then
	LAST="$3"
fi

# set test directory
TEST_DIR="${HOME_DIR}/${TEST}"

# check for other processes
if [[ -n `ps aux | grep -vE '[0-9]:[0-9][0-9] grep ' | grep -E "pencil-test .*--pencil-home=${TEST_DIR}"` ]] ; then
	echo "The same PC autotest is running: ${TEST_DIR}"
	exit 1
fi

# check for a scheduled test
icon_hash=`md5sum "${TARGET_DIR}/${TEST}/result.svg" | head -c 32`
if [[ "$icon_hash" == "980c856192f62bb123aa580d9f2a8a97" ]] ; then
	scheduled="yes"
fi

# check for new commits
if [[ ! -e "${UPDATE_DIR}/${TEST}" ]] ; then
	echo "No update infrastructure: ${TEST}"
	exit 1
fi
if [[ ! -s "${UPDATE_DIR}/${TEST}" && -z "$sheduled" ]] ; then
	# no new commits
	exit 0
fi

# prepare test directory
if [[ ! -e "${TEST_DIR}" ]] ; then
	cd "${HOME_DIR}"
	svn checkout -q ${SVN_CHECKOUT} "${TEST}"
else
	if [[ ! -d "${TEST_DIR}" ]] ; then
		echo "There is something in the way: ${TEST_DIR}"
		exit 1
	fi
	cd "${TEST_DIR}"
	updated=`svn cleanup && svn up | grep -E '^(Restored |Updated to revision [0-9]+\.|Wieder hergestellt |Aktualisiert zu Revision [0-9]+\.)'`
fi

# prepare test command
#
TEST_CMD="bin/pencil-test"
TEST_OPT="--use-pc_auto-test --auto-test-options='--level=${LEVEL},--config-files=GNU-GCC_MPI+GNU-GCC_debug,--auto-clean,--bisect,--local-lock' --nice=15 --short --local --html"
# --mail=me@myself.org

# actual code
if [[ -n "$updated" || -n "$scheduled" ]] ; then
	cd "${TEST_DIR}"
	rm -f .run_directories.log

	: > "${UPDATE_DIR}/${TEST}"
	${TEST_CMD} ${TEST_OPT} --log-dir="${TARGET_DIR}/${TEST}" --previous-dir="${TARGET_DIR}/${TEST}-${LAST}" --pencil-home="${TEST_DIR}" > "${TARGET_DIR}/${TEST}/${LOGFILE}"

#else
#	echo "Triggered while code is up-to-date."
#	: > "${UPDATE_DIR}/${TEST}"
fi
