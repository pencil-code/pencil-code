#!/bin/bash

set -e

pc_cleansrc
if test -e "src/astaroth"; then
	rm -r src/astaroth
fi

pc_build

if test -e data; then
		rm -r data
fi
mkdir data
pc_run
