#!/bin/bash

# Sample /etc/crontab entry:
# 59 *		* * *   pencil	/home/pencil/run-test.sh

# kill all processes that ran for 3 hours or more:
killall -9 -q --older-than 3h -u pencil
