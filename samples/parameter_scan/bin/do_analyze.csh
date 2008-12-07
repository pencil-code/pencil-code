#!/bin/sh
#
#  $Id: do_analyze.csh,v 1.1 2008-12-07 17:46:20 brandenb Exp $
#
./bin/pc_analyze 0.01
./bin/pc_analyze 0.1
./bin/pc_analyze 0.2
./bin/pc_analyze 0.5
./bin/pc_analyze 1.0
./bin/pc_analyze 2.0
./bin/pc_analyze 5.0
