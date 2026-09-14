#!/bin/bash
for i in $(seq 1 999); do
	ln -s host-arrhenius1.hpc.arrhenius.naiss.se-GNU_Linux-AlmaLinux.conf "host-n${i}-GNU_Linux-AlmaLinux.conf" &> /dev/null
done
