#!/bin/bash

set -euo pipefail

is_first=1
ls hap.py/*.summary.csv | while read x
do
	id="${x%.summary.csv}"
	id="${id##*/}"
	
	if [ 1 -eq ${is_first} ]
	then
		echo -n "ID,"
		head -n 1 "${x}"
	fi

	tail -n +2 "${x}" | while read y
	do
		echo -n "${id},"
		echo "${y}"
	done

	is_first=0
done
