#!/bin/bash

find -type f -name families_core.txt | while read family; do
	echo "$(dirname $family)"
	grep -c 'fam' $family
done
