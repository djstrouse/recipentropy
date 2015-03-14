#!/bin/sh

# test2.sh
# 
#
# Created by Ari Strandburg-Peshkin on 8/23/13.
# Copyright 2013 __MyCompanyName__. All rights reserved.

# declare an array called array and define 3 vales
array=( one two three )
for i in "${array[@]}"
do
	echo $i
done