#!/bin/bash
if [ -z "$KSW2COV" ]; then
    KSW2COV=./ksw2-test
fi

# Test code coverage of ksw2
u=$1
t=$2

# Normal run, no flags enabled
$KSW2COV $u $t ${@:2}

# Enables KSW_EZ_RIGHT
$KSW2COV $u $t -r ${@:2}

## Enables KSW_EZ_GENERIC_SC
$KSW2COV $u $t -c ${@:2}

# Enables KSW_EZ_REV_CIGAR
$KSW2COV $u $t -v ${@:2}

# Enables KSW_EZ_EXTZ_ONLY
$KSW2COV $u $t -x ${@:2}

# Enables KSW_EZ_APPROX_MAX
$KSW2COV $u $t -g ${@:2}

# Enables KSW_EZ_EQX
# This option is disabled because the baseline kernel segmentation faults with this option enabled
# $KSW2COV $u $t -q ${@:2}

# Sets end_bonus and zdrop
$KSW2COV $u $t ${@:2} -x -e 1 -z 1

# Sets "score only" mode
$KSW2COV $u $t -s ${@:2}

# Modify the gap penalty to trigger long thres
$KSW2COV $u $t -x -E 2 ${@:2}
