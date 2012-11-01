#!/usr/bin/env bash

export OMP_NUM_THREADS=1

# load shared functions and variables
. `dirname $0`/run-tests.functions

inform "Running Velvet Test Suite"

# check we have our binaries and test files
for FILE in $VH $VG $SEQ $ROADMAP $FQL $FQR $FQI $FAL $FAR $FAI ; do
  if [ ! -r $FILE ]; then
    problem "required testing file '$FILE' not found"
  else
    inform "Found $FILE, ok" 
  fi
done

# create a temp dir
export DIR=$(mktemp -d velvet.test.XXXXXX)
inform "making test folder: $DIR" 

# run each XXXX.t test script
NUMTESTS=$(ls -1 *.t | wc -l)
COUNT=0
for TEST in *.t ; do 
  COUNT=$((COUNT + 1))
  inform "running test $COUNT of $NUMTESTS: $TEST"
  . ./$TEST
  inform "passed test $COUNT"
done

# clean up
inform "removing test folder: $DIR"
rm -fr ./$DIR

# all done
inform "passed all $NUMTESTS tests"
inform "hooray!"
