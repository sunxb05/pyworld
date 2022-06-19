#!/bin/bash
#
# simple script for CTest

# checks if test_gen1int exists
if [ -e test_gen1int ]; then
  # we should find "number of failed tests:     0" once if the test passed
  PASSED=`./test_gen1int $1 | grep "number of failed tests:     0" | wc -l`
  # for CTest
  if [ $PASSED -eq 1 ]; then
    exit 0
  else
    exit 1
  fi
else
  exit 1
fi
