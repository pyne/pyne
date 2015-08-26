#!/bin/bash
for x in $(ls test*.py */test*.py); do 
  echo
  echo
  echo "Testing $x:" 
  echo 
  nosetests "$x"
  status=$?
  if [ $status -ne 0 ]; then 
    echo "Testing failed on $x" 
    exit $status
  fi
done
