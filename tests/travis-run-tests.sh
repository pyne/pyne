#!/bin/bash
    
for x in $(ls test*.py */test*.py); do 
  echo
  echo
  echo "Testing $x:" 
  echo 
  if [ $x == "ensdf_processing.py" ]; then
    nosetests test_ensdf_processing.py --process-timeout=120
  else
    nosetests "$x"
  fi
  status=$?
  if [ $status -ne 0 ]; then 
    echo "Testing failed on $x" 
    exit $status
  fi
done
