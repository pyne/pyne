#!/bin/bash
    
for x in $(ls test*.py */test*.py); do 
  echo
  echo
  echo "Testing $x:" 
  echo 
  if [ $x == "ensdf_processing.py" ]; then
    pytest test_ensdf_processing.py --process-timeout=120
  else
    pytest "$x"
  fi
  status=$?
  if [ $status -ne 0 ] && [ $status -ne 5 ]; then 
    echo "Testing failed on $x" 
    exit $status
  fi
done
