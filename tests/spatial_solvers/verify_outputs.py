# Open a file


with open("fort.8") as infile:
  for line in infile:
    print line
    if ("@##@##@##@SCTSTEP" in line):
      print "Found SCTSTEP"
      with open("outs/sct_simple") as checksct:
        for sct_line in checksct:
          print sct_line
      break;

#str = fo.read(10);
print "END"

