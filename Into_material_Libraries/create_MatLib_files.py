from pyne.material import Material

searchfile = open("Material_Description.txt", "r")
for line in searchfile:
    dictionary={}
    if 'MCNP Form' in line:
        matName_line = searchfile.next()
	mat_name= matName_line[matName_line.find(" ")+1: matName_line.find("rho")-2]
	print '\n', "Material Name: ", mat_name
        for line in searchfile:
	    if "-----------" in line:
		break
	for line in searchfile:
	    if "___________" in line:
		break
	    else:
	        print "Isotopes: ", line
		data= line.split()
		dictionary.update({str(data[0]): float(data[1])})
	mat = Material(dictionary, 1.0, mat_name, 1.0)
	filename= "Material Library/" +mat_name + ".txt"
    	mat.write_text(filename)
