from pyne.material import Material

searchfile = open("Material_Description.txt", "r")
for line in searchfile:
    dictionary={}
    if 'MCNP Form' in line:
        matName_line = searchfile.next()
	mat_name= matName_line[matName_line.find(" ")+1: matName_line.find("rho")-1]
	density= 0
	if "rho = " in matName_line:
	    density= matName_line[matName_line.find("rho = ")+6: matName_line.find(" ", matName_line.find("rho = ")+7)]
	mat_name= mat_name.lower()
	mat_name= mat_name.replace("(","_")
	mat_name= mat_name.replace(" ","_")
	mat_name= mat_name.replace("-","_")
	mat_name= mat_name.replace(",","")
	mat_name= mat_name.replace("__","_")
	mat_name= mat_name.replace(")","")
	mat_name= mat_name.replace(".","")
	if mat_name[-1]== '_':
	    mat_name= mat_name[0:-1]
	print '\n', "Material Name: ", mat_name
	print "DENSITY: ", density
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
	mat = Material(dictionary, float(density), mat_name, 1.0)
	filename= "Material Library/" +mat_name + ".txt"
    	mat.write_text(filename)
