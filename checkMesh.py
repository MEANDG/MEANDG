'''
This file checks the following files 
points, faces, owner, neighbour
for the structure anomaly. If the data is printed on a single line, this program splits it over multiple lines the way DG solver needs it.
To run this program,

$python checkMesh.py dirName

e.g. if the directory is stored with the name ./app/tmp/
note: each dir should have OpenFOAM file structure inside of it. i.e. the geometry files for the tmp directory 
should be inside of app/tmp/constant/polymesh 

then type

$python checkMesh.py app/tmp

or

$python checkMesh.py app/tmp/

To restore the directory to the previous state, run restoreDir.py 
'''

import numpy as np
import os 
import sys


# This function checks whether the file (given as an argument) contains the data in a sigle line or not
# returns 1 if data is over multiple lines, 0 if data is in the single line
# Code currently is configured to read the data over multiple lines only. So the function makeFile (described next) is run
# to re-arrange the file
def checkFile(fileName):
	a = 0;
	count = 0;
	with open (fileName) as file:
		for line in file:
			if (line[0] == "("): # This will happen only if the data is already given over multiple lines
				a = 1
				break;
			count += 1;
			if (count > 50): # This condition in order to not spend too much time going through large files. 
				break;
			pass;
		pass;
	return a;


# Once it is confirmed by the checkFile function that the data is presented only on a single line, 
# as done by OpenFOAM if the dataFile is small, 
# This function changes the file so that it appears as if the data is on multiple lines
def makeFile(fileName):
	fileCopy = fileName + "_copy"; # a copy is created of the original file
	os.system(" touch " + fileCopy); 
	data = [];

	# All lines of the file are stored in an empty list 'data'
	with open (fileName) as file:
		for line in file:
			data.append(line);
			pass;
		pass;

	file.close()

	# The copy file is opened for writing the data
	with open (fileCopy, "w") as file:
		# Because if there is only 1 cell in the domain (lets say with 6 faces), then the owner file will look as follows:

		# OpenFOAM header
		#//***************//
		# 6{0}
		#//***************//
	
		#What we want instead is
		# OpenFOAM header
		#//***************//
		# 6
		#(
		#0
		#0
		#0
		#0
		#0
		#0
		#)
		#//***************//

		#That is exactly done here
		if ("owner" in fileName):
			for line in data:
				if ("{" in line and "}" in line):
					line = line.replace(" ","");
					noOfFaces = line[0]; # The first letter gives total number of faces
					print >> file, noOfFaces;
					for word in line[1:]:
						if (word == "{"):
							word = "("
							print >> file, word;
						elif (word == "}"):
							word = ")"
							print >> file, word;
						else:
							for k in range(int(noOfFaces)):
								print >> file, word;
								pass;
							pass
						pass;
				else:
					print >> file, line,
					pass;

		# Similarly if the total number of points are less, then blockMesh utitility of OpenFOAM prints everything on a single line
		#e.g.

		# OpenFOAM header
		#//***************//
		#2((0 1 1) (2 1 2)) 
		#//***************//
	
		#What we want instead is
		# OpenFOAM header
		#//***************//
		# 2
		#(
		#(0 1 1)
		#(2 1 2)
		#)
		#//***************//

		elif ("point" in fileName):
			for line in data:
				if ("(" in line):
					a = "("
					noOfPoints = line.split("(",1)[0];
					print >> file, noOfPoints;
					print >> file, "(";
					line1 = line.split("(",1)[1];
					line2 = line1.rsplit(")",1)[0]
					data = [a+b for b in line2.split("(") if b];
					for d in data:
						print >> file, d;
						pass;
					print >> file , ")";
				else:
					print >> file, line,
					pass;
				pass;
				
		else:
			for line in data:
				if ("(" in line):
					line = line.replace(" ","");
					for word in line:
						print >> file, word;
						pass;
				else:
					print >> file, line,
					pass;
		pass;
	file.close();
	os.system("cp "+fileName + " " + fileName + "_old");
	os.system("cp "+fileCopy + " " + fileName);
	os.system("rm "+fileCopy );



dirName = sys.argv[1];
if (dirName[-1] == "/"):
	pointFile = dirName + "constant/polyMesh/points";
	faceFile = dirName + "constant/polyMesh/faces";
	ownerFile = dirName + "constant/polyMesh/owner";
	neighbourFile = dirName + "constant/polyMesh/neighbour";
else:
	pointFile = dirName + "/constant/polyMesh/points";
	faceFile = dirName + "/constant/polyMesh/faces";
	ownerFile = dirName + "/constant/polyMesh/owner";
	neighbourFile = dirName + "/constant/polyMesh/neighbour";
	pass;

files = [pointFile, faceFile, ownerFile, neighbourFile];

for fileName in files:
	print "Processing ", fileName, "...";
	flag = checkFile(fileName); # Check is performed
	if (flag == 0):
		# i.e. the file contains the single line data
		print "File ", fileName, " is a single-line data file. Converting into a multi-line data file. The current file is copied with the name ", fileName + "_old";
		makeFile(fileName);
		print "Done.\n"
	else:
		print "File ", fileName, " looks good. No processing required.\n";
	pass;





