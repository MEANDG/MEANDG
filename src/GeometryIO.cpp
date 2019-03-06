#include "GeometryIO.h"

// Class constructor 
Geometry::Geometry(){ 
};

Geometry::~Geometry(){ 
};

int Geometry::setNoOfPoints(string &dir){
	int noOfPoints;
	string pointFile;
	pointFile = dir + "/points";

	// assert that the file exists
	assert(fileExists(pointFile) && "Geometry::setNoOfPoints() returns error. Make sure that the file exists.\n");
	// Open points file to know total number of points
	
	ifstream file(pointFile.c_str());
	string line;
	int numberVar = 0; // running variable for indexing purpose
	if(file.good()){
		while(getline(file,line)){
			if (line[0] != '('){
				numberVar ++; // numberVar stops at the first occurance of '('
			}
			else{
				break;
			}
		};
	};
	file.close();
	
	// At this point, we know the line number where '(' starts. The previous line contains total number of points. 
	// If the OpenFOAM file does not contain the total number of points at the begining, a routine needs to compute the total number
	// points.


	// open the file again
	ifstream file1(pointFile.c_str());
	if(file1.good()){
		int localVar = 0;
		while(getline(file1,line)){
			if (localVar == numberVar-1){ //get the line number where total number of points are mentioned
				istringstream input(line);
				input >> noOfPoints;  // transfer that number to the int variable
				break;
			}
			else{
				localVar ++;
			};
		};
	};
	file1.close();
	// information automatically stored for this object
	this->noOfPoints = noOfPoints;
	return noOfPoints;
};


int Geometry::getNoOfPoints()const{
	return this->noOfPoints;
};

int Geometry::setNoOfFaces(string& dir){
	// The working is very similar to this->getNoOfPoints(). Please refer to that code.
	int noOfFaces;

	string faceFile;
	faceFile = dir + "/faces";

	assert(fileExists(faceFile) && "Geometry::setNoOfFaces() returns error. Make sure that the file exists.\n");
	//open points file to know total number of points
	
	ifstream file(faceFile.c_str());
	string line;
	int numberVar = 0;
	if(file.good()){
		while(getline(file,line)){
			if (line[0] != '('){
				numberVar ++;
			}
			else{
				break;
			}
		};
	};
	file.close();
	
	ifstream file1(faceFile.c_str());
	if(file1.good()){
		int localVar = 0;
		while(getline(file1,line)){
			if (localVar == numberVar-1){
				istringstream input(line);
				input >> noOfFaces;
				break;
			}
			else{
				localVar ++;
			};
		};
	};
	file1.close();
	this->noOfFaces = noOfFaces;
	return noOfFaces;
}

int Geometry::getNoOfFaces()const{
	return this->noOfFaces;
};

int Geometry::setNoOfBoundaryConditions(string& dir){
	int noOfBoundaryConditions;
	// The boundary file contains all releveant details regarding boundary conditions
	string boundaryFile;
	boundaryFile = dir + "/boundary";

	assert(fileExists(boundaryFile) && "Geometry::setNoOfBoundaryConditions() returns error. Make sure that the file exists.\n");
	ifstream file(boundaryFile.c_str());
	string line;
	int i,j,noOfPoints,flag=0,jmax=0;
	
	if(file.good()){
		while(getline(file,line)){
			if(isdigit(line[0])){
				flag = 1;
				}
			if(line[0]=='(' && flag==1){
				flag = 0;
				}
				
			if(flag==1){
				istringstream input(line);
				input>>noOfBoundaryConditions;
			}
			
		}
	}
	file.close();

	this->noOfBoundaryConditions = noOfBoundaryConditions;
	return noOfBoundaryConditions;
	
}

int Geometry::getNoOfBoundaryConditions()const{
	return this->noOfBoundaryConditions;
};

int Geometry::setNoOfBoundaryFaces(string& dir){
	int noOfFaces = this->getNoOfFaces(); // total number of faces
	int noOfInternalFaces;
	int noOfBoundaryFaces;

	// Next we find the number of internal faces
	// The working is very similar to this->getNoOfPoints(). Please refer to that code.

	//open neighbour file to know total number of internal faces. This is because only the internal faces have neighbours.
	string neighbourFile;
	neighbourFile = dir + "/neighbour";

	assert(fileExists(neighbourFile) && "Geometry::setNoOfBoundaryFaces() returns error. Make sure that the file exists.\n");
	
	ifstream file(neighbourFile.c_str());
	string line;
	int numberVar = 0;
	if(file.good()){
		while(getline(file,line)){
			if (line[0] != '('){
				numberVar ++;
			}
			else{
				break;
			}
		};
	};
	file.close();
	
	ifstream file1(neighbourFile.c_str());
	if(file1.good()){
		int localVar = 0;
		while(getline(file1,line)){
			if (localVar == numberVar-1){
				istringstream input(line);
				input >> noOfInternalFaces;
				break;
			}
			else{
				localVar ++;
			};
		};
	};
	file1.close();
	noOfBoundaryFaces = noOfFaces - noOfInternalFaces;
	this->noOfInternalFaces = noOfInternalFaces;
	this->noOfBoundaryFaces = noOfBoundaryFaces;
	return noOfBoundaryFaces;
}

int Geometry::getNoOfBoundaryFaces()const{
	return this->noOfBoundaryFaces;
};

int Geometry::getNoOfInternalFaces()const{
	return this->noOfInternalFaces;
};

int Geometry::setNoOfCells(string& dir){
	
	int noOfCells;
	//The owner file contains all possible cells. Since each cell is the owner cell of atleast one face.
	string ownerFile;
	ownerFile = dir + "/owner";
	assert(fileExists(ownerFile) && "Geometry::setNoOfCells() returns error. Make sure that the file exists.\n");

	ifstream file(ownerFile.c_str());
	string line;
	int i,j,noOfPoints,flag=0,jmax=0;
	
	if(file.good()){
		while(getline(file,line)){
			if(line[0]=='('){
				flag = 1;
				}
			if(line[0]==')' && flag==1){
				flag = 0;
				}
				
			if(flag==1){
				istringstream input(line);
				input>>j;
				if(j>jmax){
					jmax = j;
				}
			}
			
		}
	}
	jmax = jmax + 1; //since jmax starts from 0. 
	noOfCells = jmax;
	file.close();
	this->noOfCells = noOfCells;
	return noOfCells;
};

int Geometry::getNoOfCells()const{
	return this->noOfCells;
};

// This function returns the total number of points, cells, faces and boundary faces 
// This sets the noOfCells, noOfFaces, noOfBoundaryFaces and noOfPoints variables corresponding to class Geometry
void Geometry::getProblemData(string& dir){
	int p,c,f,b,bf;
	p = this->setNoOfPoints(dir);
	f = this->setNoOfFaces(dir);
	c = this->setNoOfCells(dir);
	b = this->setNoOfBoundaryFaces(dir); //also sets number of internal faces implicitely
	bf = this->setNoOfBoundaryConditions(dir); 
}



// Function to populate the list of points. Each member of this list is an object of class Point.
// note that when we pass an array to the function as an argument, it gets passed automatically by the pointer to the first element
// We can however, use it as a regular array inside the function as done here.
// Also, the array already has right amount of space since we first find out the number of points before creating the array.
void Geometry::readPointFile(string& pointFile, Point* points){
	
	// assert that the file exists
	assert(fileExists(pointFile) && "Geometry::readPointFile() returns error. Make sure that the file exists.\n");
	ifstream file(pointFile.c_str());
	string line;
	double x,y,z;
	
	int id = 0;
	int count = 0;
	if(file.good()){
		while(getline(file,line)){
			if(line[0]=='(' && line[1]!='\0'){
				line.erase(line.find("("), 1);
				line.erase(line.find(")"),-1);
				istringstream input(line);
				input>>x>>y>>z;
				// p[id] = Point(x,y,z,id); // This is also a valid way of doing this. It creates a new object and copies it into the array p. The created object then gets destroyed as soon as the scope ends. However, the copy remains. This is not recommended for big classes, since creating and destroying is costly, but for a small class such as the Point, it should be ok.
				// Instead we adopt the following way. We know that the array has been created for the class Point.
				// thus, each member is already of the type Point. Here we just set the required information.
				points[id].setPoint(id,x,y,z);   
				id++;
				count ++;
			}
		}
	}
	assert(count == this->getNoOfPoints() && "Geometry::readPointFile() returns error. Make sure that the file is read correctly.\n");
};



// Function to read the face file in OpenFOAM format 
// takes the name of face file and pointer to the array of faces (i.e. each member an object of class Face)
// The function reads the points and then creates a new object of type Face. 
// This object then gets 'copied' (not pass by reference) into the array at the required index number. 
// Thus, the created object gets destroyed as soon as the scope ends, however the copy remains inside the array.
void Geometry::readFaceFile(string& faceFile, Face *faces, Point* points){
	assert(fileExists(faceFile) && "Geometry::readFaceFile() returns error. Make sure that the file exists.\n");

	ifstream file(faceFile.c_str());
	string line;
	int i,j,noOfPoints;
	
	int id = 0;
	int count = 0;
	if(file.good()){
		while(getline(file,line)){
			int noOfPoints;
			if(line[1]=='(' && line[0]!=')'){
				//faces[id] = Face();
				faces[id].setId(id);
				istringstream input0(line); 
				input0 >> noOfPoints;
				line.erase(0,1);
				line.erase(line.find("("), 1);
				line.erase(line.find(")"),-1);
				istringstream input(line);
				for(int i=0;i<noOfPoints;i++){
					input>>j;
					faces[id].addDefiningPoint(&points[j]);
				}
				id++;
				count ++;
				
			}
		}
	}
	
	assert(count == this->getNoOfFaces() && "Geometry::readFaceFile() returns error. Make sure that the file is read correctly.\n");
};

// Function to read the boundary conditions files in OpenFOAM format
// takes the name of the boundary file and appends the array of boundary condition 
void Geometry::readBoundaryFile(string& boundaryFile, BoundaryConditions *bcond){
	assert(fileExists(boundaryFile) && "Geometry::readBoundaryFile() returns error. Make sure that the file exists.\n");

	ifstream file(boundaryFile.c_str());
	string line,word1,word2,word3;
	int i,j,noOfBoundaryConditions;
	
	int id = 0;
	int flag = 0;

	int count = 0;
	if(file.good()){
		while(getline(file,line)){
			
			// The text file ends with ")" everything after 
			// this is ignored
			if(line[0]==')'){
				flag = 0;
			}

			if(line[0]=='(' && flag==0){
				flag = 1;
			}
			else if(line[0]!='(' && flag == 1){
				stringstream input(line);
				input>>word1;
				bcond[id].setId(id);
				bcond[id].setName(word1);
				flag = 2;
			}		
			else if(flag==2){
				stringstream input(line);
				input>>word1;
				if(word1[0]=='{')
					flag=3;
			}
			else if(flag==3){
				stringstream input(line);
				input>>word1;
				if(word1=="type"){
					input>>word2;
					word2.erase(word2.find(";"),-1);
					bcond[id].setType(word2);
				}
				else if(word1=="nFaces"){
					input>>word2;
					word2.erase(word2.find(";"),-1);
					stringstream input1(word2);
					input1>>i;
					bcond[id].setNoOfFaces(i);

				}
				else if(word1=="startFace"){
					input>>word2;
					word2.erase(word2.find(";"),-1);
					stringstream input1(word2);
					input1>>i;
					bcond[id].setStartFace(i);
				}
				else if(word1[0]=='}'){
					flag=1;
					id++;
					count ++;
				}
				else if(word1[0]==')'){
					flag=0;
				}
			}		

			
		}
	}

	assert(count == this->getNoOfBoundaryConditions() && "Geometry::readBoundaryFile() returns error. Make sure that the file is read correctly.\n");
}


// Function to read the owner file in OpenFOAM format 
// The function reads the owners of the faces and simultaneously builds the cells using the faces
void Geometry::readOwnerFile(string& ownerFile, Cell* cells, Face* faces){
	assert(fileExists(ownerFile) && "Geometry::readOwnerFile() returns error. Make sure that the file exists.\n");

	ifstream file(ownerFile.c_str());
	string line;
	int i,j,noOfPoints,flag=0,jmax=0;
	int id = 0;
	int count = 0;
	if(file.good()){
		while(getline(file,line)){

			if(line[0]=='('){
				flag = 1;
				}
			if(line[0]==')' && flag==1){
				flag = 0;
				}
				
			if(flag==1 && line[0]!='('){
				
				istringstream input(line);
				input>>j;
				cells[j].setId(j);
				cells[j].addDefiningFace(&faces[id]);
				id++;
				count ++;

			}
			
		}
	}

	assert(count == this->getNoOfFaces() && "Geometry::readOwnerFile() returns error. Make sure that the file is read correctly.\n");
};

// Function to read the neighbour file in OpenFOAM format 
// takes the name of face file and arrays of cells and faces (pre-initialized) after the owner faces are linked as an input
// The function reads the neighbour of the faces and adds to the list of existing faces of each cell
void Geometry::appendNeighbourFile(string& neighbourFile, Cell* cells, Face* faces){
	assert(fileExists(neighbourFile) && "Geometry::appendNeighbourFile() returns error. Make sure that the file exists.\n");

	ifstream file(neighbourFile.c_str());
	string line;
	int i,j,noOfPoints,flag=0,jmax=0;
	
	int id = 0;
	int count = 0;
	if(file.good()){
		while(getline(file,line)){

			if(line[0]=='('){
				flag = 1;
				}
			if(line[0]==')' && flag==1){
				flag = 0;
				}
				
			if(flag==1 && line[0]!='('){
				istringstream input(line);
				input>>j;
				cells[j].addDefiningFace(&faces[id]);
				id++;
				count ++;
			}
		}
	}
	assert(count == this->getNoOfInternalFaces() && "Geometry::readNeighbourFile() returns error. Make sure that the file is read correctly.\n");
};



// Assign owner cells to faces. All faces have an owner cell. This is read from the owner file generated by blockMesh (or any other OpenFOAM utility
// The function takes an array of pointers to faces and another array of pointers to cells. The implementation is straightforward. 
// the function returns void because we pass arrays by reference
void Geometry::assignOwnerCellsToFaces(string& ownerFile, Cell* cells, Face* faces){
	int i,j,flag=0;
	string line;
	ifstream file(ownerFile.c_str());
	int id = 0;
	int count = 0;
	if(file.good()){
		while(getline(file,line)){
			if(line[0]=='('){
				flag = 1;
				}
			if(line[0]==')' && flag==1){
				flag = 0;
				}
				
			if(flag==1 && line[0]!='('){
				istringstream input(line);
				input>>j;
				faces[id].setOwnerCell(&cells[j]);
				id++;
				count ++;
			}
			
		}
	}
	assert(count == this->getNoOfFaces() && "Geometry::assignOwnerCellsToFaces() returns error. Make sure that the file is read correctly.\n");
};

// Assign neighbour cells to faces. All faces have an owner cell. This is read from the owner file generated by blockMesh (or any other OpenFOAM utility
// The function takes an array of pointers to faces and another array of pointers to cells. The implementation is straightforward. 
// the function returns void because we pass arrays by reference
void Geometry::assignNeighbourCellsToFaces(string& neighbourFile, Cell* cells, Face* faces){
	int i,j,flag=0;
	string line;
	ifstream file(neighbourFile.c_str());
	int id = 0;
	int count = 0;
	if(file.good()){
		while(getline(file,line)){
			if(line[0]=='('){
				flag = 1;
				}
			if(line[0]==')' && flag==1){
				flag = 0;
				}
				
			if(flag==1 && line[0]!='('){
				istringstream input(line);
				input>>j;
				faces[id].setNeighbourCell(&cells[j]);
				id++;
				count ++;
			}
			
		}
	};


	assert(count == this->getNoOfInternalFaces() && "Geometry::assignNeighbourCellsToFaces() returns error. Make sure that the file is read correctly.\n");

	// faces having a neighbour cell have been processed
	// rest of the faces are assigned a 'NULL' cell as the GHOST cell
	// note, we start from 'id' (which is nonzero at this point), i.e. after all the internal faces have been already processed
	for (i=id; i< this->noOfFaces; i++){
		faces[i].setNeighbourCell(NULL);
		count ++;
	};

	assert(count == this->getNoOfFaces() && "Geometry::assignNeighbourCellsToFaces() returns error. Make sure that the file is read correctly.\n");
};


// The following function fills an array inside each cell, such that a member in that array tells whether the same index
// member of definingFaces has this cell as owner or neighbour.
// Array name: Cell::faceRelationshipArray[]
// 
// Example, 
// let cells[i].faceRelationshipArray = [0,0,0,1,1,0]
// This means that definingFaces have [owner, owner, owner, neigh, neigh, owner] relationship with ith cell

void Geometry::fillFaceRelationshipArray(Cell *cells, Face* faces){

	for (int i=0; i<this->noOfCells; i++){
		int noOfFaces = cells[i].getNoOfFaces();
		for (int j=0; j< noOfFaces; j++){
			if (cells[i].getDefiningFace(j)->getOwnerCell()->getId() == cells[i].getId()){
				cells[i].faceRelationshipArray.setValue(j, 0);
			}
			else{
				cells[i].faceRelationshipArray.setValue(j, 1);
			};
		};

		double area = 0;

		for (int j=0; j< noOfFaces; j++){
			if (area < cells[i].getDefiningFace(j)->getArea()){
				area = cells[i].getDefiningFace(j)->getArea();
				cells[i].setLargestFace(cells[i].getDefiningFace(j));
			};
		};
	};
};


// Some of the faces may have the normal vector pointing inside the owner cell. OpenFOAM mesh generation takes care that,
// the face always points away from the owner cell. However, this may not be the case with the third-party mesh generation softwares.
// thus, it is important to find whether face has the normal pointing into the owner or away from it.
// This is indicated by a flag called Face::normalDirFlag;
void Geometry::findFaceNormalsOrientation(Cell *cells,Face *faces){
	for (int f=0; f<this->noOfFaces; f++){
		Face *face = &faces[f];

		// CtoC is the vector starting at cell center, terminating at face center
		TensorO1<double> CtoC(3);

		CtoC.setValue(0, face->getCenter()->getX() - face->getOwnerCell()->getCenter()->getX());
		CtoC.setValue(1, face->getCenter()->getY() - face->getOwnerCell()->getCenter()->getY());
		CtoC.setValue(2, face->getCenter()->getZ() - face->getOwnerCell()->getCenter()->getZ());

		// dot product of CtoC with face->normal

		double dotProd = Math::dot(&CtoC, face->getNormal());

		if (dotProd > 0){
			// face-normal is pointing away from the owner cell
			face->normalDirFlag = 1.0;
		}
		else{
			face->normalDirFlag = -1.0;
		};

	};
};

// Following function computes the shortest distance between two points in the grid. 
// This will be used for computing time-step size in the simulation.

double Geometry::findShortestDistance(Cell *cells){
	this->shortestDistance = LARGE10; //arbitrary large number

	double x1, y1, z1;
	double x2, y2, z2;

	for (int c=0; c<noOfCells; c++){
		Cell *cell = &cells[c];

		for (int v1=0; v1<cell->getNoOfPoints(); v1++){
			// vertex v1 of the cell, it takes value of each of the vertex once
			x1 = cell->getDefiningPoint(v1)->getX();
			y1 = cell->getDefiningPoint(v1)->getY();
			z1 = cell->getDefiningPoint(v1)->getZ();

			for (int v2=v1+1; v2<cell->getNoOfPoints(); v2++){
				// the next vertex, it takes value of only the remaining vertices for distance computing
				x2 = cell->getDefiningPoint(v2)->getX();
				y2 = cell->getDefiningPoint(v2)->getY();
				z2 = cell->getDefiningPoint(v2)->getZ();

				double distance = sqrt(sqr(x2-x1) + sqr(y2-y1) + sqr(z2-z1));

				if (this->shortestDistance > distance){
					this->shortestDistance = distance;
				};
			};
		};
	};

	return this->shortestDistance;
};

double Geometry::getShortestDistance(){
	return this->shortestDistance;
};


double Geometry::findCharacteristicLength(Cell *cells){
	// Find characteristic length 'h' to be used for time-step computation.
	// DT = CFL * h/v , h = characteristic length, v = max(Lambda) i.e. fastest wave-speed, CFL: Courant number, DT: timestep

	double h = LARGE10; // arbitrary large number

	for (int c=0; c<noOfCells; c++){
		Cell *cell = &cells[c];

		double V = 0;
		// cell volume is characterized by the cell Jacobian
		// find max cell J
		for (int QP =0; QP<cell->getNQuad(); QP++){
			if (V< abs(cell->getJacobian()->getValue(QP))){
				V = abs(cell->getJacobian()->getValue(QP));
			};
		};

		double S = NANO; 

		for (int f=0; f<cell->getNoOfFaces(); f++){
			Face *face = cell->getDefiningFace(f);
		 	for (int QP=0; QP<face->getNQuad(); QP++){
				if (S < abs(face->getJacobian()->getValue(QP))){
					S = abs(face->getJacobian()->getValue(QP));
				};
			};
		};

		if (h > V/S){
			h = V/S;
		};
	};

	this->charLength = h;
	return h;
};

double Geometry::getCharacteristicLength(){
	return this->charLength;
};

void Geometry::assignVertexSigns(Cell *cells, Face *faces){

	//* Cell Vertex Signs
	/// The defining points (Vertices) are arranged such that the first four are the part of the base face
	/// And the next four are directly opposite to the first four.
	/// Refer Cell::orderDefiningPoints() for details

	// Note: The vertexSign array has been configured (ie. memory allocation) in Cell::orderDefiningPoints() function, which in turn is called from Cell::setCellType() function.

	for (int c=0; c<noOfCells; c++){
		Cell *cell = &cells[c];
		if (cell->getCellType() == CellType::Hex){
			assert(cell->getVertexSignArray()->getSize0() == cell->getNoOfPoints() and 
			       cell->getVertexSignArray()->getSize1() == 3);

			enum coord{x, y, z};

			// x sign: +ve for vertices 3, 2, 6, 7; else negative	
			cell->setVertexSign(3, x, 1.0);
			cell->setVertexSign(2, x, 1.0);
			cell->setVertexSign(6, x, 1.0);
			cell->setVertexSign(7, x, 1.0);

			cell->setVertexSign(0, x, -1.0);
			cell->setVertexSign(1, x, -1.0);
			cell->setVertexSign(5, x, -1.0);
			cell->setVertexSign(4, x, -1.0);

			// y sign: +ve for vertices 1, 5, 6, 2; else negative
			cell->setVertexSign(1, y, 1.0);
			cell->setVertexSign(5, y, 1.0);
			cell->setVertexSign(6, y, 1.0);
			cell->setVertexSign(2, y, 1.0);

			cell->setVertexSign(0, y, -1.0);
			cell->setVertexSign(4, y, -1.0);
			cell->setVertexSign(7, y, -1.0);
			cell->setVertexSign(3, y, -1.0);

			// z sign: +ve for vertices 4, 5, 6, 7; else negative
			cell->setVertexSign(4, z, 1.0);
			cell->setVertexSign(5, z, 1.0);
			cell->setVertexSign(6, z, 1.0);
			cell->setVertexSign(7, z, 1.0);

			cell->setVertexSign(0, z, -1.0);
			cell->setVertexSign(1, z, -1.0);
			cell->setVertexSign(2, z, -1.0);
			cell->setVertexSign(3, z, -1.0);
		};
	};


	//* Face Vertex signs

	for (int f=0; f<noOfFaces; f++){
		Face *face = &faces[f];
		if (face->getFaceType() == FaceType::Quad){

			// Refer Face::orderDefiningPoints for the sequence of the vertices.
			// point 0: lies in the 3rd quadrant (-ve, -ve)
			// point 1: lies in the 4th quadrant (+ve, -ve)
			// point 2: lies in the 1st quadrant (+ve, +ve)
			// point 3: lies in the 2nd quadrant (-ve, +ve)

			assert(face->getVertexSignArray()->getSize0() == face->getNoOfPoints() and 
			       face->getVertexSignArray()->getSize1() == 2);

			enum coord{x, y};
			// x sign: +ve for vertices 1, 2; else negative	
			face->setVertexSign(0, x, -1.0);
			face->setVertexSign(1, x,  1.0);
			face->setVertexSign(2, x,  1.0);
			face->setVertexSign(3, x, -1.0);

			// y sign: +ve for vertices 2, 3; else negative
			face->setVertexSign(0, y, -1.0);
			face->setVertexSign(1, y, -1.0);
			face->setVertexSign(2, y,  1.0);
			face->setVertexSign(3, y,  1.0);
		};
	};
};



void Geometry::fillDataArrays(string& dir, Point* points, Face* faces, Cell* cells, BoundaryConditions* bcond){

	string pointFile = dir+"/points";
	string faceFile = dir+"/faces";
	string ownerFile = dir+"/owner";
	string neighbourFile = dir+"/neighbour";
	string boundaryFile = dir+"/boundary";

	this->readPointFile(pointFile,points); //populate points array
	this->readFaceFile(faceFile,faces,points); //populate faces array

	for (int i=0; i<this->noOfFaces; i++){
		faces[i].setFaceType(); //refer include/faceTypes.h and Face.cpp
		faces[i].calculateNormal(); //calculate normal, also implicitly calls orderDefiningPoints()
		faces[i].calculateArea();   // and area
		faces[i].calculateCenter();
		faces[i].calculateRotationMatrixParallelToXY(); //Calculates the rotation matrix required to rotate the given face such as to make it parallel to the xy axis, i.e. make normal vector || to z axis. Stores this matrix as Face::rotationMatrixPrallelToXY;
		faces[i].calculateInverseRotationMatrixParallelToXY(); //Calculates the inverse of the rotation matrix calculated above 
		faces[i].calculateRotationMatrixNormal();
		faces[i].calculateInverseRotationMatrixNormal();
	};

	// Addition for boundary conditions (new)
	this->readBoundaryFile(boundaryFile,bcond); //populate boundary array
	
	this->readOwnerFile(ownerFile,cells,faces); //populate cells array and partially fill with definingFaces
	this->appendNeighbourFile(neighbourFile,cells,faces); //complete the list of definingFaces per cell
	for (int i=0; i<this->noOfCells; i++){
		cells[i].findDefiningPoints(); //find defining Points (from associated faces) 
		cells[i].calculateCenter(); // calculate center coordinates
		cells[i].setCellType(); // Set defining point as well as call orderDefiningPoints()  
	};

	this->assignOwnerCellsToFaces(ownerFile, cells, faces); //after this, each face has an owner cell
	this->assignNeighbourCellsToFaces(neighbourFile, cells, faces); //after this each face has a neighbour cell (or GHOST cell)

	this->fillFaceRelationshipArray(cells,faces);		// Largest face is implicitly set for each cell 
	this->findFaceNormalsOrientation(cells,faces);

	this->assignVertexSigns(cells, faces);
};


void Geometry::print(){
	cout << "The total number of points: " << this->noOfPoints << endl;
	cout << "The total number of cells: " << this->noOfCells << endl;
	cout << "The total number of faces: " << this->noOfFaces << endl;
	cout << "The total number of internal faces: " << this->noOfInternalFaces << endl;
	cout << "The total number of boundary faces: " << this->noOfBoundaryFaces << endl;
	
	// Addition after boundary conditions 
	cout << "The total number of boundary conditions: " << this->noOfBoundaryConditions << endl;
};

/*******************************************************************************************************/

namespace Test{

	void TestGeometryIO(){
		cout << "\nRunning tests on test cases stored in app/Tests/  ...\n";



		string dirName = "app/Tests/1x1x1/constant/polyMesh";
		Geometry geo;
		geo.getProblemData(dirName);

		assert(geo.getNoOfPoints() == 8 and geo.getNoOfCells() == 1 and geo.getNoOfBoundaryFaces() == 6 and geo.getNoOfFaces() == 6 and geo.getNoOfInternalFaces() == 0 and geo.getNoOfBoundaryConditions() == 1);

		dirName = "app/Tests/2x1x1/constant/polyMesh";
		geo.getProblemData(dirName);

		assert(geo.getNoOfPoints() == 12 and geo.getNoOfCells() == 2 and geo.getNoOfBoundaryFaces() == 10 and geo.getNoOfFaces() == 11 and geo.getNoOfInternalFaces() == 1 and geo.getNoOfBoundaryConditions() == 3);

		dirName = "app/Tests/5x5x5/constant/polyMesh";
		geo.getProblemData(dirName);

		assert(geo.getNoOfPoints() == 216 and geo.getNoOfCells() == 125 and geo.getNoOfBoundaryFaces() == 150 and geo.getNoOfFaces() == 450 and geo.getNoOfInternalFaces() == 300 and geo.getNoOfBoundaryConditions() == 1);

		dirName = "app/Tests/pitzDaily/constant/polyMesh";
		geo.getProblemData(dirName);

		assert(geo.getNoOfPoints() ==25012 and geo.getNoOfCells() == 12225 and geo.getNoOfBoundaryFaces() == 49180-24170 and geo.getNoOfFaces() == 49180 and geo.getNoOfInternalFaces() == 24170 and geo.getNoOfBoundaryConditions() == 5);

		Point *points;
		points = new Point[25012];

		Cell *cells;
		cells = new Cell[12225];

		Face *faces;
		faces = new Face[49180];

		BoundaryConditions* bcond;
		bcond = new BoundaryConditions[5];

		geo.fillDataArrays(dirName, points, faces, cells, bcond);

		enum coord{x,y,z};
		for (Index c = 0; c< geo.getNoOfCells(); c++){
			assert(cells[c].getCellType() == CellType::Hex);
			assert(cells[c].getVertexSign(0,x)==-1.0 && cells[c].getVertexSign(0,y)==-1.0 && cells[c].getVertexSign(0,z)==-1.0);
			assert(cells[c].getVertexSign(1,x)==-1.0 && cells[c].getVertexSign(1,y)== 1.0 && cells[c].getVertexSign(1,z)==-1.0);
			assert(cells[c].getVertexSign(2,x)== 1.0 && cells[c].getVertexSign(2,y)== 1.0 && cells[c].getVertexSign(2,z)==-1.0);
			assert(cells[c].getVertexSign(3,x)== 1.0 && cells[c].getVertexSign(3,y)==-1.0 && cells[c].getVertexSign(3,z)==-1.0);
			assert(cells[c].getVertexSign(4,x)==-1.0 && cells[c].getVertexSign(4,y)==-1.0 && cells[c].getVertexSign(4,z)== 1.0);
			assert(cells[c].getVertexSign(5,x)==-1.0 && cells[c].getVertexSign(5,y)== 1.0 && cells[c].getVertexSign(5,z)== 1.0);
			assert(cells[c].getVertexSign(6,x)== 1.0 && cells[c].getVertexSign(6,y)== 1.0 && cells[c].getVertexSign(6,z)== 1.0);
			assert(cells[c].getVertexSign(7,x)== 1.0 && cells[c].getVertexSign(7,y)==-1.0 && cells[c].getVertexSign(7,z)== 1.0);
		};

		delete[] bcond;
		delete[] faces;
		delete[] cells;
		delete[] points;





		cout <<"Test::TestGeometryIO() passed.\n";

	};

};
