#include "DG.h"

// This (large) file is divided in several sections as follows:
// 1. constructor, destructor etc.
// 2. Data reading and writing as well as memory allocations, domain related functions.
// 3. Cell related functions
// 4. Face related functions
// 5. Flux and residual related functions
// 6. Time integrator related functions
// Kindly navigate to access the corresponding functions accordingly.


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 1. Construction ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

// Class constructor 
DG::DG(string path, int order, IntFlag::intflag intType){ 
	//intType: 0:inexact integration, 1:exact integration

	this->path = path;
	this->order = order;
	this->intType = intType;

	noOfVariable = 0;
	totNoOfVariable = 0;

	// set initial flags to False, since no memory is explicitly allocated so far.
	// These flags are set to true when memory is allocated using the 'new' method.
	this->VariableNameFlag = false;
	this->VariableSizeFlag = false;
	this->cummulativeVariableSizeFlag = false;
	this->pointsFlag = false;
	this->facesFlag = false;
	this->cellsFlag = false;
	this->bcondFlag = false;
	this->laplaceFlag = false;
	this->variableFlag = false;
	this->refVariableFlag = false;
};

// Class destructor 
DG::~DG(void){
	// Memory needs to be freed when the object gets destroyed. 
	// This will happen when the code ends, the object gets out of the local scope or object is manually deleted. 
	if (VariableNameFlag == true){
		delete[] VariableName;
	};
	if (pointsFlag == true){
		delete[] points;
	};
	if (facesFlag == true){
		delete[] faces;
	};
	if (cellsFlag == true){
		delete[] cells;
	};
	if (bcondFlag == true){
		delete[] bcond;
	};
}; 



//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 2. Memory and domain related functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//



// Setter function for the order of interpolation scheme
void DG::setOrder(int order){
	this->order = order;
};

// Getter function for the order of interpolation scheme
int DG::getOrder()const{
	return this->order;	
};

void DG::setPath(string path){
	this->path = path;
};

string DG::getPath()const{
	return this->path;
};

void DG::setupReferenceCellsArray(){
	// refer RefCell.cpp and cellTypes.h for details
	referenceCells[0].init(order, intType, CellType::Hex );
	referenceCells[1].init(order, intType, CellType::Tet );
	referenceCells[2].init(order, intType, CellType::Prism );
	referenceCells[3].init(order, intType, CellType::Pyramid );
};

RefCell* DG::getReferenceCellsArray(){
	return this->referenceCells;
};

void DG::setupReferenceFacesArray(){
	// refer refCell.cpp and cellTypes.h for details

	referenceFaces[0].init(order, intType, FaceType::Quad );
	referenceFaces[1].init(order, intType, FaceType::Tri );
};

RefFace* DG::getReferenceFacesArray(){
	return this->referenceFaces;
};

void DG::setNoOfVariable(int noOfVariable){
	this->noOfVariable = noOfVariable;

	// stores a list of variable names, such as p, U, rho ...
	VariableName = new string[noOfVariable];
	this->VariableNameFlag = true;

	// Stores a list of what dim each variable is. So 1 for p, 3 for U etc.
	// Example:
	// If VariableName = [rho, U, p], then VariableSize = [1, 3, 1]
	VariableSize.setSize(noOfVariable);
	this->VariableSizeFlag = true;

	// Stores the start position of the next variable in the global array "variable"
	// so cummulativeVariableSize for above example will be :  [0,1,4] 
	// i.e. rho starts at 0, U starts at 1 (because dim of rho is 1); and p starts at 4 (because dim of U is 3) 
	// refer setVariableSize() function for more details
	cummulativeVariableSize.setSize(noOfVariable);
	this->cummulativeVariableSizeFlag = true;
};
	

// Getter function for number of Variables
int DG::getNoOfVariable(){
	return noOfVariable;	
};


// Setter function for names of Variables
void DG::setVariableName(string VariableName[]){
	for(int i=0;i<noOfVariable;i++)
		this->VariableName[i] = VariableName[i];
};
	
// Setter function for the sizes of Variables
// Takes the array of Variable Sizes as an input
void DG::setVariableSize(int VariableSize[]){
	this->cummulativeVariableSize.setValue(0, 0);
	for(int i=0;i<noOfVariable;i++){
		this->VariableSize.setValue(i, VariableSize[i]);
	};

	for(int i=1;i<noOfVariable;i++){
		this->cummulativeVariableSize.setValue(i, this->cummulativeVariableSize.getValue(i-1) + this->VariableSize.getValue(i-1));
	};
};


// Total Setter function for all functional details
void DG::setFunctionalDetails(int order, IntFlag::intflag intType, Solver::solver sol, System::system sys){
	this->order = order;
	this->intType = intType;
	this->solverType = sol;
	this->systemType = sys;
};

void DG::assignSystemOfEquations(){
	// refer src/internalFluxSolver.cpp for implementation of each of the fluxes
	switch (this->systemType)
	{
		case System::scalarLinear:
			{
				this->findInternalFlux = scalarTransportFlux;
				break;
			};
		case System::scalarBurger:
			{
				this->findInternalFlux = BurgerFlux;
				break;
			};
		case System::Heat:
			{
				this->findInternalFlux = HeatEquationFlux;
				break;
			};
		case System::Euler:
			{	
				this->findInternalFlux = EulerFlux;
				break;
			};
		case System::NavierStokes:
			{	
				this->findInternalFlux = NavierStokesFlux;
				break;
			};
	};
};

void DG::assignRiemannSolver(){
	// refer src/RiemannSolver.cpp for implementation details of each of the Riemann solvers
	switch (this->solverType)
	{
		case Solver::scalarTransportSolver:
			{
				this->solveRiemannProblem = scalarRiemannSolver;
				break;
			};
		case Solver::BurgerSolver:
			{
				this->solveRiemannProblem = BurgerRiemannSolver;
				break;
			};
		case Solver::LaplaceSolver:
			{	
				this->solveRiemannProblem = centralDifferenceHeatEquation;
				break;
			};
		case Solver::Roe:
			{	
				this->solveRiemannProblem = RoeSolver;
				break;
			};
		case Solver::LLF:
			{	
				this->solveRiemannProblem = LLFSolver;
				break;
			};
	};
};

// Total Printer function for all functional details
void DG::printFunctionalDetails(){
	cout<<"~~~~~~~~~~~~~~~~~~Functional Details~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
	cout<<"Order of the polynomial: "<<getOrder()<<endl;
	cout<<"Integration type (exact / inexact): " << this->getQuadratureType() << endl; 
}

// Setter function for the time step
void DG::setDeltaT(double deltaT){
	this->deltaT = deltaT;
};

// Getter function for the time step
double DG::getDeltaT(){
	return deltaT;	
};

// Setter function for the start time
void DG::setStartTime(double startTime){
	this->startTime = startTime;
};

// Getter function for the start time
double DG::getStartTime(){
	return startTime;	
};

// Setter function for the end time
void DG::setEndTime(double endTime){
	this->endTime = endTime;
};

// Getter function for the end time
double DG::getEndTime(){
	return endTime;	
};

// Setter function for the Print time
void DG::setPrintTime(double printTime){
	this->printTime = printTime;
};

// Getter function for the Print time
double DG::getPrintTime(){
	return printTime;	
};

// Setter function for integrator Type
void DG::setIntegratorType(int integratorType){
	this->integratorType = integratorType;
};

// Getter function for integrator Type
string DG::getIntegratorType(){
	string integratorTypeName[5] = {"Euler", "1st Order Rk Method", "2nd Order Rk Method", "3rd Order Rk Method", "4th Order Rk Method"};
	return integratorTypeName[integratorType];
};

// Getter function for quadrature Type
string DG::getQuadratureType(){
	string quadratureTypeName[2] = {"Inexact", "Exact"};
	return quadratureTypeName[intType];
};


void DG::setTotNoOfVariableTimeInt(int totNoOfVariableTimeInt){
	this->totNoOfVariableTimeInt = totNoOfVariableTimeInt;
};

int DG::getTotNoOfVariables(){
	return this->totNoOfVariable;
};

// Total Setter function for all temporal details
void DG::setTemporalDetails(double deltaT, double startTime, double endTime, double printTime, int integratorType){
	this->deltaT = deltaT;
	this->startTime = startTime;
	this->endTime = endTime;
	this->printTime = printTime;
	this->integratorType = integratorType;
}

// Total Printer function for all temporal details
void DG::printTemporalDetails(){
	cout<<"~~~~~~~~~~~~~~~~~~~Temporal Details~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
	cout<<"Time Step: "<<getDeltaT()<<endl;
	cout<<"Start Time: "<<getStartTime()<<endl;
	cout<<"End Time: "<<getEndTime()<<endl;
	cout<<"Data Printing Time: "<<getPrintTime()<<endl;
	cout<<"Integrator: "<<getIntegratorType()<<endl;
}
						

// Adds all the variable sizes
// Scalar has a single value where as vector may have more than
// Used subsequently for allocating space to the property Matrix
// Number of variables assuming all are scalar only
void DG::addTotNoOfVariable(){
	for(int i=0;i<noOfVariable;i++)
		this->totNoOfVariable += VariableSize.getValue(i);
	// By default all variables need time integration
	this->totNoOfVariableTimeInt = totNoOfVariable;	
}
	
// Allocates size to the variable matrix
// Row denotes : number of total property assuming all are scalar
// Coloumn denotes: number of cells
void DG::allocateVariableArraySize(){
	// Domain variable array
	this->variable.setSize(totNoOfVariable, Domain.getNoOfCells());
	assert(this->variable.getSize0() == this->totNoOfVariable and this->variable.getSize1() == Domain.getNoOfCells());

	// reference variable array
	this->variableRef.setSize(totNoOfVariable, Domain.getNoOfCells());
	assert(this->variableRef.getSize0() == this->totNoOfVariable and this->variableRef.getSize1() == Domain.getNoOfCells());

	// Initialize both the arrays with zero value
	for (Index var=0; var<totNoOfVariable; var++){
		for (Index cell=0; cell<Domain.getNoOfCells(); cell++){
			variable.setValue(var, cell, 0.0);
			variableRef.setValue(var, cell, 0.0);
		};
	};
	// Set the flag to true
	this->variableFlag = true;

	// Set the size of the Sod's solution
	// 8 variables (rho, u, v, w, p, T, a, M) for 4 locations 
	this->SodAnalytical.setSize(8, 4);
};
	
/// Allocates size to the variable type and value matrix for boundary conditions
void DG::allocateBoundaryVariableArraySize(){
	for(int i = 0; i < noOfBoundaryConditions; i++){
		this->bcond[i].getVariableTypeArray()->setSize(totNoOfVariable);
		this->bcond[i].getVariableValueArray()->setSize(totNoOfVariable);
		// Initialize all the value to zero
		this->bcond[i].getVariableValueArray()->setAll(0.0);
		this->bcond[i].func2bcond.resize(totNoOfVariable);
		this->bcond[i].getFunc2bcondTypeArray()->setSize(totNoOfVariable);
		this->bcond[i].getFunc2bcondTypeArray()->setAll("NoValue");
	};
}


// Reads all the files in folder 0
// Uniform and Non-Uniform Input as well
// Non-uniform example in tutorials>> incompressible >> shallowWaterFoam >> squareBump >> 0
void DG::readVariableArray(double Time){
	string timePath,fileName,word,input;
	string time;
	double value;
	time = varToString(Time); //double Time converted into a string. Refer mathAndOtherFunctions.cpp 
	timePath = path + "/"+time+"/";
	cout << "Reading variable array from folder " << timePath << "\n\n" <<  endl;

	
	for (int i=0;i<noOfVariable;i++){
		int vectorFlag = 0;
		fileName = timePath + VariableName[i];	// Readin the variable file in timePath folder
		ifstream file(fileName.c_str());
		cout << "Reading data from file named: " << fileName << endl;
		//Making sure the file exists
		assert (file.good() && "The above mentioned file (i.e. Input file) not found. Please check again.");
		string line;
		int flagBCond = 0;
		int bcondId = 0;
		if(file.good()){	
			// Looping over lines
			while(getline(file,line)){
				istringstream stringline(line);		// line input to stringline (a stringstream type object)
				stringline >> word;			// word-by-word input to string word
				if(word == "class"){
					stringline >> word;
					if(word=="volVectorField;"){
						vectorFlag = 1;
					}
				}
				if(word == "internalField"){
					stringline >> word;
					if(word=="uniform"){
						if(vectorFlag == 0){
							// means a scalar field
							stringline >> value;
							for (int j=0;j<Domain.getNoOfCells();j++){
								variable.setValue(cummulativeVariableSize.getValue(i),j, value);
							}
						}
						else if(vectorFlag == 1){
							// Meaning a vector in the following format
							// (number1 number2 ... numberN)
							// This requires erasing brackets for the first and the last numbers.
							//
							// Processing the first number
							stringline >> word;
							word = word.erase(0,1); //erase first bracket
							istringstream stringword(word);
							stringword>>value;
							for (int j=0;j<Domain.getNoOfCells();j++){
								variable.setValue(cummulativeVariableSize.getValue(i),j, value);
							}
							// Processing intermediate components of the vector. 
							// notice the -1 in the upper limit of the loop. Indicates that it doesn't process
							// the last member
							for(int k=1;k<VariableSize.getValue(i)-1;k++){
								stringline >> word;
								istringstream stringword1(word);
								stringword1>>value;
								for (int j=0;j<Domain.getNoOfCells();j++){
									variable.setValue(cummulativeVariableSize.getValue(i)+k, j, value);
								}
							}	
							// Processing the last member of the array.
							stringline >> word;
							word = word.erase(word.size()-2); //-2 because of ';' and ')' 
							istringstream stringword2(word);
							stringword2>>value;
							for (int j=0;j<Domain.getNoOfCells();j++){
								//variable[cummulativeVariableSize[i]+VariableSize[i]-1][j] = value;
								int VAR = cummulativeVariableSize.getValue(i)+VariableSize.getValue(i)-1;
								variable.setValue(VAR, j, value);
							}
						}
					}
					else if(word == "nonuniform"){
						if(vectorFlag == 0){
							getline(file,line);
							getline(file,line);
							for (int j=0;j<Domain.getNoOfCells();j++){
								getline(file,line);
								istringstream stringline(line);
								stringline>>value;
								variable.setValue(cummulativeVariableSize.getValue(i), j, value);
							}
						}
						else if(vectorFlag == 1){
							getline(file,line);
							getline(file,line);
							for (int j=0;j<Domain.getNoOfCells();j++){
								getline(file,line);
								// First member of the vector
								istringstream stringline(line);
								stringline >> word;
								word = word.erase(0,1);
								istringstream stringword(word);
								stringword>>value;
								variable.setValue(cummulativeVariableSize.getValue(i), j, value);
								// intermediate vectors
								for(int k=1;k<VariableSize.getValue(i)-1;k++){
									stringline >> word;
									istringstream stringword1(word);
									stringword1>>value;
									variable.setValue(cummulativeVariableSize.getValue(i)+k, j, value);
								}	
								// Last member of the vector
								stringline >> word;
								word = word.erase(word.size()-1); //-1 because no ';' in this case. Only ')'
								istringstream stringword2(word);
								stringword2 >> value;
								//variable[cummulativeVariableSize[i]+VariableSize[i]-1][j] = value;
								int VAR1 = cummulativeVariableSize.getValue(i)+VariableSize.getValue(i)-1;
								variable.setValue(VAR1, j, value);
							};
						};
					};
				};
				
				if(word == "boundaryField"){
					flagBCond = 1;
				}
				// Line has already reached Boudnary Field in files
				else if(flagBCond==1 && word[0]=='{'){
					flagBCond = 2;
					// "{" character starts the Boundary Fields
				}	
				else if(flagBCond==2 && word[0]!='{' && word[0]!='}' && word[0]!='/'){
					// the word denotes the boundary name that has been defined in polyMesh/boundary file
					for (int bcondi=0;bcondi<noOfBoundaryConditions;bcondi++){
						// to check if the this word matches the "name" keyword stored for each boundary array element
						if(word==bcond[bcondi].getName()){
							bcondId = bcondi;				// bcondId denotes the Id of boundary array element
						}
					}
					flagBCond = 3;
				}
				else if(flagBCond==3 && word[0]=='{'){		// details of the boundary conditions for each boundary element starts from here
					flagBCond = 4;
				}
				else if(flagBCond == 4){
					if(word=="type"){
						stringline >> word;			
						word.erase(word.find(";"),-1);
						// Loop in case the variable is a vector
						for(int k=0;k<VariableSize.getValue(i);k++){
							bcond[bcondId].setVariableType(cummulativeVariableSize.getValue(i)+k, word);		
						}
						// No separate boundary conditions given for each variable.
						// The overall boundary condition is given based on the boundary type specified in constant/polymesh/boundary file
						// This routine just sets the flags for fixedValue boundary condition, for which we need to read the value from a file
						if(word=="zeroGradient"){
							flagBCond = 5;		
							for(int k=0;k<VariableSize.getValue(i);k++){
								if (this->systemType == System::scalarLinear){
									bcond[bcondId].func2bcond[cummulativeVariableSize.getValue(i)+k] = &BoundaryConditions::zeroGradientScalar;
									bcond[bcondId].getFunc2bcondTypeArray()->setValue(cummulativeVariableSize.getValue(i)+k, "zeroGradientScalar");
								}
								else if (this->systemType == System::Euler || System::NavierStokes){
									bcond[bcondId].func2bcond[cummulativeVariableSize.getValue(i)+k] = &BoundaryConditions::zeroGradientEuler;
									bcond[bcondId].getFunc2bcondTypeArray()->setValue(cummulativeVariableSize.getValue(i)+k, "zeroGradientEuler");
								}
								else if (this->systemType == System::Heat){
									bcond[bcondId].func2bcond[cummulativeVariableSize.getValue(i)+k] = &BoundaryConditions::zeroGradientHeat;
									bcond[bcondId].getFunc2bcondTypeArray()->setValue(cummulativeVariableSize.getValue(i)+k, "zeroGradientHeat");
								}
							}
						}
						else if(word=="empty"){
							flagBCond = 5;		
							for(int k=0;k<VariableSize.getValue(i);k++){
								if (this->systemType == System::scalarLinear){
									bcond[bcondId].func2bcond[cummulativeVariableSize.getValue(i)+k] = &BoundaryConditions::zeroGradientScalar;
									bcond[bcondId].getFunc2bcondTypeArray()->setValue(cummulativeVariableSize.getValue(i)+k, "zeroGradientScalar");
								}
								else if (this->systemType == System::Euler || System::NavierStokes){
									bcond[bcondId].func2bcond[cummulativeVariableSize.getValue(i)+k] = &BoundaryConditions::emptyEuler;
									bcond[bcondId].getFunc2bcondTypeArray()->setValue(cummulativeVariableSize.getValue(i)+k, "emptyEuler");
								}
							}
						}	
						else if(word=="fixedValue"){
							flagBCond = 6; 	// flag set for reading the value from a file

							for(int k=0;k<VariableSize.getValue(i);k++){
								if (this->systemType == System::scalarLinear){
									bcond[bcondId].func2bcond[cummulativeVariableSize.getValue(i)+k] = &BoundaryConditions::fixedValueScalar;
									bcond[bcondId].getFunc2bcondTypeArray()->setValue(cummulativeVariableSize.getValue(i)+k, "fixedValueScalar");
								}
								else if (this->systemType == System::Euler || System::NavierStokes){
									bcond[bcondId].func2bcond[cummulativeVariableSize.getValue(i)+k] = &BoundaryConditions::fixedValueEuler;
									bcond[bcondId].getFunc2bcondTypeArray()->setValue(cummulativeVariableSize.getValue(i)+k, "fixedValueEuler");
								}
								else if (this->systemType == System::Heat){
									bcond[bcondId].func2bcond[cummulativeVariableSize.getValue(i)+k] = &BoundaryConditions::fixedValueHeat;
									bcond[bcondId].getFunc2bcondTypeArray()->setValue(cummulativeVariableSize.getValue(i)+k, "fixedValueHeat");
								}

							}
						}	
						else if(word=="noSlip"){
							flagBCond = 5; 		
							for(int k=0;k<VariableSize.getValue(i);k++){
								if (this->systemType == System::scalarLinear){
									bcond[bcondId].func2bcond[cummulativeVariableSize.getValue(i)+k] = &BoundaryConditions::zeroGradientScalar;
									bcond[bcondId].getFunc2bcondTypeArray()->setValue(cummulativeVariableSize.getValue(i)+k, "zeroGradientScalar");
								}
								else if (this->systemType == System::Euler || System::NavierStokes){
									bcond[bcondId].func2bcond[cummulativeVariableSize.getValue(i)+k] = &BoundaryConditions::noSlipEuler;
									bcond[bcondId].getFunc2bcondTypeArray()->setValue(cummulativeVariableSize.getValue(i)+k, "noSlipEuler");
								}
							}
						}	
					}
				}
				// Fixed Value input or Dirichlet Boundary Conditions
				else if(flagBCond==6 && word=="value"){
					stringline >> word;
					if(word=="uniform"){
						if(vectorFlag==0){
							stringline >> value;
							// Set value of the variable
							bcond[bcondId].setVariableValue(cummulativeVariableSize.getValue(i), value);
						}
						else if(vectorFlag==1){
							
							// Meaning a vector in the following format
							// (number1 number2 ... numberN)
							// This requires erasing brackets for the first and the last numbers.
							//
							// Processing the first number
							stringline >> word;
							word = word.erase(0,1); //erase first bracket
							istringstream stringword(word);
							stringword>>value;
							bcond[bcondId].setVariableValue(cummulativeVariableSize.getValue(i), value);
							// Processing intermediate components of the vector. 
							// notice the -1 in the upper limit of the loop. Indicates that it doesn't process
							// the last member
							for(int k=1;k<VariableSize.getValue(i)-1;k++){
								stringline >> word;
								istringstream stringword1(word);
								stringword1>>value;
								//bcond[bcondId].variableValue[cummulativeVariableSize[i]+k] = value;
								bcond[bcondId].setVariableValue(cummulativeVariableSize.getValue(i) + k, value);
							}	
							// Processing the last member of the array.
							stringline >> word;
							word = word.erase(word.size()-2); //-2 because of ';' and ')' 
							istringstream stringword2(word);
							stringword2>>value;
							//bcond[bcondId].variableValue[cummulativeVariableSize[i]+VariableSize[i]-1] = value;
							int VAR2 = cummulativeVariableSize.getValue(i)+VariableSize.getValue(i)-1;
							bcond[bcondId].setVariableValue(VAR2, value);

						};
					
					}
					else if (word=="nonuniform"){
						// Add relevant file input code segment in
						// case the boundary condition is also non-unform
						// will have to add face element in BoundaryConditions class
						// to store face wise non uniform value for different variable
						// Not needed in short time and for most applications
					};
					
					flagBCond = 5;
				}
				else if(flagBCond==5 && word[0]=='}'){
					flagBCond = 2;
				};
			};
		};
		file.close();
	};
};
	

// Assigns the variable values of the different properties 
// to all the DOF points of the cell
void DG::assignVariabletoCell(){
	int nDOF,Np;

	// these arrays are used for Euler and NS equations
	TensorO1<double> primitiveVariableVector(totNoOfVariable);
	TensorO1<double> conservedVariableVector(totNoOfVariable);

	for(Index ncell=0;ncell<this->noOfCells;ncell++){

		//first reference variable value is stored in variableRef
		for(int nvar = 0; nvar < totNoOfVariable; ++nvar){
			//this->variableRef[j][i] = this->variable[j][i]; // only for the scalar solver
			this->variableRef.setValue(nvar, ncell, this->variable.getValue(nvar, ncell)); // only for the scalar solver
		};
		
		// Number of DOF points for different cells
		if(cells[ncell].getCellType()==CellType::Hex){
			Np = this->order+1;
			nDOF = pow((Np),3);  //number of DOF points
		}
		else if(cells[ncell].getCellType()==CellType::Tet){
			nDOF = (order+1)*(order+2)*(order+3)/6;  //number of DOF points
		}
		else if(cells[ncell].getCellType()==CellType::Prism){
			nDOF = (order+1)*(order+2)*(order+1)/2;  //number of DOF points
		}
		else if(cells[ncell].getCellType()==CellType::Pyramid){
		}
	
		// Allocating space to variable 
		cells[ncell].getVariableArray()->setSize(integratorType, totNoOfVariable, nDOF);

		// Allocating space to variableRef
		cells[ncell].getVariableRefArray()->setSize(totNoOfVariable, nDOF);
		cells[ncell].getVariableRefArray()->setAll(0.0);
		
		// Assigning the same value at nDOF points as the centre value
		// Preferable because interpolated values might smoothen the shock structures

		for(int nvar = 0; nvar < totNoOfVariable; ++nvar){
			for(int DOF = 0; DOF < nDOF; ++DOF){
				if (abs(this->variable.getValue(nvar,ncell)) > SMALL10){
					cells[ncell].getVariableArray()->setValue(0,nvar,DOF, this->variable.getValue(nvar,ncell));
				}
				else{
					cells[ncell].getVariableArray()->setValue(0,nvar,DOF, 0.0);
				};
			}
		}

		// Following part converts the primitive variable vector to the conserved variable vector 
		// This is then overwritten on the cell.variable values (which are set above)
		if (this->systemType == System::Euler or this->systemType == System::NavierStokes ){
			// first extract the primitive variable vector from the IC 
			for(int nvar = 0; nvar < totNoOfVariable; ++nvar){  // 5 variables for Euler equations
				primitiveVariableVector.setValue(nvar, this->variable.getValue(nvar,ncell));
			};

			// then get a conserved variable vector from this
			// refer src/gasdynamics.cpp for implementation
			getConservedVariableVector(&primitiveVariableVector, &conservedVariableVector);
			
			// If the system == NavierStokes
			// Quantities from 5 onwards are derivative of u, v and z
			// the getConservedVariableVector() function only deals with the first 5 variables
			// rest of the quantities are same for both conserved and primitive form
			if(this->systemType == System::NavierStokes){
				for(int nvar = 5; nvar < totNoOfVariable; ++nvar){
					conservedVariableVector.setValue(nvar, primitiveVariableVector.getValue(nvar));
				}
			}


			// then copy this value to each DOF location
			for(int DOF = 0; DOF < nDOF; ++DOF){
				for(int nvar = 0; nvar < totNoOfVariable; ++nvar){  // 5 variables for Euler equations
					cells[ncell].getVariableArray()->setValue(0,nvar,DOF, conservedVariableVector.getValue(nvar));
				};
			};
		};

		// also allocate space for flux vectors in cells 
		cells[ncell].getFluxVectorx()->setSize(totNoOfVariable, nDOF);
		cells[ncell].getFluxVectory()->setSize(totNoOfVariable, nDOF);
		cells[ncell].getFluxVectorz()->setSize(totNoOfVariable, nDOF);

		// initialize to zero
		for (Index nvar=0; nvar<totNoOfVariable; nvar++){
			for (Index DOF=0; DOF< nDOF; DOF++){
				cells[ncell].getFluxVectorx()->setValue(nvar, DOF, 0.0);
				cells[ncell].getFluxVectory()->setValue(nvar, DOF, 0.0);
				cells[ncell].getFluxVectorz()->setValue(nvar, DOF, 0.0);
			};
		};


		// As well as Riemann flux vector from all the faces in cell
		cells[ncell].getFluxStarVector()->setSize(6, totNoOfVariable, nDOF); 
		//so cells[i].fluxStarVector[f][k][j] indiacates flux from fth face for kth variable at jth DOF location

		// initialization
		for (Index f=0; f<6; f++){
			for (Index nvar=0; nvar<totNoOfVariable; nvar++){
				for (Index DOF=0; DOF<nDOF; DOF++){
					cells[ncell].getFluxStarVector()->setValue(f,nvar,DOF,  0.0);
				};
			};
		};

		// In addition, we need to assign the variable Residual vector to the cell
		cells[ncell].getVariableResidualArray()->setSize(integratorType, totNoOfVariable, nDOF);

		// initialized to zero
		for (Index nIntType=0; nIntType< integratorType ; nIntType++){
			for (Index nvar=0; nvar<totNoOfVariable; nvar++){
				for (Index DOF=0; DOF<nDOF; DOF++){
					cells[ncell].getVariableResidualArray()->setValue(nIntType,nvar,DOF, 0.0);
				};
			};
		};
	};
};


// Assigns the variable values of the different properties 
// to all the DOF points of the cell
void DG::integrateCellVariable(){
	int nDOF,Np;
	FunctionalSpace F(this->order, this->intType);

	// these arrays are used for Euler and NS equations
	TensorO1<double> primitiveVariableVector(5);
	TensorO1<double> conservedVariableVector(5);
	TensorO1<double> localVariableVector(5);

	TensorO1<double> x;
	TensorO1<double> y;
	TensorO1<double> z;
	TensorO1<double> w;
	for(Index ncell=0; ncell<this->noOfCells; ncell++){
		
		// Number of DOF points for different cells
		if(cells[ncell].getCellType() == CellType::Hex){
			Np = this->order+1;
			nDOF = pow((Np),3);  //number of DOF points
			x.setSize(nDOF);
			y.setSize(nDOF);
			z.setSize(nDOF);
			w.setSize(nDOF);
			F.LGLRootsAndWeights3D(Np, Np, Np, &x, &y, &z, &w);
		}
		else if(cells[ncell].getCellType() == CellType::Tet){
		}
		else if(cells[ncell].getCellType() == CellType::Prism){
		}
		else if(cells[ncell].getCellType() == CellType::Pyramid){
		}
	
		double sum_w=0.0;
		for (int DOF = 0; DOF < nDOF; DOF++){
			sum_w += w.getValue(DOF);
		};
			
		// Assigning the same value at nDOF points as the centre value
		// Preferable because interpolated values might smoothen the shock structures

		// first get the cell averaged values of the conserved variable vector.
		for(int nvar = 0; nvar < totNoOfVariable; ++nvar){
			this->variable.setValue(nvar, ncell, 0.0);
			for(int DOF = 0; DOF < nDOF; ++DOF){
				//this->variable[nvar][ncell] += cells[ncell].variable[0][nvar][DOF]*w[DOF];
				double value =  cells[ncell].getVariable(0,nvar,DOF) * w.getValue(DOF);
				this->variable.addValue(nvar,ncell, value);
			}
			double value = this->variable.getValue(nvar, ncell);
			this->variable.setValue(nvar,ncell, value/ sum_w);
		}
		if (this->systemType == System::Euler or this->systemType == System::NavierStokes){
			// Euler and Navier Stokes equations
			// Extract the primitive variable vector
			for(int nvar = 0; nvar < totNoOfVariable; ++nvar){
				conservedVariableVector.setValue(nvar, this->variable.getValue(nvar,ncell));
			};

			getPrimitiveVariableVector(&conservedVariableVector, &primitiveVariableVector);

			// copy it to variable vector
			for(int nvar = 0; nvar < totNoOfVariable; ++nvar){
				this->variable.setValue(nvar,ncell, primitiveVariableVector.getValue(nvar));
			};
		};
	};
};

// Assigns the flux values of the different properties 
// to all the DOF points of the Face
void DG::assignFluxVectortoFace(){
	int nDOF,Np;

	for(int i=0;i<this->noOfFaces;i++){
		
		// Number of DOF points for different faces
		if(faces[i].getFaceType() == FaceType::Quad){
			Np = this->order+1;
			nDOF = pow((Np),2);  //number of DOF points
		}

		else if(faces[i].getFaceType() == FaceType::Tri){
			nDOF = (order+1)*(order+2)/2.0;  //number of DOF points
		}



	
		// Allocating space to variable 
		faces[i].getFluxArray()->setSize(totNoOfVariable, nDOF);
		faces[i].getPrimitiveVariableArray()->setSize(totNoOfVariable, nDOF);

		// initializing this values to be zero. This will be re-written in the Riemann solver section
		// based on upwinding.
		for(Index nvar = 0; nvar < totNoOfVariable; ++nvar){
			for(Index DOF = 0; DOF < nDOF; ++DOF){
				faces[i].setFlux(nvar,DOF, 0.0);
				faces[i].setPrimitiveVariable(nvar,DOF, 0.0);
			};
		};
	};
};


// write Reference variable array.
// This function writes the data (reference variables) in the openFoam style. 
// Run paraFoam to visualize using paraview
void DG::writeVariableRefArray(double Time){
	this->computeReferenceSolution(Time);
	string time;
	time = varToString(Time);
	string timePath,fileName,word; //for the output files
	string headerFile; //for a sample file for copying the header
	timePath = path + "/" + time + "/";
	mkdir(timePath.c_str(),0777);
	double value;
	for (Index nvar=0; nvar < noOfVariable; nvar++){
		int vectorFlag = 0;
		fileName = timePath + VariableName[nvar] + "_Ref";	        // Writing the variable to the timePath folder
		headerFile = path + "/0/" + VariableName[nvar];    // Assuming that there will always be a 0 folder. 
		                                        // For finer control, even headerFile address can be accepted from the user

		// Open 0 dir files for copying the OpenFOAM header;
		ifstream readFile(headerFile.c_str()); // Reading headerFile for copying the OpenFOAM header (for plotting purpose)
		ofstream writeFile(fileName.c_str());  // Actual file where writing will be done
		string line;

		int Flag = 0;    // Flag 0 indicates the header and boundaryField information. This will be copied from the headerFile
		// Flag 1 indiacates the part where actual data is written down. This will be written from Variable array.
 
		while(getline(readFile,line)){ 		//looping over the readFile 
			istringstream stringline(line);
			stringline >> word;


			if (word == "internalField"){
				writeFile << word;
				stringline >> word; // This is just so that the word is not stuck at "internalField" and inadvertently traces the loop twice.
				Flag = 1;   //this is where data writing takes place
				if (this->VariableSize.getValue(nvar) == 1){
					// Scalar Field
					writeFile << " nonuniform List<scalar>" << endl;
					vectorFlag = 0;
				}
				else{
					// Vector Field
					writeFile << " nonuniform List<vector>" << endl;
					vectorFlag = 1;
				};
				writeFile << this->Domain.getNoOfCells() << endl;
				writeFile << "(" << endl;
				if (vectorFlag == 0){
					for (int ncell=0; ncell<Domain.getNoOfCells(); ncell++){
						writeFile << variableRef.getValue(cummulativeVariableSize.getValue(nvar),ncell) << endl;
					}
				}
				else{
					// For vectors, we need to enter in the following format
					// (number1 number2 ... numberN)
					for (int ncell=0;ncell<Domain.getNoOfCells();ncell++){
						writeFile << "(" ;
						for(int k=0;k<VariableSize.getValue(nvar);k++){
							if (k == VariableSize.getValue(nvar) - 1){
								writeFile  << variableRef.getValue(cummulativeVariableSize.getValue(nvar)+k, ncell) ; // no tab required after the final member of the vector
							}
							else{
								writeFile  << variableRef.getValue(cummulativeVariableSize.getValue(nvar)+k, ncell) << " "; //spacing required after each member till the final member.
							};
						}
						writeFile << ")" << endl;
					}	
				};
				writeFile << ")" << endl;
				writeFile << ";" << endl;
				writeFile << endl;
			}
			else if (word == "boundaryField"){
				// Boundary data needs to be implemented
				Flag = 0;
			};

			if (Flag == 0 ){
				if (word == "class"){
					// identifier for the "location" line.
					// If we are reading from the 0 dir, the location information in the OpenFOAM header is absent sometimes!
					// However, for all Time > 0, the location needs to be specified.
					writeFile << line << endl; // The line starting with "class" written as it is 

					// next line streamed
					// First word of this (next) line. Ideally it should read "location"
					getline(readFile,line);
					istringstream stringline(line);
					stringline >> word; 
					if (word == "location"){
						// Indicates that the headerFile contains a keyword "location" (default in 0 dir)
						// We replace this file with our location
						writeFile << "    location    \"" << Time << "\";" << endl; 
						// then we move on.
					}
					else{
						// Now the location information needs to be added explicitely.
						writeFile << "    location    \"" << Time << "\";" << endl; 
						// Also remember that, if this line doesn't contain the location information, then
						// it contains the next field (which is 'object' info).
						// That info needs to be printed here. Otherwise, the loop will proceed naturally
						// to the next line skipping object info altogether.
						writeFile << line << endl;
					};

				}
				else{
					// Rest of the lines are written as it is
					writeFile << line << endl;
				};

			};
		};

		writeFile.close();
		readFile.close();
	}
};



// write variable array.
// This function writes the data (field data or variables) in the openFoam style. 
// Run paraFoam to visualize using paraview
void DG::writeVariableArray(double Time){
	string time;
	time = varToString(Time);
	string timePath,fileName,word; //for the output files
	string headerFile; //for a sample file for copying the header
	timePath = path + "/" + time + "/";
	mkdir(timePath.c_str(),0777);
	double value;
	for (Index nvar=0; nvar < noOfVariable; nvar++){
		int vectorFlag = 0;
		fileName = timePath + VariableName[nvar];	        // Writing the variable to the timePath folder
		headerFile = path + "/0/" + VariableName[nvar];    // Assuming that there will always be a 0 folder. 
		                                        // For finer control, even headerFile address can be accepted from the user

		// Open 0 dir files for copying the OpenFOAM header;
		ifstream readFile(headerFile.c_str()); // Reading headerFile for copying the OpenFOAM header (for plotting purpose)
		ofstream writeFile(fileName.c_str());  // Actual file where writing will be done
		string line;

		int Flag = 0;    // Flag 0 indicates the header and boundaryField information. This will be copied from the headerFile
		// Flag 1 indiacates the part where actual data is written down. This will be written from Variable array.
 
		while(getline(readFile,line)){ 		//looping over the readFile 
			istringstream stringline(line);
			stringline >> word;


			if (word == "internalField"){
				writeFile << word;
				stringline >> word; // This is just so that the word is not stuck at "internalField" and inadvertently traces the loop twice.
				Flag = 1;   //this is where data writing takes place
				if (this->VariableSize.getValue(nvar) == 1){
					// Scalar Field
					writeFile << " nonuniform List<scalar>" << endl;
					vectorFlag = 0;
				}
				else{
					// Vector Field
					writeFile << " nonuniform List<vector>" << endl;
					vectorFlag = 1;
				};
				writeFile << this->Domain.getNoOfCells() << endl;
				writeFile << "(" << endl;
				if (vectorFlag == 0){
					for (int ncell=0; ncell<Domain.getNoOfCells(); ncell++){
						writeFile << variable.getValue(cummulativeVariableSize.getValue(nvar),ncell) << endl;
					}
				}
				else{
					// For vectors, we need to enter in the following format
					// (number1 number2 ... numberN)
					for (int ncell=0;ncell<Domain.getNoOfCells();ncell++){
						writeFile << "(" ;
						for(int k=0;k<VariableSize.getValue(nvar);k++){
							if (k == VariableSize.getValue(nvar) - 1){
								writeFile  << variable.getValue(cummulativeVariableSize.getValue(nvar)+k, ncell) ; // no tab required after the final member of the vector
							}
							else{
								writeFile  << variable.getValue(cummulativeVariableSize.getValue(nvar)+k, ncell) << " "; //spacing required after each member till the final member.
							};
						}
						writeFile << ")" << endl;
					}	
				};
				writeFile << ")" << endl;
				writeFile << ";" << endl;
				writeFile << endl;
			}
			else if (word == "boundaryField"){
				// Boundary data needs to be implemented
				Flag = 0;
			};

			if (Flag == 0 ){
				if (word == "class"){
					// identifier for the "location" line.
					// If we are reading from the 0 dir, the location information in the OpenFOAM header is absent sometimes!
					// However, for all Time > 0, the location needs to be specified.
					writeFile << line << endl; // The line starting with "class" written as it is 

					// next line streamed
					// First word of this (next) line. Ideally it should read "location"
					getline(readFile,line);
					istringstream stringline(line);
					stringline >> word; 
					if (word == "location"){
						// Indicates that the headerFile contains a keyword "location" (default in 0 dir)
						// We replace this file with our location
						writeFile << "    location    \"" << Time << "\";" << endl; 
						// then we move on.
					}
					else{
						// Now the location information needs to be added explicitely.
						writeFile << "    location    \"" << Time << "\";" << endl; 
						// Also remember that, if this line doesn't contain the location information, then
						// it contains the next field (which is 'object' info).
						// That info needs to be printed here. Otherwise, the loop will proceed naturally
						// to the next line skipping object info altogether.
						writeFile << line << endl;
					};

				}
				else{
					// Rest of the lines are written as it is
					writeFile << line << endl;
				};

			};
		};

		writeFile.close();
		readFile.close();
	}
};


	
// This function initializes arrays and reads geometry files
// and populates the points, faces and cells arrays.
void DG::createDomain(){
	string geomPath;
	geomPath = path + "/constant/polyMesh"; // The files are in polymesh folder
	Domain.getProblemData(geomPath);
	//this->noOfPoints
	
	// Allocating space to points, faces and cells array
	this->noOfPoints = Domain.getNoOfPoints();
	points = new Point[this->noOfPoints];
	this->pointsFlag = true;

	this->noOfFaces = Domain.getNoOfFaces();
	faces = new Face[this->noOfFaces];
	this->facesFlag = true;

	this->noOfCells = Domain.getNoOfCells();
	cells = new Cell[this->noOfCells];
	this->cellsFlag = true;
	
	// Addition for boundary conditions
	this->noOfBoundaryConditions = Domain.getNoOfBoundaryConditions();
	bcond = new BoundaryConditions[this->noOfBoundaryConditions];
	this->bcondFlag = true;
	
	// Populating the data for the arrays
	// bcond for boundary conditions
	Domain.fillDataArrays(geomPath, points, faces, cells, bcond); //setCellType is called internally

	// Get cell information for numerical integration 
	this->setupReferenceCellsArray();
	this->setupReferenceFacesArray();

	// Assign reference cell to all the cells
	for (Index ncell=0; ncell<this->noOfCells; ncell++){
		if (cells[ncell].getCellType() == CellType::Hex){
			cells[ncell].setReferenceCell(&referenceCells[0]);
			assert(cells[ncell].getReferenceCell()->getCellType() == CellType::Hex);
		}
		else if (cells[ncell].getCellType() == CellType::Tet){
		}
		else if (cells[ncell].getCellType() == CellType::Prism){
		}
		else if (cells[ncell].getCellType() == CellType::Pyramid){
		};
	};

	// Assign reference face to all the faces
	for (Index nface =0; nface<this->noOfFaces; nface++){
		if (faces[nface].getFaceType() == FaceType::Quad){
			faces[nface].setReferenceFace(&referenceFaces[0]);
			assert(faces[nface].getReferenceFace()->getFaceType() == FaceType::Quad);
		}
		else if (faces[nface].getFaceType() == FaceType::Tri){
		};
	};


	this->getQuadPointsGlobalLocation(); // Also sets nQuad number for each cell
	this->getDOFPointsGlobalLocation(); // Also sets nDOF number for each cell

	this->calculateJacobian(); //Jacobian of transformation for each cell computed at all quadrature points
	this->calculateInverseJacobian(); // Inverse Jacobian matrix (metric of transformation). To be used for derivatives.
					  // dRST_by_dXYZ filled for all the nQuad points for each cell.


	// next we need to find these quantities for all the faces as well
	this->getFaceQuadPointsGlobalLocation(); // Gets (x,y,z) location of all the quadrature points for the faces. Also sets nQuad numberfor each face
	this->getFaceDOFPointsGlobalLocation(); // Gets (x,y,z) location of all the quadrature points for the faces. Also sets nDOF number for each face

	this->mapFaceDOFPointsToCellDOFPoints(); // Maps face quadrature points (2D) to the cell quadrature points (3D)

	this->mapFaceQuadPointsToCellQuadPoints(); // Maps face quadrature points (2D) to the cell quadrature points (3D)
	this->calculateFaceJacobian(); // This is actually ||dX/dr x dX/ds||. Refer Face::calculateJacobianQuadFace();


	this->shortestDistance = Domain.findShortestDistance(cells);
	this->charLength = Domain.findCharacteristicLength(cells);


	this->assignBoundaryConditions();
};

// Printing the Domain information
// Simply calls print function of the Geometry class.
void DG::printDomain(){
	cout<<"~~~~~~~~~~~~~~~~~~~Domain Details~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
	Domain.print();
};


double DG::getShortestDistance(){
	return this->shortestDistance;
};

double DG::getCharacteristicLength(){
	return this->charLength;
};

void DG::assignBoundaryConditions(){
};

void DG::generateVariableArray(int noOfVariables, string VariableName[],int VariableSize[]){
	setNoOfVariable(noOfVariables);
	setVariableName(VariableName);
	setVariableSize(VariableSize);
	addTotNoOfVariable();
	allocateVariableArraySize();
	allocateBoundaryVariableArraySize();
	
};

void DG::assignFaceName(){
	// in this function, each face is given a name
	// if the face is an internal face, the name is "internal"
	// else the face shares the same name as the boundary it is the part of.
	for (int f=0; f<noOfFaces; f++){
		if (faces[f].getNeighbourCell() != NULL){
			faces[f].setName("internal");
		};
	};

	for (int bcondi=0; bcondi < noOfBoundaryConditions; bcondi++){
		int startFace = bcond[bcondi].getStartFace();
		int endFace = bcond[bcondi].getStartFace() + bcond[bcondi].getNoOfFaces();
		for (int fi=startFace; fi < endFace; fi++){
			faces[fi].setName(bcond[bcondi].getName());
		};
	};
};


void DG::assignNoOfBoundaryFacesToCells(){
	// in this function, a count 'noOfBoundaryFaces' is set inside each cell.
	// if the cell is such that all the definingFaces are "internal" faces, then the cell is obviously not near to the boundary
	// then noOfBoundaryFaces = 0
	// if the cell has one of its faces as the boundary face, then noOfBoundaryFaces = 1, and so on
	// for the cell at the 3D corner, this number will be 3

	// initialize with zero
	for (int c=0; c<noOfCells; c++){
		cells[c].setNoOfBoundaryFaces(0);
	};

	for (int bcondi=0; bcondi < noOfBoundaryConditions; bcondi++){
		int startFace = bcond[bcondi].getStartFace();
		int endFace = bcond[bcondi].getStartFace() + bcond[bcondi].getNoOfFaces();
		for (int fi=startFace; fi < endFace; fi++){
			faces[fi].getOwnerCell()->addNoOfBoundaryFaces(1);
		};
	};

};


// readData() function basically calls different function from DG class
// takes variable names and size as input
// the entire data for the dependent variables is stored in the variable matrix
void DG::readData(double Time){
    	readVariableArray(Time);
   	assignVariabletoCell();
   	assignFluxVectortoFace();
	assignFaceName();
	assignNoOfBoundaryFacesToCells();
};

void DG::applyIC(){
	if (this->systemType == System::scalarLinear){
		// scalar Advection problem implemented.
		// make sure that generateVariableArray() is called first.
		double rho, u, v, w;
		double x, y, z;

		double xc = -0.5;
		double yc = 0.0;
		double zc = 0.0;
		double sigmac = 1.0/8.0;

		for(int i=0;i<this->noOfCells;i++){
			x = cells[i].getCenter()->getX();
			y = cells[i].getCenter()->getY();
			z = cells[i].getCenter()->getZ();

			u = y;
			v = -x;
			w = 0;

			rho = exp(-( sqr(x-xc) + sqr(y-yc) + sqr(z-zc) )/ (2.0 * sigmac *sigmac));

			this->variable.setValue(0, i , rho);
			this->variable.setValue(1, i , u );
			this->variable.setValue(2, i , v );
			this->variable.setValue(3, i , w );
		};
	}
	else if (this->systemType == System::Euler){

		// Explosion (Euler equations).
		// make sure that generateVariableArray() is called first.
		/*
		double rho, u, v, w, p;
		double x, y, z;

		double xc = 0.5;
		double yc = 0.5;
		double zc = 0.5;
		double sigmac = 1.0/8.0;

		for(int i=0;i<this->noOfCells;i++){
			x = cells[i].getCenter()->getX();
			y = cells[i].getCenter()->getY();
			z = cells[i].getCenter()->getZ();

			u = 0.0;
			v = 0.0;
			w = 0.0;

			rho = 1.0 + exp(-( sqr(x-xc) + sqr(y-yc) + sqr(z-zc) )/ (2.0 * sigmac *sigmac));
			p = 1.0 + 10.0*exp(-( sqr(x-xc) + sqr(y-yc) + sqr(z-zc) )/ (2.0 * sigmac *sigmac));

			this->variable.setValue(0, i , rho);
			this->variable.setValue(1, i , u );
			this->variable.setValue(2, i , v );
			this->variable.setValue(3, i , w );
			this->variable.setValue(4, i , p );
		};
		*/

		// Sod's shock tube problem
		double rho, u, v, w, P;
		double x, y, z;

		double xc = 5.0;
		double R = 287;

		for(int i=0;i<this->noOfCells;i++){
			x = cells[i].getCenter()->getX();
			y = cells[i].getCenter()->getY();
			z = cells[i].getCenter()->getZ();

			u = 0.0;
			v = 0.0;
			w = 0.0;

			if (x < xc){
				P = 50000;
				rho = P/ (R * 300);
			}
			else{
				P = 10000;
				rho = P/ (R * 300);
			};

			this->variable.setValue(0, i , rho);
			this->variable.setValue(1, i , u );
			this->variable.setValue(2, i , v );
			this->variable.setValue(3, i , w );
			this->variable.setValue(4, i , P );
		};
	};

	this->writeVariableArray(1);
	// later copy this folder manually as 0 folder and change the 'location' field inside of all files.
};


void DG::computeReferenceSolution(double time){
	//Depends on the problem.
	//The variableRef array already stores the IC values. 
	//For problems with an analytical solution available (such as scalarTransportSolver, this value is set here.

	if (this->refVariableFlag == false){
		getAnalyticalShockTube(&this->SodAnalytical);
		this->refVariableFlag = true;
	};

	double rho2 = SodAnalytical.getValue(0,0);
	double u2 = SodAnalytical.getValue(1,0);
	double v2 = SodAnalytical.getValue(2,0);
	double w2 = SodAnalytical.getValue(3,0);
	double p2 = SodAnalytical.getValue(4,0);
	double T2 = SodAnalytical.getValue(5,0);
	double a2 = SodAnalytical.getValue(6,0);
	double M2 = SodAnalytical.getValue(7,0);

	double rho4 = SodAnalytical.getValue(0,1);
	double u4 = SodAnalytical.getValue(1,1);
	double v4 = SodAnalytical.getValue(2,1);
	double w4 = SodAnalytical.getValue(3,1);
	double p4 = SodAnalytical.getValue(4,1);
	double T4 = SodAnalytical.getValue(5,1);
	double a4 = SodAnalytical.getValue(6,1);
	double M4 = SodAnalytical.getValue(7,1);

	double rho3 = SodAnalytical.getValue(0,2);
	double u3 = SodAnalytical.getValue(1,2);
	double v3 = SodAnalytical.getValue(2,2);
	double w3 = SodAnalytical.getValue(3,2);
	double p3 = SodAnalytical.getValue(4,2);
	double T3 = SodAnalytical.getValue(5,2);
	double a3 = SodAnalytical.getValue(6,2);
	double M3 = SodAnalytical.getValue(7,2);

	double rho1 = SodAnalytical.getValue(0,3);
	double u1 = SodAnalytical.getValue(1,3);
	double v1 = SodAnalytical.getValue(2,3);
	double w1 = SodAnalytical.getValue(3,3);
	double p1 = SodAnalytical.getValue(4,3);
	double T1 = SodAnalytical.getValue(5,3);
	double a1 = SodAnalytical.getValue(6,3);
	double M1 = SodAnalytical.getValue(7,3);

	double Ws = M1*a1;

	for(int i=0;i<this->noOfCells;i++){
		double x = cells[i].getCenter()->getX();
		double L = 10.0; 
		double x_shock=L/2+Ws*time;
		double x_exp_h=L/2-a2*time;
		double x_exp_t=L/2+(u4-a4)*time;
		double x_cd=L/2+u3*time;

		double rho, u, v, w, p;
		if (x >= x_shock){
			rho = rho1;
			u = u1;
			v = v1;
			w = w1;
			p = p1;
		}
		else if (x < x_shock and x >= x_cd){
			rho = rho3;
			u = u3;
			v = v3;
			w = w3;
			p = p3;
		}
		else if (x < x_cd and x >= x_exp_t){
			rho = rho4;
			u = u4;
			v = v4;
			w = w4;
			p = p4;
		}
		else if (x < x_exp_h){
			rho = rho2;
			u = u2;
			v = v2;
			w = w2;
			p = p2;
		}
		else{
			double rat=(x_exp_t-x_exp_h)/(x-x_exp_h);
			u=u2+(u4-u2)/rat;
			v = 0.0;
			w = 0.0;
			double T2 = p2/(0.287*rho2);
			double T=T2*(pow(expp(u,a2),2));
			p=p2*(pow(expp(u,a2),2*1.4/0.4));
			rho=rho2*(pow(expp(u,a2),2/0.4));

		};

		this->variableRef.setValue(0, i, rho); 
		this->variableRef.setValue(1, i, u); 
		this->variableRef.setValue(2, i, v); 
		this->variableRef.setValue(3, i, w); 
		this->variableRef.setValue(4, i, p*1000.0); 

		for(Index DOF=0; DOF < cells[i].getNDOF(); DOF++){
			double x = cells[i].getDOFPointsGlobalLocation()->getValue(DOF,0);
			double L = 10.0; 
			double x_shock=L/2+Ws*time;
			double x_exp_h=L/2-a2*time;
			double x_exp_t=L/2+(u4-a4)*time;
			double x_cd=L/2+u3*time;

			double rho, u, v, w, p;
			if (x >= x_shock){
				rho = rho1;
				u = u1;
				v = v1;
				w = w1;
				p = p1;
			}
			else if (x < x_shock and x >= x_cd){
				rho = rho3;
				u = u3;
				v = v3;
				w = w3;
				p = p3;
			}
			else if (x < x_cd and x >= x_exp_t){
				rho = rho4;
				u = u4;
				v = v4;
				w = w4;
				p = p4;
			}
			else if (x < x_exp_h){
				rho = rho2;
				u = u2;
				v = v2;
				w = w2;
				p = p2;
			}
			else{
				double rat=(x_exp_t-x_exp_h)/(x-x_exp_h);
				u=u2+(u4-u2)/rat;
				v = 0.0;
				w = 0.0;
				double T2 = p2/(0.287*rho2);
				double T=T2*(pow(expp(u,a2),2));
				p=p2*(pow(expp(u,a2),2*1.4/0.4));
				rho=rho2*(pow(expp(u,a2),2/0.4));

			};

			cells[i].getVariableRefArray()->setValue(0, DOF, rho);
			cells[i].getVariableRefArray()->setValue(1, DOF, u);
			cells[i].getVariableRefArray()->setValue(2, DOF, v);
			cells[i].getVariableRefArray()->setValue(3, DOF, w);
			cells[i].getVariableRefArray()->setValue(4, DOF, p*1000.0);
		};
	};
};


void DG::computeError(){
	//Computes L_infty, L_1 and L_2 errors between the solution (variable array) and the reference solution (variableRef)

	// L_infinity error:
	/*

	   // Error based on cell averages
	double max = 0.0;
	for(int i=0;i<this->noOfCells;i++){
		if (abs(this->variable.getValue(0,i) - this->variableRef.getValue(0,i)) > max){
			max = abs(this->variable.getValue(0,i) - this->variableRef.getValue(0,i));
		};

	};
	this->error[0] = max;   // infinity norm error
	*/

	
	/* 
	*/
	   // Error based on DOF Value
	double max = 0.0;
	for(int i=0;i<this->noOfCells;i++){
		for (Index DOF =0; DOF < cells[i].getNDOF(); DOF++){
			if (abs(cells[i].getVariableRefArray()->getValue(0,DOF) - cells[i].getVariableArray()->getValue(0,0,DOF)) > max){
				max = abs(cells[i].getVariableRefArray()->getValue(0,DOF) - cells[i].getVariableArray()->getValue(0,0,DOF));
			};
		};
	};
	this->error[0] = max;   // infinity norm error





	// L_1 error 
	/*
	   // Error based on cell averages
	this->error[1] = 0.0;
	double num = 0.0;
	double den = 0.0;
	for(int i=0;i<this->noOfCells;i++){
		num += abs(this->variableRef.getValue(0,i) - this->variable.getValue(0,i));
		den += abs(this->variableRef.getValue(0,i)); 
	};
	this->error[1] = num/den;
        */

	   // Error based on DOF Value

	this->error[1] = 0.0;
	for(int i=0;i<this->noOfCells;i++){
		for (Index DOF =0; DOF < cells[i].getNDOF(); DOF++){
			double E = abs(cells[i].getVariableArray()->getValue(0,0,DOF) - cells[i].getVariableRefArray()->getValue(0,DOF));
			error[1]+= E ;
		};
	};
	error[1] = error[1]/(this->noOfCells * cells[0].getNDOF());

	/*
	*/






	// L_2 error
	/*
	   // Error based on cell averages
	this->error[2] = 0.0;
	num = 0.0;
	den = 0.0;
	for(int i=0;i<this->noOfCells;i++){
		num += sqr(this->variableRef.getValue(0, i) - this->variable.getValue(0,i));
		den += sqr(this->variableRef.getValue(0, i)); 
	};
	this->error[2] = sqrt(num/den);
	*/

	   // Error based on DOF value
	this->error[2] = 0.0;
	for(int i=0;i<this->noOfCells;i++){
		for (Index DOF =0; DOF < cells[i].getNDOF(); DOF++){
			double E = sqr(cells[i].getVariableArray()->getValue(0,0,DOF) - cells[i].getVariableRefArray()->getValue(0,DOF));
			error[2]+= E; 
		};
	};
	error[2] = sqrt(error[2])/(this->noOfCells * cells[0].getNDOF());
	/*
	*/

};

double* DG::getError(){
	return this->error;
};







//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 3. Cell Functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//


// Get global (x,y,z) location of  every quadrature point (dim: nQuad x 3)
// for all cells
// This functions in turn calls mapRSTtoXYZ functions in the Cell class (e.g. Cell:mapRSTtoXYZ3DTensor() for a hex cell)
void DG::getQuadPointsGlobalLocation(){
	for(Index ncell=0; ncell<this->noOfCells; ncell++){
		if(cells[ncell].getCellType() == CellType::Hex){
		
			FunctionalSpace F(this->order, this->intType);
			int Np = this->order+1 + this->intType; //intType: 0->inexact, 1->exact
			int nQuad = pow((Np),3);  //number of quadrature points
			cells[ncell].setNQuad(nQuad);   

			TensorO1<double> r(nQuad);
			TensorO1<double> s(nQuad);
			TensorO1<double> t(nQuad);
			TensorO1<double> w(nQuad);
			
			F.LGLRootsAndWeights3D(Np,Np,Np, &r, &s, &t, &w);
			
			// Assign size of the quadPointsGlobalLocation array
			cells[ncell].getQuadPointsGlobalLocation()->setSize(nQuad, 3);

			Point cellQuadPoint; 	// global coordinate system
			Point refCellQuadPoint; // local coordinate system
			// Main loop for calculating the mapping at each quadrature point
			for (Index QP=0; QP<nQuad; QP++){	

				refCellQuadPoint.setX(r.getValue(QP));
				refCellQuadPoint.setY(s.getValue(QP));
				refCellQuadPoint.setZ(t.getValue(QP));

				cells[ncell].mapRSTtoXYZ3DTensor(&refCellQuadPoint, &cellQuadPoint);

				cells[ncell].getQuadPointsGlobalLocation()->setValue(QP, 0, cellQuadPoint.getX());
				cells[ncell].getQuadPointsGlobalLocation()->setValue(QP, 1, cellQuadPoint.getY());
				cells[ncell].getQuadPointsGlobalLocation()->setValue(QP, 2, cellQuadPoint.getZ());
			};
	
		}
		else if(cells[ncell].getCellType()==CellType::Tet){
		}
		else if(cells[ncell].getCellType()==CellType::Prism){
		}
		else if(cells[ncell].getCellType()==CellType::Pyramid){
		}
	}
}


// Get global (x,y,z) location of  every DOF point (dim: nDOF x 3)
// for all cells
// This functions in turn calls mapRSTtoXYZ functions in the Cell class (e.g. Cell:mapRSTtoXYZ3DTensor() for a hex cell)
void DG::getDOFPointsGlobalLocation(){
	for(Index ncell=0; ncell<this->noOfCells; ncell++){
		if(cells[ncell].getCellType() == CellType::Hex){
		
			FunctionalSpace F(this->order, this->intType);

			int Np = this->order+1; 
			int nDOF = pow((Np),3);  //number of DOF points
			cells[ncell].setNDOF(nDOF);   

			TensorO1<double> r(nDOF);
			TensorO1<double> s(nDOF);
			TensorO1<double> t(nDOF);
			TensorO1<double> w(nDOF);
			
			F.LGLRootsAndWeights3D(Np,Np,Np, &r, &s, &t, &w);
			
			// Assign size of the quadPointsGlobalLocation array
			cells[ncell].getDOFPointsGlobalLocation()->setSize(nDOF, 3);

			Point cellDOFPoint; 	// global coordinate system
			Point refCellDOFPoint; // local coordinate system
			// Main loop for calculating the mapping at each quadrature point
			for (Index DOF=0; DOF<nDOF; DOF++){	

				refCellDOFPoint.setX(r.getValue(DOF));
				refCellDOFPoint.setY(s.getValue(DOF));
				refCellDOFPoint.setZ(t.getValue(DOF));

				cells[ncell].mapRSTtoXYZ3DTensor(&refCellDOFPoint, &cellDOFPoint);

				cells[ncell].getDOFPointsGlobalLocation()->setValue(DOF, 0, cellDOFPoint.getX());
				cells[ncell].getDOFPointsGlobalLocation()->setValue(DOF, 1, cellDOFPoint.getY());
				cells[ncell].getDOFPointsGlobalLocation()->setValue(DOF, 2, cellDOFPoint.getZ());
			};
		
		}
		else if(cells[ncell].getCellType()==CellType::Tet){
		}
		else if(cells[ncell].getCellType()==CellType::Prism){
		}
		else if(cells[ncell].getCellType()==CellType::Pyramid){
		};
	};
};



// Constructs Jacobian at each and every quadrature point (nQuad in each dim = nDOF in each dim + intType)
// for all cells
// This functions calls computeJacobian functions in the Cell class
void DG::calculateJacobian(){
	for(Index ncell=0; ncell<this->noOfCells; ncell++){
		if(cells[ncell].getCellType() == CellType::Hex){

			// Declare the functional space
			FunctionalSpace F(this->order, this->intType);

			// 1D quadrature points number 
			int Np = this->order+1 + this->intType; //intType: 0->inexact, 1->exact

			// Total quadrature points number 
			int nQuad = pow((Np),3);  //number of quadrature points
			assert(cells[ncell].getNQuad() == nQuad);


			TensorO1<double> r(nQuad);
			TensorO1<double> s(nQuad);
			TensorO1<double> t(nQuad);
			TensorO1<double> w(nQuad);
			
			F.LGLRootsAndWeights3D(Np,Np,Np, &r, &s, &t, &w);
			
			// Set the corresponding flag to true
			cells[ncell].JacobianFlag = true;

			// Assign space for the cell Jacobian array 
			cells[ncell].getJacobian()->setSize(nQuad);	

			// Initialize with zero
			for (Index QP=0; QP<nQuad; QP++){
			        cells[ncell].getJacobian()->setValue(QP, 0.0);
			};

			// Main loop for calculating the Jacobian at each interpolating points
			for (Index QP=0; QP<nQuad; QP++){	
				cells[ncell].getJacobian()->setValue(QP, cells[ncell].calculateJacobian3DTensor(r.getValue(QP),
							                                                        s.getValue(QP),
														t.getValue(QP)));
			};
		
		}
		else if(cells[ncell].getCellType()==CellType::Tet){
		}
		else if(cells[ncell].getCellType()==CellType::Prism){
		}
		else if(cells[ncell].getCellType()==CellType::Pyramid){
		}
	};
};



// Constructs Inverse Jacobian at each and every quadrature point (nQuad x 3 x 3)
// for all cells
// This functions calls Cell::computeInverseJacobian functions in the Cell class
void DG::calculateInverseJacobian(){
	for(Index ncell=0; ncell<this->noOfCells; ncell++){
		if(cells[ncell].getCellType()==CellType::Hex){
		
			FunctionalSpace F(this->order, this->intType);

			// 1D quadrature points number 
			int Np = this->order+1 + this->intType; //intType: 0->inexact, 1->exact

			// Total quadrature points number 
			int nQuad = pow((Np),3);  //number of quadrature points
			assert(cells[ncell].getNQuad() == nQuad);


			TensorO1<double> r(nQuad);
			TensorO1<double> s(nQuad);
			TensorO1<double> t(nQuad);
			TensorO1<double> w(nQuad);
			
			F.LGLRootsAndWeights3D(Np,Np,Np, &r, &s, &t, &w);
			
			// Assign memory for dRST_by_dXYZ array
			cells[ncell].getInverseJacobian()->setSize(nQuad, 3, 3);

			// Main loop for calculating the Jacobian at each interpolating points
			for (int QP=0; QP<nQuad; QP++){	
				// Local matrix of the Inverse Jacobian at quadrature point QP :
				Matrix<double> drstbydxyz_at_QP(3);
				cells[ncell].calculateInverseJacobianMatrix3DTensor(r.getValue(QP),
						                                    s.getValue(QP),
										    t.getValue(QP), &drstbydxyz_at_QP);

				for (Index i=0; i<3; i++){
					for (Index j=0; j<3; j++){
						cells[ncell].getInverseJacobian()->setValue(QP, i, j, drstbydxyz_at_QP.getValue(i,j));
					};
				};
			};
		
		}
		else if(cells[ncell].getCellType()==CellType::Tet){
		}
		else if(cells[ncell].getCellType()==CellType::Prism){
		}
		else if(cells[ncell].getCellType()==CellType::Pyramid){
		};
	};
};




//Cell matrices are created
void DG::computeCellMatrix(){
	FunctionalSpace F(this->order, this->intType);

	cout <<"Computing cell matrices...\n";
	for(Index ncell=0; ncell<this->noOfCells; ncell++){

		// printing progress bar. Refer mathAndOtherFunctions.cpp
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		double percentComplete = 1.0*ncell / (1.0*this->noOfCells-1.0) * 100;
		if (this->noOfCells == 1){
			percentComplete = 100;
		};
		printProgressBar(percentComplete);
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

		Cell *cell;
		cell = &this->cells[ncell]; //c is a pointer to ith cell

		assert (cell->getCellType() == CellType::Hex or cell->getCellType() == CellType::Tet or cell->getCellType() == CellType::Prism or cell->getCellType() == CellType::Pyramid);
		assert (cell->refCellFlag == true and cell->JacobianFlag == true);

		// Re-size all the matrices 
		cell->resizeMatrices(cell->getNDOF());

		if (cell->getCellType() == CellType::Hex){
			// Accumulating necessary data before computing matrices (common data for cell matrices)
			assert(cell->getCellType() == cell->getReferenceCell()->getCellType() and cell->getCellType() == CellType::Hex);

			// Generate mass matrix from FunctionalSpace class.
			F.generateMassMatrix3DTensor(cell);

			// Get inverse mass matrix. Refere Tensor/Matrix for getInverse() method description.
			cell->getMassMatrix()->getInverse(cell->getInverseMassMatrix());

			//generate Vandermonde matrix
			F.generateVandermondeMatrix3DTensor(cell);


			F.generateDiffMatrix3DTensor(cell); //it's ok even if the formFlag argument is not passed. 0 is default.

			// generation of shortcut matrices to reduce time in every step
			Math::dot(cell->getInverseMassMatrix(), cell->getDxMatrix(), cell->getDtildexMatrix());
			Math::dot(cell->getInverseMassMatrix(), cell->getDyMatrix(), cell->getDtildeyMatrix());
			Math::dot(cell->getInverseMassMatrix(), cell->getDzMatrix(), cell->getDtildezMatrix());

			//generate Flux matrix
			F.generateFluxMatrix3DTensor(cell);

			// flux times m_inverse needs to change according to the flux structures
			for (int f=0; f<6; f++){
				Matrix<double> FtildeLocal(cell->getNDOF());
				Matrix<double> F(cell->getNDOF());
				for (Index i=0; i<cell->getNDOF(); i++){
					for (Index j=0; j<cell->getNDOF(); j++){
						F.setValue(i,j, cell->getFMatrix()->getValue(f, i,j));
					};
				};
				Math::dot(cell->getInverseMassMatrix(), &F, &FtildeLocal);

				for (Index i=0; i<cell->getNDOF(); i++){
					for (Index j=0; j<cell->getNDOF(); j++){
						cell->getFtildeMatrix()->setValue(f, i,j, FtildeLocal.getValue(i,j));
					};
				};
			};

			//generateOtherMatrices next


		}
		else if (cell->getCellType() == CellType::Tet){
		}
		else if (cell->getCellType() == CellType::Prism){
		}
		else if (cell->getCellType() == CellType::Pyramid){
		}
		else{
			cout << "Error! check the cell type in function DG::computeCellMatrix()\n";
		}
	};
}










//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 4. Face Functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//


// Get global (x,y,z) location of  every quadrature point on 2D face(dim: nQuadFace x 3)
// for all faces
// This functions in turn calls mapRSTtoXYZ functions in the Face class (e.g. Face:mapRSTtoXYZQuad() for a quadrilateral face)
void DG::getFaceQuadPointsGlobalLocation(){
	for(Index nface=0; nface<this->noOfFaces; nface++){
		if(faces[nface].getFaceType() == FaceType::Quad){
			FunctionalSpace F(this->order, this->intType);

			int Np = this->order+ 1 + this->intType; //intType: 0->inexact, 1->exact
			int nQuad = pow((Np),2);  //number of quadrature points
			faces[nface].setNQuad (nQuad);   

			// Assign size of the quadPointsGlobalLocation array
			faces[nface].getQuadPointsGlobalLocation()->setSize(nQuad, 3);

			TensorO1<double> r(nQuad);
			TensorO1<double> s(nQuad);
			TensorO1<double> w(nQuad);
			
		
			F.LGLRootsAndWeights2D(Np,Np, &r, &s, &w);
			
			Point faceQuadPoint; 	// global coordinate system
			Point refFaceQuadPoint; // local coordinate system
			// Main loop for calculating the mapping at each quadrature point
			for (Index QP=0; QP<nQuad; QP++){	

				refFaceQuadPoint.setX(r.getValue(QP));
				refFaceQuadPoint.setY(s.getValue(QP));
				refFaceQuadPoint.setZ(0.0);	// doesn't matter, since it doesn't affect the shape functions

				faces[nface].mapRSTtoXYZQuad(&refFaceQuadPoint, &faceQuadPoint);

				faces[nface].getQuadPointsGlobalLocation()->setValue(QP, 0, faceQuadPoint.getX());
				faces[nface].getQuadPointsGlobalLocation()->setValue(QP, 1, faceQuadPoint.getY());
				faces[nface].getQuadPointsGlobalLocation()->setValue(QP, 2, faceQuadPoint.getZ());
			};
		
		}
		else if(faces[nface].getFaceType() == FaceType::Tri){
		};
	};
};


// Get global (x,y,z) location of  every DOFrature point on 2D face(dim: nDOFFace x 3)
// for all faces
// This functions in turn calls mapRSTtoXYZ functions in the Face class (e.g. Face:mapRSTtoXYZDOF() for a DOFrilateral face)
void DG::getFaceDOFPointsGlobalLocation(){
	for(Index nface=0; nface<this->noOfFaces; nface++){
		if(faces[nface].getFaceType() == FaceType::Quad){
			FunctionalSpace F(this->order, this->intType);

			int Np = this->order+ 1; //intType: 0->inexact, 1->exact
			int nDOF = pow((Np),2);  //number of DOFrature points
			faces[nface].setNDOF (nDOF);   

			// Set array size.
			faces[nface].getDOFPointsGlobalLocation()->setSize(nDOF, 3);

			TensorO1<double> r(nDOF);
			TensorO1<double> s(nDOF);
			TensorO1<double> w(nDOF);
			
		
			F.LGLRootsAndWeights2D(Np,Np, &r, &s, &w);
			
			Point faceDOFPoint; 	// global coordinate system
			Point refFaceDOFPoint; // local coordinate system
			// Main loop for calculating the mapping at each DOFrature point
			for (Index DOF=0; DOF<nDOF; DOF++){	

				refFaceDOFPoint.setX(r.getValue(DOF));
				refFaceDOFPoint.setY(s.getValue(DOF));
				refFaceDOFPoint.setZ(0.0);	// doesn't matter, since it doesn't affect the shape functions

				faces[nface].mapRSTtoXYZQuad(&refFaceDOFPoint, &faceDOFPoint);

				faces[nface].getDOFPointsGlobalLocation()->setValue(DOF, 0, faceDOFPoint.getX());
				faces[nface].getDOFPointsGlobalLocation()->setValue(DOF, 1, faceDOFPoint.getY());
				faces[nface].getDOFPointsGlobalLocation()->setValue(DOF, 2, faceDOFPoint.getZ());
			};
		
		}
		else if(faces[nface].getFaceType() == FaceType::Tri){
		};
	};
};


void DG::mapFaceDOFPointsToCellDOFPoints(){
	for (Index nface=0; nface<this->noOfFaces; nface++){
		int nDOFFace = faces[nface].getNDOF();  //number of DOF points for face
		int nDOFOwnerCell = faces[nface].getOwnerCell()->getNDOF();  //number of DOF points for owner cell

		faces[nface].getOwnerDOFPointsArray()->setSize(nDOFFace);
		faces[nface].getNeighbourDOFPointsArray()->setSize(nDOFFace);
		faces[nface].mapFlag = true;


		unsigned int numberOfPointsOr = 0;
		unsigned int numberOfPointsNr = 0;
		//first owner cells
		for (int DOF=0; DOF<nDOFFace; DOF++){
			double xf =faces[nface].getDOFPointsGlobalLocation()->getValue(DOF, 0);
			double yf =faces[nface].getDOFPointsGlobalLocation()->getValue(DOF, 1);
			double zf =faces[nface].getDOFPointsGlobalLocation()->getValue(DOF, 2);

			for (int c =0; c<nDOFOwnerCell; c++){
				double xc =faces[nface].getOwnerCell()->getDOFPointsGlobalLocation()->getValue(c, 0);
				double yc =faces[nface].getOwnerCell()->getDOFPointsGlobalLocation()->getValue(c, 1);
				double zc =faces[nface].getOwnerCell()->getDOFPointsGlobalLocation()->getValue(c, 2);

				if (abs(xf-xc) < SMALL and abs(yf-yc) < SMALL and abs(zf-zc) < SMALL){
					faces[nface].getOwnerDOFPointsArray()->setValue(DOF, c);
					numberOfPointsOr++;
					break;
				};
			};

			if (faces[nface].getNeighbourCell() != NULL){
				int nDOFNeighbourCell = faces[nface].getNeighbourCell()->getNDOF();  //number of DOF points for neighbour cell
				for (int c =0; c<nDOFNeighbourCell; c++){
					double xc =faces[nface].getNeighbourCell()->getDOFPointsGlobalLocation()->getValue(c, 0);
					double yc =faces[nface].getNeighbourCell()->getDOFPointsGlobalLocation()->getValue(c, 1);
					double zc =faces[nface].getNeighbourCell()->getDOFPointsGlobalLocation()->getValue(c, 2);

					if (abs(xf-xc) < SMALL and abs(yf-yc) < SMALL and abs(zf-zc) < SMALL){
						faces[nface].getNeighbourDOFPointsArray()->setValue(DOF, c);
						numberOfPointsNr++;
						break;
					};
				};
			}
			else{
				faces[nface].getNeighbourDOFPointsArray()->setValue(DOF, 0);
			};
		};

		assert(numberOfPointsOr == nDOFFace);
		if (faces[nface].getNeighbourCell() != NULL){
			assert(numberOfPointsNr == nDOFFace);
		}
		else{
			assert(numberOfPointsNr == 0);
		};
	};
};



void DG::mapFaceQuadPointsToCellQuadPoints(){
	for (Index nface=0; nface<this->noOfFaces; nface++){
		int nQuadFace = faces[nface].getNQuad();  //number of Quad points for face
		int nQuadOwnerCell = faces[nface].getOwnerCell()->getNQuad();  //number of Quad points for owner cell

		faces[nface].getOwnerQuadPointsArray()->setSize(nQuadFace);
		faces[nface].getNeighbourQuadPointsArray()->setSize(nQuadFace);
		faces[nface].mapFlag = true;


		unsigned int numberOfPointsOr = 0;
		unsigned int numberOfPointsNr = 0;
		//first owner cells
		for (int Quad=0; Quad<nQuadFace; Quad++){
			double xf =faces[nface].getQuadPointsGlobalLocation()->getValue(Quad, 0);
			double yf =faces[nface].getQuadPointsGlobalLocation()->getValue(Quad, 1);
			double zf =faces[nface].getQuadPointsGlobalLocation()->getValue(Quad, 2);

			for (int c =0; c<nQuadOwnerCell; c++){
				double xc =faces[nface].getOwnerCell()->getQuadPointsGlobalLocation()->getValue(c, 0);
				double yc =faces[nface].getOwnerCell()->getQuadPointsGlobalLocation()->getValue(c, 1);
				double zc =faces[nface].getOwnerCell()->getQuadPointsGlobalLocation()->getValue(c, 2);

				if (abs(xf-xc) < SMALL and abs(yf-yc) < SMALL and abs(zf-zc) < SMALL){
					faces[nface].getOwnerQuadPointsArray()->setValue(Quad, c);
					numberOfPointsOr++;
					break;
				};
			};

			if (faces[nface].getNeighbourCell() != NULL){
				int nQuadNeighbourCell = faces[nface].getNeighbourCell()->getNQuad();  //number of Quad points for neighbour cell
				for (int c =0; c<nQuadNeighbourCell; c++){
					double xc =faces[nface].getNeighbourCell()->getQuadPointsGlobalLocation()->getValue(c, 0);
					double yc =faces[nface].getNeighbourCell()->getQuadPointsGlobalLocation()->getValue(c, 1);
					double zc =faces[nface].getNeighbourCell()->getQuadPointsGlobalLocation()->getValue(c, 2);

					if (abs(xf-xc) < SMALL and abs(yf-yc) < SMALL and abs(zf-zc) < SMALL){
						faces[nface].getNeighbourQuadPointsArray()->setValue(Quad, c);
						numberOfPointsNr++;
						break;
					};
				};
			}
			else{
				faces[nface].getNeighbourQuadPointsArray()->setValue(Quad, 0);
			};
		};

		assert(numberOfPointsOr == nQuadFace);
		if (faces[nface].getNeighbourCell() != NULL){
			assert(numberOfPointsNr == nQuadFace);
		}
		else{
			assert(numberOfPointsNr == 0);
		};
	};
};



// Calculate face Jacobian (det(J)) array
void DG::calculateFaceJacobian(){
	for(Index nface=0; nface<this->noOfFaces; nface++){
		if(faces[nface].getFaceType() == FaceType::Quad){
			FunctionalSpace F(this->order, this->intType);
			int Np = this->order+1 + this->intType; //intType: 0->inexact, 1->exact
			int nQuad = pow((Np),2);  //number of quadrature points

			TensorO1<double> r(nQuad);
			TensorO1<double> s(nQuad);
			TensorO1<double> w(nQuad);
			
			faces[nface].JacobianFlag = true; 
			faces[nface].getJacobian()->setSize(nQuad);

			F.LGLRootsAndWeights2D(Np,Np, &r, &s, &w);
			
			// Main loop for calculating the Jacobian at each interpolating points
			for (int QP=0; QP<nQuad; QP++){	
				faces[nface].getJacobian()->setValue(QP, faces[nface].calculateJacobianQuadFace(r.getValue(QP),
							                                                        s.getValue(QP)));
			};

			faces[nface].JacobianFlag = true;
		}
		//Triangular Faces
		else if(faces[nface].getFaceType() == FaceType::Tri){
		};
	};
};


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 5. Flux and RES related ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//


void DG::computeFlux(int rkstep){

//#pragma omp parallel for
	// 1. Compute the flux for faces based on RK substep data. Riemann problem solution.
	//    This fills face.Flux array (called as face->getFlux(varNo, DOF));
	for(Index nface=0; nface<this->noOfFaces; nface++){

		(* solveRiemannProblem)(&this->faces[nface],  rkstep); 

		// solveRiemannSolver is a function pointer which points to an appropriate solver from src/RiemannSolver.cpp
		// The value of the function pointer is set in the function DG::assignRiemannSolver() 
		// We have to pass an argument as to which solver we want while calling DG::setFunctionalDetails(). 
		// All available solvers are given in include/Solvers.h

		// Refer src/RiemannSolver.cpp for implementation of the Riemann Solver (depending on which solver is called)
	};


	// 2. Compute the flux for Cell DOF locations. This is analytical value of flux at the DOF location.
	//    This fills the Cell.fluxVectorx, Cell.fluxVectory and Cell.fluxVectorz

//#pragma omp parallel for
	for (Index ncell=0; ncell< this->noOfCells; ncell++){
		(* findInternalFlux)(&this->cells[ncell], rkstep);
	};

	// 3. Compute the boundary fluxes
	this->computeBoundaryConditions(rkstep);
};

void DG::computeBoundaryConditions(int rkstep){

	TensorO1<double> primitiveVariableVector(5);
	TensorO1<double> fluxVectorX(5);
	TensorO1<double> fluxVectorY(5);
	TensorO1<double> fluxVectorZ(5);


	// Start loop for different boundary conditions
	for (Index bcondi=0; bcondi < noOfBoundaryConditions; bcondi++){
		int startFace = bcond[bcondi].getStartFace();
		int endFace = bcond[bcondi].getStartFace() + bcond[bcondi].getNoOfFaces();

		for (Index fi=startFace; fi < endFace; fi++){
			for(int vi = 0; vi < totNoOfVariable; vi++){
				// Explanation:
				// bcond[bcondi].func2bcond[vi] stores a function pointer to one of the boundary functions.
				// Boundary functions are described in BoundaryConditions.cpp
				// The following line calls that respective boundary condition function 

				(bcond[bcondi].*(bcond[bcondi].func2bcond[vi]))(&this->faces[fi],vi,rkstep);
			};


			if (systemType == System::Euler){
				for (int DOF=0; DOF< this->faces[fi].getNDOF(); DOF++){
					for (int nVar=0; nVar<5; nVar++){
						primitiveVariableVector.setValue(nVar, this->faces[fi].getPrimitiveVariable(nVar,DOF));
					};

					getEulerFluxFromPrimitiveVariables(&primitiveVariableVector, &fluxVectorX, &fluxVectorY, &fluxVectorZ);

					for (int nVar=0; nVar<5; nVar++){
						double FLUX = this->faces[fi].getNormal()->getValue(0)*fluxVectorX.getValue(nVar) + 
							      this->faces[fi].getNormal()->getValue(1)*fluxVectorY.getValue(nVar) + 
							      this->faces[fi].getNormal()->getValue(2)*fluxVectorZ.getValue(nVar);
						this->faces[fi].setFlux(nVar,DOF, FLUX); 
					};
				};
			};
		};
	};
};




void DG::computeRES(int rkstep){

	// _______________________________________________________________________________________
	// step 1: add the Riemann (and boundary) fluxes to the cell at appropriate DOF locations
	// _______________________________________________________________________________________

//#pragma omp parallel for 
	for(Index ncell=0; ncell<this->noOfCells; ncell++){
		// we are at cell ncell

		int noOfFaces = cells[ncell].getNoOfFaces();


		// Set initial flux for all DOF locations = 0.0
		cells[ncell].getFluxStarVector()->setAll(0.0);


		// add Riemann flux to appropriate location
		for (int f=0; f< noOfFaces; f++){
			// we are at face f of cell ncell
			

			Face *face = cells[ncell].getDefiningFace(f);

			double addFlag;

			// copy the mapQuadPoints array based on owner or neighbour cell
			TensorO1<int> * FCArray; //stores mapOwnerDOFPoints or neighbourDOFPoints depending on faceRelationshipArray.
			if (cells[ncell].faceRelationshipArray.getValue(f) == 0){
				// 0 means owner i.e. cells[ncell] is an owner cell of cells[ncell].faces[f]
				FCArray = face->getOwnerDOFPointsArray();
				addFlag = 1.0 * face->normalDirFlag; 
			}
			else{
				//1 means neighbour
				FCArray = face->getNeighbourDOFPointsArray();
				addFlag = -1.0 * face->normalDirFlag;
			}


			// copying the flux_star from face to cell array
			for (int DOF=0; DOF< face->getNDOF(); DOF++){
				// Get cell DOF value corresponding to the Face DOF point
				int DOFc = FCArray->getValue(DOF); 

				// Setting flux for the scalar problem
				if (systemType == System::scalarLinear or systemType == System::scalarBurger){
					cells[ncell].getFluxStarVector()->setValue(f, 0, DOFc,  addFlag*face->getFlux(0, DOF));  
					// 0 indicates scalar flux
				} 
				// Setting flux for the system of equations
				else if (systemType == System::Euler or systemType == System::Heat or systemType == System::NavierStokes){
					for (int nVar=0; nVar<totNoOfVariable; nVar++){
						double fluxValue = addFlag*face->getFlux(nVar, DOF);
						cells[ncell].getFluxStarVector()->addValue(f, nVar, DOFc, fluxValue); 
					};
				};
				// ADD other equations here e.g. heat equation
			};
		};
	};

	// _______________________________________________________________________________________
	//Step 2: Compute the Right hand side containing all the spatial data.
	// _______________________________________________________________________________________

	// actual RES computation for each cell

//#pragma omp parallel for 
	for (Index ncell=0; ncell< this->noOfCells; ncell++){
		int nDOF = cells[ncell].getNDOF();
		//we are at cell i

		// Creating local variable arrays for processing declared inside the loop for openmp's sake			
		// Initialize res with zero
		TensorO2<double> rhs1x(totNoOfVariable, nDOF);        rhs1x.setAll(0.0);
		TensorO2<double> rhs1y(totNoOfVariable, nDOF);        rhs1y.setAll(0.0);
		TensorO2<double> rhs1z(totNoOfVariable, nDOF);        rhs1z.setAll(0.0);
		TensorO2<double> rhs1(totNoOfVariable, nDOF);         rhs1.setAll(0.0); 
		TensorO2<double> rhs2(totNoOfVariable, nDOF);         rhs2.setAll(0.0); 

		// actual loops over variables and DOFs start here:

		for (int varNo=0; varNo<totNoOfVariable; varNo++){
			//we are dealing with variable varNo 

			// First analytical flux is multiplied by D matrix 
			for (int DOF=0; DOF<cells[ncell].getNDOF(); DOF++){

				// Confirm initialization with zero
				rhs1x.setValue(varNo, DOF, 0.0);
				rhs1y.setValue(varNo, DOF, 0.0);
				rhs1z.setValue(varNo, DOF, 0.0);

				// Matrix-vector multiplication
				for (int k=0; k<cells[ncell].getNDOF(); k++){
					// Dot product of Dtilde and internalFluxVector
					double valuex = cells[ncell].getDtildexMatrix()->getValue(DOF, k) * cells[ncell].getFluxVectorx()->getValue(varNo, k);
					double valuey = cells[ncell].getDtildeyMatrix()->getValue(DOF, k) * cells[ncell].getFluxVectory()->getValue(varNo, k);
					double valuez = cells[ncell].getDtildezMatrix()->getValue(DOF, k) * cells[ncell].getFluxVectorz()->getValue(varNo, k);
					rhs1x.addValue(varNo, DOF, valuex); 
					rhs1y.addValue(varNo, DOF, valuey) ; 
					rhs1z.addValue(varNo, DOF, valuez) ; 
				};
				rhs1.setValue(varNo, DOF, rhs1x.getValue(varNo,DOF) + 
						          rhs1y.getValue(varNo,DOF) + 
							  rhs1z.getValue(varNo,DOF));
			};
		};



		// Next, Riemann fluxes are multiplied by the Flux matrix and added (remember surface integration is for a closed cell)

		int noOfFaces = cells[ncell].getNoOfFaces();
		rhs2.setAll(0.0);

		for (int varNo=0; varNo<totNoOfVariable; varNo++){
			// for a particular variable
			for (int DOF=0; DOF<cells[ncell].getNDOF(); DOF++){
				// Initialize with zero
				rhs2.setValue(varNo,DOF, 0.0);
				for (int f=0; f< noOfFaces ; f++){
					// at face f of cell ncell
					double dummy = 0;
					for (int k=0; k<cells[ncell].getNDOF(); k++){
						// Summation of (Dot product of Ftilde and F_star)
						dummy += cells[ncell].getFtildeMatrix()->getValue(f,DOF,k) * 
							 cells[ncell].getFluxStarVector()->getValue(f,varNo,k) ; 
					};

					rhs2.addValue(varNo,DOF, dummy);
				};
			};
		};

		// assembly of the final residual (it should be Dtilde*F - Ftilde*FStar)

		cells[ncell].getVariableResidualArray()->setAll(0.0);
		for (int varNo=0; varNo<totNoOfVariable; varNo++){
			// final RHS computation 
			for (int DOF=0; DOF<cells[ncell].getNDOF(); DOF++){
				cells[ncell].getVariableResidualArray()->setValue(rkstep, varNo, DOF, rhs1.getValue(varNo,DOF) - 
						                                                 rhs2.getValue(varNo,DOF)) ;
			};
		};
	};
};










//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 6. Time integrator related ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//


void DG:: calculateCFL(){
	this->CFL = 0.0;
	this->charSpeed = 0.0;

	TensorO1<double> consVarVector(5);
	TensorO1<double> primVarVector(5);

	for (Index ncell =0; ncell<noOfCells; ncell++){
		Cell *cell = &cells[ncell];
		
		for (Index DOF =0; DOF < cell->getNDOF(); DOF++){
			consVarVector.setAll(0.0);
			primVarVector.setAll(0.0);
			for (Index nVar=0; nVar<5; nVar++){
				consVarVector.setValue(nVar, cell->getVariable(0, nVar, DOF));
			};
			
			getPrimitiveVariableVector(&consVarVector, &primVarVector);

			double rho = primVarVector.getValue(0);
			double u = primVarVector.getValue(1);
			double v = primVarVector.getValue(2);
			double w = primVarVector.getValue(3);
			double p = primVarVector.getValue(4);
			double c   = sqrt(GAMMA* p/rho);

			double speedMax = max(abs(u), max(abs(v), abs(w)));

			if (this->charSpeed < speedMax + c){
				this->charSpeed = speedMax+c;
			};
		};
	};

	this->CFL = this->charSpeed * this->deltaT / this->charLength; 
};



void DG:: runApplication(){
	
	// Main function to run the simulation
	// includes integrator and Riemann solver calls
	
	// Pre-processing
	computeCellMatrix();			// Calculating all important matrices

	assignSystemOfEquations();		//assign system of equations such as scalar, euler etc.
						//refer include/systems.h for available systems

	assignRiemannSolver(); 			// Assigns the Riemann solver from file src/RiemannSolver.cpp 
						// first we have to give argument of which Riemann solver in DG::setFunctionalDetails() 
	
	//applyFilterVariable(0);
	// Initializing temporary time data
	double nowTime = this->startTime;
	double nowTimePrint = this->deltaT/2.0;

	// Main time loop
	int IterationNumber = 0;

	cout << "\nRunning Time Loop: " << endl;
	double tstart = omp_get_wtime();

	string fileName = path + "/" + "Error.dat";
	ofstream file(fileName.c_str());

	file << this->cells[0].getNDOF()*this->noOfCells << " " << this->cells[0].getNQuad()*this->noOfCells << " " << 0  << " " << 0 << " " << 0 << endl;   

	while (nowTime < this->endTime){
	//for (int i = 0; i< 1; i++){
		
		if (abs(nowTime - this->endTime) < this->deltaT){
			this->deltaT = abs(nowTime - this->endTime);
		};
		
		// Call Integrator
		
		this->calculateCFL();
		assert(this->CFL <= 1 && "DG::runApplication returns error. CFL number exceeding 1. Solution becoming unstable.\n");

		this->integratorRK3();
		// update time
		nowTimePrint += this->deltaT;

		nowTime += this->deltaT;
		IterationNumber ++;

		string message = "  CFL: " + varToString(this->CFL);
		printProgressBar(nowTime / this->endTime * 100.0, message);

		// Data writing if condition
		if(nowTimePrint >= this->printTime){ 
			// if nowTime is greater than printTime then print data
			this->integrateCellVariable();				// get center value first
			this->writeVariableArray(nowTime);		
			this->writeVariableRefArray(nowTime);		
			this->computeError();
			nowTimePrint = this->deltaT/2.0;

			double runTime = omp_get_wtime() - tstart;
			file << nowTime << " " << this->error[0] << " " << this->error[1] << " " << this->error[2] << " " << runTime << endl;
		}
		else{
		};
	};
	double tend = omp_get_wtime();
	double exTime = tend - tstart;
	cout << "\nTime (RK3) iterations complete. Results written in " << this->path << endl;
	cout << "exTime: " << exTime << endl;
	cout << endl;

	file.close();
};



// Integrator RK3 
// As of now not added any defitions, as face is already initialized in DG and point to face list
void DG:: integratorRK3(){
	
	// First Step
	//applyFilterVariable(0);
	this->computeFlux(0);					// Uncomment 
	computeRES(0);
	//applyFilterResidual(0);

	if (this->systemType == System::scalarLinear or this->systemType == System::scalarBurger){
		totNoOfVariable = 1;
	};


//#pragma omp parallel for collapse(3) 
	for(Index ncell=0; ncell<this->noOfCells; ncell++){	
		for(int varNo = 0; varNo < totNoOfVariable; ++varNo){
			for(int DOF = 0; DOF < cells[ncell].getNDOF(); ++DOF){
				double oldValue =         cells[ncell].getVariable(0, varNo, DOF);
				double resValue = cells[ncell].getVariableResidualArray()->getValue(0, varNo, DOF);
				double newValue = oldValue + this->deltaT * resValue;
				cells[ncell].getVariableArray()->setValue(1, varNo, DOF, newValue );
			};
		};
	};

	// Second Step
	//applyFilterVariable(1);
	this->computeFlux(1);					// Uncomment 
	computeRES(1);

//#pragma omp parallel for collapse(3) 
	for(Index ncell=0; ncell<this->noOfCells; ncell++){	
		for(int varNo = 0; varNo < totNoOfVariable; ++varNo){
			for(int DOF = 0; DOF < cells[ncell].getNDOF(); ++DOF){
				double oldValue0 =         cells[ncell].getVariable(0, varNo, DOF);
				double oldValue1 =         cells[ncell].getVariable(1, varNo, DOF);
				double resValue  = cells[ncell].getVariableResidualArray()->getValue(1, varNo, DOF);
				double newValue = 0.75*oldValue0 + 0.25*oldValue1 + 0.25*this->deltaT*resValue;

				cells[ncell].getVariableArray()->setValue(2, varNo, DOF, newValue);
				//cells[ncell].variable[2][varNo][DOF] = 0.75*cells[ncell].variable[0][varNo][DOF] + 0.25*cells[ncell].variable[1][varNo][DOF] + 0.25*this->deltaT*cells[ncell].variableResidual[1][varNo][DOF];
			}
		}
	}

	// Third and Final Step
	//applyFilterVariable(2);
	this->computeFlux(2);					// Uncomment 
	computeRES(2);
	//applyFilterResidual(2);

//#pragma omp parallel for collapse(3) 
	for(Index ncell=0; ncell<this->noOfCells; ncell++){	
		for(int varNo = 0; varNo < totNoOfVariable; ++varNo){
			for(int DOF = 0; DOF < cells[ncell].getNDOF(); ++DOF){
				double oldValue0 =         cells[ncell].getVariable(0, varNo, DOF);
				double oldValue2 =         cells[ncell].getVariable(2, varNo, DOF);
				double resValue  = cells[ncell].getVariableResidualArray()->getValue(2, varNo, DOF);
				double newValue = 1.0/3.0 * oldValue0  + 2.0/3.0 * oldValue2 + 2.0/3.0 * this->deltaT * resValue;

				cells[ncell].getVariableArray()->setValue(0, varNo, DOF, newValue);
				
				//cells[ncell].variable[0][varNo][DOF] = cells[ncell].variable[0][varNo][DOF]/3. + 2.*cells[ncell].variable[2][varNo][DOF]/3. + 2.*this->deltaT*cells[ncell].variableResidual[2][varNo][DOF]/3.;
			}
		}
	}

};

