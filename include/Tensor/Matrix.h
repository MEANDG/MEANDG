template<typename numType>
class Matrix:public TensorO2<numType>{
	public:
		// Constructor
		Matrix(unsigned int size=1);
		~Matrix();
		// Data members
		static const int order = 2;	///< Order of the tensor
		unsigned int size[2];	        ///< Size of the Matrix. firstIndex x secondIndex
		numType **value;		///< Value of the Matrix (of datatype numType)
		double  **inverse;		///< Inverse of the Matrix.
		bool inverseFlag;		///< Sets to 1 if inverse if computed. 



		// Methods related to size and accessing functions
		unsigned int getSize0()const; 	
		unsigned int getSize1()const; 
		void setSize(unsigned int matsize);
		numType getValue(unsigned int n, unsigned int m)const;
		void setValue(unsigned int n, unsigned int m, numType N);
		void addValue(unsigned int n, unsigned int m, numType N);
		void setAll(numType N);

		// Methods related to algebra of matrix
		void invert();			///< Inverts the matrix and stores the entries in inverse variable
		void getInverse(Matrix<numType> *MI);
	private:
		// Prevent copy constructor
		Matrix(const Matrix&);
		// copy assignment operator 
		Matrix& operator=(const Matrix& tmp_obj){ 
		        return *this; 
	  	};
};


// Constructor
template<typename numType>
Matrix<numType>::Matrix(unsigned int matsize):TensorO2<numType>(matsize, matsize){
	this->size[0] = matsize; 
	this->size[1] = matsize;

	// assigning space for the variable array 'value'
	this->value = new(nothrow) numType*[size[0]];
	if (! this->value) cout << "Allocation of memory failed in Matrix." << endl;

	for (Index i=0; i<size[0]; i++){
		this->value[i] = new(nothrow) numType[size[1]];
		if (! this->value[i]) cout << "Allocation of memory failed in Matrix." << endl;
	};
	this->inverseFlag = false;
};


// Destructor
template<typename numType>
Matrix<numType>::~Matrix(){
	for (Index i=0; i<size[0]; i++){
		delete[] this->value[i];
	};
	delete[] this->value;

	if (inverseFlag == true){
		// free memory of the inverse matrix
		for (Index i=0; i<size[0]; i++){
			delete[] this->inverse[i];
		};
		delete[] this->inverse;
	};
};



/// Returns number of rows of the matrix
template<typename numType>
unsigned int Matrix<numType>::getSize0()const{
	assert(this->size[0] == this->size[1] && "Matrix does not have same dimensions for rows and columns.");
	return this->size[0];
};

/// Returns number of columns of the matrix
template<typename numType>
unsigned int Matrix<numType>::getSize1()const{
	assert(this->size[1] == this->size[0] && "Matrix does not have same dimensions for rows and columns.");
	return this->size[1];
};

/// Sets number of rows and columns of the matrix
template<typename numType>
void Matrix<numType>::setSize(unsigned int matsize){
	// first delete the value array 
	for (Index i=0; i<size[0]; i++){
		delete[] this->value[i];
	};
	delete[] this->value;


	// then create a new array
	this->size[0] = matsize;
	this->size[1] = matsize;

	this->value = new(nothrow) numType*[size[0]];
	if (! this->value) cout << "Allocation of memory failed in Matrix." << endl;
	for (Index i=0; i<size[0]; i++){
		this->value[i] = new(nothrow) numType[size[1]];
		if (! this->value[i]) cout << "Allocation of memory failed in Matrix." << endl;
	};
};

/// Returns the value at requested coordinates
template<typename numType>
numType Matrix<numType>::getValue(unsigned int n, unsigned int m) const{
	assert(n<this->size[0] and m<this->size[1]  && "Index exceeds the size of the vector.");
	return this->value[n][m];
};

/// Sets value at the requested coordinates (value passed as the last argument)
template<typename numType>
void Matrix<numType>::setValue(unsigned int n, unsigned int m, numType N){
	assert(n<this->size[0]  and m<this->size[1] && "Index exceeds the size of the vector.");
	this->value[n][m] = N;
};

/// Adds value at the requested coordinates (value passed as the last argument)
template<typename numType>
void Matrix<numType>::addValue(unsigned int n, unsigned int m, numType N){
	assert(n<this->size[0]  and m<this->size[1] && "Index exceeds the size of the vector.");
	this->value[n][m] += N;
};

/// Sets value at all coordinates (value passed as the argument)
template<typename numType>
void Matrix<numType>::setAll(numType N){
	for (Index n=0; n<this->getSize0(); n++){
		for (Index m=0; m<this->getSize1(); m++){
			this->setValue(n, m, N);
		};
	};
};


/// Stores inverse of the matrix in the inverse variable
template<typename numType>
void Matrix<numType>::invert(){
	assert(this->inverseFlag == false && "Matrix::invert() returns error. Invert called more than once.\n");
	// first make sure that the matrix is initiated
	assert(this->size[0] == this->size[1] and this->size[0]>0 && "Matrix::invert() can not invert a matrix. Make sure that the matrix is initialized properly and has the same number of rows and columns.");

	// next, allocate space for the inverse variable
	this->inverse = new(nothrow) numType*[size[0]];
	if (! this->inverse) cout << "Allocation of memory failed in Matrix::invert()." << endl;
	for (Index i=0; i<size[0]; i++){
		this->inverse[i] = new(nothrow) numType[size[1]];
		if (! this->inverse[i]) cout << "Allocation of memory failed in Matrix::invert()." << endl;
	};
	this->inverseFlag = true; // set the inverse flag to true, which will be used to free the memory in the destructor

	// now initialize an Eigen matrix
	Eigen::MatrixXd EValue(this->size[0], this->size[1]);
	Eigen::MatrixXd EValueI;
	for (Index i=0; i<this->size[0]; i++){
		for (Index j=0; j<this->size[1]; j++){
			EValue(i,j) = this->getValue(i,j);
		};
	};

	EValueI = EValue.inverse();

	//copy this value into the inverse variable
	for (Index i=0; i<this->size[0]; i++){
		for (Index j=0; j<this->size[1]; j++){
			this->inverse[i][j] = EValueI(i,j);
		};
	};
};

/// Returns the inverse of the matrix 
template<typename numType>
void Matrix<numType>::getInverse(Matrix<numType> *MI){
	// perform inversion only once
	if (this->inverseFlag == false){
		// First make sure that the result matrix is of appropriate dimensions
		assert(this->size[0] == MI->size[0] and this->size[1] == MI->size[1] && "Matrix::getInverse() can not invert this matrix. Make sure that the matrix is initialized properly and has the same number of rows and columns.");

		// Next invert the matrix
		this->invert();


	}


	// Next copy the values into the result matrix

	for (Index i=0; i<this->size[0]; i++){
		for (Index j=0; j<this->size[0]; j++){
			MI->setValue(i,j, this->inverse[i][j]);
		};
	};
};
