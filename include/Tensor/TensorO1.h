template<typename numType>
class TensorO1{
	public:
		/// Constructor
		TensorO1(unsigned int size = 0);

		/// Destructor
		~TensorO1();

		// Data members

		/// Order of the tensor 
		static const int order = 1;	

		/// Size of the O1 tensor (i.e. dimension)
		unsigned int size;	        

		/// Value of the vector (of datatype numType)
		numType *value;			

		/// Returns 'true' if space has been allocated for value
		bool sizeFlag;			

		/** norm of a vector **/
		struct Norm{
			double L1, L2, Linf;
		};
		Norm norm;

		// Methods related to size and accessing functions

		/// Get method for size 
		unsigned int getSize()const; 

		/// Set the size of the vector 
		void setSize(unsigned int size);

		/// Get value at nth location
		numType getValue(unsigned int n)const;

		/// Set value at nth location
		void setValue(unsigned int n, numType N);

		/// Add a given value to previously existing value at nth location
		void addValue(unsigned int n, numType N);

		/// Set all the values to an argument value
		void setAll(numType N);



		// Methods related to algebra

		/// Computes L1, L2 and Linf norm of the data
		void computeNorm();   	

		/// Get L1 norm 
		double getL1Norm()const;	

		/// Get L2 norm 
		double getL2Norm()const;	

		/// Get Linf norm 
		double getLinfNorm()const;	
		
		void randomize();

	private:
		// Prevent copy constructor
		TensorO1(const TensorO1&);
		// copy assignment operator 
		TensorO1& operator=(const TensorO1& tmp_obj){ 
		        return *this; 
	  	};
};


// Constructor
template<typename numType>
TensorO1<numType>::TensorO1(unsigned int size){
	this->size = size;
	this->value = new(nothrow) numType[size];
	if (! this->value){
		cout << "Memory allocation failed for TensorO1." << endl;
	};
	sizeFlag = true;
	this->norm.L1=0.0;
	this->norm.L2=0.0;
	this->norm.Linf=0.0;
};


// Destructor
template<typename numType>
TensorO1<numType>::~TensorO1(){
	delete[] this->value;
};



template<typename numType>
unsigned int TensorO1<numType>::getSize()const{
	return this->size;
};

template<typename numType>
void TensorO1<numType>::setSize(unsigned int size){
	// first delete the existing value array
	delete[] this->value;
	sizeFlag = false;

	// then set the new size and array
	this->size = size;
	this->value = new(nothrow) numType[size];
	if (! this->value){
		cout << "Memory allocation failed for TensorO1." << endl;
	};
	sizeFlag = true;
};

template<typename numType>
numType TensorO1<numType>::getValue(unsigned int n) const{
	assert(sizeFlag == true);
	assert(n<this->size && "Index exceeds the size of the vector.");
	return this->value[n];
};

template<typename numType>
void TensorO1<numType>::setValue(unsigned int n, numType N){
	assert(sizeFlag == true);
	assert(n<this->size && "Index exceeds the size of the vector.");
	this->value[n] = N;
};


template<typename numType>
void TensorO1<numType>::addValue(unsigned int n, numType N){
	assert(sizeFlag == true);
	assert(n<this->size && "Index exceeds the size of the vector.");
	this->value[n] += N;
};

template<typename numType>
void TensorO1<numType>::setAll(numType N){
	assert(sizeFlag == true);
	for (Index i=0; i<this->getSize(); i++){
		this->setValue(i, N);
	};
};

/// L1 norm of the function: sum(abs(members))
template<typename numType>
double TensorO1<numType>::getL1Norm()const{
	assert(this->size > 0 && "TensorO1::getL1Norm() can not process an empty array.");
	double norm = 0.0;
	for(Index i=0; i< this->size ; i++){
		norm += abs(this->value[i]);
	};
	return norm;
};

/// L2 norm of the function: sqrt(sum(sqr(members)))
template<typename numType>
double TensorO1<numType>::getL2Norm()const{
	assert(this->size > 0 && "TensorO1::getL1Norm() can not process an empty array.");
	double norm = 0.0;
	for(Index i=0; i< this->size ; i++){
		norm += abs(this->value[i] * this->value[i]);
	};
	return sqrt(norm);
};

/// Linf norm of the function: max(abs(members))
template<typename numType>
double TensorO1<numType>::getLinfNorm()const{
	assert(this->size > 0 && "TensorO1::getL1Norm() can not process an empty array.");
	double norm = 0.0;
	for(Index i=0; i< this->size ; i++){
		if (norm < abs(this->value[i])){ 
			norm = abs(this->value[i]);
		};
	};
	return 1.0*norm;
};


/// Randomize (shuffle) entries of the O1 Tensor object
template<typename Type>
void TensorO1<Type> :: randomize(){
	assert(this->size > 0 && "TensorO1::randomize() can not process an empty array.");
	unsigned int Size = this->getSize();

	for (Index iter=0; iter <100; iter ++){
		for (Index i=Size-1; i>0; i--){
			int j = rand()%(i+1);
			Type dummy = this->getValue(i);

			this->setValue(i, this->getValue(j));
			this->setValue(j, dummy);
		};
	};
};


/// Computes all three norms of the functions
template<typename numType>
void TensorO1<numType>::computeNorm(){
	this->norm.L1 = getL1Norm();
	this->norm.L2 = getL2Norm();
	this->norm.Linf = getLinfNorm();
};












