template<typename numType>
class TensorO3{
	public:
		/// Constructor
		TensorO3(unsigned int r=0, unsigned int s=0, unsigned int t=0);
			
		/// Destructor
		~TensorO3();

		// Data members


		/// Order of the tensor
		static const int order = 3;	

		/// Size of the 3D Tensor. firstIndex x secondIndex x thirdIndex 
		unsigned int size[3];	        

		/// Value of the 3D Tensor (of datatype numType)
		numType ***value;

		// Methods related to size and accessing functions

		/// Get size of 0 dim
		unsigned int getSize0()const; 

		/// Get size of 1 dim
		unsigned int getSize1()const; 

		/// Get size of 2 dim
		unsigned int getSize2()const; 

		/// Set size of all the dims
		void setSize(unsigned int r, unsigned int s, unsigned int t);

		/// Get value at r,s,t coordinates
		numType getValue(unsigned int r, unsigned int s, unsigned int t)const;

		/// Set value of r,s,t coordinate
		void setValue(unsigned int r, unsigned int s, unsigned int t, numType N);

		/// Add value to previously existing value at r,s,t coordinate
		void addValue(unsigned int r, unsigned int s, unsigned int t, numType N);
		
		/// Set all the values same as the argument value
		void setAll(numType N);

	private:
		// Prevent copy constructor
		TensorO3(const TensorO3&);
		// copy assignment operator 
		TensorO3& operator=(const TensorO3& tmp_obj){ 
		        return *this; 
	  	};
};


// Constructor
template<typename numType>
TensorO3<numType>::TensorO3(unsigned int r, unsigned int s, unsigned int t){
	this->size[0] = r; 
	this->size[1] = s;
	this->size[2] = t;

	this->value = new(nothrow) numType**[size[0]];
	if (! this->value) cout << "Allocation of memory failed in TensorO2." << endl;

	for (Index i=0; i<size[0]; i++){
		this->value[i] = new(nothrow) numType*[size[1]];
		if (! this->value[i]) cout << "Allocation of memory failed in TensorO2." << endl;

		for (Index j=0; j<size[1]; j++){
			this->value[i][j] = new(nothrow) numType[size[2]];
			if (! this->value[i][j]) cout << "Allocation of memory failed in TensorO2." << endl;
		};
	};
};


// Destructor
template<typename numType>
TensorO3<numType>::~TensorO3(){
	for (Index i=0; i<size[0]; i++){
		for (Index j=0; j<size[1]; j++){
			delete[] this->value[i][j];
		};
		delete[] this->value[i];
	};
	delete[] this->value;
};



template<typename numType>
unsigned int TensorO3<numType>::getSize0()const{
	return this->size[0];
};

template<typename numType>
unsigned int TensorO3<numType>::getSize1()const{
	return this->size[1];
};

template<typename numType>
unsigned int TensorO3<numType>::getSize2()const{
	return this->size[2];
};

template<typename numType>
void TensorO3<numType>::setSize(unsigned int r, unsigned int s, unsigned int t){
	// first delete the value array
	for (Index i=0; i<size[0]; i++){
		for (Index j=0; j<size[1]; j++){
			delete[] this->value[i][j];
		};
		delete[] this->value[i];
	};
	delete[] this->value;


	// then set the new array
	this->size[0] = r;
	this->size[1] = s;
	this->size[2] = t;

	this->value = new(nothrow) numType**[size[0]];
	if (! this->value) cout << "Allocation of memory failed in TensorO2." << endl;

	for (Index i=0; i<size[0]; i++){
		this->value[i] = new(nothrow) numType*[size[1]];
		if (! this->value[i]) cout << "Allocation of memory failed in TensorO2." << endl;

		for (Index j=0; j<size[1]; j++){
			this->value[i][j] = new(nothrow) numType[size[2]];
			if (! this->value[i][j]) cout << "Allocation of memory failed in TensorO2." << endl; 
		};
	};
};

template<typename numType>
numType TensorO3<numType>::getValue(unsigned int r, unsigned int s, unsigned int t) const{
	assert(r<this->size[0] and s<this->size[1] and t<this->size[2]  && "Index exceeds the size of the vector.");
	return this->value[r][s][t];
};

template<typename numType>
void TensorO3<numType>::setValue(unsigned int r, unsigned int s, unsigned int t, numType N){
	assert(r<this->size[0]  and s<this->size[1]  and t<this->size[2] && "Index exceeds the size of the vector.");
	this->value[r][s][t] = N;
};

template<typename numType>
void TensorO3<numType>::addValue(unsigned int r, unsigned int s, unsigned int t, numType N){
	assert(r<this->size[0]  and s<this->size[1]  and t<this->size[2] && "Index exceeds the size of the vector.");
	this->value[r][s][t] += N;
};


template<typename numType>
void TensorO3<numType>::setAll(numType N){
	for (Index r=0; r<this->getSize0(); r++){
		for (Index s=0; s<this->getSize1(); s++){
			for (Index t=0; t<this->getSize2(); t++){
				this->setValue(r, s, t, N);
			};
		};
	};
};

