template<typename numType>
class TensorO2{
	public:
		/// Constructor
		TensorO2(unsigned int n=0, unsigned int m=0);

		/// Destructor
		~TensorO2();




		// Data members

		/// Order of the tensor 
		static const int order = 2;	

		/// Size of the 2D Tensor. firstIndex x secondIndex (i.d. dimension of the product space)
		unsigned int size[2];	        

		/// Value of the 2D Tensor (of datatype numType)
		numType **value;		




		// Methods related to size and accessing functions

		/// Get size in 0 dim 
		unsigned int getSize0()const; 

		/// Get size in 1 dim
		unsigned int getSize1()const; 

		/// Set size of the vector (row, column)
		void setSize(unsigned int n, unsigned int m);

		/// Get value at the n,m coordinate
		numType getValue(unsigned int n, unsigned int m)const;

		/// Set value at the n,m coordinate 
		void setValue(unsigned int n, unsigned int m, numType N);

		/// Add value to previously existing value at n,m location
		void addValue(unsigned int n, unsigned int m, numType N);

		/// Set all the values equal to the argument value
		void setAll(numType N);

	private:
		// Prevent copy constructor
		TensorO2(const TensorO2&);
		// copy assignment operator 
		TensorO2& operator=(const TensorO2& tmp_obj){ 
		        return *this; 
	  	};

};


// Constructor
template<typename numType>
TensorO2<numType>::TensorO2(unsigned int n, unsigned int m){
	this->size[0] = n; 
	this->size[1] = m;

	this->value = new(nothrow) numType*[size[0]];
	if (! this->value) cout << "Allocation of memory failed in TensorO2." << endl;
	for (Index i=0; i<size[0]; i++){
		this->value[i] = new(nothrow) numType[size[1]];
		if (! this->value[i]) cout << "Allocation of memory failed in TensorO2." << endl;
	};
};


// Destructor
template<typename numType>
TensorO2<numType>::~TensorO2(){
	for (Index i=0; i<size[0]; i++){
		delete[] this->value[i];
	};
	delete[] this->value;
};



template<typename numType>
unsigned int TensorO2<numType>::getSize0()const{
	return this->size[0];
};

template<typename numType>
unsigned int TensorO2<numType>::getSize1()const{
	return this->size[1];
};

template<typename numType>
void TensorO2<numType>::setSize(unsigned int n, unsigned int m){
	// first delete the value array 
	for (Index i=0; i<size[0]; i++){
		delete[] this->value[i];
	};
	delete[] this->value;


	// then create a new array
	this->size[0] = n;
	this->size[1] = m;

	this->value = new(nothrow) numType*[size[0]];
	if (! this->value) cout << "Allocation of memory failed in TensorO2." << endl;
	for (Index i=0; i<size[0]; i++){
		this->value[i] = new(nothrow) numType[size[1]];
		if (! this->value[i]) cout << "Allocation of memory failed in TensorO2." << endl;
	};
};

template<typename numType>
numType TensorO2<numType>::getValue(unsigned int n, unsigned int m) const{
	assert(n<this->size[0] and m<this->size[1]  && "Index exceeds the size of the vector.");
	return this->value[n][m];
};

template<typename numType>
void TensorO2<numType>::setValue(unsigned int n, unsigned int m, numType N){
	assert(n<this->size[0]  and m<this->size[1] && "Index exceeds the size of the vector.");
	this->value[n][m] = N;
};

template<typename numType>
void TensorO2<numType>::addValue(unsigned int n, unsigned int m, numType N){
	assert(n<this->size[0]  and m<this->size[1] && "Index exceeds the size of the vector.");
	this->value[n][m] += N;
};

template<typename numType>
void TensorO2<numType>::setAll(numType N){
	for (Index n=0; n<this->getSize0(); n++){
		for (Index m=0; m<this->getSize1(); m++){
			this->setValue(n,m, N);
		};
	};
};

