/** @file */
#include"Tensor.h"
using namespace std;

// Note: 
// The class Tensor is a template class, hence declared completely in the header file include/Tensor.h
// It further includes the template classes include/Tensor/TensorO1.h, include/Tensor/TensorO2.h, include/Tensor/TensorO3.h, include/Tensor/Matrix.h 
// This file only contains the unit tests for tensor class (declared in include/Tensor.h and called in Test.cpp)

namespace Test{
	void TestTensor(){
		/// Testing TensorO1 class
		TestTensorO1();
		TestTensorO2();
		TestMatrix();  		//special type of TensorO2 with nRow = nCol
		TestTensorO3();
	};

	void TestTensorO1(){
		cout << "\nRunning test on class TensorO1...  " << endl;
		// testing with positive intigers

		// initialize a random sequence
		// take a random array size between [1, 100]
		int size = getRandom(10,100);
		// create a new TensorO1 (of dtype int) of this size
		TensorO1<int> A(size);

		unt N = A.getSize();  //test getSize() function
		assert (N == size && "TensorO1::getSize() does not return correct size");

		// create a new random variable for setSize() trial
		int newSize = getRandom(10,100);
		A.setSize(newSize);	//set new size
		// test the new size
		N = A.getSize();
		assert (N == newSize && "TensorO1::setSize() does not set correct size");

		// loop for testing setValue and getValue functions
		for (Index i=0; i<N; i++){
			A.value[i] = i;
			int value = A.getValue(i);
			assert (value==i && "TensorO1::getValue() does not return the correct value");

			int newValue = getRandom(1,10000001);
			A.setValue(i, newValue); //set ith variable to a new value i.e. each variable in A.value

			value = A.getValue(i);
			assert (value==newValue && "TensorO1::setValue() does not set the correct value");



		};



		// same steps followed for testing double
		// testing with positive doubles
		size = getRandom(10,100);
		TensorO1<double> B(size);

		N = B.getSize();
		assert (N == size && "TensorO1::getSize() does not return correct size");

		newSize = getRandom(10,100); 
		B.setSize(newSize);
		N = B.getSize();
		assert (N == newSize && "TensorO1::setSize() does not set correct size");

		for (Index i=0; i<N; i++){
			double value1 = getRandom(0.0000001, 99999.999999);
			B.value[i] = value1;
			double value = B.getValue(i);
			assert (value== value1 && "TensorO1::getValue() does not return the correct value");

			double value2 =  getRandom(0.0000001, 99999.999999);
			B.setValue(i, value2);

			value = B.getValue(i);
			assert (value== value2 && "TensorO1::setValue() does not set the correct value");

			double rand = getRandom(-100.0, 100.0);
			B.addValue(i, rand);

			value = B.getValue(i);
			assert (value== value2+rand && "TensorO1::addValue() does not add the correct value");
		};

		double value4 = getRandom(-1000.11, 10000.11);
		B.setAll(value4);
		for (Index i=0; i<B.getSize(); i++){
			assert(B.getValue(i) == value4 && "TensorO1::setAll() does not set correct value.\n" );
		};


		//similar process is used for testing with negative integers
		//testing with negative integers
		size = getRandom(10,100);
		TensorO1<int> C(size);

		N = C.getSize();
		assert (N == size && "TensorO1::getSize() does not return correct size");

		newSize = getRandom(10,100);
		C.setSize(newSize);
		N = C.getSize();
		assert (N == newSize && "TensorO1::setSize() does not set correct size");


		for (Index i=0; i<N; i++){
			int value1 = getRandom(-999,-1);
			C.value[i] = value1;
			int value = C.getValue(i);
			assert (value== value1 && "TensorO1::getValue() does not return the correct value");

			int value2 = getRandom(-999,-1);
			C.setValue(i, value2);

			value = C.getValue(i);
			assert (value== value2 && "TensorO1::setValue() does not set the correct value");
		};

		// testing with negative doubles

		size = getRandom(10,100);
		TensorO1<double> D(size);

		N = D.getSize();
		assert (N == size && "TensorO1::getSize() does not return correct size");

		newSize = getRandom(10,100); 
		D.setSize(newSize);
		N = D.getSize();
		assert (N == newSize && "TensorO1::setSize() does not set correct size");

		for (Index i=0; i<N; i++){
			double value1 = getRandom(-999.9999, -0.001);
			D.value[i] = value1;
			double value = D.getValue(i);

			assert (value== value1 && "TensorO1::getValue() does not return the correct value");

			double value2 = getRandom(-999.999, -0.001);
			D.setValue(i, value2);

			value = D.getValue(i);
			assert (value== value2 && "TensorO1::setValue() does not set the correct value");
		};

		// Testing norm function
		size = 10000;
		A.setSize(size);
		// Testing L1 norm for positive numbers
		for (Index i=0; i<size; i++){
			A.setValue(i,i);
		};
		assert(A.getL1Norm() == 49995000.0 && "TensorO1::getL1Norm() does not compute L1 norm correctly.");
		assert(A.getL2Norm() == 577306.967739001 && "TensorO1::getL2Norm() does not compute L2 norm correctly.");
		assert(A.getLinfNorm() == 9999.0 && "TensorO1::getLinfNorm() does not compute Linf norm correctly.");
		// Testing L1 norm for negative numbers
		for (Index i=0; i<size; i++){
			A.setValue(i,-1*i);
		};
		assert(A.getL1Norm() == 49995000.0 && "TensorO1::getL1Norm() does not compute L1 norm correctly.");
		assert(A.getL2Norm() == 577306.967739001 && "TensorO1::getL2Norm() does not compute L2 norm correctly.");
		assert(A.getLinfNorm() == 9999.0 && "TensorO1::getLinfNorm() does not compute Linf norm correctly.");

		A.computeNorm();
		assert(A.norm.L1 == 49995000.0 and A.norm.L2 == 577306.967739001 and A.norm.Linf ==  9999.0 && "TensorO1::computeNorm() does not return correct value.");



		cout << "Test::TestTensorO1() passed." << endl;
	};






	void TestTensorO2(){
		cout << "\nRunning test on class TensorO2...  " << endl;
		// testing with intigers (positive and negative)

		// initialize a random sequence
		// take a random array size between [1, 100]
		int n = getRandom(10,100);
		int m = getRandom(10,100);
		// create a new TensorO2 of this size
		TensorO2<int> A(n,m);

		unt N = A.getSize0();  //test getSize0() function
		unt M = A.getSize1();  //test getSize1() function
		assert (N == n && "TensorO2::getSize0() does not return correct size");
		assert (M == m && "TensorO2::getSize1() does not return correct size");

		// create a new random variable for setSize() trial
		int newN = getRandom(10,100);
		int newM = getRandom(10,100);
		A.setSize(newN, newM);	//set new size
		// test the new size
		N = A.getSize0();
		M = A.getSize1();
		assert (N == newN and M == newM && "TensorO2::setSize() does not set correct size");

		// loop for testing setValue and getValue functions
		for (Index i=0; i<N; i++){
			for (Index j=0; j<M; j++){
				A.value[i][j] = i+j;
				int value = A.getValue(i,j);
				assert (value==i+j && "TensorO2::getValue() does not return the correct value");

				// test setValue with positive integers
				int value1 = getRandom(1,100000);
				A.setValue(i,j, value1); //set ith variable to a new value i.e. each variable in A.value

				value = A.getValue(i,j);
				assert (value== value1 && "TensorO2::setValue() does not set the correct value");

				// test setValue with negative integers
				int value2 = getRandom(-9999999, -1);
				A.setValue(i,j, value2); //set ith variable to a new value i.e. each variable in A.value

				value = A.getValue(i,j);
				assert (value== value2 && "TensorO2::setValue() does not set the correct value");

			};
		};



		// testing with double (positive and negative)

		// initialize a random sequence
		// take a random array size between [1, 100]
		N = getRandom(10,100);
		M = getRandom(10,100);
		// create a new TensorO2 of this size
		TensorO2<double> C(N,M);

		// loop for testing setValue and getValue functions
		for (Index i=0; i<N; i++){
			for (Index j=0; j<M; j++){
				// test getValue function
				double value1 = getRandom(0.000001, 99999.9999999);
				C.value[i][j] = value1;
				double value = C.getValue(i,j);
				assert (value== value1 && "TensorO2::getValue() does not return the correct value");

				// test setValue with positive double
				double value2 = getRandom(0.0000001, 99999.9999999);
				C.setValue(i,j, value2); //set ith variable to a new value i.e. each variable in C.value

				value = C.getValue(i,j);
				assert (value== value2 && "TensorO2::setValue() does not set the correct value");

				// test setValue with negative double
				double value3 = getRandom(-9999.99999, -0.000001); 
				C.setValue(i,j, value3); //set ith variable to a new value i.e. each variable in C.value

				value = C.getValue(i,j);
				assert (value== value3 && "TensorO2::setValue() does not set the correct value");

				double rand = getRandom(-100.0, 100.0);
				C.addValue(i, j, rand);

				value = C.getValue(i, j);
				assert (value== value3+rand && "TensorO2::addValue() does not add the correct value");
			};
		};


		double Value4 = getRandom(-1000.00, 1000.00);
		C.setAll(Value4);

		for (Index n=0; n<C.getSize0(); n++){
			for (Index m=0; m<C.getSize1(); m++){
				assert(C.getValue(n,m) == Value4);
			};
		};



		cout << "Test::TestTensorO2() passed." << endl;
	};






	void TestMatrix(){
		cout << "\nRunning test on class Matrix...  " << endl;
		// testing with intigers (positive and negative)

		// initialize a random sequence
		// take a random array size between [1, 100]
		int n = getRandom(10,100);
		Matrix<int> A(n);

		unt N = A.getSize0();  //test getSize0() function
		unt M = A.getSize1();  //test getSize1() function
		assert (N == n && "Matrix::getSize0() does not return correct size");
		assert (M == n && "Matrix::getSize1() does not return correct size");

		// create a new random variable for setSize() trial
		int newN = getRandom(10,100);
		A.setSize(newN);	//set new size
		// test the new size
		N = A.getSize0();
		M = A.getSize1();
		assert (N == newN and M == newN && "Matrix::setSize() does not set correct size");

		// loop for testing setValue and getValue functions
		for (Index i=0; i<N; i++){
			for (Index j=0; j<M; j++){
				// testing getValue function
				A.value[i][j] = i+j;
				int value = A.getValue(i,j);
				assert (value==i+j && "Matrix::getValue() does not return the correct value");

				// test setValue with positive integers
				int value1 = getRandom(10,100);
				A.setValue(i,j, value1); //set ith variable to a new value i.e. each variable in A.value

				value = A.getValue(i,j);
				assert (value== value1 && "Matrix::setValue() does not set the correct value");

				// test setValue with negative integers
				int value2 = getRandom(-100, -10);
				A.setValue(i,j, value2); //set ith variable to a new value i.e. each variable in A.value

				value = A.getValue(i,j);
				assert (value== value2 && "Matrix::setValue() does not set the correct value");
			};
		};



		// testing with double (positive and negative)

		// initialize a random sequence
		// take a random array size between [1, 100]
		N = getRandom(10,100);
		Matrix<double> C(N);

		// loop for testing setValue and getValue functions
		for (Index i=0; i<N; i++){
			for (Index j=0; j<N; j++){
				double value1 = getRandom(0.00001, 99999.9999);
				C.value[i][j] = value1;
				double value = C.getValue(i,j);
				assert (value== value1 && "Matrix::getValue() does not return the correct value");

				// test setValue with positive double
				double value2 =  getRandom(0.00001, 99999.9999);
				C.setValue(i,j, value2); //set ith variable to a new value i.e. each variable in C.value

				value = C.getValue(i,j);
				assert (value== value2 && "Matrix::setValue() does not set the correct value");

				// test setValue with negative double
				double value3 =  getRandom(-99999.9999, -0.0000001);
				C.setValue(i,j, value3); //set ith variable to a new value i.e. each variable in C.value

				value = C.getValue(i,j);
				assert (value== value3 && "Matrix::setValue() does not set the correct value");

				double rand = getRandom(-100.0, 100.0);
				C.addValue(i, j, rand);

				value = C.getValue(i, j);
				assert (value== value3+rand && "Matrix::addValue() does not add the correct value");
			};
		};

		double Value4 = getRandom(-1000.00, 1000.00);
		C.setAll(Value4);

		for (Index n=0; n<C.getSize0(); n++){
			for (Index m=0; m<C.getSize1(); m++){
				assert(C.getValue(n,m) == Value4);
			};
		};

  


		// Testing inverse of the matrix

		// initialize a random sequence
		// take a random array size between [1, 100]
		N = getRandom(5,35);  //size between 5 and 15
		Matrix<double> Mat(N);
		
		// Also, open an Eigen matrix for comparing
		Eigen::MatrixXd EMat(N,N);

		// fill both the matrices
		for (Index i=0; i<N; i++){
			for (Index j=0; j<N; j++){
				// generate a random value
				double value = getRandom(-500.0, 500.0);

				// set this value at i,j location for both the matrices
				Mat.setValue(i,j,value);
				EMat(i,j) = value;
			};
		};

		// invert the matrix
		Mat.invert();
		Eigen::MatrixXd EI(N,N);
		EI = EMat.inverse();

		// compare the two values
		for (Index i=0; i<N; i++){
			for (Index j=0; j<N; j++){
				// extract value from i,j location
				double eigVal = EI(i,j); 		//Eigen matrix value
				double meanVal = Mat.inverse[i][j];	//MeanDG matrix value

				assert(eigVal == meanVal && "Matrix::inverse() returns error. Please check the code again.");
			};
		};


		// Test with a diagonal matrix
		// take a random array size between [1, 100]
		N = getRandom(10,1000);		//set a larger matrix between 10 and 110
		Matrix<double> NewMat(N);
		// Also, open an Eigen matrix for comparing
		Eigen::MatrixXd ENewMat(N,N);
		// fill the matrix
		for (Index i=0; i<N; i++){
			for (Index j=0; j<N; j++){
				// generate a random number
				double value = getRandom(-999.999, 999.999);
				if (i == j){
					NewMat.setValue(i,j,value);
					ENewMat(i,j) = value;
				}
				else{
					NewMat.setValue(i,j,0.0);
					ENewMat(i,j) = 0.0;
				};
			};
		};

		// invert the matrix
		NewMat.invert();
		Eigen::MatrixXd EI1(N,N);
		EI1 = ENewMat.inverse();

		// compare the two values
		for (Index i=0; i<N; i++){
			for (Index j=0; j<N; j++){
				// extract the values
				double eigVal = EI1(i,j);
				double meanVal = NewMat.inverse[i][j];
				// compare values between the two matrices
				assert(eigVal == meanVal && "Matrix::inverse() returns error. Please check the code again.");

				// compare with the analytical value
				if (i == j){
					assert(meanVal == 1.0/NewMat.getValue(i,j) && "Matrix::inverse() returns error. Please check the code again.");
				}
				else{
					meanVal = abs(meanVal);
					assert(meanVal < pow(10,-15) && "Matrix::inverse() returns error. Please check the code again.");
				};
			};
		};



		// testing getInverse() method

		int Size = getRandom(5,35);  //size between 5 and 15
		Matrix<double> Mat1(Size);
		
		// Also, open an Eigen matrix for comparing
		Eigen::MatrixXd EMat1(Size,Size);

		// fill both the matrices
		for (Index i=0; i<Size; i++){
			for (Index j=0; j<Size; j++){
				// generate a random value
				double value = getRandom(-500.0, 500.0);

				// set this value at i,j location for both the matrices
				Mat1.setValue(i,j,value);
				EMat1(i,j) = value;
			};
		};

		Matrix<double> inverse(Size);
		Eigen::MatrixXd inverseEig(Size, Size);

		// getInverse method called
		Mat1.getInverse(&inverse);
		inverseEig = EMat1.inverse();
		for (Index i=0; i<Size; i++){
			for (Index j=0; j<Size; j++){
				double error = abs(inverse.getValue(i,j) - inverseEig(i,j));
				assert(error <  PICO && "Mat1h::TestDot() returns error. Kindly check the dot product formulation. Alternately, check the order of the error term, it may be larger than the suggested tolerance (Refer sysInclude.h for EPS, MICRO, NANO etc.)");
			};
		};


		// getInverse method called again, this time only the pre-stored inverse gets copied
		Mat1.getInverse(&inverse);
		for (Index i=0; i<Size; i++){
			for (Index j=0; j<Size; j++){
				double error = abs(inverse.getValue(i,j) - inverseEig(i,j));
				assert(error <  PICO && "Mat1h::TestDot() returns error. Kindly check the dot product formulation. Alternately, check the order of the error term, it may be larger than the suggested tolerance (Refer sysInclude.h for EPS, MICRO, NANO etc.)");
			};
		};




		cout << "Test::TestMatrix() passed." << endl;
	};





	void TestTensorO3(){
		cout << "\nRunning test on class TensorO3...  " << endl;

		// initialize a random sequence
		// take a random array size between [1, 100]
		int r = getRandom(10,100);
		int s = getRandom(10,100);
		int t = getRandom(10,100);
		// create a new TensorO3 of this size
		TensorO3<int> A(r,s,t);

		unt R = A.getSize0();  //test getSize0() function
		unt S = A.getSize1();  //test getSize1() function
		unt T = A.getSize2();  //test getSize2() function
		assert (R == r && "TensorO3::getSize0() does not return correct size");
		assert (S == s && "TensorO3::getSize1() does not return correct size");
		assert (T == t && "TensorO3::getSize2() does not return correct size");

		// create a new random variable for setSize() trial
		int newR = getRandom(10,100);
		int newS = getRandom(10,100);
		int newT = getRandom(10,100);
		A.setSize(newR, newS, newT);	//set new size
		// test the new size
		R = A.getSize0();
		S = A.getSize1();
		T = A.getSize2();
		assert (R == newR and S == newS and T == newT && "TensorO3::setSize() does not set correct size");

		// loop for testing setValue and getValue functions
		for (Index i=0; i<R; i++){
			for (Index j=0; j<S; j++){
				for (Index k=0; k<T; k++){
					// testing getValue function
					A.value[i][j][k] = i+j+k;
					int value = A.getValue(i,j,k);
					assert (value==i+j+k && "TensorO3::getValue() does not return the correct value");

					// testing with positive integers
					int value1 = getRandom(10,100);
					A.setValue(i,j,k, value1); //set ith variable to a new value i.e. each variable in A.value

					value = A.getValue(i,j,k);
					assert (value== value1 && "TensorO3::setValue() does not set the correct value");


					// testing with negative integers
					int value2 = getRandom(-100,-10);
					A.setValue(i,j,k, value2); //set ith variable to a new value i.e. each variable in A.value

					value = A.getValue(i,j,k);
					assert (value== value2 && "TensorO3::setValue() does not set the correct value");

				};
			};
		};


		// Testing for double datatype
		// initialize a random sequence
		// take a random array size between [1, 100]
		r = getRandom(10,100);
		s = getRandom(10,100);
		t = getRandom(10,100);
		// create a new TensorO3 of this size (note the datatype double)
		TensorO3<double> B(r,s,t);

		R = B.getSize0();  //test getSize0() function
		S = B.getSize1();  //test getSize1() function
		T = B.getSize2();  //test getSize2() function
		assert (R == r && "TensorO3::getSize0() does not return correct size");
		assert (S == s && "TensorO3::getSize1() does not return correct size");
		assert (T == t && "TensorO3::getSize2() does not return correct size");

		// create a new random variable for setSize() trial
		newR = getRandom(10,100);
		newS = getRandom(10,100);
		newT = getRandom(10,100);
		B.setSize(newR, newS, newT);	//set new size
		// test the new size
		R = B.getSize0();
		S = B.getSize1();
		T = B.getSize2();
		assert (R == newR and S == newS and T == newT && "TensorO3::setSize() does not set correct size");

		// loop for testing setValue and getValue functions
		for (Index i=0; i<R; i++){
			for (Index j=0; j<S; j++){
				for (Index k=0; k<T; k++){
					// testing getValue() function
					double value0 = getRandom(0.0001, 999.999);
					B.value[i][j][k] = value0;
					double value = B.getValue(i,j,k);
					assert (value== value0 && "TensorO3::getValue() does not return the correct value");

					// testing setValue() with positive double
					double value1 = getRandom(0.0001, 999.999);
					B.setValue(i,j,k, value1); //set ith variable to a new value i.e. each variable in B.value

					value = B.getValue(i,j,k);
					assert (value== value1 && "TensorO3::setValue() does not set the correct value");


					// testing with negative double
					double value2 = getRandom(0.0001, 999.999);
					B.setValue(i,j,k, value2); //set ith variable to a new value i.e. each variable in B.value

					value = B.getValue(i,j,k);
					assert (value== value2 && "TensorO3::setValue() does not set the correct value");

					double rand = getRandom(-100.0, 100.0);
					B.addValue(i, j, k, rand);

					value = B.getValue(i, j, k);
					assert (value== value2+rand && "TensorO3::addValue() does not add the correct value");
				};
			};
		};
		double Value4 = getRandom(-1000.00, 1000.00);
		B.setAll(Value4);

		for (Index n=0; n<B.getSize0(); n++){
			for (Index m=0; m<B.getSize1(); m++){
				for (Index h=0; h<B.getSize2(); h++){
					assert(B.getValue(n,m,h) == Value4);
				};
			};
		};
		cout << "Test::TestTensorO3() passed." << endl;
	};
};
