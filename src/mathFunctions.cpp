#include "mathFunctions.h"

using namespace std;

namespace Math{

	double getDistance(Point &A, Point &B){
		double distance;

		TensorO1<double> v(3);

		v.setValue(0, B.getX()-A.getX());
		v.setValue(1, B.getY()-A.getY());
		v.setValue(2, B.getZ()-A.getZ());

		return v.getL2Norm();
	};

	double getDistance(Point *A, Point *B){
		double distance;

		TensorO1<double> v(3);

		v.setValue(0, B->getX()-A->getX());
		v.setValue(1, B->getY()-A->getY());
		v.setValue(2, B->getZ()-A->getZ());

		return v.getL2Norm();
	};


	// Calculate area of a triangle defined by three vertices
	double calculateAreaOfTriangle(Point &A, Point &B, Point&C){
		// Using Heron's formula for computation of the area

		double S;	//semi-perimeter
		double L1, L2, L3; // Lengths of the three sides


		L1 = Math::getDistance(A, B);
		L2 = Math::getDistance(B, C);
		L3 = Math::getDistance(A, C);

		S = (L1 + L2 + L3) / 2.0;

		double Area = sqrt(S * (S-L1) * (S-L2) * (S-L3));

		return Area;
	};


	// Calculate area of a triangle defined by three vertices
	double calculateAreaOfTriangle(Point *A, Point *B, Point*C){
		// Using Heron's formula for computation of the area

		double S;	//semi-perimeter
		double L1, L2, L3; // Lengths of the three sides


		L1 = Math::getDistance(A, B);
		L2 = Math::getDistance(B, C);
		L3 = Math::getDistance(A, C);

		S = (L1 + L2 + L3) / 2.0;

		double Area = sqrt(S * (S-L1) * (S-L2) * (S-L3));

		return Area;
	};

};

// Note:
// All the other functions are described in the header file. Cpp file only contains tests

namespace Test{
	void TestDot(){
		TestDistance();
		TestArea();
		TestDotT1T1();
		TestDotT2T1();
		TestDotMT1();
		TestDotMM();
		TestDotT2T2();

		TestCrossT1T1();
	};

	void TestDistance(){
		cout << "\nRunning test for 'getDistance' operation..." << endl;
		Point A(0);
		Point B(0);

		for (Index i=0; i<10000; i++){
			double randNo;

			randNo = getRandom(-999.99, 999.99);
			A.setX(randNo);
			randNo = getRandom(-999.99, 999.99);
			A.setY(randNo);
			randNo = getRandom(-999.99, 999.99);
			A.setZ(randNo);

			randNo = getRandom(-999.99, 999.99);
			B.setX(randNo);
			randNo = getRandom(-999.99, 999.99);
			B.setY(randNo);
			randNo = getRandom(-999.99, 999.99);
			B.setZ(randNo);

			double distanceA, distanceN;

			distanceA = sqrt(sqr(B.getX() - A.getX()) + sqr(B.getY() - A.getY()) + sqr(B.getZ() - A.getZ()));
			distanceN = Math::getDistance(A, B);

			assert(abs(distanceA - distanceN) < FEMTO && "Test::TestDistance in mathFunctions.cpp returns error. Test the distance formula.\n");

			double distanceN2 = Math::getDistance(A,B);
			assert(abs(distanceA - distanceN2) < FEMTO && "Test::TestDistance in mathFunctions.cpp returns error. Test the distance formula.\n");
		};


		A.setX(-10.0);
		A.setY(0.0);
		A.setZ(0.0);

		B.setX(10.0);
		B.setY(0.0);
		B.setZ(0.0);
		double distanceN = Math::getDistance(A, B);

		assert(distanceN == 20.0 && "Test::TestDistance in mathFunctions.cpp returns error. Test the distance formula.\n");

		A.setX(0.0);
		A.setY(-2.0);
		A.setZ(0.0);

		B.setX(0.0);
		B.setY(-10.0);
		B.setZ(0.0);
		distanceN = Math::getDistance(A, B);

		assert(distanceN == 8.0 && "Test::TestDistance in mathFunctions.cpp returns error. Test the distance formula.\n");

		cout << "Test::TestDistance() passed.\n" << endl;
	};

	void TestArea(){
		cout << "\nRunning test for 'getArea' operation..." << endl;
		Point A(0);
		Point B(1);
		Point C(2);

		for (Index i=0; i<1000; i++){
			double height = getRandom(-10000.00, 10000.00);
			double length = getRandom(-10000.00, 10000.00);

			A.setX(getRandom(-10000.00, 100000.00));
			A.setY(0.0);
			A.setZ(0.0);

			B.setX(A.getX() + length);
			B.setY(0.0);
			B.setZ(0.0);

			C.setX(0.0);
			C.setY(height);
			C.setZ(0.0);

			double areaA = 1.0/2.0 * (abs(length) * abs(height));
			double areaN = Math::calculateAreaOfTriangle(A, B, C);

			double Tolerance = areaA * pow(10,-3);


			assert(abs(areaA - areaN) < Tolerance && "Test::TestArea() returns error in mathFunctions.cpp. Check the area calculations again.\n");
			
			height = getRandom(-10000.00, 10000.00);	
			length = getRandom(-10000.00, 10000.00);

			A.setX(0.0);
			A.setY(getRandom(-10000.00, 10000.00));
			A.setZ(0.0);

			B.setX(0.0);
			B.setY(A.getY() + length);
			B.setZ(0.0);

			C.setX(0.0);
			C.setY(A.getY());
			C.setZ(height);
			
			
			areaA = 1.0/2.0 * (abs(length) * abs(height));
                        areaN = Math::calculateAreaOfTriangle(A, B, C);
                                                                        
                        Tolerance = areaA * pow(10,-3);

			assert(abs(areaA - areaN) < Tolerance && "Test::TestArea() returns error in mathFunctions.cpp. Check the area calculations again.\n");
		};

		cout << "Test::TestArea() passed.\n" << endl;
	};

	void TestDotT1T1(){
		// Tests dot product between two vectors (objects of class TensorO1)

		cout << "\nRunning test for TensorO1.dot.TensorO1 operation..." << endl;

		int Size = getRandom(10,100);
		TensorO1<double> A(Size);
		TensorO1<double> B(Size);

		// Initiate Eigen vectors for comparison of the results
		Eigen::VectorXd AEig(Size);
		Eigen::VectorXd BEig(Size);

		for (Index i=0; i<Size; i++){
			// Initialize all the vectors to random values
			double valueA = getRandom(-99.9999, 99.9999);
			double valueB = getRandom(-99.9999, 99.9999);

			A.setValue(i,valueA);	//MEANDG vector
			AEig(i) = valueA;	//Eigen vector

			B.setValue(i,valueB);   //MEANDG vector
			BEig(i) = valueB;	//Eigen vector

		};

		double resultMEANDG = Math::dot(A,B);
		double resultEIGEN = AEig.transpose()*BEig;

		double error = abs(resultMEANDG - resultEIGEN);
		assert(error <  NANO && "Math::TestDot() returns error. Kindly check the dot product formulation. Alternately, check the order of the error term, it may be larger than the suggested tolerance (Refer sysInclude.h for EPS, MICRO, NANO etc.)");

		// Testing pass by pointer
		resultMEANDG = Math::dot(&A,&B);

		error = abs(resultMEANDG - resultEIGEN);
		assert(error <  NANO && "Math::TestDot() returns error. Kindly check the dot product formulation. Alternately, check the order of the error term, it may be larger than the suggested tolerance (Refer sysInclude.h for EPS, MICRO, NANO etc.)");
		cout << "Test::TestDot() passed for TensorO1.dot.TensorO1 operation.\n" << endl;
	};

	void TestDotT2T1(){
		// Testing dot product between objects of TensorO2 and TensorO1

		cout << "\nRunning test for TensorO2.dot.TensorO1 operation..." << endl;

		int SizeRow = getRandom(10,100);
		int SizeCol = getRandom(10,100);
		TensorO2<double> A(SizeRow, SizeCol);
		TensorO1<double> B(SizeCol);

		// Initiate Eigen vectors for comparison of the results
		Eigen::MatrixXd AEig(SizeRow, SizeCol);
		Eigen::VectorXd BEig(SizeCol);

		for (Index i=0; i<SizeRow; i++){
			for (Index j=0; j<SizeCol; j++){
				// Initialize all the vectors to random values
				double valueA = getRandom(-99.9999, 99.9999);

				A.setValue(i,j,valueA);	//MEANDG vector
				AEig(i,j) = valueA;	//Eigen vector

			};
		};
		for (Index i=0; i<SizeCol; i++){
			double valueB = getRandom(-99.9999, 99.9999);
			B.setValue(i,valueB);   //MEANDG vector
			BEig(i) = valueB;	//Eigen vector
		};

		TensorO1<double> resultMEANDG(SizeRow);
		Eigen::VectorXd resultEIGEN(SizeRow);

		Math::dot(A,B, &resultMEANDG);
		resultEIGEN = AEig*BEig;

		for (Index i=0; i<SizeRow; i++){
			double error = abs(resultMEANDG.getValue(i) - resultEIGEN(i));
			assert(error <  NANO && "Math::TestDot() returns error. Kindly check the dot product formulation. Alternately, check the order of the error term, it may be larger than the suggested tolerance (Refer sysInclude.h for EPS, MICRO, NANO etc.)");
		};

			
		// Pass by pointer

		Math::dot(&A, &B, &resultMEANDG);
		cout << "Test::TestDot() passed for TensorO2.dot.TensorO1 operation.\n" << endl;
		for (Index i=0; i<SizeRow; i++){
			double error = abs(resultMEANDG.getValue(i) - resultEIGEN(i));
			assert(error <  NANO && "Math::TestDot() returns error. Kindly check the dot product formulation. Alternately, check the order of the error term, it may be larger than the suggested tolerance (Refer sysInclude.h for EPS, MICRO, NANO etc.)");
		};
	};

	void TestDotMT1(){
		// Testing dot product between objects of Matrix and TensorO1

		cout << "\nRunning test for Matrix.dot.TensorO1 operation..." << endl;

		int SizeRow = getRandom(10,100);
		int SizeCol = SizeRow;
		Matrix<double> A(SizeRow);
		TensorO1<double> B(SizeCol);

		// Initiate Eigen vectors for comparison of the results
		Eigen::MatrixXd AEig(SizeRow, SizeCol);
		Eigen::VectorXd BEig(SizeCol);

		for (Index i=0; i<SizeRow; i++){
			for (Index j=0; j<SizeCol; j++){
				// Initialize all the vectors to random values
				double valueA = getRandom(-99.9999, 99.9999);

				A.setValue(i,j,valueA);	//MEANDG vector
				AEig(i,j) = valueA;	//Eigen vector

			};
		};
		for (Index i=0; i<SizeCol; i++){
			double valueB = getRandom(-99.9999, 99.9999);
			B.setValue(i,valueB);   //MEANDG vector
			BEig(i) = valueB;	//Eigen vector
		};

		TensorO1<double> resultMEANDG(SizeRow);
		Eigen::VectorXd resultEIGEN(SizeRow);

		Math::dot(A,B, &resultMEANDG);
		resultEIGEN = AEig*BEig;

		for (Index i=0; i<SizeRow; i++){
			double error = abs(resultMEANDG.getValue(i) - resultEIGEN(i));
			assert(error <  NANO && "Math::TestDot() returns error. Kindly check the dot product formulation. Alternately, check the order of the error term, it may be larger than the suggested tolerance (Refer sysInclude.h for EPS, MICRO, NANO etc.)");
		};

		// Pass by pointer
		Math::dot(&A, &B, &resultMEANDG);
		for (Index i=0; i<SizeRow; i++){
			double error = abs(resultMEANDG.getValue(i) - resultEIGEN(i));
			assert(error <  NANO && "Math::TestDot() returns error. Kindly check the dot product formulation. Alternately, check the order of the error term, it may be larger than the suggested tolerance (Refer sysInclude.h for EPS, MICRO, NANO etc.)");
		};

			
		cout << "Test::TestDot() passed for Matrix.dot.TensorO1 operation.\n" << endl;
	};

	void TestDotMM(){
		// Testing dot product between objects of Matrix and Matrix

		cout << "\nRunning test for Matrix.dot.Matrix operation..." << endl;

		int Size = getRandom(10,200);
		Matrix<double> A(Size);
		Matrix<double> B(Size);

		// Initiate Eigen vectors for comparison of the results
		Eigen::MatrixXd AEig(Size, Size);
		Eigen::MatrixXd BEig(Size, Size);

		for (Index i=0; i<Size; i++){
			for (Index j=0; j<Size; j++){
				// Initialize all the vectors to random values
				double valueA = getRandom(-99.9999, 99.9999);

				A.setValue(i,j,valueA);	//MEANDG vector
				AEig(i,j) = valueA;	//Eigen vector

				double valueB = getRandom(-99.9999, 99.9999);
				B.setValue(i,j, valueB);   //MEANDG vector
				BEig(i,j) = valueB;	//Eigen vector
			};
		};

		Matrix<double> resultMEANDG(Size);
		Eigen::MatrixXd resultEIGEN(Size, Size);

		// Dot product performed 
		Math::dot(A,B, &resultMEANDG);		//MEANDG vector
		resultEIGEN = AEig*BEig;		//Eigen vector

		for (Index i=0; i<Size; i++){
			for (Index j=0; j<Size; j++){
				double error = abs(resultMEANDG.getValue(i,j) - resultEIGEN(i,j));
				assert(error <  PICO && "Math::TestDot() returns error. Kindly check the dot product formulation. Alternately, check the order of the error term, it may be larger than the suggested tolerance (Refer sysInclude.h for EPS, MICRO, NANO etc.)");
			};
		};

		
		// Pass by pointer
		Math::dot(&A, &B, &resultMEANDG);		//MEANDG vector
		for (Index i=0; i<Size; i++){
			for (Index j=0; j<Size; j++){
				double error = abs(resultMEANDG.getValue(i,j) - resultEIGEN(i,j));
				assert(error <  PICO && "Math::TestDot() returns error. Kindly check the dot product formulation. Alternately, check the order of the error term, it may be larger than the suggested tolerance (Refer sysInclude.h for EPS, MICRO, NANO etc.)");
			};
		};

			
		cout << "Test::TestDot() passed for Matrix.dot.Matrix operation.\n" << endl;
	};

	void TestDotT2T2(){
		// Testing dot product between objects of TensorO2 and TensorO2

		cout << "\nRunning test for TensorO2.dot.TensorO2 operation..." << endl;

		int SizeR1 = getRandom(10,200);
		int SizeCR = getRandom(10,200);
		int SizeC2 = getRandom(10,200);
		TensorO2<double> A(SizeR1, SizeCR);
		TensorO2<double> B(SizeCR, SizeC2);

		// Initiate Eigen vectors for comparison of the results
		Eigen::MatrixXd AEig(SizeR1, SizeCR);
		Eigen::MatrixXd BEig(SizeCR, SizeC2);

		for (Index i=0; i<SizeR1; i++){
			for (Index j=0; j<SizeCR; j++){
				// Initialize all the vectors to random values
				double valueA = getRandom(-99.9999, 99.9999);

				A.setValue(i,j,valueA);	//MEANDG vector
				AEig(i,j) = valueA;	//Eigen vector
			};
		};


		for (Index i=0; i<SizeCR; i++){
			for (Index j=0; j<SizeC2; j++){
				// Initialize all the vectors to random values
				double valueB = getRandom(-99.9999, 99.9999);
				B.setValue(i,j, valueB);   //MEANDG vector
				BEig(i,j) = valueB;	//Eigen vector
			};
		};

		TensorO2<double> resultMEANDG(SizeR1, SizeC2);
		Eigen::MatrixXd resultEIGEN(SizeR1, SizeC2);

		// Dot product performed 
		Math::dot(A,B, &resultMEANDG);		//MEANDG vector
		resultEIGEN = AEig*BEig;		//Eigen vector

		for (Index i=0; i<SizeR1; i++){
			for (Index j=0; j<SizeC2; j++){
				double error = abs(resultMEANDG.getValue(i,j) - resultEIGEN(i,j));
				assert(error <  PICO && "Math::TestDot() returns error. Kindly check the dot product formulation. Alternately, check the order of the error term, it may be larger than the suggested tolerance (Refer sysInclude.h for EPS, MICRO, NANO etc.)");
			};
		};

		// Pass by pointer
		Math::dot(&A, &B, &resultMEANDG);		//MEANDG vector
		for (Index i=0; i<SizeR1; i++){
			for (Index j=0; j<SizeC2; j++){
				double error = abs(resultMEANDG.getValue(i,j) - resultEIGEN(i,j));
				assert(error <  PICO && "Math::TestDot() returns error. Kindly check the dot product formulation. Alternately, check the order of the error term, it may be larger than the suggested tolerance (Refer sysInclude.h for EPS, MICRO, NANO etc.)");
			};
		};

			
		cout << "Test::TestDot() passed for TensorO2.dot.TensorO2 operation.\n" << endl;
	};


	void TestCrossT1T1(){
		cout << "\nRunning test for TensorO1.cross.TensorO1 operation..." << endl;

		for (Index i=0; i<1000; i++){
			TensorO1<double> A(3);
			TensorO1<double> B(3);
			TensorO1<double> C(3);


			double a1 = getRandom(-1000.0, 1000.0); 	double b1 = getRandom(-1000.0, 1000.0); 
			double a2 = getRandom(-1000.0, 1000.0);         double b2 = getRandom(-1000.0, 1000.0); 
			double a3 = getRandom(-1000.0, 1000.0);         double b3 = getRandom(-1000.0, 1000.0); 


			A.setValue(0, a1);	B.setValue(0, b1);
			A.setValue(1, a2);      B.setValue(1, b2);
			A.setValue(2, a3);      B.setValue(2, b3);


			Math::cross(A, B, &C);

			double c1, c2, c3;

			c1 = C.getValue(0);
			c2 = C.getValue(1);
			c3 = C.getValue(2);


			assert(c1 == a2*b3 - a3*b2 && "Test::TestCrossT1T1() in mathFunctions.cpp returns error. Make sure that the cross product is computed correctly.\n");
			assert(c2 == -(a1*b3 - a3*b1) && "Test::TestCrossT1T1() in mathFunctions.cpp returns error. Make sure that the cross product is computed correctly.\n");
			assert(c3 == a1*b2 - a2*b1 && "Test::TestCrossT1T1() in mathFunctions.cpp returns error. Make sure that the cross product is computed correctly.\n");

			// Test for pass by pointer routine 

			a1 = getRandom(-1000.0, 1000.0); 	 b1 = getRandom(-1000.0, 1000.0); 
			a2 = getRandom(-1000.0, 1000.0);         b2 = getRandom(-1000.0, 1000.0); 
			a3 = getRandom(-1000.0, 1000.0);         b3 = getRandom(-1000.0, 1000.0); 


			A.setValue(0, a1);	B.setValue(0, b1);
			A.setValue(1, a2);      B.setValue(1, b2);
			A.setValue(2, a3);      B.setValue(2, b3);


			Math::cross(&A, &B, &C);


			c1 = C.getValue(0);
			c2 = C.getValue(1);
			c3 = C.getValue(2);


			assert(c1 == a2*b3 - a3*b2 && "Test::TestCrossT1T1() in mathFunctions.cpp returns error. Make sure that the cross product is computed correctly.\n");
			assert(c2 == -(a1*b3 - a3*b1) && "Test::TestCrossT1T1() in mathFunctions.cpp returns error. Make sure that the cross product is computed correctly.\n");
			assert(c3 == a1*b2 - a2*b1 && "Test::TestCrossT1T1() in mathFunctions.cpp returns error. Make sure that the cross product is computed correctly.\n");
		};

		cout << "Test::TestCross() passed for TensorO1.cross.TensorO1 operation.\n" << endl;
	};
};


