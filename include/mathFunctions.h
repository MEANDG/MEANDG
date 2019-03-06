/** @file **/
#ifndef MATHFUNCTIONS_H
#define MATHFUNCTIONS_H

#include "sysInclude.h"
#include "Point.h"

using namespace std;

namespace Math{

	/// Find distance between two points
	double getDistance(Point &A, Point &B);
	double getDistance(Point *A, Point *B);

	double calculateAreaOfTriangle(Point &A, Point &B, Point&C);
	double calculateAreaOfTriangle(Point *A, Point *B, Point*C);

	// Algebra operations
	//______________________________________________________________________________________________________________________________


	/// Returns true if the two arguments are approximately equal (absolute value)
	template<typename T>
	bool approx(T arg1, T arg2){
		if (abs(arg1 - arg2) < MICRO){
			return true;
		}
		else{
			return false;
		};
	};


	/// Finds dot product between two objects of the type TensorO1
	template<typename T>
	double dot(TensorO1<T> &A, TensorO1<T> &B){
		// First make sure that both vectors have the same size (and not empty)
		assert(A.size == B.size and A.size > 0 && "Math::dot() returns error. Kindly check sizes of the two argument vectors.");

		double result = 0;
		for (Index i=0; i<A.size; i++){
			result += A.getValue(i)*B.getValue(i);
		};

		return result;
	};

	/// Finds dot product between two objects of the type TensorO1 (pass by pointer)
	template<typename T>
	double dot(TensorO1<T> *A, TensorO1<T> *B){
		// First make sure that both vectors have the same size (and not empty)
		assert(A->size == B->size and A->size > 0 && "Math::dot() returns error. Kindly check sizes of the two argument vectors.");

		double result = 0;
		for (Index i=0; i<A->size; i++){
			result += A->getValue(i)*B->getValue(i);
		};

		return result;
	};

	/// Finds dot product between a TensorO2 object (on the left) and TensorO1 object on the right
	template<typename T>
	void dot(TensorO2<T> &A, TensorO1<T> &B, TensorO1<T> *result){
		// First make sure that both the objects are compatible for a dot product
		assert(A.size[1] == B.size and A.size[0] > 0  and A.size[1] > 0 && "Math::dot() returns error. Kindly check sizes of the two argument objects (TensorO2 and TensorO1).");

		// Make sure that the result vector has the correct size
		int resultSize = result->getSize();
		int requiredSize = A.getSize0();
		assert(resultSize == requiredSize && "Math::dot() returns error. In the function performing dot product between TensorO2 and TensorO1 objects, the result size does not match the required size.");
		for (Index i=0; i<A.getSize0(); i++){
			double result_at_i = 0.0;
			for (Index j=0; j<B.getSize(); j++){
				result_at_i += A.getValue(i,j) * B.getValue(j);
			};
			result->setValue(i, result_at_i);
		};
	};

	/// Finds dot product between a TensorO2 object (on the left) and TensorO1 object on the right (pass by pointer)
	template<typename T>
	void dot(TensorO2<T> *A, TensorO1<T> *B, TensorO1<T> *result){
		// First make sure that both the objects are compatible for a dot product
		assert(A->size[1] == B->size and A->size[0] > 0  and A->size[1] > 0 && "Math::dot() returns error. Kindly check sizes of the two argument objects (TensorO2 and TensorO1).");

		// Make sure that the result vector has the correct size
		int resultSize = result->getSize();
		int requiredSize = A->getSize0();
		assert(resultSize == requiredSize && "Math::dot() returns error. In the function performing dot product between TensorO2 and TensorO1 objects, the result size does not match the required size.");
		for (Index i=0; i<A->getSize0(); i++){
			double result_at_i = 0.0;
			for (Index j=0; j<B->getSize(); j++){
				result_at_i += A->getValue(i,j) * B->getValue(j);
			};
			result->setValue(i, result_at_i);
		};
	};


	/// Finds dot product between a Matrix object (on the left) and TensorO1 object on the right (pass by pointer)
	template<typename T>
	void dot(Matrix<T> &A, TensorO1<T> &B, TensorO1<T> *result){
		// First make sure that both the objects are compatible for a dot product
		assert(A.size[1] == B.size and A.size[0] > 0  and A.size[1] == A.size[0] && "Math::dot() returns error. Kindly check sizes of the two argument objects (Matrix and TensorO1).");

		// Make sure that the result vector has the correct size
		int resultSize = result->getSize();
		int requiredSize = A.getSize0();
		assert(resultSize == requiredSize && "Math::dot() returns error. In the function performing dot product between Matrix and TensorO1 objects, the result size does not match the required size.");

		for (Index i=0; i<A.getSize0(); i++){
			double result_at_i = 0.0;
			for (Index j=0; j<B.getSize(); j++){
				result_at_i += A.getValue(i,j) * B.getValue(j);
			};
			result->setValue(i, result_at_i);
		};
	};

	/// Finds dot product between a Matrix object (on the left) and TensorO1 object on the right (pass by pointer)
	template<typename T>
	void dot(Matrix<T> *A, TensorO1<T> *B, TensorO1<T> *result){
		// First make sure that both the objects are compatible for a dot product
		assert(A->size[1] == B->size and A->size[0] > 0  and A->size[1] == A->size[0] && "Math::dot() returns error. Kindly check sizes of the two argument objects (Matrix and TensorO1).");

		// Make sure that the result vector has the correct size
		int resultSize = result->getSize();
		int requiredSize = A->getSize0();
		assert(resultSize == requiredSize && "Math::dot() returns error. In the function performing dot product between Matrix and TensorO1 objects, the result size does not match the required size.");

		for (Index i=0; i<A->getSize0(); i++){
			double result_at_i = 0.0;
			for (Index j=0; j<B->getSize(); j++){
				result_at_i += A->getValue(i,j) * B->getValue(j);
			};
			result->setValue(i, result_at_i);
		};
	};

	/// Finds dot product between a Matrix object (on the left) and Matrix object on the right
	template<typename T>
	void dot(Matrix<T> &A, Matrix<T> &B, Matrix<T> *result){
		// First make sure that both the objects are compatible for a dot product
		assert(A.size[0] == B.size[1] and A.size[0] > 0  and A.size[1] == A.size[0] && "Math::dot() returns error. Kindly check sizes of the two argument objects (Matrix and Matrix).");

		// Make sure that the result vector has the correct size
		int resultSize = result->getSize0();
		int requiredSize = A.getSize0();
		assert(resultSize == requiredSize && "Math::dot() returns error. In the function performing dot product between Matrix and Matrix objects, the result size does not match the required size.");

		for (Index i=0; i<result->getSize0(); i++){
			for (Index j=0; j<result->getSize1(); j++){
				double result_at_ij = 0.0;
				for (Index k=0; k<A.getSize1(); k++){
					result_at_ij += A.getValue(i,k) * B.getValue(k,j);
				};
				result->setValue(i,j, result_at_ij);
			};
		};
	};

	/// Finds dot product between a Matrix object (on the left) and Matrix object on the right (pass by pointer)
	template<typename T>
	void dot(Matrix<T> *A, Matrix<T> *B, Matrix<T> *result){
		// First make sure that both the objects are compatible for a dot product
		assert(A->size[0] == B->size[1] and A->size[0] > 0  and A->size[1] == A->size[0] && "Math::dot() returns error. Kindly check sizes of the two argument objects (Matrix and Matrix).");

		// Make sure that the result vector has the correct size
		int resultSize = result->getSize0();
		int requiredSize = A->getSize0();
		assert(resultSize == requiredSize && "Math::dot() returns error. In the function performing dot product between Matrix and Matrix objects, the result size does not match the required size.");

		for (Index i=0; i<result->getSize0(); i++){
			for (Index j=0; j<result->getSize1(); j++){
				double result_at_ij = 0.0;
				for (Index k=0; k<A->getSize1(); k++){
					result_at_ij += A->getValue(i,k) * B->getValue(k,j);
				};
				result->setValue(i,j, result_at_ij);
			};
		};
	};

	/// Finds dot product between a TensorO2 object (on the left) and TensorO2 object on the right
	template<typename T>
	void dot(TensorO2<T> &A, TensorO2<T> &B, TensorO2<T> *result){
		// First make sure that both the objects are compatible for a dot product
		assert(A.size[1] == B.size[0] and A.size[0] > 0  and A.size[1] >0 && "Math::dot() returns error. Kindly check sizes of the two argument objects (TensorO2 and TensorO2).");

		// Make sure that the result vector has the correct size
		int resultSize0 = result->getSize0();
		int resultSize1 = result->getSize1();
		int requiredSize0 = A.getSize0();
		int requiredSize1 = B.getSize1();
		assert(resultSize0 == requiredSize0 and resultSize1 == requiredSize1 && "Math::dot() returns error. In the function performing dot product between TensorO2 and TensorO2 objects, the result size does not match the required size.");

		for (Index i=0; i<result->getSize0(); i++){
			for (Index j=0; j<result->getSize1(); j++){
				double result_at_ij = 0.0;
				for (Index k=0; k<A.getSize1(); k++){
					result_at_ij += A.getValue(i,k) * B.getValue(k,j);
				};
				result->setValue(i,j, result_at_ij);
			};
		};
	};

	/// Finds dot product between a TensorO2 object (on the left) and TensorO2 object on the right (pass by pointer)
	template<typename T>
	void dot(TensorO2<T> *A, TensorO2<T> *B, TensorO2<T> *result){
		// First make sure that both the objects are compatible for a dot product
		assert(A->size[1] == B->size[0] and A->size[0] > 0  and A->size[1] >0 && "Math::dot() returns error. Kindly check sizes of the two argument objects (TensorO2 and TensorO2).");

		// Make sure that the result vector has the correct size
		int resultSize0 = result->getSize0();
		int resultSize1 = result->getSize1();
		int requiredSize0 = A->getSize0();
		int requiredSize1 = B->getSize1();
		assert(resultSize0 == requiredSize0 and resultSize1 == requiredSize1 && "Math::dot() returns error. In the function performing dot product between TensorO2 and TensorO2 objects, the result size does not match the required size.");

		for (Index i=0; i<result->getSize0(); i++){
			for (Index j=0; j<result->getSize1(); j++){
				double result_at_ij = 0.0;
				for (Index k=0; k<A->getSize1(); k++){
					result_at_ij += A->getValue(i,k) * B->getValue(k,j);
				};
				result->setValue(i,j, result_at_ij);
			};
		};
	};


	/// Finds cross product between TensorO1 objects
	template<typename T>
	void cross(TensorO1<T> &A, TensorO1<T> &B, TensorO1<T> *result){
		assert(A.size == B.size and A.size == 3 && "Math::cross() returns error in mathFunctions.h. Kindly check sizes of the two argument vectors. Currently supported for 3 only.");

		double First, Second, Third;
		
		First = A.getValue(1) * B.getValue(2) - A.getValue(2) * B.getValue(1);
		
		Second = -(A.getValue(0)* B.getValue(2) - A.getValue(2) * B.getValue(0));

		Third = A.getValue(0) * B.getValue(1) - A.getValue(1) * B.getValue(0);

		result->setValue(0, First);
		result->setValue(1, Second);
		result->setValue(2, Third);
	};

	/// Finds cross product between TensorO1 objects (pass by pointer)
	template<typename T>
	void cross(TensorO1<T> *A, TensorO1<T> *B, TensorO1<T> *result){
		assert(A->size == B->size and A->size == 3 && "Math::cross() returns error in mathFunctions.h. Kindly check sizes of the two argument vectors. Currently supported for 3 only.");

		double First, Second, Third;
		
		First = A->getValue(1) * B->getValue(2) - A->getValue(2) * B->getValue(1);
		
		Second = -(A->getValue(0)* B->getValue(2) - A->getValue(2) * B->getValue(0));

		Third = A->getValue(0) * B->getValue(1) - A->getValue(1) * B->getValue(0);

		result->setValue(0, First);
		result->setValue(1, Second);
		result->setValue(2, Third);
	};

};


namespace Test{
	void TestDistance();
	void TestArea();
	void TestDot();
	void TestDotT1T1();   // dot product of TensorO1 and TensorO1 objects
	void TestDotT2T1();   // dot product of TensorO2 and TensorO1 objects
	void TestDotMT1();   // dot product of Matrix and TensorO1 objects
	void TestDotMM();   // dot product of Matrix and Matrix objects
	void TestDotT2T2();   // dot product of TensorO2 and TensorO2 objects

	void TestCrossT1T1();

};

#endif
