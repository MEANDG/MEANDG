/** @file */
#include "Point.h"

using namespace std;


	/// Constructor definition 
	Point::Point(double x, double y, double z, int id){
		this->x = x;
		this->y = y;
		this->z = z;
		this->id = id;
	}


	/// Setter function for point. Takes the id and coordinates as an input.
	void Point::setPoint(int id, double x, double y, double z){
		this->id = id;
		this->x = x;
		this->y = y;
		this->z = z;
	}

	/// Getter function for obtaining the ID of a point
	int Point::getId() const{
		return id;	
	}

	/// Setter function for obtaining the ID of a point
	void Point::setId(int id){
		this->id = id;
	}

	/// Getter function for obtaining the x-coordinate of a point
	double Point::getX() const{
		return x;	
	}

	/// Setter function for obtaining the x-coordinate of a point
	void Point::setX(double x){
		this->x = x;
	}

	/// Getter function for obtaining the y-coordinate of a point
	double Point::getY() const{
		return y;
	}

	/// Setter function for obtaining the y-coordinate of a point
	void Point::setY(double y){
		this->y = y;
	}

	/// Getter function for obtaining the z-coordinate of a point
	double Point::getZ() const{
		return z;
	}

	/// Setter function for obtaining the z-coordinate of a point
	void Point::setZ(double z){
		this->z = z;
	}


	/// Print function. Prints information about the Point object.
	void Point::print(){
		cout << "______________________________" << endl;
		cout << "Point id: " << id << endl;
		cout << "x:" << x << " y:" << y << " z:" << z << endl; 
		cout << "______________________________" << endl;
		cout << endl;
	};




/// Test for class Point. Takes dummy values for initialization and tests all the member functions.
namespace Test{

	void TestPoint(){
		cout << "\nRunning test on class Point... " << endl;

		Point A(0.0, 1.0, 2.0,0);
		int id = A.getId();
		assert(id==0 && "Point::getId() does not return correct id");
		A.setId(1);
		id = A.getId();
		assert(id==1 && "Point::getId() does not return correct id");

		double x = A.getX();
		double y = A.getY();
		double z = A.getZ();
		assert(x == 0.0 && "Point::getX() does not return the correct value.");
		assert(y == 1.0 && "Point::getY() does not return the correct value.");
		assert(z == 2.0 && "Point::getZ() does not return the correct value.");

		A.setX(10.0);
		A.setY(20.0);
		A.setZ(30.0);

		x = A.getX();
                y = A.getY();
	        z = A.getZ();

		assert(x == 10.0 && "Point::setX() does not set the correct value.");
		assert(y == 20.0 && "Point::setY() does not set the correct value.");
		assert(z == 30.0 && "Point::setZ() does not set the correct value.");

		A.setPoint(100, -2.4, -4.3, 1000);
		id = A.getId();
		x = A.getX();
                y = A.getY();
	        z = A.getZ();

		assert(id==100 && "Point::setPoint() does not set the correct correct id");
		assert(x == -2.4 &&   "Point::setPoint() does not set the correct x value.");
		assert(y == -4.3 &&   "Point::setPoint() does not set the correct y value.");
		assert(z == 1000.0 && "Point::setPoint() does not set the correct z value.");

		cout << "Test::TestPoint() passed.\n";
	};

};
