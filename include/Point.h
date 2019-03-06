/** @file */
#ifndef POINT_H
#define POINT_H

#include "sysInclude.h"

class Point {

	
public:

	/// Global Variables denoting the three coordinates of a point
	double x;
	double y;
	double z;
	int id;

	/// Constructor with dfault value pointing at the origin
	Point(double x = 0.0 , double y = 0.0, double z = 0.0, int id = 0);

	/// Getter and Setter functions
	void setPoint(int id, double x, double y, double z);

	int getId() const;
	void setId(int Id);

	double getX() const;
	void setX(double x);

	double getY() const;
	void setY(double y);

	double getZ() const;
	void setZ(double z);

	void print();


};

namespace Test{
	void TestPoint();
};
#endif
