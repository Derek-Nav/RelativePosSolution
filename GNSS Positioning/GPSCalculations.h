#pragma once

#include "GPSTime.h"
#include "GPSStructures.h"
#include <Array>

#define MU 3.986005e14
#define ROTATION 7.2921151467e-5
#define PI 3.14159265358979323846264338327950288419716939937

using namespace std;
using namespace Eigen;

double ionopcode(double, double);

MatrixXd topocent(MatrixXd x, MatrixXd satpos);
double ionoL1(double, ephemeris, MatrixXd, MatrixXd, double);


MatrixXd getGeocentric(MatrixXd);

double tropo(double, double, double, double);

bool needsAdjustment(MatrixXd, double);

class SatellitePositions
{
public:
	//	Define fields
	vector<vector<ephemeris>> ephem;

	//	Constructor
	SatellitePositions() {}
	SatellitePositions(vector<vector<ephemeris>>);

	MatrixXd getSatellitePosition(int satellite_number, double t, double P1);

	bool hasNav(int satnum);

};
