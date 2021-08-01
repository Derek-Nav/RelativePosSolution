#pragma once

#include "GPSTime.h"
#include "GPSStructures.h"
#include "FileReaders.h"
#include <Eigen/Dense>
#include "GPSCalculations.h"

using namespace std;
using namespace Eigen;
using namespace readers;

class SPPStat
{
public:

	//	Field definition
	ObservationFileReader obs_file;
	NavigationFileReader nav_file;
	SatellitePositions satellite_positions;
	string outPath;

	//	Vector containing navigation data
	//vector<vector<ephemeris>> nav_data;
	////	Current epoch under consideration
	//epoch current_epoch;
	////	Solution vector (x y z t0)'
	//MatrixXd X;
	//MatrixXd Covar;
	//double varFac;

	////	Constructor definition
	//SPPStat(string, string, string);
	//~SPPStat();

	////	Solve all epochs
	//void solveAllEpochs();
	////	Add a single epoch to the current solution
	//void addEpoch();
	////	Get covariance matrix for a single epoch
	//MatrixXd getCovar(MatrixXd, MatrixXd);
	////	Get solution for a single epoch

};

