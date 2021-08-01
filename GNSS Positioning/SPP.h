
#include <iostream>
#include <iomanip>
#include "GPSTime.h"
#include "GPSStructures.h"
#include "FileReaders.h"
#include <Eigen/Dense>
#include "GPSCalculations.h"

using namespace std;
using namespace Eigen;
using namespace readers;

#pragma once
class SPP
{
public:
	
	//	Cut-off angle
	double cutoff = 15;

	//	File readers for single point positioning
	ObservationFileReader obs_file;
	NavigationFileReader nav_file;
	SatellitePositions satellite_positions;
	string outPath;

	//	Vector containing navigation data
	vector<vector<ephemeris>> nav_data;
	//	Current epoch under consideration
	epoch current_epoch;

	//	Solution vector (x y z t0)' and its covariance matrix
	MatrixXd X;
	MatrixXd Covar;

	//	Variance factor
	double varFac;

	//	Residual vector
	MatrixXd r;

	//	Constructor definition
	SPP(string, string, string, double);	//	Initialize using strings of the GPS files
	~SPP();

	//	Solve all epochs
	void solveAllEpochs();
	//	Solve a single epoch
	void solveEpoch();

	//	Iterate to the next epoch
	void nextEpoch();

	void PrintVector(const VectorXd& vec, const string& vecName);

	//	Get covariance matrix for a single epoch
	MatrixXd getCovar(MatrixXd);
	//	Get solution for a single epoch
	MatrixXd getdx(MatrixXd, MatrixXd);
	//	Get design matrix for a single epoch
	MatrixXd getA();
	//	Get misclosure vector for a single epoch
	MatrixXd getW();

private:
	SPP() {};
};

