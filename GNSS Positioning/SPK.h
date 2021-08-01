/*
	***Single-Point Kinematic Positioning Module***

	Description: Performs kinematic single-point positioning for GPS observations, using a PVA Kalman filter.
		This method only uses pseudo-ranges, not carrier phase observations.

	TODO: Implementation

	Last Edited: 
		February 22, 2020 - Class definition, definition of required functions, etc.
*/

#pragma once
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

struct AMB_STRUCT
{
	double value;
	double rms;
	bool isFixed;
	bool isValid;
};

class SPK
{
public:

	/*
		DEFINE OBJECT PARAMETERS
	*/

	/*
		Define cutoff angle and scaling parameter
	*/

	//	Define cutoff angle for GPS observations
	double cutoff = 15;

	//	Define the scaling parameter for the process noise in the Kalman filter
	double scaling_parameter_x = 5.0;// 1e-10;
	double scaling_parameter_y = 5.0;//1e-10;
	double scaling_parameter_z = 5.0;//1e-10;
	double scaling_parameter_t = 5.0;//1e-10;

	//	Define variance components for the process noise, system state, and observations in the Kalman filter
	double sigma_wx = 1;
	double sigma_wy = 1;
	double sigma_wz = 1;
	double sigma_wt = 1;

	double sigma_x = 1;

	double sigma_z = 1;

	/*
		Define input/output files, and data structure for organizing GNSS data
	*/

	//	File readers for single point positioning
	ObservationFileReader obs_file, obs_file2;
	NavigationFileReader nav_file, nav_file2;
	SatellitePositions satellite_positions;
	string outPath;
	 
	//	Vector containing navigation data
	vector<vector<ephemeris>> nav_data;

	//	Current epoch under consideration
	epoch current_epoch, current_epoch2;

	/*
		Define Kalman filter variables
	*/

	//	Variance factor
	double varFac;

	//	Measurement residual vector
	//MatrixXd r;

	//	Time of the current epoch
	double tcurr;
	//	Time of the previous epoch
	double tprev;

	//	Solution vector (x x' x'' y y' y'' z z' z'' t0 t0')'
	//  Solution vector(x, x', x'', y, y', y'', z, z', z'')'
	MatrixXd X;
	//	Covariance matrix for the solution vector
	MatrixXd Pk;
	MatrixXd Q;
	//	Predicted value of the solution vector after the time update
	MatrixXd Xp;
	//	Covariance matrix for the predicted solution vector
	MatrixXd PXp;

	//	Double differential pseudornage Observation design matrix
	MatrixXd H_rou;
	MatrixXd H_phi;
	MatrixXd H_rp;// H includes both pseudornage and phase observables 
	MatrixXd H_lsq;

	//	Double differential Pseudorange Observation vector
	MatrixXd Z_rou;
	MatrixXd Z_phi;
	MatrixXd Z_rp;// Z includes both pseudornage and phase observables 
	MatrixXd Z_lsq;

	//	Covariance matrix for the observation vector
	MatrixXd R_rou;
	MatrixXd R_phi;
	MatrixXd R_rp;//R of pseudorange and carrierphase
	MatrixXd R_lsq;

	//	Gain matrix
	MatrixXd G;
	//	Change in solution due to system innovation
	MatrixXd dX;
	//	RESIDUALS
	//	State vector residuals
	MatrixXd vx;
	//	Observation residuals
	MatrixXd vz;
	//	Process noise residuals
	MatrixXd vw;

	/*
		DEFINE OBJECT FUNCTIONS
	*/

	/*
		Constructors/Destructors
	*/

	//	Constructor definition
	SPK(string, string, string, double);	//	Initialize using strings of the GPS files
	SPK(string, string, string, string, string, int, double);
	~SPK();

	/*
		Broad-Scale calculations (epoch-wise solutions, solutions for all epochs, etc.)
	*/

	//	Solve all epochs
	void solveAllEpochs();
	//	Solve a single epoch
	void solveEpoch(bool& initEstimate);

	bool NavigationFloatSolution();

	//	Get the next epoch
	void nextEpoch();

	bool NavsolutionByPseudorange();
	/*
		Calculating specific vector/matrix quantities
	*/

	//	Return the transition matrix
	MatrixXd getTransition();

	//	Return the process noise covariance matrix
	MatrixXd getQ();

	//	Return the gamma matrix, used to derive the process noise covariance matrix
	MatrixXd getGamma();

	//	Predict solution for a single epoch
	void setXp();
	//	Set covariance matrix for the predicted state in a single epoch
	void setPXp();
	//	Get observation design matrix for a single epoch
	void setPseudorangeDesignMatrix();
	//	Get the gain matrix
	void setG();
	//	Calculate the change in estimated state due to the system innovation
	void setdX();

	//	Get solution for a single epoch
	void ToUpdate();
	double GetMaxAbsOfVector(const VectorXd& vec);
	//	Set covariance matrix for the estimated state in a single epoch
	//void setPk();

	
	//	Get the state residuals
	//void setVx();
	//	Get observation residuals
	//void setVz();
	//	Get process noise residuals
	//void setVw();

	//	Get the variance component estimates
	//void setVarianceComponents();

public: 
	// set the double differential pseudorange innovation: (psueodrnage - range)_rover - (psueodrnage - range)_base
	void SetPseudorangeInnovation();
	// set covariance of the double differencital pseudorange innovation  
	void SetPseudorangeCov();
	void SetCarrierphaseCov();// unit in meter^2
	//-----------------------------------------------------------
	void CheckValidObs();
	long long SetSatelltieBitmask();
	MatrixXd GetFloatProcessNoise(const int& sizeOfAmb);
	MatrixXd GetQn(const int& sizeOfAmb);
	MatrixXd GetPn(const int& sizeOfAmb);
	void PrintVector(const VectorXd& vec, const string& vecName);
	void PrintMatrix(const MatrixXd& mat, const string& matName);
	void SetCarrierphaseAmbDesignAndInnovationMatrix(const int& sizeOfAmb);
	MatrixXd SetCarrierphaseDesignInnovMatrix(const int& sizeOfAmb);
	void SetKalmanGain();
	void SaveAmbValues();
	double ComputeFloatApproxAmb(int satID);
public:
	MatrixXd Xn;
	MatrixXd K;
	MatrixXd Px;
	MatrixXd TempH1, TempZ1;
	AMB_STRUCT mSatAmbArray[64];
	int mBasePrn;
	long long mSatBitMask;
	double mSigmaAmbNoise;
	double mWaveLengthL1;
	int mLightSpeed;
	double mSignalFrequencyL1;
	double GetRangeFromCoordinates(const MatrixXd& a, const MatrixXd& b);
	bool mFistSolution;
	double mSigmaPhi;//cycle
	double msigmaRou;//meter
	double mTimeOfFirstsln;
	int mNumEpoches;
private:
	SPK() {};

};

