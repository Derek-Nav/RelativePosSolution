#include "SPK.h"
#include "SPP.h"
#include "GPSTime.h"
#include "GPSStructures.h"
#include "FileReaders.h"
#include <Eigen/Dense>

//#define _IS_POSITION_STATE_ONLY_
#define _IS_PVA_STATE_
//#define _IS_DEBUG_
//#define _IS_MEM_DEBUG_
//#define _IS_BOTH_PSR_CARRIERPHASE_
#define _IS_CARRIER_PHASE_ONLY_
//#define _IS_PSR_ONLY_

using namespace readers;
extern Matrix3d x_base;
using namespace std;
SPK::SPK(string obs_path, string nav_path, string path, double cutoff)
{
	//	Set the cutoff angle for processing
	this->cutoff = cutoff;

	//	Initialize file reader objects
	this->obs_file = ObservationFileReader(obs_path);
	this->nav_file = NavigationFileReader(nav_path);

	//	Ensure that file readers have initialized properly and that files are properly formatted
	if (!obs_file.hasValidFormat())
		throw new exception("Observation file must be RINEX Version 2.10 or 2.11");
	if (!obs_file.hasValidObservationType())
		throw new exception("Observation file must specify GPS observations");
	if (!obs_file.hasValidRINEXType()) 
		throw new exception("Specified observation file does not have compatible file type");
	if (!nav_file.hasValidFormat())
		throw new exception("Navigation file must be RINEX Version 2.10 or 2.11");
	if (!nav_file.hasValidRINEXType())
		throw new exception("Specified navigation file does not have compatible file type");

	//	Read the navigation data
	if (nav_file.hasNextNav())
		this->nav_data = nav_file.getNavigationData();

	//	Get approximate XYZ data
	this->X = VectorXd(3, 1); //MatrixXd(9, 1);
	X.setZero();

	//	Initialize predicted solution vector
	this->Xp = VectorXd(3, 1);//MatrixXd(9, 1);
	Xp.setZero();

	//	Initialize covariance matrices to zero
	this->Pk = MatrixXd(3, 3);//MatrixXd(9, 9);
	Pk.setZero();

	this->PXp = MatrixXd(3, 3);//MatrixXd(9, 9);
	PXp.setZero();

	//	Set the assumed variance factor of 1
	this->varFac = 1;

	//	Initialize the satellite position information
	satellite_positions = SatellitePositions(nav_data);

	//	Set the output path for processing
	outPath = path;

	//	Process the first epoch using the standard least-squares solution
	this->current_epoch = this->obs_file.getNextEpoch();

	this->tcurr = this->current_epoch.time;

	SPP SPP_solution(obs_path, nav_path, "testfile", cutoff);
	SPP_solution.nextEpoch();

	epoch e = SPP_solution.current_epoch;
	SPP_solution.solveEpoch();

	MatrixXd SPP_X = SPP_solution.X;
	this->X(0, 0) = SPP_X(0, 0);
	this->X(1, 0) = SPP_X(1, 0);
	this->X(2, 0) = SPP_X(2, 0);
	/*this->X(0, 0) = SPP_X(0, 0);
	this->X(3, 0) = SPP_X(1, 0);
	this->X(6, 0) = SPP_X(2, 0);*/
#ifdef _IS_DEBUG_
	PrintVector(this->X, "X");
#endif
	MatrixXd SPP_PX = SPP_solution.Covar;
	this->Pk(0, 0) = SPP_PX(0, 0);
	this->Pk(1, 0) = SPP_PX(1, 0);
	this->Pk(2, 0) = SPP_PX(2, 0);
	this->Pk(0, 1) = SPP_PX(0, 1);
	this->Pk(1, 1) = SPP_PX(1, 1);
	this->Pk(2, 1) = SPP_PX(2, 1);
	this->Pk(0, 2) = SPP_PX(0, 2);
	this->Pk(1, 2) = SPP_PX(1, 2);
	this->Pk(2, 2) = SPP_PX(2, 2);
	/*this->Pk(0, 0) = SPP_PX(0, 0);
	this->Pk(3, 0) = SPP_PX(1, 0);
	this->Pk(6, 0) = SPP_PX(2, 0);
	this->Pk(0, 3) = SPP_PX(0, 1);
	this->Pk(3, 3) = SPP_PX(1, 1);
	this->Pk(6, 3) = SPP_PX(2, 1);
	this->Pk(0, 6) = SPP_PX(0, 2);
	this->Pk(3, 6) = SPP_PX(1, 2);
	this->Pk(6, 6) = SPP_PX(2, 2);*/
}

#ifdef _IS_POSITION_STATE_ONLY_
SPK::SPK(string obs_path1, string nav_path1, string obs_path2, 
	     string nav_path2, string path, int basePrn, double cutoff)
	:mBasePrn(basePrn),
	mSatBitMask(0),
	mSigmaAmbNoise(2.0),
	mWaveLengthL1(0.19),
	mLightSpeed(299792458),
	mSigmaPhi(1),// unit of cycle
	mFistSolution(true)
{
	//	Set the cutoff angle for processing
	this->cutoff = cutoff;

	//	Initialize file reader objects
	this->obs_file  = ObservationFileReader(obs_path1);
	this->nav_file  = NavigationFileReader(nav_path1);
	this->obs_file2 = ObservationFileReader(obs_path2);
	this->nav_file2 = NavigationFileReader(nav_path2);

	if (nav_file.hasNextNav())
		this->nav_data = nav_file.getNavigationData();
	//	Get approximate XYZ data
	this->X = MatrixXd(3, 1); //MatrixXd(9, 1);//this->X = MatrixXd(11, 1);
	X.setZero();

	//	Initialize predicted solution vector
	this->Xp = MatrixXd(3, 1);// MatrixXd(9, 1);//this->Xp = MatrixXd(11, 1);
	Xp.setZero();

	//	Initialize covariance matrices to zero
	this->Pk = MatrixXd(3, 3);// MatrixXd(9, 9);//this->Pk = MatrixXd(11, 11);
	Pk.setZero();

	this->PXp = MatrixXd(3, 3);//MatrixXd(9, 9);//this->PXp = MatrixXd(11, 11);
	PXp.setZero();

	//	Set the assumed variance factor of 1
	this->varFac = 1;
	//	Set the assumed variance factor of 1
	this->varFac = 1;

	//	Initialize the satellite position information
	satellite_positions = SatellitePositions(nav_data);

	//	Set the output path for processing
	outPath = path;

	//	Process the first epoch using the standard least-squares solution
	this->current_epoch  = this->obs_file.getNextEpoch();
	this->current_epoch2 = this->obs_file2.getNextEpoch();

	this->tcurr = this->current_epoch.time;

	SPP SPP_solution(obs_path1, nav_path1, "testfile", cutoff);
	SPP_solution.nextEpoch();

	epoch e = SPP_solution.current_epoch;
	SPP_solution.solveEpoch();
	PrintVector(this->X, "X");

	MatrixXd SPP_X = SPP_solution.X;
	this->X(0, 0) = SPP_X(0, 0);
	this->X(1, 0) = SPP_X(1, 0);
	this->X(2, 0) = SPP_X(2, 0);
	PrintVector(this->X, "X");

	MatrixXd SPP_PX = SPP_solution.Covar;
	this->Pk(0, 0) = SPP_PX(0, 0);
	this->Pk(1, 0) = SPP_PX(1, 0);
	this->Pk(2, 0) = SPP_PX(2, 0);
	this->Pk(0, 1) = SPP_PX(0, 1);
	this->Pk(1, 1) = SPP_PX(1, 1);
	this->Pk(2, 1) = SPP_PX(2, 1);
	this->Pk(0, 2) = SPP_PX(0, 2);
	this->Pk(1, 2) = SPP_PX(1, 2);
	this->Pk(2, 2) = SPP_PX(2, 2);
	//this->Pk(9, 9) = SPP_PX(3, 3);
}
#endif

#ifdef _IS_PVA_STATE_
SPK::SPK(string obs_path1, string nav_path1, string obs_path2,
	string nav_path2, string path, int basePrn, double cutoff)
	:mBasePrn(basePrn),
	mSatBitMask(0),
	mSigmaAmbNoise(5.0),
	mWaveLengthL1(0.190293),
	mLightSpeed(299792458),
	mSigmaPhi(1),// unit of cycle
	msigmaRou(5),
	mNumEpoches(0),
	mSignalFrequencyL1(1575.42e6),
	mFistSolution(true)
{
	mWaveLengthL1 = mLightSpeed / mSignalFrequencyL1;
	//	Set the cutoff angle for processing
	this->cutoff = cutoff;

	//	Initialize file reader objects
	this->obs_file = ObservationFileReader(obs_path1);
	this->nav_file = NavigationFileReader(nav_path1);
	this->obs_file2 = ObservationFileReader(obs_path2);
	this->nav_file2 = NavigationFileReader(nav_path2);

	if (nav_file.hasNextNav())
		this->nav_data = nav_file.getNavigationData();
	//	Get approximate XYZ data
	this->X = VectorXd(9, 1);
	X.setZero();

	//	Initialize predicted solution vector
	this->Xp = VectorXd(9, 1);
	Xp.setZero();

	//	Initialize covariance matrices to zero
	this->Pk = MatrixXd(9, 9);
	Pk.setIdentity();

	this->PXp = MatrixXd(9, 9);
	PXp.setIdentity();

	//	Set the assumed variance factor of 1
	this->varFac = 1;
	//	Set the assumed variance factor of 1
	this->varFac = 1;

	//	Initialize the satellite position information
	satellite_positions = SatellitePositions(nav_data);

	//	Set the output path for processing
	outPath = path;

	//	Process the first epoch using the standard least-squares solution
	this->current_epoch = this->obs_file.getNextEpoch();
	this->current_epoch2 = this->obs_file2.getNextEpoch();

	this->tcurr = this->current_epoch.time;

	SPP SPP_solution(obs_path1, nav_path1, "testfile", cutoff);
	SPP_solution.nextEpoch();

	epoch e = SPP_solution.current_epoch;
	SPP_solution.solveEpoch();
	
	MatrixXd SPP_X = SPP_solution.X;
	/*this->X(0, 0) = SPP_X(0, 0);
	this->X(1, 0) = SPP_X(1, 0);
	this->X(2, 0) = SPP_X(2, 0);*/
	this->X(0, 0) = SPP_X(0, 0);
	this->X(3, 0) = SPP_X(1, 0);
	this->X(6, 0) = SPP_X(2, 0);
	//X(0, 0) = -2440668;
	//X(3, 0) = -3416438;
	//X(6, 0) = 4785136;// testing.
#ifdef _IS_DEBUG_
	PrintVector(this->X, "initial X");
#endif
	MatrixXd SPP_PX = SPP_solution.Covar;
	/*this->Pk(0, 0) = SPP_PX(0, 0);
	this->Pk(1, 0) = SPP_PX(1, 0);
	this->Pk(2, 0) = SPP_PX(2, 0);
	this->Pk(0, 1) = SPP_PX(0, 1);
	this->Pk(1, 1) = SPP_PX(1, 1);
	this->Pk(2, 1) = SPP_PX(2, 1);
	this->Pk(0, 2) = SPP_PX(0, 2);
	this->Pk(1, 2) = SPP_PX(1, 2);
	this->Pk(2, 2) = SPP_PX(2, 2);*/
	this->Pk(0, 0) = SPP_PX(0, 0);
	this->Pk(3, 0) = SPP_PX(1, 0);
	this->Pk(6, 0) = SPP_PX(2, 0);
	this->Pk(0, 3) = SPP_PX(0, 1);
	this->Pk(3, 3) = SPP_PX(1, 1);
	this->Pk(6, 3) = SPP_PX(2, 1);
	this->Pk(0, 6) = SPP_PX(0, 2);
	this->Pk(3, 6) = SPP_PX(1, 2);
	this->Pk(6, 6) = SPP_PX(2, 2);

	//***
	memset(mSatAmbArray, 0, sizeof(mSatAmbArray));
}
#endif

SPK::~SPK()
{
}

#ifdef _IS_POSITION_STATE_ONLY_
MatrixXd SPK::getTransition()
{
	//	Get time difference between current and previous epoch
	double dt = 24 * 3600 * (this->tcurr - this->tprev);

	/*	
		Populate the transition matrix
	*/
	MatrixXd transitMatrix = MatrixXd(3, 3);//MatrixXd phi = MatrixXd(11, 11);
	transitMatrix.setIdentity();
	return transitMatrix;
}
#endif

#ifdef _IS_PVA_STATE_
MatrixXd SPK::getTransition()
{
	//	Get time difference between current and previous epoch
	double dt = 24 * 3600 * (this->tcurr - this->tprev);

	/*
		Populate the transition matrix
	*/
	MatrixXd phi = MatrixXd(9, 9);
	phi.setIdentity();

	//	Set x component values
	phi(0, 0) = 1;
	phi(0, 1) = dt;
	phi(0, 2) = 0.5 * dt * dt;
	phi(1, 1) = 1;
	phi(1, 2) = dt;
	phi(2, 2) = 1;

	//	Set y component values
	phi(3, 3) = 1;
	phi(3, 4) = dt;
	phi(3, 5) = 0.5 * dt * dt;
	phi(4, 4) = 1;
	phi(4, 5) = dt;
	phi(5, 5) = 1;

	//	Set z component values
	phi(6, 6) = 1;
	phi(6, 7) = dt;
	phi(6, 8) = 0.5 * dt * dt;
	phi(7, 7) = 1;
	phi(7, 8) = dt;
	phi(8, 8) = 1;

	//	Set receiver clock error values
	//phi(9, 9) = 1;
	//phi(9, 10) = dt;
	//phi(10, 10) = 1;

	return phi;
}
#endif

#ifdef _IS_POSITION_STATE_ONLY_
MatrixXd SPK::getQ()
{
	//	Get time difference between current and previous epoch
	double dt = 24 * 3600 * (this->tcurr - this->tprev);

	/*
		Calculate the process noise covariance matrix
	*/
	MatrixXd Q = MatrixXd(3, 3);//MatrixXd Q = MatrixXd(11, 11);
	Q.setZero();

	//	Block for the x component
	Q(0, 0) = pow(dt, 6) / 9. * scaling_parameter_x;
	Q(0, 1) = pow(dt, 5) / 6. * scaling_parameter_x;
	Q(0, 2) = pow(dt, 4) / 3. * scaling_parameter_x;
	Q(1, 0) = pow(dt, 6) / 9. * scaling_parameter_y;
	Q(1, 1) = pow(dt, 5) / 6. * scaling_parameter_y;
	Q(1, 2) = pow(dt, 4) / 3. * scaling_parameter_y;
	Q(2, 0) = pow(dt, 6) / 9. * scaling_parameter_z;
	Q(2, 1) = pow(dt, 5) / 6. * scaling_parameter_z;
	Q(2, 2) = pow(dt, 4) / 3. * scaling_parameter_z;

	return Q;
}
#endif

#ifdef _IS_PVA_STATE_
MatrixXd SPK::getQ()
{
	//	Get time difference between current and previous epoch
	double dt = 24 * 3600 * (this->tcurr - this->tprev);

	/*
		Calculate the process noise covariance matrix
	*/
	MatrixXd Q = MatrixXd(this->X.rows(), this->X.rows());//MatrixXd Q = MatrixXd(11, 11);
	Q.setOnes();
	Q = Q * scaling_parameter_x;

	////	Block for the x component
	//Q(0, 0) = pow(dt, 6) / 9. * scaling_parameter_x;
	//Q(0, 1) = pow(dt, 5) / 6. * scaling_parameter_x;
	//Q(0, 2) = pow(dt, 4) / 3. * scaling_parameter_x;
	//Q(1, 0) = pow(dt, 5) / 6. * scaling_parameter_x;
	//Q(1, 1) = pow(dt, 4) / 4. * scaling_parameter_x;
	//Q(1, 2) = pow(dt, 3) / 2. * scaling_parameter_x;
	//Q(2, 0) = pow(dt, 4) / 3. * scaling_parameter_x;
	//Q(2, 1) = pow(dt, 3) / 2. * scaling_parameter_x;
	//Q(2, 2) = pow(dt, 2) * scaling_parameter_x;

	////	Block for the y component
	//Q(3, 3) = pow(dt, 6) / 9. * scaling_parameter_y;
	//Q(3, 4) = pow(dt, 5) / 6. * scaling_parameter_y;
	//Q(3, 5) = pow(dt, 4) / 3. * scaling_parameter_y;
	//Q(4, 3) = pow(dt, 5) / 6. * scaling_parameter_y;
	//Q(4, 4) = pow(dt, 4) / 4. * scaling_parameter_y;
	//Q(4, 5) = pow(dt, 3) / 2. * scaling_parameter_y;
	//Q(5, 3) = pow(dt, 4) / 3. * scaling_parameter_y;
	//Q(5, 4) = pow(dt, 3) / 2. * scaling_parameter_y;
	//Q(5, 5) = pow(dt, 2) * scaling_parameter_y;

	////	Block for the z component
	//Q(6, 6) = pow(dt, 6) / 9. * scaling_parameter_z;
	//Q(6, 7) = pow(dt, 5) / 6. * scaling_parameter_z;
	//Q(6, 8) = pow(dt, 4) / 3. * scaling_parameter_z;
	//Q(7, 6) = pow(dt, 5) / 6. * scaling_parameter_z;
	//Q(7, 7) = pow(dt, 4) / 4. * scaling_parameter_z;
	//Q(7, 8) = pow(dt, 3) / 2. * scaling_parameter_z;
	//Q(8, 6) = pow(dt, 4) / 3. * scaling_parameter_z;
	//Q(8, 7) = pow(dt, 3) / 2. * scaling_parameter_z;
	//Q(8, 8) = pow(dt, 2) * scaling_parameter_z;

	//	Block for the receiver clock error
	/*Q(9, 9) = pow(dt, 4) / 4. * scaling_parameter_t;
	Q(9, 10) = pow(dt, 3) / 2. * scaling_parameter_t;
	Q(10, 9) = pow(dt, 3) / 2. * scaling_parameter_t;
	Q(10, 10) = pow(dt, 2) * scaling_parameter_t;*/

	return Q;
}
#endif

MatrixXd SPK::getGamma()
{
	//	Get time difference between current and previous epoch
	double dt = 24 * 3600 * (this->tcurr - this->tprev);

	/*
		Calculate the gamma matrix
	*/
	double dt3 = pow(dt, 3);
	double dt2 = pow(dt, 2);

	MatrixXd Gamma = MatrixXd(9, 3);//MatrixXd(11, 4);
	Gamma.setZero();

	Gamma(0, 0) = dt3 / 3.;
	Gamma(1, 0) = dt2 / 2.;
	Gamma(2, 0) = dt;
	Gamma(3, 1) = dt3 / 3.;
	Gamma(4, 1) = dt2 / 2.;
	Gamma(5, 1) = dt;
	Gamma(6, 2) = dt3 / 3.;
	Gamma(7, 2) = dt2 / 2.;
	Gamma(8, 2) = dt;
	//Gamma(9, 3) = dt2 / 2.;
	//Gamma(10, 3) = dt;
	//Gamma(10, 0) = pow(dt, 1);

	return Gamma;
}

void SPK::setXp()
{
	//	Get the transition matrix
	MatrixXd phi = this->getTransition();

	//	Propagate the previous state to the current epoch
	this->Xp = phi * this->X;
}

void SPK::setPXp()
{
	//	Get time difference between current and previous epoch
	double dt = 24 * 3600 * (this->tcurr - this->tprev);

	//	Get the process noise covariance matrix
	MatrixXd Q = this->getQ();
	//DEBUG:
	cout << "Q size: " << Q.rows() << endl;

	/*
		Propagate the error from the previous state
	*/
	//	Get the transition matrix
	MatrixXd phi = this->getTransition();

	//	Propagate error from the previous state through the transition function
	MatrixXd P;// = phi * this->Pk * phi.transpose();

	//DEBUGG:
	cout << "P size: " << P.rows() << endl;

	//	Add the process noise to the propagated covariance matrix
	this->PXp = P + Q;
}


#ifdef _IS_POSITION_STATE_ONLY_
void SPK::setPseudorangeDesignMatrix()
{
	epoch e1  = current_epoch;
	epoch e2 = current_epoch2;

	// Determine the number of L1 observations that are in the data
	int count = 0;
	for (int i = 0; i < signed(e1.epoch_observations.size()); i++)
	{
		observation curr_obs = e1.epoch_observations[i];

		if (curr_obs.C1 && satellite_positions.hasNav(curr_obs.PRN))
			count++;
	}

	//	Initialize the size of the first design matrix
	MatrixXd H1(count, this->X.rows());// A(count, 11);
	H1.setZero();
	double t = e1.time;
	vector<observation> obs1 = e1.epoch_observations;
	vector<observation> obs2 = e2.epoch_observations;

	//	For each observation in e, determine the satellite positions
	int row = 0;
	for (int i = 0; i < signed(obs1.size()); i++)
	{
		observation curr_obs1 = obs1[i];
		observation curr_obs2;
		curr_obs2.C1 = -1;
		for (int j = 0; j < obs2.size(); j++) {
			if (obs2[j].PRN == curr_obs1.PRN) {
				curr_obs2 = obs2[j];
				break;
			}
		}
		if (curr_obs2.C1 == -1) continue;

		if (curr_obs1.C1 == 0 || !satellite_positions.hasNav(curr_obs1.PRN))
			continue;


		//	Get the satellite position and dts value
		MatrixXd sat_info = satellite_positions.getSatellitePosition(curr_obs1.PRN, t, curr_obs1.P1);

		//	Initialize the first derivatives with respect to X, Y, Z, and dt
		MatrixXd dx = this->Xp.block(0, 0, 3, 1) - sat_info.block(0, 0, 3, 1);
		double range = sqrt(pow(dx(0, 0), 2) + pow(dx(1, 0), 2) + pow(dx(2, 0), 2));
		double partialX = (dx(0, 0)) / range;
		double partialY = (dx(1, 0)) / range;
		double partialZ = (dx(2, 0)) / range;

		H1(row, 0) = partialX;
		H1(row, 1) = partialY;
		H1(row, 2) = partialZ;
		row++;
	}
	MatrixXd H2(row - 1, this->X.rows());
	H2.setZero();
	// differential H
	for (int i = 1; i < row; i++) {
		H2(i-1, 0) = H1(i, 0) - H1(0, 0);
		H2(i-1, 1) = H1(i, 1) - H1(0, 1);
		H2(i-1, 2) = H1(i, 2) - H1(0, 2);
	}

	this->H_rou = H2;
	// DEBUG:
	std::cout << "------------ H ---------------\n";
	for (int i = 0; i < H2.rows(); i++) {
		for (int j = 0; j < this->X.rows(); j++) {
			std::cout << H2(i, j) << ", ";
		}
		std::cout << endl;
	}
	std::cout << "-------------------------------\n";
}
#endif

#ifdef _IS_PVA_STATE_
void SPK::setPseudorangeDesignMatrix()
{
	epoch e1 = current_epoch;
	epoch e2 = current_epoch2;

	// Determine the number of L1 observations that are in the data
	int count = 0;
	for (int i = 0; i < signed(e1.epoch_observations.size()); i++)
	{
		observation curr_obs = e1.epoch_observations[i];

		if (curr_obs.C1 && satellite_positions.hasNav(curr_obs.PRN))
			count++;
	}

	//	Initialize the size of the first design matrix
	MatrixXd H1(count, this->X.rows());// A(count, 11);
	H1.setZero();

	double t = e1.time;
	vector<observation> obs1 = e1.epoch_observations;
	vector<observation> obs2 = e2.epoch_observations;

	//	For each observation in e, determine the satellite positions
	int row = 0;
	for (int i = 0; i < signed(obs1.size()); i++)
	{
		observation curr_obs1 = obs1[i];
		observation curr_obs2;
		curr_obs2.C1 = -1;
		for (int j = 0; j < obs2.size(); j++) {
			if (obs2[j].PRN == curr_obs1.PRN) {
				curr_obs2 = obs2[j];
				break;
			}
		}
		if (curr_obs2.C1 == -1) continue;

		if (curr_obs1.C1 == 0 || !satellite_positions.hasNav(curr_obs1.PRN))
			continue;


		//	Get the satellite position and dts value
		MatrixXd sat_info = satellite_positions.getSatellitePosition(curr_obs1.PRN, t, curr_obs1.P1);

		//	Initialize the first derivatives with respect to X, Y, Z, and dt
		MatrixXd dx(3, 1);
		dx(0, 0) = this->X(0, 0) - sat_info(0, 0);
		dx(1, 0) = this->X(3, 0) - sat_info(1, 0);
		dx(2, 0) = this->X(6, 0) - sat_info(2, 0);

		double range = sqrt(pow(dx(0, 0), 2) + pow(dx(1, 0), 2) + pow(dx(2, 0), 2));
		double partialX = (dx(0, 0)) / range;
		double partialY = (dx(1, 0)) / range;
		double partialZ = (dx(2, 0)) / range;
		//double partialT = 1.;

		H1(row, 0) = partialX;
		H1(row, 3) = partialY;
		H1(row, 6) = partialZ;
		//A(row, 9) = partialT;
		row++;
	}
	// differential H
	MatrixXd H2(row - 1, this->X.rows());
	H2.setZero();
	for (int i = 1; i < row; i++) {
		H2(i - 1, 0) = H1(i, 0) - H1(0, 0);
		H2(i - 1, 3) = H1(i, 3) - H1(0, 3);
		H2(i - 1, 6) = H1(i, 6) - H1(0, 6);
	}

	this->H_rou = H2;
	// DEBUG:
	std::cout << "------------ H ---------------\n";
	for (int i = 0; i < H2.rows(); i++) {
		for (int j = 0; j < this->X.rows(); j++) {
			std::cout << H2(i, j) << ", ";
		}
		std::cout << endl;
	}
	std::cout << "-------------------------------\n";
}

#endif

//	Returns the system innovation vector
#ifdef _IS_POSITION_STATE_ONLY_
void SPK::SetPseudorangeInnovation()
{
	epoch e1 = current_epoch;
	epoch e2 = current_epoch2;// base

	//	Determine the number of L1 observations that are in the data
	int count = 0;
	for (int i = 0; i < signed(e1.epoch_observations.size()); i++)
	{
		observation curr_obs = e1.epoch_observations[i];
		if ((curr_obs.C1 != 0) && satellite_positions.hasNav(curr_obs.PRN))
			count++;
	}

	//	Initialize the size of the observation vector
	MatrixXd l(count, 1);
	l.setZero();
	double t = e1.time;
	vector<observation> obs1 = e1.epoch_observations;
	vector<observation> obs2 = e2.epoch_observations;// base

	//	Populate the observation vector
	int row = 0;
	int numOfObs = signed(obs1.size());

	for (int i = 0; i < numOfObs; i++)
	{
		observation curr_obs1 = obs1[i];
		observation curr_obs2;
		curr_obs2.C1 = -1;
		for (int j = 0; j < obs2.size(); j++) {
			if (obs2[j].PRN == curr_obs1.PRN) {
				curr_obs2 = obs2[j];
				break;
			}
		}
		if (curr_obs2.C1 == -1) continue;

		if (curr_obs1.C1 == 0 || !satellite_positions.hasNav(curr_obs1.PRN))
			continue;

		//	Get L1 values
		double pseudorange1, pseudorange2;
		pseudorange1 = curr_obs1.C1;
		pseudorange2 = curr_obs2.C1;

		//	Get the satellite position and dts value
		MatrixXd sat_info = satellite_positions.getSatellitePosition(curr_obs1.PRN, t, curr_obs1.C1);

		//	Correct for ionospheric delay
		double iono = ionoL1(pseudorange1, this->nav_data[curr_obs1.PRN-1][0], this->X, sat_info, t);

		double dts = sat_info(3, 0);
		MatrixXd svX = sat_info.block(0, 0, 3, 1);

		//	Correct the pseudorange for the dts value
		iono = iono + dts * 299792458;

		//	Get the elevation angle of the satellite
		MatrixXd topocentric = topocent(this->X, svX);
		double el = topocentric(1, 0);

		//	Get the tropospheric delay
		double trop = 0; // tropo(el, 1013., 0, 293);

		double range1 = sqrt(pow(svX(0, 0) - Xp(0, 0), 2) + pow(svX(1, 0) - Xp(1, 0), 2) + pow(svX(2, 0) - Xp(2, 0), 2));
		double range2 = sqrt(pow(svX(0, 0) - x_base(0, 0), 2) + pow(svX(1, 0) - x_base(1, 0), 2) + pow(svX(2, 0) - x_base(2, 0), 2));
		//	Correct the pseudorange for tropospheric delay
		l(row, 0) = (pseudorange1 - range1) - (pseudorange2 - range2);
		//int tmp = l(row,0);
		row++;

	}
	//DEBUG:
	std::cout << "\n-------------- l ----------------\n";
	for (int i = 0; i < row; i++) {
		cout << l(i, 0) << ", ";
	}
	std::cout << "\n-------------------------------\n";

	MatrixXd ll(row-1, 1);
	ll.setZero();
	for (int i = 1; i < row; i++) {
		ll(i - 1, 0) = l(i, 0) - l(0, 0);
	}

	this->Z_rou = ll;

	// DEBUG:
	std::cout << "-------------- Z ----------------\n";
	for (int i = 0; i < row-1; i++) {
		std::cout << ll(i, 0) << ", ";
	}
	std::cout << "\n-------------------------------\n";
}
#endif

#ifdef _IS_PVA_STATE_
long long SPK::SetSatelltieBitmask()
{
	long long satBitMask = 0;
	long long UL1 = 1;
	epoch e1 = current_epoch;
	epoch e2 = current_epoch2;// base

	//	Initialize the size of the observation vector
	vector<observation> obs1 = e1.epoch_observations;
	vector<observation> obs2 = e2.epoch_observations;// base

	//	Populate the observation vector
	int row = 0;
	int numOfObs = signed(obs1.size());
	for (int i = 0; i < numOfObs; i++)
	{
		observation curr_obs1 = obs1[i];
		observation curr_obs2;
		curr_obs2.L1 = -1;
		for (int j = 0; j < obs2.size(); j++) {
			if (obs2[j].PRN == curr_obs1.PRN) {
				curr_obs2 = obs2[j];
				break;
			}
		}
		if (curr_obs2.L1 == -1) continue;

		if (curr_obs1.L1 == 0 || curr_obs1.C1 == 0 ||
			curr_obs2.L1 == 0 || curr_obs2.C1 == 0 || 
			!satellite_positions.hasNav(curr_obs1.PRN))
			continue;

		satBitMask |= (UL1 << (curr_obs1.PRN - 1));
	}
	return satBitMask;
}
#endif
#ifdef _IS_PVA_STATE_
void SPK::SetPseudorangeInnovation()
{
	epoch e1 = current_epoch;
	epoch e2 = current_epoch2;// base

	//	Determine the number of L1 observations that are in the data
	int count = 0;
	for (int i = 0; i < signed(e1.epoch_observations.size()); i++)
	{
		observation curr_obs = e1.epoch_observations[i];
		if ((curr_obs.C1 != 0) && satellite_positions.hasNav(curr_obs.PRN))
			count++;
	}

	//	Initialize the size of the observation vector
	MatrixXd l(count, 1);
	l.setZero();
	double t = e1.time;
	vector<observation> obs1 = e1.epoch_observations;
	vector<observation> obs2 = e2.epoch_observations;// base

	//	Populate the observation vector
	int row = 0;
	int numOfObs = signed(obs1.size());
	for (int i = 0; i < numOfObs; i++)
	{
		observation curr_obs1 = obs1[i];
		observation curr_obs2;
		curr_obs2.C1 = -1;
		for (int j = 0; j < obs2.size(); j++) {
			if (obs2[j].PRN == curr_obs1.PRN) {
				curr_obs2 = obs2[j];
				break;
			}
		}
		if (curr_obs2.C1 == -1) continue;

		if (curr_obs1.C1 == 0 || !satellite_positions.hasNav(curr_obs1.PRN))
			continue;

		//	Get L1 values
		double pseudorange1, pseudorange2;
		pseudorange1 = curr_obs1.C1;
		pseudorange2 = curr_obs2.C1;

		//	Get the satellite position and dts value
		MatrixXd sat_info = satellite_positions.getSatellitePosition(curr_obs1.PRN, t, curr_obs1.C1);

		//	Correct for ionospheric delay
		double iono = ionoL1(pseudorange1, this->nav_data[curr_obs1.PRN - 1][0], this->X, sat_info, t);

		double dts = sat_info(3, 0);
		MatrixXd svX = sat_info.block(0, 0, 3, 1);

		//	Correct the pseudorange for the dts value
		iono = iono + dts * 299792458;

		//	Get the elevation angle of the satellite
		MatrixXd topocentric = topocent(this->X, svX);
		double el = topocentric(1, 0);

		//	Get the tropospheric delay
		double trop = 0; // tropo(el, 1013., 0, 293);

		double range1 = sqrt(pow(svX(0, 0) - Xp(0, 0), 2) + pow(svX(1, 0) - Xp(3, 0), 2) + pow(svX(2, 0) - Xp(6, 0), 2));
		double range2 = sqrt(pow(svX(0, 0) - x_base(0, 0), 2) + pow(svX(1, 0) - x_base(1, 0), 2) + pow(svX(2, 0) - x_base(2, 0), 2));
		//	Correct the pseudorange for tropospheric delay
		//L = L - trop;

		//	Update the observation vector
		//l(row, 0) = (iono - trop) - range - (curr_obs2.C1 - range2);//l(row, 0) = L - D -this->Xp(9, 0);
		l(row, 0) = (pseudorange1 - range1) - (pseudorange2 - range2);
		int tmp = l(row, 0);
		row++;

	}

	MatrixXd ll(row - 1, 1);
	ll.setZero();
#ifdef _IS_DEBUG_
	PrintVector(l, "l");
#endif

	for (int i = 1; i < row; i++) {
		ll(i - 1, 0) = l(i, 0) - l(0, 0);
	}

	this->Z_rou = ll;
#ifdef _IS_DEBUG_
	PrintVector(ll, "ll");
#endif
}


#endif

void SPK::SetPseudorangeCov()
{
	int count = this->Z_rou.rows();
	//	Calculate covariance matrix for the observation vector
	MatrixXd R(count, count);
	MatrixXd Temp(count, count);
	Temp.setOnes(count, count);
	R.setIdentity();
	R = R + Temp;	
	R = R * this->varFac;

	//	Calculate covariance matrix for the predicted observation vector
    MatrixXd PZp(count, count);
	
	PZp = H_rou * this->PXp * H_rou.transpose();

	//	Calculate covariance matrix for the system innovation
	this->R_rou = R + PZp;
	//DEBUG:
#ifdef _IS_DEBUG_
	//PrintMatrix(R, "R");
#endif
	//std::cout << "----------------- R ---------------\n";
	//for (int i = 0; i < count; i++) {
	//	for (int j = 0; j < count; j++) {
	//		std::cout << this->R_rou(i, j) << ", ";
	//	}
	//	std::cout << endl;
	//}
	//std::cout << "-----------------------------------\n";
}

void SPK::SetCarrierphaseCov()
{
	// R_phi
	int count = this->Z_phi.rows();
	MatrixXd R(count, count);
	MatrixXd Temp(count, count);
	R.setOnes();
	Temp.setIdentity();
	R += Temp;
	this->R_phi = R * pow(this->mSigmaPhi * mWaveLengthL1, 2);

	// R_rou
	int count2 = this->Z_rou.rows();
	MatrixXd R2(count2, count2);
	MatrixXd Temp2(count2, count2);
	R2.setOnes();
	Temp2.setIdentity();
	R2 += Temp2;
	this->R_rou = R * pow(this->msigmaRou, 2);

	this->R_rp = MatrixXd(R_rou.rows() + R_phi.rows(), R_rou.cols() + R_phi.cols());
	R_rp.setZero();
	R_rp.block(0, 0, R_rou.rows(), R_rou.cols()) = R_rou;
	R_rp.block(R_rou.rows(), R_rou.cols(), R_phi.rows(), R_phi.cols()) = R_phi;
#ifdef _IS_DEBUG_
	//PrintMatrix(R_rp, "R_rp");
#endif
}

void SPK::setG()
{
	//	Calculate the gain matrix
	//this->G = this->PXp * this->H_rou.transpose() * this->R_rou.inverse();
	this->G = PXp * H_rou.transpose() * (H_rou * PXp * H_rou.transpose() + R_rou).inverse();
#ifdef _IS_DEBUG_
	//PrintMatrix(this->G, "G");
#endif
	/*std::cout << "\n------------- G --------------\n";
	for (int i = 0; i < G.rows(); i++) {
		for (int j = 0; j < G.cols(); j++) {
			cout << G(i, j) << ", ";
		}
		cout << endl;
	}
	std::cout << "\n------------------------------\n"*/;
}

void SPK::SetKalmanGain()
{
#ifdef _IS_CARRIER_PHASE_ONLY_
	this->K = Px * H_phi.transpose() * (H_phi * Px * H_phi.transpose() + R_phi).inverse();
#endif 
#ifdef _IS_PSR_ONLY_
	this->K = Px * H_rp1.transpose() * (H_rp1 * Px * H_rp1.transpose() + R_rou).inverse();
#endif
#ifdef _IS_BOTH_PSR_CARRIERPHASE_
	//***
	this->K = Px * H_rp.transpose() * (H_rp * Px * H_rp.transpose() + R_rp).inverse();
#endif
}

void SPK::setdX()
{
	// DEBUG:
	std::cout << "G size: " << this->G.rows() << " " << this->G.cols() << endl;
	this->dX = this->G * this->Z_rou;// (this->Z_rou - this->H_rou * this->Xp);
}


/*
	This function performs the iteration to conduct the measurement update
	Pre-condition -- do the time update (i.e. have Xp and PXp estimated)
*/
void SPK::ToUpdate()
{
	//	Set the initial estimate of the state vector to the predicted state vector
	this->X = this->Xp;
#ifdef _IS_DEBUG_
	PrintVector(this->X, "X");
#endif
	int count = 0;
	//	Termination condition -- should be false if at termination point
	bool isFinished = false;

	while (!isFinished)
	{		
		//	H
		this->setPseudorangeDesignMatrix();
		//	Z
		this->SetPseudorangeInnovation();
		//	R
		this->SetPseudorangeCov();
		//	G
		this->setG();
		//  Pk
#ifdef _IS_POSITION_STATE_ONLY_
		MatrixXd I = MatrixXd(3, 3);
#endif
#ifdef _IS_PVA_STATE_
		MatrixXd I = MatrixXd(9, 9);
#endif
		I.setIdentity();
		this->Pk = (I - this->G * this->H_rou) * this->PXp;

		//	Calculate the change to the state vector
		this->dX = this->G * this->Z_rou;
#ifdef _IS_DEBUG_
		//** Debug
		PrintVector(this->dX, "dX");
#endif
		//	Update the solution vector
		this->Xp = this->Xp + this->dX;
		this->X = this->Xp;//this->Xp = this->Xp + this->dX;
#ifdef _IS_DEBUG_
		PrintVector(this->X, "X");
#endif
		//	Check to see if iteration should terminate
		if (abs(GetMaxAbsOfVector(this->dX)) < 0.1)
		{
			isFinished = true;
		}
		count++;

		if ((isFinished == false) && (count > 30))
		{
			cout << "Solution failed to converge" << endl;
			break;
		}
	}
}

void SPK::PrintMatrix(const MatrixXd& mat, const string& matName)
{
	std::cout << "\n-----------" << matName << "-------------" << endl;
	for (int i = 0; i < mat.rows(); i++) {
		for (int j = 0; j < mat.cols(); j++) {
			std::cout << std::fixed << setprecision(3) << mat(i, j) << ", ";
		}
		std::cout << endl;
	}
	std::cout << "\n-----------------------------------------" << endl;
}

void SPK::PrintVector(const VectorXd& vec, const string& vecName)
{
	std::cout << "\n-----------" << vecName << "-------------" << endl;
	for (int i = 0; i < vec.rows(); i++) {
		std::cout << std::fixed << setprecision(3) << vec(i, 0) << endl;
	}
	std::cout << "\n-----------------------------------------" << endl;
}

void SPK::solveEpoch(bool& initEstimate)
{
	//	Update the epoch times to properly refer to the current and previous epoch
	this->tprev = this->tcurr;
	this->tcurr = this->current_epoch.time;
	initEstimate = true;// set for debug 
	if (initEstimate) {
		// carrier phase udpate
		NavigationFloatSolution();
	}
	else {
		// PSR udpate
		if (NavsolutionByPseudorange())
			initEstimate = true;
	}

	/*
		Update the residual vectors
	*/
	// to do....
	//this->setVw();
	//this->setVz();// need to change
	//this->setVx();

	/*
		Calculate variance components
	*/

	//this->setVarianceComponents();
}
/*
EKF
1. Predict:
	X = Phi * X
	y = F(x) + v
	  = F(X) + H*dX + v
	=> dZ = y - F(X) = H*dX + v
	Px = Phi*Px*Phi' + Q
	dX = 0
2. Update:
	K  = Px * H' * inv(H*Px*H' + R)
	dX = K * dZ;
	X  = X + dX;
	Px = (I - K*H)*Px
*/
void SPK::solveAllEpochs()
{
	//	Initialize the output stream for positioning results
	ofstream outPoints(outPath + "Point Position.txt");
	outPoints << "SPK POINT COORDINATES (ECEF)" << endl;

	//	Initialize the output stream for process noise residuals
	ofstream outVw(outPath + "Process Noise Residuals.txt");
	outVw << "PROCESS NOISE RESIDUALS" << endl;

	//	Iniitialize output stream for state vector residuals
	ofstream outVx(outPath + "State Vector Residuals.txt");
	outVx << "STATE VECTOR RESIDUALS" << endl;

	//	Initialize output stream for measurement residuals
	ofstream outVz(outPath + "Measurement Residuals.txt");
	outVz << "MEASUREMENT RESIDUALS" << endl;

	ofstream outVarFac(outPath + "Variance Factors.txt");
	outVarFac << "VARIANCE COMPONENTS" << endl;

	IOFormat mat_format(10);

	int count = 0;

	bool initDone = false;// use pseudorange DD model to initalize firstly.
	//	Analyze the next epoch while there are epochs left in the data file
	while (obs_file.hasNextEpoch())
	{
		count++;
		
		//** Update the current epoch
		this->current_epoch  = obs_file.getNextEpoch();//rover obs
		this->current_epoch2 = obs_file2.getNextEpoch();//base obs


		//** Get the solution to the current epoch
		this->solveEpoch(initDone);

		// Output
		// X:
		long long UL1 = 1;
		int time = this->tcurr * 24 * 3600 - mTimeOfFirstsln;
		outPoints << time << ",";
		for (int i = 0; i < this->X.rows(); i++)
			outPoints<< fixed << setprecision(3) << this->X(i, 0) << ",";
		for (int i = 0; i < this->Xn.rows(); i++)
			outPoints << fixed << setprecision(3) << this->Xn(i, 0) << ",";
		for (int i = 0; i < sizeof(mSatBitMask) * 8; i++)
		{
			if ((mSatBitMask & (UL1 << i)) > 0)
			{
				outPoints << i + 1 << ", ";
			}
		}

		outPoints << endl;
	}
	outPoints.close();
}

void SPK::nextEpoch()
{
	this->current_epoch = this->obs_file.getNextEpoch();
}

double SPK::GetMaxAbsOfVector(const VectorXd& vec) {
	if (vec.rows() > 0)
	{
		double max = abs(vec(0,0));
		for (int i = 1; i < vec.rows(); i++) {
			if (max < abs(vec(i, 0))) {
				max = abs(vec(i, 0));
			}
		}
		return max;
	}
	else
		return 0;
}

#ifdef _IS_PVA_STATE_
bool SPK::NavsolutionByPseudorange()
{
	//-----------------------------------------------------------
// 1. Predict
//-----------------------------------------------------------

//this->setXp();
//this->setPXp();
// Xp
	MatrixXd transitMatrix = this->getTransition();
	this->Xp = transitMatrix * this->X;
#ifdef _IS_DEBUG_
	//PrintVector(this->Xp, "Xp");
#endif
	// Pxp
	// time difference between current and previous epoch
	double dt = 24 * 3600 * (this->tcurr - this->tprev);
	// process noise covariance matrix
	MatrixXd Q = this->getQ();
	this->PXp = transitMatrix * this->Pk * transitMatrix.transpose() + Q;

	//-----------------------------------------------------------
	// 2. Update
	//-----------------------------------------------------------
	epoch currObs1 = this->current_epoch;
	epoch currObs2 = this->current_epoch2;//base

	MatrixXd xx(3, 1);
#ifdef _IS_POSITION_STATE_ONLY_
	xx = this->Xp;
#endif
#ifdef _IS_PVA_STATE_
	xx(0, 0) = this->Xp(0, 0);
	xx(1, 0) = this->Xp(3, 0);
	xx(2, 0) = this->Xp(6, 0);
#endif

	//	Filter the observations to ensure they are above the cutoff angle, set all observations below the cutoff angle to zero (to prevent processing)
	for (int i = 0; i < signed(currObs1.epoch_observations.size()); i++)
	{
		//	Get current set of observations for the epoch
		observation obs = currObs1.epoch_observations[i];

		//	Ensure that there is an L1 signal coming from the satellite
		if (obs.C1 == 0)//if (obs.P1 == 0)
			continue;

		//	Get topocentric coordinates of the satellite
		MatrixXd topo = topocent(xx, satellite_positions.getSatellitePosition(obs.PRN, currObs1.time, obs.C1));//(obs.satnum, e.time, obs.P1));

		//	Get elevation angle from the topocentric coordinate vector
		double elevation = topo(1) * 180 / PI;

		//	Ensure that the elevation angle is less than the threshold cutoff
		if (elevation < cutoff || elevation < 0)
		{
			this->current_epoch.epoch_observations[i].C1 = 0; //this->current_epoch.epoch_observations[i].P1 = 0;
		}
	}
	for (int i = 0; i < signed(currObs2.epoch_observations.size()); i++)
	{

		//	Get current set of observations for the epoch
		observation obs = currObs2.epoch_observations[i];

		//	Ensure that there is an L1 signal coming from the satellite
		if (obs.C1 == 0) continue;

		//	Get topocentric coordinates of the satellite
		MatrixXd topo = topocent(xx, satellite_positions.getSatellitePosition(obs.PRN, currObs1.time, obs.C1));

		//	Get elevation angle from the topocentric coordinate vector
		double elevation = topo(1) * 180 / PI;

		//	Ensure that the elevation angle is less than the threshold cutoff
		if (elevation < cutoff || elevation < 0)
		{
			this->current_epoch2.epoch_observations[i].C1 = 0; //this->current_epoch.epoch_observations[i].P1 = 0;
		}
	}

	/*
		Perform the measurement update
	*/

	this->ToUpdate();
	
	//this->setPk();
	if (abs(GetMaxAbsOfVector(this->dX)) < 0.1)
	{
		return true;
	}
	else return false;
}

bool SPK::NavigationFloatSolution()
{
	//------------------------------
	// 1.Predict
	//------------------------------
	/*******************************
	 * X update:
	 * X = Phi1 * X;
	 * N = Phi2 * N
	*/
	int time = this->tcurr * 24 * 3600 - mTimeOfFirstsln;
	
#ifdef _IS_DEBUG_
	PrintVector(this->X, "X(+)");
#endif
	MatrixXd transitionMatrix1 = this->getTransition();
	this->X = transitionMatrix1 * this->X;
#ifdef _IS_DEBUG_
	PrintVector(this->X, "X(-)");
#endif
	CheckValidObs();
	// to form the ambiguities which are not fixed yet.
	long long UL1 = 1;
	int nFloatAmb = 0;
	mSatBitMask = SetSatelltieBitmask();
	MatrixXd TempXn(40, 1);
	TempXn.setZero();
	bool resetFlag = false;
	for (int i = 0; i < sizeof(mSatBitMask)*8; i++) {//8 bits for each byte			
		if ((mSatBitMask & (UL1 << i)) > 0 && i != (mBasePrn-1)) 
		{
			if (!mSatAmbArray[i].isFixed) 
			{
				if (mSatAmbArray[i].rms == 0) mSatAmbArray[i].rms = 5;//if rms=0, means this is first estimation. initialize its rms with 0.
				if (mSatAmbArray[i].value == 0 && !mFistSolution)
				{
					//mFistSolution = true;
					double value = ComputeFloatApproxAmb(i+1);
					mSatAmbArray[i].value = value;
					mSatBitMask &= ~(UL1 << i);
				}
				else 
				{
					TempXn(nFloatAmb, 0) = mSatAmbArray[i].value;
					nFloatAmb++;
				}
			}
		}
		else if (mSatAmbArray[i].value != 0 && i != (mBasePrn - 1))
		{
			mSatAmbArray[i].value = 0;
			mSatAmbArray[i].rms   = 0;
		}
	}
	this->Xn = MatrixXd(nFloatAmb, 1);
	this->Xn = TempXn.block(0, 0, nFloatAmb, 1);

	this->Q = MatrixXd(this->X.rows() + nFloatAmb, this->X.rows() + nFloatAmb);
	this->Q.setIdentity();
	this->Q = GetFloatProcessNoise(nFloatAmb);

#ifdef _IS_DEBUG_
	//PrintMatrix(Q, "float: Q");
#endif
	/********************************************************
	 * Px = [ Pk, 0 ]
	        [ 0 , Pn]
	 *
	 * Px = [ Phi, 0 ]*[ Pk, 0 ]*[ Phi, 0 ]' + [ Qx, 0 ]
	 *      [ 0 ,  I ] [ 0 , I ] [ 0  , I ]    [ 0 , Qn]
	 *
	 *********************************************************/
	int numStates = this->X.rows() + nFloatAmb;
	this->Px = MatrixXd(numStates, numStates);
	Px.setIdentity();

	this->Pk = transitionMatrix1 * this->Pk * transitionMatrix1.transpose();
	MatrixXd Pn = MatrixXd(nFloatAmb, nFloatAmb);
	Pn = GetPn(nFloatAmb);

	Px.block(0, 0, this->X.rows(), this->X.rows()) = this->Pk;
	Px.block(this->X.rows(), this->X.rows(), nFloatAmb, nFloatAmb) = Pn;
	Px = Px + Q;

#ifdef _IS_DEBUG_
	PrintMatrix(Px, "float: Px");
#endif
	//------------------------------
	// 1.Update
	//------------------------------
	//*** Prepare
	//CheckValidObs();
	/******************* first solution LSQ ************************/
	//*** use multi-epochs to calculate the inital values.
	//*** probably two epochs are enough.
	if (mFistSolution)
	{
		// Pseudorange and carrier phase are all used in obs
		// H and Z_phi
		SetCarrierphaseAmbDesignAndInnovationMatrix(nFloatAmb);
		mNumEpoches++;
		
		if (mNumEpoches == 1)
		{
			TempH1 = MatrixXd(H_phi.rows(), H_phi.cols());
			TempZ1 = MatrixXd(Z_phi.rows(), 1);
			TempH1 = H_phi;
			TempZ1 = Z_phi;
		}
		else if(mNumEpoches == 2)
		{
			H_lsq = MatrixXd(TempH1.rows() + H_phi.rows(), TempH1.cols());
			Z_lsq = MatrixXd(TempZ1.rows() + Z_phi.rows(), 1);
			H_lsq << TempH1, H_phi;
			Z_lsq << TempZ1, Z_phi;

			//this->dX = (H_rp.transpose() * H_rp).inverse() * H_rp.transpose() * Z_rp;
			//this->dX = (H_rp.transpose() * H_rp).inverse() * H_rp.transpose() * Z_rp;
			this->dX = (H_lsq.transpose() * H_lsq).inverse() * H_lsq.transpose() * Z_lsq;
			this->X(0, 0) += this->dX(0, 0);
			this->X(3, 0) += this->dX(1, 0);
			this->X(6, 0) += this->dX(2, 0);
			this->Xn      += this->dX.block(3, 0, nFloatAmb, 1);
			/*this->X(0, 0) += this->dX(0, 0);
			this->X(1, 0) += this->dX(1, 0);
			this->X(2, 0) += this->dX(2, 0);
			this->Xn += this->dX.block(3, 0, nFloatAmb, 1);*/
			//* save the float amb values.
			SaveAmbValues();
#ifdef _IS_DEBUG_
			PrintMatrix(H_lsq, "H_lsq");
			PrintMatrix(Z_lsq, "Z_lsq");
			//PrintMatrix(HtH.inverse(), "HH");
			PrintMatrix(dX, "dX");
			//cout << "H row: " << H_phi.rows() << endl
			//	<< "H col: " << H_phi.cols() << endl;
			PrintMatrix(X, "X");
			PrintMatrix(Xn, "Xn");
#endif
			//MatrixXd Hx  = H_lsq.block(0, 0, H_lsq.rows(), 3);
			//MatrixXd HtH = H_phi.transpose() * H_phi;
			mTimeOfFirstsln = this->tcurr*24*3600;
			mFistSolution   = false;//test
			mNumEpoches     = 0;
		}
 	}

	/****************************************************************/
	else
	{
		//*** 
		int count = 0;
		bool isFinished = false;
		if (time >= 0) {
			cout << "Epochs: " << time << endl;
		}

		MatrixXd I = MatrixXd(numStates, numStates);
		I.setIdentity();
		while (!isFinished)
		{
			// H and Z_phi
			SetCarrierphaseAmbDesignAndInnovationMatrix(nFloatAmb);

			// R
			SetCarrierphaseCov();

#ifdef _IS_DEBUG_
			//PrintMatrix(this->R_phi, "R_phi");
#endif
			// K
			SetKalmanGain();

#ifdef _IS_DEBUG_
			PrintMatrix(this->K, "K");
			PrintMatrix(Z_rp, "Z_rp");
#endif
			// Px
			//MatrixXd I = MatrixXd(numStates, numStates);
			//I.setIdentity();
#ifdef _IS_CARRIER_PHASE_ONLY_
			this->Px = (I - this->K * this->H_phi) * this->Px;
			this->dX = this->K * this->Z_phi;
#endif
#ifdef _IS_PSR_ONLY_
			this->Px = (I - this->K * this->H_rou) * this->Px;
			this->dX = this->K * this->Z_rou;
#endif
#ifdef _IS_BOTH_PSR_CARRIERPHASE_
			this->Px = (I - this->K * this->H_rp) * this->Px;
			this->dX = this->K * this->Z_rp;
#endif

#ifdef _IS_DEBUG_
			PrintVector(this->dX, "dX");
			//PrintMatrix(Px, "Px");
#endif
			//	Update the solution vector
			this->X = this->X + this->dX.block(0, 0, this->X.rows(), 1);
			this->Xn = this->Xn + this->dX.block(this->X.rows(), 0, nFloatAmb, 1);
			SaveAmbValues();
#ifdef _IS_DEBUG_
			PrintVector(this->X, "X");
			PrintVector(this->Xn, "Xn");
#endif
			count++;
			if (abs(GetMaxAbsOfVector(this->dX)) < 0.1 || count > 30)
			{
				isFinished = true;
				/*if (count > 30) {
					time = this->tcurr*24*3600 - mTimeOfFirstsln;
				}*/
			}
#ifdef _IS_MEM_DEBUG_
			//cout << "1455: " << count << "time: " << time << endl;
#endif
		}
		return true;
	}
	
	return true;
}

double SPK::ComputeFloatApproxAmb(int satID)
{
	epoch e1 = current_epoch;
	epoch e2 = current_epoch2;
	double t = e1.time;
	vector<observation> obs1 = e1.epoch_observations;
	vector<observation> obs2 = e2.epoch_observations;
	observation curr_obs1;
	observation curr_obs2;
	int rowH1 = 0, rowHx = 0;
	long long UL1 = 1;

	for (int i = 0; i < sizeof(mSatBitMask) * 8; i++)
	{
		int PRN = 0;
		double floatAmb = 0;
		if ((mSatBitMask & (UL1 << i)) > 0)
		{
			PRN = i + 1;
			// amb value:
			floatAmb = mSatAmbArray[i].value;
			// to find the related obs of this satellite PRN
			for (int j = 0; j < obs1.size(); j++) {
				if (obs1[j].PRN == PRN) {
					curr_obs1 = obs1[j];
					break;
				}
			}
			for (int j = 0; j < obs2.size(); j++) {
				if (obs2[j].PRN == PRN) {
					curr_obs2 = obs2[j];
					break;
				}
			}
		}
		else continue;

		//	Get the satellite position and dts value
		double pseudorange1 = curr_obs1.C1;
		double pseudorange2 = curr_obs2.C1;
		double carrierphase1 = curr_obs1.L1;
		double carrierphase2 = curr_obs2.L1;
		MatrixXd sat_info = satellite_positions.getSatellitePosition(curr_obs1.PRN, t, curr_obs1.P1);
		double iono = ionoL1(pseudorange1, this->nav_data[curr_obs1.PRN - 1][0], this->X, sat_info, t);
		double dts = sat_info(3, 0);
		MatrixXd svX = sat_info.block(0, 0, 3, 1);
		iono = iono + dts * mLightSpeed;
		MatrixXd topocentric = topocent(this->X, svX);
		double elevation = topocentric(1, 0);
		double trop = tropo(elevation, 1013., 0, 293);

		// range
		MatrixXd xPos(3, 1);//xPos = this->X;
		xPos(0, 0) = X(0, 0);
		xPos(1, 0) = X(3, 0);
		xPos(2, 0) = X(6, 0);
		double range1 = GetRangeFromCoordinates(xPos, svX);// (pow(svX(0, 0) - X(0, 0), 2) + pow(svX(1, 0) - X(3, 0), 2) + pow(svX(2, 0) - X(6, 0), 2));
		double range2 = GetRangeFromCoordinates(x_base, svX); //sqrt(pow(svX(0, 0) - x_base(0, 0), 2) + pow(svX(1, 0) - x_base(1, 0), 2) + pow(svX(2, 0) - x_base(2, 0), 2));

		//	Initialize the first derivatives with respect to X, Y, Z, and dt
		MatrixXd dx(3, 1);
		/*dx(0, 0) = this->X(0, 0) - sat_info(0, 0);
		dx(1, 0) = this->X(1, 0) - sat_info(1, 0);
		dx(2, 0) = this->X(2, 0) - sat_info(2, 0);*/
		dx(0, 0) = this->X(0, 0) - sat_info(0, 0);
		dx(1, 0) = this->X(3, 0) - sat_info(1, 0);
		dx(2, 0) = this->X(6, 0) - sat_info(2, 0);
		//double range = sqrt(pow(dx(0, 0), 2) + pow(dx(1, 0), 2) + pow(dx(2, 0), 2));
		double l = 0;
		if (curr_obs1.PRN == satID) {
			double amb = ((carrierphase1 - carrierphase2) - (range1 - range2) / mWaveLengthL1) - mSatAmbArray[mBasePrn-1].value;
			return amb;
		}
	}
}

void SPK::SaveAmbValues()
{
	long long UL1 = 1;
	int iPos = 0;
	for (int i = 0; i < sizeof(mSatBitMask) * 8; i++) {//8 bits for each byte
		if ((mSatBitMask & (UL1 << i)) > 0 && i != (mBasePrn - 1) /*&& i != 26 && i != 6 && i != 7 && i != 27*/) {//test
			if (!mSatAmbArray[i].isFixed) {
				mSatAmbArray[i].value = this->Xn(iPos, 0);
				iPos++;
			}
		}
	}
}

MatrixXd SPK::GetPn(const int& sizeOfAmb)
{
	long long UL1 = 1;
	MatrixXd Pn = MatrixXd(sizeOfAmb, sizeOfAmb);
	Pn.setIdentity();
	int nFloatAmb = 0;
	for (int i = 0; i < sizeof(mSatBitMask) * 8; i++) {
		if ((mSatBitMask & (UL1 << i)) > 0 && i != (mBasePrn - 1)/* && i != 26 && i != 6 && i != 7 && i != 27*/) {//test
			if (!mSatAmbArray[i].isFixed && i != (mBasePrn-1)) {
				Pn(nFloatAmb, nFloatAmb) = mSatAmbArray[i].rms * mSatAmbArray[i].rms;
				nFloatAmb++;
			}
		}
	}
	return Pn;
}
MatrixXd SPK::GetQn(const int& sizeOfAmb)
{
	double dt = 24 * 3600 * (this->tcurr - this->tprev);
	MatrixXd Qn = MatrixXd(sizeOfAmb, sizeOfAmb);
	Qn.setOnes();
	Qn *= mSigmaAmbNoise;
	return Qn;
}

MatrixXd SPK::GetFloatProcessNoise(const int& sizeOfAmb) 
{

	// unknown X size: 9, unknown float amb size: sizeOfAmb
	MatrixXd Qxn = MatrixXd(this->X.rows() +sizeOfAmb, this->X.rows() +sizeOfAmb);
	MatrixXd Qx  = MatrixXd(this->X.rows(), this->X.rows());
	MatrixXd Qn  = MatrixXd(sizeOfAmb, sizeOfAmb);
	Qxn.setIdentity();
	Qx.setIdentity();
	Qn.setIdentity();

	//*** Qx
	Qx = getQ();

	//*** Qn
	Qn = GetQn(sizeOfAmb);

	Qxn.block(0, 0, this->X.rows(), this->X.rows()) = Qx;
	Qxn.block(this->X.rows(), this->X.rows(), sizeOfAmb, sizeOfAmb) = Qn;

	return Q;
}

void SPK::CheckValidObs()
{
	epoch currObs1 = this->current_epoch;
	epoch currObs2 = this->current_epoch2;//base

	MatrixXd xx(3, 1);
#ifdef _IS_POSITION_STATE_ONLY_
	xx = this->Xp;
#endif
#ifdef _IS_PVA_STATE_
	xx(0, 0) = this->Xp(0, 0);
	xx(1, 0) = this->Xp(1, 0);
	xx(2, 0) = this->Xp(2, 0);
	/*xx(0, 0) = this->Xp(0, 0);
	xx(1, 0) = this->Xp(3, 0);
	xx(2, 0) = this->Xp(6, 0);*/
#endif

	//	Filter the observations to ensure they are above the cutoff angle, set all observations below the cutoff angle to zero (to prevent processing)
	for (int i = 0; i < signed(currObs1.epoch_observations.size()); i++)
	{
		//	Get current set of observations for the epoch
		observation obs = currObs1.epoch_observations[i];

		//	Ensure that there is an L1 signal coming from the satellite
		if (obs.C1 == 0)//if (obs.P1 == 0)
			continue;

		//	Get topocentric coordinates of the satellite
		MatrixXd topo = topocent(xx, satellite_positions.getSatellitePosition(obs.PRN, currObs1.time, obs.C1));//(obs.satnum, e.time, obs.P1));

		//	Get elevation angle from the topocentric coordinate vector
		double elevation = topo(1) * 180 / PI;

		//	Ensure that the elevation angle is less than the threshold cutoff
		if (elevation < cutoff || elevation < 0)
		{
			this->current_epoch.epoch_observations[i].C1 = 0; //this->current_epoch.epoch_observations[i].P1 = 0;
		}
	}
	for (int i = 0; i < signed(currObs2.epoch_observations.size()); i++)
	{

		//	Get current set of observations for the epoch
		observation obs = currObs2.epoch_observations[i];

		//	Ensure that there is an L1 signal coming from the satellite
		if (obs.C1 == 0) continue;

		//	Get topocentric coordinates of the satellite
		MatrixXd topo = topocent(xx, satellite_positions.getSatellitePosition(obs.PRN, currObs1.time, obs.C1));

		//	Get elevation angle from the topocentric coordinate vector
		double elevation = topo(1) * 180 / PI;

		//	Ensure that the elevation angle is less than the threshold cutoff
		if (elevation < cutoff || elevation < 0)
		{
			this->current_epoch2.epoch_observations[i].C1 = 0; //this->current_epoch.epoch_observations[i].P1 = 0;
		}
	}
}

MatrixXd SPK::SetCarrierphaseDesignInnovMatrix(const int& sizeOfAmb)
{
	int numObs = sizeOfAmb + 1;
	MatrixXd H1;
	if (mFistSolution)
		H1 = MatrixXd(numObs, 3);
	else
		H1 = MatrixXd(numObs, this->X.rows());
	MatrixXd L(numObs, 1);
	H1.setZero();
	L.setZero();
	MatrixXd Lr(numObs, 1);
	Lr.setZero();
	
	epoch e1 = current_epoch;
	epoch e2 = current_epoch2;
	double t = e1.time;
	vector<observation> obs1 = e1.epoch_observations;
	vector<observation> obs2 = e2.epoch_observations;
	observation curr_obs1;
	observation curr_obs2;
	int rowH1 = 0, rowHx = 0;
	long long UL1 = 1;
	
	for (int i = 0; i < sizeof(mSatBitMask) * 8; i++)
	{
		int PRN = 0;	
		double floatAmb = 0;
		if ((mSatBitMask & (UL1 << i)) > 0) 
		{
			PRN = i + 1;
			// amb value:
			floatAmb = mSatAmbArray[i].value;
			// to find the related obs of this satellite PRN
			for (int j = 0; j < obs1.size(); j++) {
				if (obs1[j].PRN == PRN) {
					curr_obs1 = obs1[j];
					break;
				}
			}
			for (int j = 0; j < obs2.size(); j++) {
				if (obs2[j].PRN == PRN) {
					curr_obs2 = obs2[j];
					break;
				}
			}
		}
		else continue;

		//	Get the satellite position and dts value
		double pseudorange1  = curr_obs1.C1;
		double pseudorange2 = curr_obs2.C1;
		double carrierphase1 = curr_obs1.L1;
		double carrierphase2 = curr_obs2.L1;
		MatrixXd sat_info   = satellite_positions.getSatellitePosition(curr_obs1.PRN, t, curr_obs1.P1);
		double iono  = ionoL1(pseudorange1, this->nav_data[curr_obs1.PRN - 1][0], this->X, sat_info, t);
		double dts   = sat_info(3, 0);
		MatrixXd svX = sat_info.block(0, 0, 3, 1);
		iono         = iono + dts * mLightSpeed;
		MatrixXd topocentric = topocent(this->X, svX);
		double elevation     = topocentric(1, 0);
		double trop = tropo(elevation, 1013., 0, 293);

		// range
		MatrixXd xPos(3, 1);//xPos = this->X;
		xPos(0, 0) = X(0, 0);
		xPos(1, 0) = X(3, 0);
		xPos(2, 0) = X(6, 0);

		double range1 = GetRangeFromCoordinates(xPos, svX);
		double range2 = GetRangeFromCoordinates(x_base, svX); 

		//	Initialize the first derivatives with respect to X, Y, Z, and dt
		MatrixXd dx(3, 1);
		if (mFistSolution) {
			dx(0, 0) = this->X(0, 0) - sat_info(0, 0);
			dx(1, 0) = this->X(1, 0) - sat_info(1, 0);
			dx(2, 0) = this->X(2, 0) - sat_info(2, 0);
		}
		else { 
			dx(0, 0) = this->X(0, 0) - sat_info(0, 0);
			dx(1, 0) = this->X(3, 0) - sat_info(1, 0);
			dx(2, 0) = this->X(6, 0) - sat_info(2, 0);		
		}

		//double range = sqrt(pow(dx(0, 0), 2) + pow(dx(1, 0), 2) + pow(dx(2, 0), 2));
		double partialX = (dx(0, 0)) / range1;
		double partialY = (dx(1, 0)) / range1;
		double partialZ = (dx(2, 0)) / range1;
		double l = 0;
		if (curr_obs1.PRN == mBasePrn) {
			if (mFistSolution) {
				H1(0, 0) = partialX;
				H1(0, 1) = partialY;
				H1(0, 2) = partialZ;
			}
			else {
				H1(0, 0) = partialX;
				H1(0, 3) = partialY;
				H1(0, 6) = partialZ;
			}

			L(0, 0)  = mWaveLengthL1 * (carrierphase1 - carrierphase2) - (range1 - range2);
			//if (mFistSolution)
			{
				Lr(0, 0) = (pseudorange1 - range1) - (pseudorange2 - range2);
			}
			l = L(0, 0);
			mSatAmbArray[mBasePrn - 1].value = l / mWaveLengthL1;
		}
		else {
			rowHx++;
			if (mFistSolution) {
				H1(rowHx, 0) = partialX;
				H1(rowHx, 1) = partialY;
				H1(rowHx, 2) = partialZ;
			}
			else {
				H1(rowHx, 0) = partialX;
				H1(rowHx, 3) = partialY;
				H1(rowHx, 6) = partialZ;
			}

			L(rowHx, 0)  = mWaveLengthL1 * (carrierphase1 - carrierphase2 - floatAmb) - (range1 - range2);
			//if (mFistSolution)
			{
				Lr(rowHx, 0) = (pseudorange1 - range1) - (pseudorange2 - range2);
			}
			l = L(rowHx, 0);
		}
#ifdef _IS_DEBUG_
		cout << "PRN: " << PRN << " ...... l: " << l << endl;
#endif
		rowH1++;
	}
	
	// phase
	MatrixXd Hx;
	if (mFistSolution) {
		Hx = MatrixXd(rowHx, 3);
	}
	else {
		Hx = MatrixXd(rowHx, this->X.rows());
	}
	Z_phi = MatrixXd(rowHx, 1);
	Hx.setZero();
	Z_phi.setZero();

	for (int i = 1; i < rowH1; i++) {
		if (mFistSolution) {
			Hx(i - 1, 0) = H1(i, 0) - H1(0, 0);
			Hx(i - 1, 1) = H1(i, 1) - H1(0, 1);
			Hx(i - 1, 2) = H1(i, 2) - H1(0, 2);
		}
		else {
			Hx(i - 1, 0) = H1(i, 0) - H1(0, 0);
			Hx(i - 1, 3) = H1(i, 3) - H1(0, 3);
			Hx(i - 1, 6) = H1(i, 6) - H1(0, 6);
		}
		Z_phi(i - 1, 0) = L(i, 0) - L(0, 0);
	}

	MatrixXd LLr(rowHx, 1);
	LLr.setZero();
	//if (mFistSolution)
	{
		for (int i = 1; i < rowH1; i++) {
			LLr(i - 1, 0) = Lr(i, 0) - Lr(0, 0);
		}
		this->Z_rp = MatrixXd(2 * rowHx, 1);// (LLr.rows() + Z_phi.rows(), LLr.cols());
		Z_rp << LLr, Z_phi;
	}
	// psr
#ifdef _IS_DEBUG_
	//PrintMatrix(L, "L");
	PrintMatrix(Lr, "Lr");
#endif
	return Hx;
}

void SPK::SetCarrierphaseAmbDesignAndInnovationMatrix(const int& sizeOfAmb)
{
	/* H_phi
	 * H_phi = [ Hx, Hn ]
	*/
	MatrixXd Hx;
	if (mFistSolution) {
		H_phi = MatrixXd(sizeOfAmb, 3 + sizeOfAmb);
		Hx = MatrixXd(sizeOfAmb, 3);
	}
	else {
		H_phi = MatrixXd(sizeOfAmb, this->X.rows() + sizeOfAmb);
		Hx = MatrixXd(sizeOfAmb, this->X.rows());
	}
	
	MatrixXd Hn = MatrixXd(sizeOfAmb, sizeOfAmb);
	H_phi.setZero();
	Hx.setZero();
	Hn.setIdentity();
	//*** Hx
	Hx = SetCarrierphaseDesignInnovMatrix(sizeOfAmb);
#ifdef _IS_DEBUG_
	PrintMatrix(Hx, "Hx");
#endif
	//*** Hn
	Hn *= mWaveLengthL1;
	
	//*** H_phi = [ Hx, Hn ]
	H_phi << Hx, Hn;
#ifdef _IS_DEBUG_
	//PrintMatrix(H_phi, "H_phi");
	//PrintMatrix(Z_phi, "Z_phi");
#endif
	if (!mFistSolution) {
		MatrixXd H_rp1 = MatrixXd(Hx.rows(), X.rows() + Hn.cols());
		MatrixXd H_rp2 = MatrixXd(Hx.rows(), X.rows() + Hn.cols());
		H_rp1 << Hx, Hn * 0;
		H_rp2 << Hx, Hn;
		this->H_rp = MatrixXd(Hx.rows() * 2, X.rows() + Hn.cols());
		this->H_rp << H_rp1, H_rp2;
	}
}

double SPK::GetRangeFromCoordinates(const MatrixXd& a, const MatrixXd& b)
{
	return sqrt(pow(a(0, 0) - b(0, 0), 2) + pow(a(1, 0) - b(1, 0), 2) + pow(a(2, 0) - b(2, 0), 2)); 
}
#endif







//void SPK::setVw()
//{
//	MatrixXd QQ = MatrixXd::Zero(3, 3);//MatrixXd::Zero(4, 4);
//	QQ(0, 0) = scaling_parameter_x;
//	QQ(1, 1) = scaling_parameter_y;
//	QQ(2, 2) = scaling_parameter_z;
//	//QQ(3, 3) = scaling_parameter_t;
//	MatrixXd Gamma = getGamma();
//
//	this->vw = QQ * (Gamma.transpose()) * (PXp.inverse()) * G * D;
//}
//
//void SPK::setVx()
//{
//	MatrixXd phi = getTransition();
//	this->vx = phi * Pk * (phi.transpose()) * (PXp.inverse()) * G * D;
//}
//
//void SPK::setVz()
//{
//	int siz = A.rows();
//	MatrixXd I = MatrixXd::Ones(siz, siz);
//	this->vz = (A * G - I) * D;
//}
//
//void SPK::setVarianceComponents()
//{
//	//	Get helping matrices
//	MatrixXd phi = getTransition();
//	MatrixXd Gamma = getGamma();
//	MatrixXd Pw = MatrixXd::Zero(4, 4);
//	Pw(0, 0) = scaling_parameter_x;// covariance of the 
//	Pw(1, 1) = scaling_parameter_y;
//	Pw(2, 2) = scaling_parameter_z;
//	Pw(3, 3) = scaling_parameter_t;
//
//	MatrixXd V(vx.rows() + vw.rows() + vz.rows(), 1);
//	V << vx, vw, vz;
//
//	//	Define the annihilator matrix
//	MatrixXd D1 = phi * Pk * phi.transpose() * PXp.inverse() * G;// vx
//	MatrixXd D2 = Pw * Gamma.transpose() * PXp.inverse() * G;// for vw
//	MatrixXd D3 = A * G - MatrixXd::Identity(D.rows(), D.rows());// for vz
//
//	MatrixXd Dk(D1.rows() + D2.rows() + D3.rows(), D1.cols());
//	Dk << D1, D2, D3;
//
//	//	Define Pi values for different quantities
//
//	//	For state vector
//	MatrixXd P_x = MatrixXd::Zero(Dk.rows(), Dk.rows());
//	for (int i = 0; i < Pk.rows(); i++)
//	{
//		for (int j = 0; j < Pk.cols(); j++)
//		{
//			P_x(i, j) = Pk(i, j);
//		}
//	}
//
//	//	For process noise components
//	MatrixXd P_wx = MatrixXd::Zero(Dk.rows(), Dk.rows());
//	P_wx(11, 11) = scaling_parameter_x;
//	MatrixXd P_wy = MatrixXd::Zero(Dk.rows(), Dk.rows());
//	P_wy(12, 12) = scaling_parameter_y;
//	MatrixXd P_wz = MatrixXd::Zero(Dk.rows(), Dk.rows());
//	P_wz(13, 13) = scaling_parameter_z;
//	MatrixXd P_wt = MatrixXd::Zero(Dk.rows(), Dk.rows());
//	P_wt(14, 14) = scaling_parameter_t;
//
//	//	For observations
//	MatrixXd P_z = MatrixXd::Zero(Dk.rows(), Dk.rows());
//	for (int i = 0; i < vz.rows(); i++)
//	{
//		P_z(i + 15, i + 15) = 1;
//	}
//
//	MatrixXd Pz = MatrixXd::Identity(D.rows(), D.rows());
//
//	//	Overall
//	MatrixXd P_temp = P_wx + P_wy + P_wz + P_wt + P_z + P_x;
//
//	//	Define the Qvv matrix// estimated covariance for the residuals
//	MatrixXd P_vw = Pw * Gamma.transpose() * A.transpose() * PD.inverse() * A * Gamma * Pw;
//	MatrixXd P_vx = Pk * A.transpose() * PD.inverse() * A * Pk;
//	MatrixXd P_vz = (MatrixXd::Identity(A.rows(), G.cols()) - A * G) * Pz;
//	MatrixXd Qvv = MatrixXd::Zero(P_vw.rows() + P_vx.rows() + P_vz.rows(), P_vw.rows() + P_vx.rows() + P_vz.rows());
//	for (int i = 0; i < P_vw.rows(); i++)
//	{
//		for (int j = 0; j < P_vw.cols(); j++)
//		{
//			Qvv(i, j) = P_vw(i, j);
//		}
//	}
//	for (int i = 0; i < P_vx.rows(); i++)
//	{
//		for (int j = 0; j < P_vx.cols(); j++)
//		{
//			Qvv(i + P_vw.rows(), j + P_vw.cols()) = P_vx(i, j);
//		}
//	}
//	for (int i = 0; i < P_vz.rows(); i++)
//	{
//		for (int j = 0; j < P_vz.cols(); j++)
//		{
//			Qvv(i + P_vw.rows() + P_vx.rows(), j + P_vw.cols() + P_vx.cols()) = P_vz(i, j);
//		}
//	}
//
//	//	Define the Wk matrix
//	//MatrixXd Wk = P_temp.inverse() * Dk;
//
//	//	Calculate variance component for observations
//	double Numerator = (V.transpose() * P_wx * V).trace();
//	double Denominator = (Qvv * P_wx).trace();
//
//	this->sigma_wx = Numerator / Denominator;
//
//	Numerator = (V.transpose() * P_wy * V).trace();
//	Denominator = (Qvv * P_wy).trace();
//
//	this->sigma_wy = Numerator / Denominator;
//
//	Numerator = (V.transpose() * P_wz * V).trace();
//	Denominator = (Qvv * P_wz).trace();
//
//	this->sigma_wz = Numerator / Denominator;
//
//	Numerator = (V.transpose() * P_wt * V).trace();
//	Denominator = (Qvv * P_wt).trace();
//
//	this->sigma_wt = Numerator / Denominator;
//
//	Numerator = (V.transpose() * P_x * V).trace();
//	Denominator = (Qvv * P_x).trace();
//
//	this->sigma_x = Numerator / Denominator;
//
//	Numerator = (V.transpose() * P_z * V).trace();
//	Denominator = (Qvv * P_z).trace();
//
//	this->sigma_z = Numerator / Denominator;
//}