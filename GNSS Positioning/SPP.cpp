#include "SPP.h"
#include "GPSTime.h"
#include "GPSStructures.h"
#include "FileReaders.h"
#include <Eigen/Dense>

using namespace std;
using namespace readers;

SPP::SPP(string obs_path, string nav_path, string path, double cutoff)
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

	//	Get the navigation data
	if (nav_file.hasNextNav())
		this->nav_data = nav_file.getNavigationData();

	//	Get approximate XYZ data
	this->X = MatrixXd(4, 1);
	X.setZero();
	X(0, 0) = obs_file.file_info.approx_XYZ(0, 0);
	X(1, 0) = obs_file.file_info.approx_XYZ(1, 0);
	X(2, 0) = obs_file.file_info.approx_XYZ(2, 0);

	this->Covar = MatrixXd(4, 4);
	Covar.setZero();

	this->varFac = 1;

	satellite_positions = SatellitePositions(nav_data);

	outPath = path;
}


SPP::~SPP()
{
}

MatrixXd SPP::getA()
{
	epoch e  = current_epoch;

	//	Determine the number of L1 observations that are in the data
	int count = 0;
	for (int i = 0; i < signed(e.epoch_observations.size()); i++)
	{
		observation curr_obs = e.epoch_observations[i];

		if (curr_obs.C1 && satellite_positions.hasNav(curr_obs.PRN))
			count++;
	}

	//	Initialize the size of the first design matrix
	MatrixXd A(count, 4);
	
	double t = e.time;
	vector<observation> obs = e.epoch_observations;
	
	//	For each observation in e, determine the satellite positions
	int row = 0;
	for (int i = 0; i < signed(obs.size()); i++)
	{
		observation curr_obs = obs[i];
		if (curr_obs.C1 == 0 || !satellite_positions.hasNav(curr_obs.PRN))
			continue;

		//	Get the satellite position and dts value
		MatrixXd sat_info = satellite_positions.getSatellitePosition(curr_obs.PRN, t, curr_obs.C1);

		//	Initialize the first derivatives with respect to X, Y, Z, and dt
		MatrixXd dx = this->X.block(0, 0, 3, 1) - sat_info.block(0, 0, 3, 1);
		double D = sqrt(pow(dx(0, 0), 2) + pow(dx(1, 0), 2) + pow(dx(2, 0), 2));
		double partialX = (dx(0, 0)) / D;
		double partialY = (dx(1, 0)) / D;
		double partialZ = (dx(2, 0)) / D;
		double partialT = 1.;

		A(row, 0) = partialX;
		A(row, 1) = partialY;
		A(row, 2) = partialZ;
		A(row, 3) = partialT;
		row++;
	}
	
	return A;
}

MatrixXd SPP::getW()
{
	epoch e = current_epoch;

	//	Determine the number of L1 observations that are in the data
	int count = 0;
	for (int i = 0; i < signed(e.epoch_observations.size()); i++)
	{
		observation curr_obs = e.epoch_observations[i];
		if ((curr_obs.C1 != 0) && satellite_positions.hasNav(curr_obs.PRN))
			count++;
	}

	//	Initialize the size of the observation vector
	MatrixXd l(count, 1);

	double t = e.time;
	vector<observation> obs = e.epoch_observations;

	//	Populate the observation vector
	int row = 0;
	for (int i = 0; i < signed(obs.size()); i++)
	{
		observation curr_obs = obs[i];
		if (curr_obs.C1 == 0 || !satellite_positions.hasNav(curr_obs.PRN))
			continue;

		//	Get L1 values
		double L1;
		L1 = curr_obs.C1;
		
		//	Get the satellite position and dts value
		MatrixXd sat_info = satellite_positions.getSatellitePosition(curr_obs.PRN, t, curr_obs.C1);

		//	Correct for ionospheric delay
		double L = ionoL1(L1, this->nav_data[curr_obs.PRN -1][0], this->X, sat_info, t);
		//double L = L1;

		double dts = sat_info(3, 0);
		MatrixXd S = sat_info.block(0, 0, 3, 1);

		//	Correct the pseudorange for the dts value
		L = L + dts * 299792458;

		//	Get the elevation angle of the satellite
		MatrixXd topocentric = topocent(this->X, S);
		double el = topocentric(1, 0);

		//	Get the tropospheric delay
		double trop = 0; //tropo(el, 1013., 0, 293);

		double D = sqrt(pow(S(0, 0) - X(0, 0), 2) + pow(S(1, 0) - X(1, 0), 2) + pow(S(2, 0) - X(2, 0), 2));

		//	Correct the pseudorange for tropospheric delay
		L = L - trop;


		/*if (abs(trop) > 100)
		{
			cout << trop << endl;
			cout << el << endl;
			cout << L1 << endl;
			cin.ignore();
		}*/

		/*cout.precision(10);
		cout << "Tropospheric effect: ";
		cout << trop << " m" << endl;
		cin.ignore();
		*/

		//	Update the observation vector
		l(row, 0) = L - D - this->X(3, 0);
		row++;
	
	}

	cout << l << endl;

	return l;
}

MatrixXd SPP::getdx(MatrixXd A, MatrixXd w)
{
	return (A.transpose() * A).inverse() * A.transpose() * w;
}

MatrixXd SPP::getCovar(MatrixXd A)
{
	return this->varFac * (A.transpose() * A).inverse();
}

void SPP::solveAllEpochs()
{
	//	Initialize the output stream
	ofstream outPoints(outPath + "Point Position.txt");
	outPoints << "SPP POINT COORDINATES" << endl;

	IOFormat mat_format(10);

	//	Analyze the next epoch while there are epochs left in the data file
	while (obs_file.hasNextEpoch())
	{

		//	Update the current epoch
		this->current_epoch = obs_file.getNextEpoch();

		//	Get the solution to the current epoch
		solveEpoch();

		//	Output the point coordinates
		outPoints << this->X.transpose().format(mat_format) << endl;
	}

	outPoints.close();
}

void SPP::nextEpoch()
{
	this->current_epoch = obs_file.getNextEpoch();
}

void SPP::PrintVector(const VectorXd& vec, const string& vecName)
{
	std::cout << vecName << endl;
	for (int i = 0; i < vec.rows(); i++) {
		std::cout << std::fixed << setprecision(3) << vec(i, 0) << endl;
	}
}

void SPP::solveEpoch()
{
	epoch e = current_epoch;

	//	Initialize approximate solution
	MatrixXd init(4, 1);
	init(0, 0) = obs_file.file_info.approx_XYZ(0, 0);
	init(1, 0) = obs_file.file_info.approx_XYZ(1, 0);
	init(2, 0) = obs_file.file_info.approx_XYZ(2, 0);
	init(3, 0) = 0;

	this->X(0, 0) = init(0, 0);
	this->X(1, 0) = init(1, 0);
	this->X(2, 0) = init(2, 0);

	PrintVector(this->X, "X");

	MatrixXd A;
	MatrixXd w;
	MatrixXd dx;

	//	Separating the coordinates of the station from the t0 parameter for determining the elevation angles to the satellites
	MatrixXd coord = this->X.block(0, 0, 3, 1);

	//	Filter the observations to ensure they are above the cutoff angle, set all observations below the cutoff angle to zero (to prevent processing)
	for (int i = 0; i < signed(e.epoch_observations.size()); i++)
	{
		//	Get current set of observations for the epoch
		observation obs = e.epoch_observations[i];

		//	Ensure that there is an L1 signal coming from the satellite
		if(obs.C1 == 0) //if (obs.P1 == 0)
			continue;

		//	Get topocentric coordinates of the satellite
		MatrixXd topo = topocent(coord, satellite_positions.getSatellitePosition(obs.PRN, e.time, obs.C1));//obs.P1));

		//	Get elevation angle from the topocentric coordinate vector
		double elevation = topo(1) * 180 / PI;

		//	Ensure that the elevation angle is less than the threshold cutoff
		if (elevation < cutoff || elevation < 0)
		{
			//this->current_epoch.epoch_observations[i].P1 = 0; 
			this->current_epoch.epoch_observations[i].C1 = 0;
		}
	}

	//do
	//{
	//	Get the design matrix and the misclosure vector
	A = getA();
	w = getW();

	//	Calculate dx
	dx = getdx(A, w);

	// debug:
	cout << "X size: " << this->X.size() << endl
			<< "dx size: " << dx.size() << endl;
	this->X(0, 0) += dx(0, 0);
	this->X(1, 0) += dx(1, 0);
	this->X(2, 0) += dx(2, 0);
	//} while (needsAdjustment(dx, 0.1));


	this->Covar = getCovar(A);
}
