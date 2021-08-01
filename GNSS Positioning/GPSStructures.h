#pragma once
#include <vector>
#include <string>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

struct ephemeris {
	//	Epoch time
	double time = 0;
	//	Clock correction coefficients
	double clock_bias = 0;
	double clock_drift = 0;
	double clock_drift_rate = 0;
	//	Issue of data ephemeris
	double IODE = 0;
	//	Radius correction sinus and cosinus components
	double Crs = 0;
	double Crc = 0;
	//	Latitude correction sinus and cosinus components
	double Cus = 0;
	double Cuc = 0;
	//	Inclination correction sinus and cosinus components
	double Cis = 0;
	double Cic = 0;
	//	Delta N value
	double delta_N = 0;
	//	Mo angle
	double Mo = 0;
	//	Orbit eccentricity
	double e = 0;
	//	Square root of orbital semi-major axis
	double root_a = 0;
	//	Time of ephemeris
	double TOE = 0;
	//	Right ascension of the ascending node for the orbit
	double OMEGA = 0;
	double OMEGA_DOT = 0; // Time rate of change of OMEGA
	//	Inclination of the orbital plane
	double I = 0;
	double I_DOT = 0; //	Time rate of change of I
	//	Total group delay
	double TGD = 0;
	//	Argument of perigee
	double omega = 0;

	//	Alpha parameters for ionospheric delay calculation
	double alpha0 = 0;
	double alpha1 = 0;
	double alpha2 = 0;
	double alpha3 = 0;

	//	Beta parameters for ionospheric delay calculation
	double beta0 = 0;
	double beta1 = 0;
	double beta2 = 0;
	double beta3 = 0;
};

struct observation {
	//	L1 pseudo-range measurement
	double C1 = 0;
	double C2 = 0;
	//	P-code measurements
	double P1 = 0;
	double P2 = 0;
	//	L1 and L2 carrier phase measurements
	double L1 = 0;
	double L2 = 0;
	//	SNR information
	double S1 = 0;
	double S2 = 0;
	//	Satellite PRN
	int PRN;
};

struct epoch {
	//	Time of observation
	double time;
	//	Vector containing all observations
	vector<observation> epoch_observations;
};

struct GPSfile {
	vector<string> observation_types;
	string marker_name = "";
	string marker_number = "";
	string file_version = "";
	char file_observation_type = 'Z';
	char RINEX_type = 'Z';
	Vector3d approx_XYZ;
	Vector3d delta_HEN;
	double start_time = 0;

	//	Alpha parameters for ionospheric delay calculation
	double alpha0 = 0;
	double alpha1 = 0;
	double alpha2 = 0;
	double alpha3 = 0;

	//	Beta parameters for ionospheric delay calculation
	double beta0 = 0;
	double beta1 = 0;
	double beta2 = 0;
	double beta3 = 0;

};