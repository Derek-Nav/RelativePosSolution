#include "GPSCalculations.h"
#include <iostream>

//	If both P code signals are available, produces the iono-free combination
double ionopcode(double P1, double P2)
{
	double gamma = pow(1575.42 / 1227.6, 2.);
	return (P2 - gamma * P1) / (1. - gamma);
}

//	Corrects L1 pseudorange for ionospheric delay using 
double ionoL1(double L1, ephemeris e, MatrixXd pos, MatrixXd satpos, double time)
{
	//	Elevation angle and azimuth from user to satellite
	MatrixXd rel = topocent(pos, satpos);
	double E = rel(1); //	Elevation angle
	double A = rel(0); //	Azimuth angle

	//	Get the ECEF coordinates of the point
	double xx = pos(0, 0);
	double yy = pos(1, 0);
	double zz = pos(2, 0);

	//	Define auxiliary parameters
	double psi = 0.0137 / (E + 0.11) - 0.022;
	double lat = asin(zz / sqrt(xx * xx + yy * yy + zz * zz)); //	TODO: Replace with actual calculation later
	double lon = atan2(xx, yy); //	TODO: Replace with actual calculation later

	double phi_i = lat + psi * cos(A);
	if (phi_i > 0.416)
		phi_i = 0.416;
	else if (phi_i < -0.416)
		phi_i = -0.416;

	double lambda_i = lon + psi * sin(A) / cos(phi_i);

	double phi_m = phi_i + 0.064 * cos(lambda_i - 1.617);

	double t = 4.32e4 * lambda_i + time;

	if (t > 86400)
		t -= 86400;
	else if (t < 0)
		t += 86400;

	//	Calculate amplitude of ionospheric delay
	double AMP = e.alpha0 + e.alpha1 * phi_m + e.alpha2 * pow(phi_m, 2) + e.alpha3 * pow(phi_m, 3);
	//	Calculate period of ionospheric delay
	double PER = e.beta0 + e.beta1 * phi_m + e.beta2 * pow(phi_m, 2) + e.beta3 * pow(phi_m, 3);

	double F = 1.0 + 16.0 * pow(0.53 - E, 3);

	double x = 2 * PI * (t - 50400) / PER;

	//	Ionospheric delay
	double Tiono = 0;
	
	if (abs(x) >= 1.57)
		Tiono = F * 5e-9;
	else
		Tiono = F * (5e-9 + AMP * (1 - pow(x, 2) / 2 + pow(x, 4) / 24));

	//	Remove effects of ionospheric delay from the L1 signal
	double corr = L1 - Tiono * 299792458;

	return corr;
}

bool needsAdjustment(MatrixXd dx, double thresh)
{
	for (int i = 0; i < dx.rows(); i++)
	{
		if (abs(dx(i, 0)) > thresh)
			return true;
	}

	return false;
}

//	Returns (az, el, D)'
MatrixXd topocent(MatrixXd x, MatrixXd satpos)
{
	//	Make sure to only use the 3D coordinates in case there are extra entries
	MatrixXd X = x.block(0, 0, 3, 1);
	MatrixXd S = satpos.block(0, 0, 3, 1);

	MatrixXd dx = S - X;

	//	Get the latitude and longitude for the topocentric point
	MatrixXd latlon = getGeocentric(X);
	double latitude = latlon(0, 0);
	double longitude = latlon(1, 0);

	//	Initialize elements of the rotation matrix that will be used
	double cl = cos(longitude);
	double sl = sin(longitude);
	double cb = cos(latitude);
	double sb = sin(latitude);

	MatrixXd F(3, 3);
	F(0, 0) = -sl;
	F(0, 1) = -sb * cl;
	F(0, 2) = cb * cl;
	F(1, 0) = cl;
	F(1, 1) = -sb * sl;
	F(1, 2) = cb * sl;
	F(2, 0) = 0;
	F(2, 1) = cb;
	F(2, 2) = sb;

	//	Apply the rotation matrix to the difference between the coordinates
	MatrixXd topocentric = F * dx;

	double az = atan2(topocentric(0, 0), topocentric(1, 0));
	double D = sqrt(pow(dx(0, 0), 2) + pow(dx(1, 0), 2) + pow(dx(2, 0), 2));
	double el = asin(dx(2, 0) / D);

	MatrixXd res(3, 1);
	res(0, 0) = az;
	res(1, 0) = el;
	res(2, 0) = D;

	return res;
}

//	Returns (phi, lambda, r)'
MatrixXd getGeocentric(MatrixXd X)
{
	double R = sqrt(pow(X(0, 0), 2) + pow(X(1, 0), 2) + pow(X(2, 0), 2));

	double phi = asin(X(2, 0) / R);

	double lambda = atan2(X(1, 0), X(0, 0));

	MatrixXd res(3, 1);
	res(0, 0) = phi;
	res(1, 0) = lambda;
	res(2, 0) = R;

	return res;
}

//	Tropospheric correction
double tropo(double el, double Pd, double Pwv, double T)
{
	//	Calculate zenith angle
	double z = PI / 2 - el;

	//	Calculate tropospheric delay
	double trop = Pd + (1255. / T + 0.05) * Pwv - pow(tan(z), 2.);
	trop = trop * 0.002277 / cos(z);

	return trop;
}
//double tropo(double sinel, double hsta, double p, double tkel, double hum, double hp, double htkel, double hhum)
//{
//	double a_e = 6378.137;
//	double b0 = 7.839257e-5;
//	double tlapse = -6.5;
//	double tkhum = tkel + tlapse * (hhum - htkel);
//	double atkel = 7.5 * (tkhum - 273.15) / (237.3 + tkhum - 273.15);
//	double e0 = 0.0611 * hum * pow(10, atkel);
//	double tksea = tkel - tlapse * htkel;
//	double em = -978.77 / (2.8704e6 * tlapse * 1.0e-5);
//	double tkelh = tksea + tlapse + hhum;
//	double e0sea = e0 * pow(tksea / tkelh, 4. * em);
//	double tkelp = tksea + tlapse * hp;
//	double psea = p * pow(tksea / tkelp, em);
//
//	if (sinel < 0)
//		sinel = 0;
//
//	double tropo = 0;
//	bool flag = false;
//
//	double refsea = 77.624e-6 / tksea;
//	double htop = 1.1385e-5 / refsea;
//	refsea = refsea * psea;
//	double ref = refsea * pow((htop - hsta) / htop, 4.);
//
//	while (true)
//	{
//		double rtop = pow(a_e + htop, 2.) - pow(a_e + hsta, 2) * (1 - pow(sinel, 2.));
//		if (rtop < 0)
//			rtop = 0;
//
//		rtop = sqrt(rtop) - (a_e + hsta) * sinel;
//		double a = -sinel / (htop - hsta);
//		double b = -b0 * (1 - pow(sinel, 2.)) / (htop - hsta);
//		
//		double alpharn;
//		
//		alpharn = pow(rtop, 2) * 2 * a;
//		alpharn = alpharn + pow(rtop, 3) * (2 * pow(a, 2) + 4 * b / 3);
//		alpharn = alpharn + pow(rtop, 4) * (a * (pow(a, 2) + 3 * b));
//		alpharn = alpharn + pow(rtop, 5) * (pow(a, 4) / 5 + 2.4 * pow(a, 2) * b + 1.2 * pow(b, 2));
//		alpharn = alpharn + pow(rtop, 6) * (2 * a * b * (pow(a, 2) + 3 * b));
//		alpharn = alpharn + pow(rtop, 7) * (pow(b, 2) * (6 * pow(a, 2) + 4 * b)) * 1.428571e-1;
//		
//		if (b > 1.0e-35)
//		{
//			alpharn = alpharn + pow(rtop, 8) * a * pow(b, 3) / 2.;
//			alpharn = alpharn + pow(rtop, 9) * pow(b, 4) / 9.;
//		}
//
//		double dr = rtop + alpharn;
//		tropo = tropo + dr * ref * 1000;
//
//		if (flag)
//			break;
//
//		flag = true;
//		refsea = (371900.0e-6 / tksea - 12.92e-6) / tksea;
//		htop = 1.1385e-5 * (1255 / tksea + 0.05) / refsea;
//		ref = refsea * e0sea * pow((htop - hsta) / htop, 4);
//	}
//
//	return tropo;
//}

SatellitePositions::SatellitePositions(vector<vector<ephemeris>> info)
{
	this->ephem = info;
}
bool SatellitePositions::hasNav(int satnum)
{
	return (this->ephem[satnum - 1].size() > 0);
}

MatrixXd SatellitePositions::getSatellitePosition(int satellite_number, double t, double P)
{
	MatrixXd res(4, 1);

	//	Get the ephemeris for the current satellite
	vector<ephemeris> curr_ephemeris = this->ephem[satellite_number - 1];

	if (curr_ephemeris.size() > 0) {
		//	Search through the ephemerides to find the reference epoch
		int epochnum = -1;
		double mindiff(1e14);

		for (int i = 0; i < signed(curr_ephemeris.size()); i++)
		{
			if (abs(curr_ephemeris[i].time - t) < mindiff)
			{
				mindiff = abs(curr_ephemeris[i].time - t);
				epochnum = i;
			}
		}

		//	Retrieve the relevent ephemeris
		ephemeris eph = curr_ephemeris[epochnum];

		//	Convert julian day to GPS time
		double tk = GPST(t);

		//	Get time since reference epoch and correct for travel time
		tk = tk - eph.TOE - P / 2.99792458e8;

		//	Correct for GPS week discrepancies
		if (tk < -302400)
			tk = tk + 604800;
		else if (tk > 302400)
			tk = tk - 604800;

		//	Get semi-major axis
		double A = pow(eph.root_a, 2.);

		//	Calculate the mean motion
		double n0 = sqrt(MU / pow(A, 3));

		//	Calculate the corrected mean motion
		double n = n0 + eph.delta_N;

		//	Calculate the mean anomaly
		double Mk = eph.Mo + n * tk;

		//	Get the eccentricity
		double e = eph.e;
		//	Calculate the eccentric anomaly with Kepler's iterative formula
		double Ek = Mk;
		double old_e = 0;
		do
		{
			old_e = Ek;
			Ek = Mk + e * sin(Ek);
		} while (abs(Ek - old_e) > 1e-6);

		//	Calculate the true anomaly
		double vk = atan2(sqrt(1 - pow(e, 2)) * sin(Ek), cos(Ek) - e);
		//	Get argument of latitude
		double phik = vk + eph.omega;
		//	Calculate argument of latitude correction
		double duk = eph.Cus * sin(2 * phik) + eph.Cuc * cos(2 * phik);
		//	Calculate radius correction
		double drk = eph.Crs * sin(2 * phik) + eph.Crc * cos(2 * phik);
		//	Calculate inclination correction
		double dik = eph.Cis * sin(2 * phik) + eph.Cic * cos(2 * phik);
		//	Calculate corrected argument of latitude
		double uk = phik + duk;
		//	Calculate corrected radius
		double rk = A * (1 - e * cos(Ek)) + drk;
		//	Calculate corrected inclination
		double ik = eph.I + dik + eph.I_DOT * tk;

		//	Calculate positions in orbital plane
		double xkp = rk * cos(uk);
		double ykp = rk * sin(uk);

		//	Calculate corrected longitude of ascending node
		double omk = eph.OMEGA + (eph.OMEGA_DOT - ROTATION) * tk - ROTATION * eph.TOE;

		//	Calculate ECEF coordinates
		double xk = xkp * cos(omk) - ykp * cos(ik) * sin(omk);
		double yk = xkp * sin(omk) + ykp * cos(ik) * cos(omk);
		double zk = ykp * sin(ik);

		//	Get dts value
		double dts = eph.clock_bias + eph.clock_drift * tk + eph.clock_drift_rate * pow(tk, 2);
		dts = dts - 4.442807633e-10 * e * eph.root_a * sin(Ek) - eph.TGD;

		//	Initialize a 4-dimensional vector containing (Xs, Ys, Zs, dts)'
		//MatrixXd res(4, 1);
		res(0, 0) = xk;
		res(1, 0) = yk;
		res(2, 0) = zk;
		res(3, 0) = dts;
	}
	else {
		res(0, 0) = 0;
		res(1, 0) = 0;
		res(2, 0) = 0;
		res(3, 0) = 0;
	}
	
	return res;
}