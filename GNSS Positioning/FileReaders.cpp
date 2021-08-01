#include "FileReaders.h"
#include <iostream>
#include <fstream>
#include "GPSStructures.h"
#include <Eigen/Dense>
#include <sstream>
#include "GPSTime.h"

#define LINEWIDTH 500

using namespace std;
using namespace Eigen;
using namespace readers;

void readers::trim(string& str)
{
	bool start = true;
	int count = 0;
	for (int i = 0; i < str.length(); i++)
	{
		start = (str[i] == ' ') || (str[i] == '\t');
		if (!start)
			break;
		count++;
	}

	str.erase(0, count);

	bool end = true;
	count = 0;
	for (int i = str.length() - 1; i >= 0; i--)
	{
		end = (str[i] == ' ') || (str[i] == '\t');
			if (!end)
				break;
		count++;
	}

	str.erase(str.length() - count, str.length() - 1);
}

string readers::headerLineType(string header)
{
	return header.substr(60);
}

string readers::headerLineValue(string header)
{
	return header.substr(0, 60);
}

bool readers::ObservationFileReader::hasNextEpoch()
{
	if (in.peek() == EOF)
		return false;
	return true;
}

GPSfile readers::getFileInfo(ifstream& in)
{
	//	File information is stored in this variable
	GPSfile info;

	//	Iterate through the header lines to get all header information
	while (true)
	{
		//	Get the line and cast it as a string
		string header_line;
		getline(in, header_line);

		//	Separate the line into header type and header value
		string header_type = headerLineType(header_line);
		string header_value = headerLineValue(header_line);
		trim(header_type);

		if (header_type == "END OF HEADER")
			break;

		//	Use header line types to choose which field of info to update,
		//	and update the field with the contents of the header line value
		if (header_type == "COMMENT") //	Skip comment lines
			continue;
		else if (header_type == "ION ALPHA")
		{
			convertD(header_line);
			stringstream str;
			str = stringstream(header_line);
			str.width(18);
			str >> info.alpha0;
			str >> info.alpha1;
			str >> info.alpha2;
			str >> info.alpha3;
		}
		else if (header_type == "ION BETA")
		{
			convertD(header_line);
			stringstream str;
			str = stringstream(header_line);
			str.width(18);
			str >> info.beta0;
			str >> info.beta1;
			str >> info.beta2;
			str >> info.beta3;
		}
		else if (header_type == "RINEX VERSION / TYPE") //	Lines describing the file type, and RINEX version
		{
			string version, type, obs_type;
			version = header_value.substr(0, 20);
			type = header_value.substr(20, 20);
			obs_type = header_value.substr(40, 20);

			stringstream stringer(version);
			stringer >> version;

			info.file_version = version;
			info.file_observation_type = obs_type[0];
			info.RINEX_type = type[0];
		}
		else if (header_type == "MARKER NAME")	//	Name of the GPS station
		{
			trim(header_value);
			info.marker_name = header_value;
		}
		else if (header_type == "MARKER NUMBER") // Number of the GPS station
		{
			trim(header_value);
			info.marker_number = header_value;
		}
		else if (header_type == "APPROX POSITION XYZ")
		{
			trim(header_value);
			stringstream stringer(header_value);
			double X, Y, Z;
			stringer >> X >> Y >> Z;
			Vector3d res;
			res(0) = X;
			res(1) = Y;
			res(2) = Z;

			info.approx_XYZ = res;
		}
		else if (header_type == "ANTENNA: DELTA H/E/N")
		{
			trim(header_value);
			stringstream stringer(header_value);
			double H, E, N;

			stringer >> H >> E >> N;
			Vector3d res;
			res(0) = H;
			res(1) = E;
			res(2) = N;

			info.delta_HEN = res;
		}
		else if (header_type == "# / TYPES OF OBSERV")
		{
			trim(header_value);
			stringstream stringer(header_value);
			int num_obs;
			stringer >> num_obs;

			info.observation_types = vector<string>();

			for (int i = 0; i < num_obs; i++)
			{
				string name;
				stringer >> name;
				info.observation_types.push_back(name);
			}
		}
		else if (header_type == "TIME OF FIRST OBS")
		{
			trim(header_value);
			stringstream stringer(header_value);
			double year, month, day, hour, minute, second;

			stringer >> year >> month >> day >> hour >> minute >> second;

			info.start_time = JulianDay(year, month, day, hour, minute, second);
		}
		else
			continue;
	}

	return info;
}

readers::ObservationFileReader::ObservationFileReader(string file)
{
	this->in = ifstream(file);
	this->file_info = getFileInfo(in);
}

epoch readers::ObservationFileReader::getNextEpoch()
{
	epoch res;

	//	Get line containing all information re: the epoch
	string info_line;
	getline(in, info_line);

	stringstream stringer(info_line);

	//	Get the date and time of observations
	double year, month, day, hour, minute, second;
	stringer >> year >> month >> day >> hour >> minute >> second;
	year = year + 2000;

	res.time = JulianDay(year, month, day, hour, minute, second);

	//	Get the flag for the epoch
	int flag;
	stringer >> flag;

	//	Gets the remainder of the line - use for debugging only
	//string GPSSats = stringer.str();

	//	Number of satellites detected for this epoch
	int numsats;
	stringer.width(3);
	stringer >> numsats;

	//	Get the satellite numbers of each observation
	vector<observation> observ;
	for (int i = 0; i < numsats; i++)
	{
		//	Stores the current observation, to be appended to the vector of observations in this epoch
		observation curr;
		char nuisance;	//	Nuisance variable
		//	Read first of three characters - should be 'G'
		stringer.width(1);
		stringer >> nuisance;
		//	Read last two of three characters - should be the satellite number
		stringer.width(2);
		stringer >> curr.PRN;

		observ.push_back(curr);
	}


	for (int j = 0; j < numsats; j++)
	{
		//	Get the first epoch line
		string epoch_line;
		getline(in, epoch_line);

		//	Get the current observation
		observation curr = observ[j];

		//	Keep a count of the number of rows read
		int numrows = 0;

		for (int i = 0; i < file_info.observation_types.size(); i++)
		{
			//	Go to the next line if there are too many observations to fit on one line
			if ((16 * i - 80 * numrows) >= 80)
			{
				getline(in, epoch_line);
				numrows++;
			}

			//	Read the next observation
			string val_string;
			if ((16 * (i+1) - 80 * numrows) > (epoch_line.length()+2))
			{
				val_string = "";
			}
			else
			{
				val_string = epoch_line.substr(16 * i - 80 * numrows, 16);
				trim(val_string);
			}

			double value;	//	Observation value
			char nuisance;	//	Nuisance variable

			//	Read the observation value
			if (val_string.size() == 0)
				value = 0;
			else
				value = stod(val_string);

			//	Read the nuisance variable
			stringer.width(1);
			stringer >> nuisance;

			//	Get the type of observation i from the file information
			string curr_type = file_info.observation_types[i];

			//	Add the observation value to the current observation set
			if (curr_type == "C1")
				curr.C1 = value;
			else if (curr_type == "C2")
				curr.C2 = value;
			else if (curr_type == "P1")
				curr.P1 = value;
			else if (curr_type == "P2")
				curr.P2 = value;
			else if (curr_type == "L1")
				curr.L1 = value;
			else if (curr_type == "L2")
				curr.L2 = value;
			else if (curr_type == "S1")
				curr.S1 = value;
			else if (curr_type == "S2")
				curr.S2 = value;
		}

		//	Update the observation vector
		observ[j] = curr;
	}

	res.epoch_observations = observ;
	return res;
}

bool readers::ObservationFileReader::hasValidFormat()
{
	if ((file_info.file_version == "2.10") || (file_info.file_version == "2.11"))
		return true;
	return false;
}

bool readers::NavigationFileReader::hasValidFormat()
{
	if ((file_info.file_version == "2.10") || (file_info.file_version == "2.11"))
		return true;
	return false;
}

bool readers::ObservationFileReader::hasValidObservationType()
{
	if (file_info.file_observation_type == 'G' || file_info.file_observation_type == 'g' ||
		file_info.file_observation_type == 'M' || file_info.file_observation_type == 'm')
		return true;
	return false;
}

bool readers::ObservationFileReader::hasValidRINEXType()
{
	if (file_info.RINEX_type == 'O' || file_info.RINEX_type == 'o')
		return true;
	return false;
}

bool readers::NavigationFileReader::hasValidRINEXType()
{
	if (file_info.RINEX_type == 'N' || file_info.RINEX_type == 'n')
		return true;
	return false;
}

readers::NavigationFileReader::NavigationFileReader(string file)
{
	this->in = ifstream(file);
	this->file_info = getFileInfo(in);
}

void readers::convertD(string& line)
{
	for (int i = 0; i < line.length(); i++)
	{
		if ((line[i] == 'D') || (line[i] == 'd'))
			line[i] = 'E';
	}
}

vector<vector<ephemeris>> readers::NavigationFileReader::getNavigationData()
{
	vector<vector<ephemeris>> res;
	res.resize(32);	//	Assuming 32 GPS satellites

	//	Read through all the navigation data
	while (this->hasNextNav())
	{

		//	The current ephemeris
		ephemeris curr;

		//	Set the ionospheric coefficients to be representative values (values taken from
		// "Single frequency ionospheric error correction using coefficients generated from
		//	regional ionospheric data for IRNSS"
		curr.alpha0 = this->file_info.alpha0;//0.93129e-8;
		curr.alpha1 = this->file_info.alpha1; //0.22350e-7;
		curr.alpha2 = this->file_info.alpha2; //-0.59601e-7;
		curr.alpha3 = this->file_info.alpha3; //-0.11919e-6;

		curr.beta0 = this->file_info.beta0;//0.92160e5;
		curr.beta1 = this->file_info.beta1;//0.11469e6;
		curr.beta2 = this->file_info.beta2;//-0.13109e6;
		curr.beta3 = this->file_info.beta3;//-0.58983e6;

		//	Get the first ephemeris line
		string ephemeris_line;
		getline(this->in, ephemeris_line);
		convertD(ephemeris_line);
		stringstream stringy(ephemeris_line);

		//	Get the satellite number, the index in the ephemeris vector
		int prn;
		stringy >> prn;

		//	Get epoch information
		double year, month, day, hour, minute, second;
		stringy >> year >> month >> day >> hour >> minute >> second;
		year = year + 2000;
		curr.time = JulianDay(year, month, day, hour, minute, second);

		//	Read the clock correction coefficients
		stringy.width(19);
		stringy >> curr.clock_bias;
		stringy >> curr.clock_drift;
		stringy >> curr.clock_drift_rate;

		//	Get second ephemeris line
		getline(this->in, ephemeris_line);
		convertD(ephemeris_line);
		stringy = stringstream(ephemeris_line);

		//	Get IODE
		stringy.width(19);
		stringy >> curr.IODE;

		//	Get Crs
		stringy >> curr.Crs;

		//	Get delta n
		stringy >> curr.delta_N;

		//	Get Mo
		stringy >> curr.Mo;

		//	Get third ephemeris line
		getline(this->in, ephemeris_line);
		convertD(ephemeris_line);
		stringy = stringstream(ephemeris_line);

		//	Get Cuc
		stringy.width(19);
		stringy >> curr.Cuc;

		//	Get eccentricity
		stringy >> curr.e;

		//	Get Cus
		stringy >> curr.Cus;

		//	Get the square root of the semi major axis
		stringy >> curr.root_a;

		//	Get fourth ephemeris line
		getline(this->in, ephemeris_line);
		convertD(ephemeris_line);
		stringy = stringstream(ephemeris_line);

		//	Get TOE
		stringy.width(19);
		stringy >> curr.TOE;

		//	Get Cic
		stringy >> curr.Cic;

		//	Get OMEGA
		stringy >> curr.OMEGA;

		//	Get Cis
		stringy >> curr.Cis;

		//	Get fifth ephemeris line
		getline(this->in, ephemeris_line);
		convertD(ephemeris_line);
		stringy = stringstream(ephemeris_line);

		//	Get I
		stringy.width(19);
		stringy >> curr.I;

		//	Get Crc
		stringy >> curr.Crc;

		//	Get omega
		stringy >> curr.omega;

		//	Get OMEGA_DOT
		stringy >> curr.OMEGA_DOT;

		//	Get sixth ephemeris line
		getline(this->in, ephemeris_line);
		convertD(ephemeris_line);
		stringy = stringstream(ephemeris_line);

		//	Get I_DOT
		stringy.width(19);
		stringy >> curr.I_DOT;

		//	Get seventh ephemeris line
		getline(this->in, ephemeris_line);
		convertD(ephemeris_line);
		stringy = stringstream(ephemeris_line);

		//	Get TGD
		stringy.width(19);
		double nuisance;
		stringy >> nuisance;
		stringy >> nuisance;
		stringy >> curr.TGD;

		//	Get eigth ephemeris line
		getline(this->in, ephemeris_line);

		res[prn - 1].push_back(curr);
	}

	return res;
}

bool readers::NavigationFileReader::hasNextNav()
{
	if (in.peek() == EOF)
		return false;
	return true;
}

void readers::NavigationFileReader::close()
{
	this->in.close();
}

void readers::ObservationFileReader::close()
{
	this->in.close();
}