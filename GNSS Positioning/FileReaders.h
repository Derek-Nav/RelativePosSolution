#pragma once
#include <iostream>
#include <fstream>
#include "GPSStructures.h"
#include <Eigen/Dense>
/**

	FILE NAME: FileReaders.h
	DESCRIPTION: Classes and functions for reading GPS navigation and observation files.  Contains
				 two classes: ObservationFileReader, for reading observation files, and
				 NavigationFileReader, for reading navigation files.
	AUTHOR: Benjamin Brunson
	LAST UPDATE: August 30, 2017 - file creation

	TODO: Implementation

*/

namespace readers {

	/*
		SHARED FUNCTIONS
	*/

	//	Removes white space from start and end of a line
	void trim(string&);
	//	Converts all 'D's in a string to 'E's, for parsing numbers
	void convertD(string&);

	//	Returns the header line type
	string headerLineType(string);
	//	Returns the header line contents
	string headerLineValue(string);
	//	Returns the file information contained in the RINEX file
	GPSfile getFileInfo(ifstream&);

	/*
		ObservationFileReader contains a set of functions for reading observation files.
		The class can read one epoch of measurements at a time (for memory conservation)
		or read a whole file.
	*/
	class ObservationFileReader
	{
	public:
		//	GPS file information
		GPSfile file_info;
		ifstream in;

		ObservationFileReader() {}
		//	Default constructor
		ObservationFileReader(string);

		/*
			FUNCTION DEFINITIONS
		*/

		epoch getNextEpoch();

		void close();

		//	Checks if observation file is in the correct format
		bool hasValidFormat();
		//	Checks whether or not the observation type is for GPS-based observations
		bool hasValidObservationType();
		//	Checks whether or not the specified RINEX file is an observation file
		bool hasValidRINEXType();
		//	Checks whether or not there is another epoch listed in the file
		bool hasNextEpoch();

	private:
		
	};

	/*
		NavigationFileReader contains a set of functions for reading GPS navigation files.
		The class reads satellite epoch information
	*/
	class NavigationFileReader
	{
	public:
		GPSfile file_info;
		ifstream in;

		NavigationFileReader() {}
		NavigationFileReader(string);

		vector<vector<ephemeris>> getNavigationData();

		bool hasNextNav();

		void close();

		//	Checks if navigation file is in the correct format
		bool hasValidFormat();
		//	Checks whether or not the specified RINEX file is a navigation file
		bool hasValidRINEXType();

	private:
	};
}