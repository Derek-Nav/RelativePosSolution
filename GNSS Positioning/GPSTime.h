#pragma once
#include <math.h>

using namespace std;

//	Converts the passed date into its equivalent Julian date
double JulianDay(double year, double month, double day, double hour, double minute, double second);

//	Converts the passed julian date into GPS time
double GPST(double julian_date);