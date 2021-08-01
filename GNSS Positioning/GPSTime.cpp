#include "GPSTime.h"

double JulianDay(double year, double month, double day, double hour, double minute, double second)
{
	if (month <= 2)
	{
		year = year - 1;
		month = month + 12;
	}

	hour = hour + minute / 60. + second / 3600.;

	double jd = floor(365.25 * (year + 4716)) + floor(30.6001 * (month + 1)) + day + hour / 24.;
	jd = jd - 1537.5;

	return jd;
}

double GPST(double julian_date)
{
	double a = floor(julian_date + 0.5);
	double b = a + 1537;
	double c = floor((b - 122.1) / 365.25);
	double e = floor(365.25 * c);
	double f = floor((b - e) / 30.6001);
	double d = b - e - floor(30.6001 * f) + fmod(julian_date + 0.5, 1.);

	double day_of_week = fmod(floor(julian_date + 0.5), 7);

	double time = (fmod(d, 1.) + day_of_week + 1) * 86400.;

	return time;
}