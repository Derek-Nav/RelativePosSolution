#include "GPSStructures.h"
#include "FileReaders.h"
#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include "GPSTime.h"
#include "SPP.h"
#include "SPK.h"

using namespace std;
using namespace readers;
using namespace Eigen;
Matrix3d x_base(3, 3);

int main()
{
	// rover
	string nav_path1 = "DyncRover1.nav";
	string obs_path1 = "DyncRover1.obs";
	// base
	string nav_path2 = "DyncBase1.nav";
	string obs_path2 = "DyncBase1.obs";

	// ecef:   
	/*x_base(0, 0) = -2420720.52;
	x_base(1, 0) = -3439857.26;
	x_base(2, 0) = 4778532.67;*/
	x_base(0, 0) = 8.568606040579976e5;
	x_base(1, 0) = -4.528420217245726e6;
	x_base(2, 0) = 4.394516244444999e6;
	int basePrn  = 8;

	cout << "---------- Positioning Started ---------" << endl;
	//SPP Positioner(obs_path, nav_path, "", 15);
	SPK Positioner(obs_path1, nav_path1, obs_path2, nav_path2, "", basePrn, 10);//cut off ele angle.

	/*for (int i = 0; i < 150; i++)
	{
		Positioner.nextEpoch();
	}*/

	Positioner.solveAllEpochs();

	cout << "Positioning Complete" << endl;
	

	//cin.ignore();
	return 0;
}
