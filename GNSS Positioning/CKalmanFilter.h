#pragma once
//#include "/usr/local/Cellar/eigen/3.3.1/include/eigen3/Eigen/Dense" //Include Eigen Library
#include <Eigen/Dense>
using namespace Eigen;

/*
*	Kalman Filter Class Definition.
*
*	by Fabio Carbone, 23/12/2016
*   www.fabiocarbone.org
*
*	Matrix Dimension must be:
*
*	A: n x n
*	B: n x m
*	H: n x n
*	Q: n x n
*	R: n x n
*	I: n x n
*	X: n x 1
*	U: m x 1
*	Z: n x 1
*	P: n x n
*	K: n x n
*
*/

class CKalmanFilter {

public:

	/* Problem Dimension */
	int n; //State vector dimension
	int m; //Control vector (input) dimension (if there is not input, set to zero)

	/* Fixed Matrix */
	MatrixXf A; //System dynamics matrix
	MatrixXf B; //Control matrix 
	MatrixXf H; //Mesaurement Adaptation matrix
	MatrixXf Q; //Process Noise Covariance matrix
	MatrixXf R; //Measurement Noise Covariance matrix
	MatrixXf I; //Identity matrix

	/* Variable Matrix */
	VectorXf X; //(Current) State vector
	MatrixXf P; //State Covariance
	MatrixXf K; //Kalman Gain matrix

	/* Inizial Value */
	VectorXf X0; //Initial State vector
	MatrixXf P0; //Initial State Covariance matrix

	/*
	* Constructor
	* _n: state vector dimension
	* _m: control vector dimension (if there is not input, set to zero)
	*/
	CKalmanFilter(int _n, int _m);

	/* Set Fixed Matrix (NO INPUT) */
	void SetFixed(MatrixXf _A, MatrixXf _H, MatrixXf _Q, MatrixXf _R);

	/* Set Fixed Matrix (WITH INPUT) */
	void SetFixed(MatrixXf _A, MatrixXf _B, MatrixXf _H, MatrixXf _Q, MatrixXf _R);

	/* Set Initial Value */
	void SetInitial(VectorXf _X0, MatrixXf _P0);

	/* Do prediction (NO INPUT) */
	void Predict(void);

	/* Do prediction (INPUT) */
	void Predict(VectorXf U);

	/* Do correction */
	void Update(VectorXf Z);

};