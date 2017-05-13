#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.85; 

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  
  // Will be set to true once 1st measurement is taken
  is_initialized_ = false;

  // No. of state variables
  n_x_ = 5;

  // No. of state + augmented variables
  n_aug_ = 7;

  // Spread factor
  lambda_ = 3 - n_x_;

  // Predicted sigma points matrix
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  // Weights
  weights_ = VectorXd(2 * n_aug_ + 1);
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
	if (is_initialized_ == false)
	{
		if (meas_package.sensor_type_ == MeasurementPackage::LASER)
		{
			// Initialize the state
			double px = meas_package.raw_measurements_(0);
			double py = meas_package.raw_measurements_(1);
			x_ << px, py, 1, 1, 1;

			// Initialize the state covariance matrix
			P_ = MatrixXd::Identity(5, 5);
			P_(0, 0) = std_laspx_ * std_laspx_;
			P_(1, 1) = std_laspy_ * std_laspy_;

			// Initialize time
			time_us_ = meas_package.timestamp_;
		}
		if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
		{
			// Initialize the state
			double ro = meas_package.raw_measurements_(0);
			double phi = meas_package.raw_measurements_(1);
			//double ro_dot = meas_package.raw_measurements_(2);
			double px = ro * cos(phi);
			double py = ro * sin(phi);
			x_ << px, py, 1, 1, 1;

			// Initialize the state covariance matrix
			P_ = MatrixXd::Identity(5, 5);
			// Variance of product of two Gaussians
			double var_radar = sqrtf((powf(std_radr_, 2)*powf(std_radphi_, 2)) / (powf(std_radr_, 2) + powf(std_radphi_, 2)));
			P_(0, 0) = var_radar;
			P_(1, 1) = var_radar;

			// Initialize time
			time_us_ = meas_package.timestamp_;
		}
	}
		
	// Compute delta_t in seconds and initialize previous time
	double delta_t = (meas_package.timestamp_ - time_us_) / pow(10, 6);
	time_us_ = meas_package.timestamp_;

	// Predict step
	if (delta_t > 0.001)
	{
		Prediction(delta_t);
	}

	// Update step based on measurement
	if (meas_package.sensor_type_ == MeasurementPackage::LASER && is_initialized_ == true)
	{
		UpdateLidar(meas_package);
	}
	else if (meas_package.sensor_type_ == MeasurementPackage::RADAR && is_initialized_ == true)
	{
		UpdateRadar(meas_package);
	}

	is_initialized_ = true;
	return;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

    //create augmented mean vector
	VectorXd x_aug = VectorXd(7);

	//create augmented state covariance
	MatrixXd P_aug = MatrixXd(7, 7);

	//create sigma point matrix
	MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

	//create augmented mean state
	x_aug << x_, 0, 0;

	//create process noise matrix
	MatrixXd Qv = MatrixXd(2, 2);
	Qv << std_a_*std_a_, 0,
		  0, std_yawdd_*std_yawdd_;

	// create augmented state covariance matrix
	P_aug << P_, MatrixXd::Zero(5, 2),
			 MatrixXd::Zero(2, 5), Qv;

	//create square root matrix
	MatrixXd A = MatrixXd(7, 7);
	A = P_aug.llt().matrixL();

	//create augmented sigma points
	Xsig_aug.col(0) = x_aug;
	for (int i = 0; i<n_aug_; i++) {
		Xsig_aug.col(i + 1) = x_aug + sqrt(lambda_ + n_aug_) * A.col(i);
		Xsig_aug.col(i + 1 + n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * A.col(i);
	}

	//create process model
	//VectorXd F = VectorXd(5);
	
	//create total process noise function of time and state
	//VectorXd Q = VectorXd(5);

	//place holder x vector
	VectorXd x_place = VectorXd(7);

	for (int i = 0; i<2 * n_aug_ + 1; i++) {
		// place holder for augmented state variables
		x_place = Xsig_aug.col(i);
		double px = x_place(0);
		double py = x_place(1);
		double v = x_place(2);
		double yaw = x_place(3);
		double yaw_dot = x_place(4);
		double nu_a = x_place(5);
		double nu_yawdd = x_place(6);

		// compute Q
		double q0 = 0.5*delta_t*delta_t*cos(yaw)*nu_a;
		double q1 = 0.5*delta_t*delta_t*sin(yaw)*nu_a;
		double q2 = delta_t*nu_a;
		double q3 = 0.5*delta_t*delta_t*nu_yawdd;
		double q4 = delta_t*nu_yawdd;
		//Q << q0, q1, q2, q3, q4;

		if (fabs(yaw_dot) > 0.0001) {
			Xsig_pred_(0, i) = px + v * (sin(yaw + yaw_dot*delta_t) - sin(yaw)) / yaw_dot + q0;
			Xsig_pred_(1, i) = py + v * (-cos(yaw + yaw_dot*delta_t) + cos(yaw)) / yaw_dot + q1;
			Xsig_pred_(2, i) = v + 0.0 + q2;
			Xsig_pred_(3, i) = yaw + yaw_dot * delta_t + q3;
			Xsig_pred_(4, i) = yaw_dot + 0.0 + q4;
			//F << f0, f1, f2, f3, f4;
		}
		else {
			Xsig_pred_(0, i) = px + v * cos(yaw) * delta_t + q0;
			Xsig_pred_(1, i) = py + v * sin(yaw) * delta_t + q1;
			Xsig_pred_(2, i) = v + 0.0 + q2;
			Xsig_pred_(3, i) = yaw + yaw_dot * delta_t + q3;
			Xsig_pred_(4, i) = yaw_dot + 0.0 + q4;
			//F << f0, f1, f2, f3, f4;
		}
		//Xsig_pred_.col(i) = x_place.head(5) + F + Q;
	}

	// set weights
	double weight_0 = lambda_ / (lambda_ + n_aug_);
	weights_(0) = weight_0;
	for (int i = 1; i<2 * n_aug_ + 1; i++) {  //2n+1 weights
		double weight = 0.5 / (n_aug_ + lambda_);
		weights_(i) = weight;
	}

	//predicted state mean
	x_.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
		x_ = x_ + weights_(i) * Xsig_pred_.col(i);
	}

	//predicted state covariance matrix
	P_.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
        // state difference
		VectorXd x_diff = Xsig_pred_.col(i) - x_;
		//angle normalization
		x_diff(3) = x_diff(3) - 2 * M_PI * floor((x_diff(3) + M_PI) / (2 * M_PI));
		P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
	}
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
	// read actual lidar measurements
	VectorXd z = VectorXd(2);
	z = meas_package.raw_measurements_;

	// lidar measurement matrix
	MatrixXd H = MatrixXd(2, 5);
	H << 1, 0, 0, 0, 0,
		0, 1, 0, 0, 0;

	// error in measurement
	VectorXd y = VectorXd(2);
	VectorXd z_pred = VectorXd(2);
	z_pred = H*x_;
	y = z - z_pred;

	// measurement noise
	MatrixXd R = MatrixXd(2, 2);
	R << std_laspx_*std_laspx_, 0,
		0, std_laspy_*std_laspy_;

	// compute S
	MatrixXd PHt = P_*H.transpose();
	MatrixXd S = MatrixXd(2, 2);
	S = H*PHt + R;

	// compute Kalman gain
	MatrixXd K = MatrixXd(5, 2);
	K = PHt*S.inverse();

	//update state and state covariance matrix
	x_ = x_ + K*y;
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K*H)*P_;

	// NIS LIDAR
	NIS_laser_ = y.transpose()*S.inverse()*y;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
   // set measurement dimension
	int n_z = 3;

  //create matrix for sigma points in measurement space
	MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
	Zsig.fill(0.0);

	//mean predicted measurement
	VectorXd z_pred = VectorXd(n_z);
	z_pred.fill(0.0);

	//measurement covariance matrix S
	MatrixXd S = MatrixXd(n_z, n_z);
	S.fill(0.0);

	// place holder for state
	VectorXd xsp = VectorXd(n_x_);
	xsp.fill(0.0);

	// place holder for difference
	VectorXd diff_z = VectorXd(n_x_);
	diff_z.fill(0.0);

	// measurement covariance matrix
	MatrixXd R = MatrixXd(3, 3);
	R << std_radr_*std_radr_, 0, 0,
		0, std_radphi_*std_radphi_, 0,
		0, 0, std_radrd_*std_radrd_;

	for (int i = 0; i<2 * n_aug_ + 1; i++) {
		//transform sigma points into measurement space
		xsp = Xsig_pred_.col(i);
		double px = xsp(0);
		double py = xsp(1);
		double v = xsp(2);
		double yaw = xsp(3);
		double yaw_dot = xsp(4);

		// check for 0
		double range = sqrt(pow(px, 2) + pow(py, 2));
		if (range<0.01) {
			range = 0.01;
		}
		double bearing = atan2(py, px);
		double range_rate = (px*v*cos(yaw) + py*v*sin(yaw)) / range;
		Zsig.col(i) << range, bearing, range_rate;
		//calculate mean predicted measurement
		z_pred = z_pred + weights_(i)*Zsig.col(i);
	}

	for (int i = 0; i<2 * n_aug_ + 1; i++) {
		//calculate measurement covariance matrix S
		diff_z = Zsig.col(i) - z_pred;
		diff_z(1) = diff_z(1) - 2 * M_PI * floor((diff_z(1) + M_PI) / (2 * M_PI));
		S = S + weights_(i)*diff_z*diff_z.transpose();
	}

	S = S + R;

	//create matrix for cross correlation Tc
	MatrixXd Tc = MatrixXd(n_x_, n_z);
	Tc.fill(0.0);

	//calculate cross correlation matrix
	VectorXd x_diff = VectorXd(5);
	x_diff.fill(0.0);
	VectorXd z_diff = VectorXd(3);
	z_diff.fill(0.0);

	for (int i = 0; i<2 * n_aug_ + 1; i++) {
		x_diff = Xsig_pred_.col(i) - x_;
		z_diff = Zsig.col(i) - z_pred;
		//angle normalization
		x_diff(3) = x_diff(3) - 2 * M_PI * floor((x_diff(3) + M_PI) / (2 * M_PI));
		z_diff(1) = z_diff(1) - 2 * M_PI * floor((z_diff(1) + M_PI) / (2 * M_PI));
		Tc = Tc + weights_(i)*x_diff*z_diff.transpose();
	}

	//calculate Kalman gain K;
	MatrixXd K = MatrixXd(5, 3);
	K = Tc*S.inverse();

	//read measurement
	VectorXd z = VectorXd(3);
	z = meas_package.raw_measurements_;

	//update state mean and covariance matrix
	VectorXd diff = VectorXd(3);
	diff.fill(0.0);
	diff = z - z_pred;
	diff(1) = diff(1) - 2 * M_PI * floor((diff(1) + M_PI) / (2 * M_PI));
	x_ = x_ + K*diff;
	P_ = P_ - K*S*K.transpose();

	// NIS RADAR
	NIS_radar_ = diff.transpose()*S.inverse()*diff;
	return;
}
