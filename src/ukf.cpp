#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;
  //use_laser_ = false;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;
  //use_radar_ = false;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.05;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
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
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
  
  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  is_initialized_ = false;
  n_x_ = 5;
  //set augmented dimension
  n_aug_ = 7;
  //define spreading parameter
  lambda_ = 3 - n_aug_;
  weights_ = VectorXd(2 * n_aug_ + 1);
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
  //NIS 

  nis_up_lim_rad_ = 7.815; 
  nis_up_lim_lid_ = 5.991; 
  //number of times nis is above upper boundry
  num_above_lim_radar_ = 0;
  num_above_lim_lidar_ = 0;
  //total observations
  num_obs_rad_ = 0;
  num_obs_lid_ = 0;   
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
  
  double d_t;
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {    
    if (!is_initialized_){      
     //initialize
      double ro = meas_package.raw_measurements_[0];      
      double phi = meas_package.raw_measurements_[1];
      double ro_dot = meas_package.raw_measurements_[2];
      x_ << ro * cos(phi),
            ro * sin(phi),
            ro_dot,
            phi,
            0;
      //!!!TODO:temp P_ update later - take process noise into account??? 
      double a = std_radr_*std_radr_;
      double b = std_radr_*std_radr_;
      double c = std_radrd_*std_radrd_;
      double d = std_radphi_*std_radphi_;
      double e = 1;
      P_ << a,0,0,0,0,
            0,b,0,0,0,
            0,0,c,0,0,
            0,0,0,d,0,
            0,0,0,0,e;
      time_us_ = meas_package.timestamp_;
      is_initialized_ = true;
    }    
    else{
      
      if (use_radar_ == true){
        d_t = (meas_package.timestamp_ - time_us_)/1000000.0;
        time_us_ = meas_package.timestamp_;          
        Prediction(d_t);
        UpdateRadar(meas_package);
      }
    }      
  }
  else { //lidar case
    if (!is_initialized_){
      //initialize
      double px = meas_package.raw_measurements_[0];
      double py = meas_package.raw_measurements_[1];          
      x_ << px,
            py,
            0,
            0, //atan2(py,px),
            0;
      //!!!TODO:temp P_ update later - take process noise into account??? 
      double a = 1; //std_laspx_*std_laspx_;
      double b = 1;//std_laspy_*std_laspy_;
      double c = 1;
      double d = 1;//std_laspx_*std_laspy_;
      double e = 1;
      P_ << a,0,0,0,0,
            0,b,0,0,0,
            0,0,c,0,0,
            0,0,0,d,0,
            0,0,0,0,e; 
      time_us_ = meas_package.timestamp_;
      is_initialized_ = true;
    }
    else {
      if(use_laser_ == true){
        
        d_t = (meas_package.timestamp_ - time_us_)/1000000.0;        
        time_us_ = meas_package.timestamp_;
        Prediction(d_t);        
        UpdateLidar(meas_package);
      }        
    }
  }
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
  VectorXd x_aug = VectorXd(n_aug_);
  x_aug << 0,0,0,0,0,0,0;
  x_aug.block(0, 0, n_x_, 1) << x_;
  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_,n_aug_);
  P_aug << 0,0,0,0,0,0,0,
           0,0,0,0,0,0,0,
           0,0,0,0,0,0,0,
           0,0,0,0,0,0,0,
           0,0,0,0,0,0,0,
           0,0,0,0,0,0,0,
           0,0,0,0,0,0,0;
  P_aug.block(0, 0, n_x_, n_x_) << P_;
  P_aug(n_x_,n_x_) = std_a_ * std_a_;
  P_aug(n_x_+1, n_x_+1) = std_yawdd_  * std_yawdd_;
  //create square root matrix
  MatrixXd A = P_aug.llt().matrixL(); 
  //create augmented sigma points matrix
  MatrixXd s1 = A * sqrt(lambda_ + n_aug_);
  s1.colwise() += x_aug;
  MatrixXd s2 = -1 * A * sqrt(lambda_ + n_aug_);
  s2.colwise() += x_aug;
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  Xsig_aug.block(0, 0, n_aug_, 1) << x_aug;
  Xsig_aug.block(0, 1, n_aug_, n_aug_) << s1;
  Xsig_aug.block(0, 1+n_aug_, n_aug_, n_aug_) << s2;  
  
  //predict sigma points
  double v;
  double yaw;
  double yaw_dot;
  double v_a;
  double v_yaw_ddot;

  for(int i = 0; i < 2 * n_aug_ + 1; i++){
    v =          Xsig_aug(2,i);
    yaw =        Xsig_aug(3,i);
    yaw_dot =    Xsig_aug(4,i);
    v_a =        Xsig_aug(5,i);
    v_yaw_ddot = Xsig_aug(6,i);        
    if (fabs(yaw_dot) >= 0.0001){      
      Xsig_pred_(0, i) = Xsig_aug(0, i) + 
                        v / yaw_dot * (sin(yaw + yaw_dot * delta_t) - sin(yaw)) +
                        1./ 2. * delta_t * delta_t * cos(yaw) * v_a;      
      Xsig_pred_(1, i) = Xsig_aug(1, i) + 
                        v / yaw_dot * (-cos(yaw + yaw_dot * delta_t) + cos(yaw)) +
                        1./ 2. * delta_t * delta_t * sin(yaw) * v_a;
    }
    else{
      Xsig_pred_(0, i) = Xsig_aug(0, i) + 
                        v * cos(yaw) * delta_t +
                        1./ 2. * delta_t * delta_t * cos(yaw) * v_a;
      Xsig_pred_(1, i) = Xsig_aug(1, i) + 
                        v * sin(yaw) * delta_t +
                        1./ 2. * delta_t * delta_t * sin(yaw) * v_a;
    }
    Xsig_pred_(2, i) = Xsig_aug(2, i) + delta_t * v_a;
    Xsig_pred_(3, i) = Xsig_aug(3, i) + yaw_dot * delta_t +
                      1./ 2. * delta_t * delta_t * v_yaw_ddot;
    Xsig_pred_(4, i) = Xsig_aug(4, i) + delta_t * v_yaw_ddot;

  }
  //create vector of weights  
  weights_(0) = lambda_/(lambda_ + n_aug_);
  for (int i = 1; i < 2 * n_aug_ + 1; i++){
    weights_(i) = 1./(2. * (lambda_ + n_aug_));
  }
  //predict state
  MatrixXd C = Xsig_pred_.array().rowwise() * weights_.array().transpose();
  x_ = C.rowwise().sum();
  //predict state covariance matrix
  MatrixXd D = Xsig_pred_.colwise() - x_;
  MatrixXd E = D.array().rowwise() * weights_.array().transpose();
  P_ = E * D.transpose();
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
  int n_z = 2;
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  for (int i = 0; i < 2 * n_aug_ + 1; i++){
    Zsig(0,i) = Xsig_pred_(0,i);
    Zsig(1,i) = Xsig_pred_(1,i);
  }
  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
      z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  //innovation covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 sigma points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z,n_z);
  R <<    std_laspx_*std_laspx_, 0,
          0, std_laspy_*std_laspy_;
  S = S + R;

  VectorXd z_lidar(2);
	z_lidar <<  meas_package.raw_measurements_[0],
	      	    meas_package.raw_measurements_[1];

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 sigma points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
    while (x_diff(4)> M_PI) x_diff(4)-=2.*M_PI;
    while (x_diff(4)<-M_PI) x_diff(4)+=2.*M_PI;
    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }
  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();
  //residual
  VectorXd z_diff = z_lidar - z_pred;
  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  while (x_(3)> M_PI) x_(3)-=2.*M_PI;
  while (x_(3)<-M_PI) x_(3)+=2.*M_PI;
  while (x_(4)> M_PI) x_(4)-=2.*M_PI;
  while (x_(4)<-M_PI) x_(4)+=2.*M_PI;
  P_ = P_ - K*S*K.transpose();

  //NIS
  double e = (z_lidar - z_pred).transpose() * S.inverse() * (z_lidar - z_pred);
  num_obs_lid_ += 1;
  if (e > nis_up_lim_lid_){
    num_above_lim_lidar_ +=1;    
  }
  std::cout << "lidar: e = " << e << std::endl;
  std::cout << num_above_lim_lidar_ / (num_obs_lid_ * 1.) << std::endl;
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
  // Radar updates
  
  int n_z = 3;
  //transform sigma points into measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 sigma points
  // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v  = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // measurement model
    Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
    Zsig(1,i) = atan2(p_y,p_x);                                 //phi
    Zsig(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
  }

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
      z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  //innovation covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z,n_z);
  R <<    std_radr_*std_radr_, 0, 0,
          0, std_radphi_*std_radphi_, 0,
          0, 0,std_radrd_*std_radrd_;
  S = S + R;
  
  VectorXd z_radar(3);

	z_radar <<  meas_package.raw_measurements_[0],
	      	    meas_package.raw_measurements_[1],
		          meas_package.raw_measurements_[2];
  
  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
    while (x_diff(4)> M_PI) x_diff(4)-=2.*M_PI;
    while (x_diff(4)<-M_PI) x_diff(4)+=2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = z_radar - z_pred;

  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  while (x_(3)> M_PI) x_(3)-=2.*M_PI;
  while (x_(3)<-M_PI) x_(3)+=2.*M_PI;
  while (x_(4)> M_PI) x_(4)-=2.*M_PI;
  while (x_(4)<-M_PI) x_(4)+=2.*M_PI;
  P_ = P_ - K*S*K.transpose();
  
  //NIS
  double e = (z_radar - z_pred).transpose() * S.inverse() * (z_radar - z_pred);
  num_obs_rad_ += 1;
  if (e > nis_up_lim_rad_){
    num_above_lim_radar_ +=1;    
  }
  std::cout << "radar: e = " << e << std::endl;
  std::cout << num_above_lim_radar_ / (num_obs_rad_ * 1.) << std::endl;
}
