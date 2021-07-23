#include "ukf.h"
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;

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
  std_a_ = 0.8;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.6;
  
  /**
   * DO NOT MODIFY measurement noise values below.
   * These are provided by the sensor manufacturer.
   */

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
   * End DO NOT MODIFY section for measurement noise values 
   */
  
  is_initialized_ = false;

  //state dimension [pos1, pos2, vel_abs, yaw_angle, yaw_rate]
  n_x_ = 5;

  //augmented state dimension
  n_aug_ = 7;

  //Spreading parameter
  lambda_ = 0;

  //Sigma point count
  n_aug_sigma_ = 2 * n_aug_ + 1;

  //Matrix for sigma points
  Xsig_pred_ = MatrixXd(n_x_, n_aug_sigma_);

  //Vector for weights
  weights_ = VectorXd(n_aug_sigma_);
  double t = lambda_ + n_aug_;
  weights_(0) = lambda_ / t;
  weights_.tail(n_aug_sigma_ - 1).fill(0.5 / t);
 
  //Start time
  time_us_ = 0;

  //NIS
  NIS_radar_ = 0;
  NIS_laser_ = 0;

  //state vector 
  x_ = VectorXd(n_x_);
  x_.fill(0);

  //state covariance matrix
  P_ = MatrixXd(n_x_, n_x_);

  Xsig_pred_ = MatrixXd(n_x_, n_aug_sigma_);
  Xsig_aug_ =  MatrixXd(n_aug_, n_aug_sigma_); 

}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  //Make sure you switch between lidar and radar
  if(!is_initialized_)  {
    if(meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      double rho = meas_package.raw_measurements_(0);
      double phi = meas_package.raw_measurements_(1);

      double px = cos(phi) * rho;
      double py = sin(phi) * rho;

      x_ << px, py, 0, 0, 0;
    
    }
    else if(meas_package.sensor_type_ == MeasurementPackage::LASER) {
      x_ << meas_package.raw_measurements_(0), meas_package.raw_measurements_(1), 0, 0, 0.0; 
    } else {
      return;
    }
    //State covariance matrix
    P_ << 1, 0, 0, 0, 0,
          0, 1, 0, 0, 0,
          0, 0, 1, 0, 0,
          0, 0, 0, 1, 0,
          0, 0, 0, 0, 1;

    is_initialized_ = true;
    time_us_ = meas_package.timestamp_;
    return;
  }

  //Calculate delta_t and store current time 
  double delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0;
  time_us_ = meas_package.timestamp_;

  while(delta_t > 0.1){
    const double dt = 0.05;
    Prediction(dt);
    delta_t -= dt;
  }
  Prediction(delta_t);

  //Update measurements
  if(meas_package.sensor_type_ == MeasurementPackage::RADAR){
    UpdateRadar(meas_package);    
  } else {
    UpdateLidar(meas_package);
  }
}

/**
  * Estimate the object's location. 
  * Modify the state vector, x_. Predict sigma points, the state,    
  * and the state covariance matrix.
*/
void UKF::Prediction(double delta_t) {
  //augmented mean vector
  VectorXd x_aug_ = VectorXd(n_aug_);
  x_aug_.fill(0.0);  

  //augmented state covariance
  MatrixXd P_aug_ = MatrixXd(n_aug_, n_aug_);

  //augmented mean state
  x_aug_.head(n_x_) = x_;  

  //augmented covariance matrix  
  P_aug_.fill(0.0);
  P_aug_.topLeftCorner(n_x_, n_x_) = P_;
  P_aug_(n_aug_ - 2, n_aug_ - 2) = std_a_ * std_a_;
  P_aug_(n_aug_ - 1, n_aug_ - 1) = std_yawdd_ * std_yawdd_;

  //square root matrix
  MatrixXd A_aug = P_aug_.llt().matrixL();

  //augmented sigma points
  Xsig_aug_.col(0) = x_aug_;
  for(int i = 0; i < n_aug_; i++) {
    Xsig_aug_.col(i + 1) = x_aug_ + std::sqrt(lambda_ + n_aug_) * A_aug.col(i);
    Xsig_aug_.col(i + 1 + n_aug_) = x_aug_ - std::sqrt(lambda_ + n_aug_) * A_aug.col(i);
  }

  //predict sigma points
  for(int i = 0; i < n_aug_sigma_; i++) {
    double px = Xsig_aug_(0, i);
    double py = Xsig_aug_(1, i);
    double v = Xsig_aug_(2, i);
    double yaw = Xsig_aug_(3, i);
    double yawd = Xsig_aug_(4, i);
    double v_aug = Xsig_aug_(5, i);
    double v_yawdd = Xsig_aug_(6, i);  

    double px_p;
    double py_p;

    //avoid division by zero
    if (fabs(yawd) > 0.001) {
      px_p = px + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
      py_p = py + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    }
    else {
      px_p = px + v*delta_t*cos(yaw);
      py_p = py + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5*v_aug*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*v_aug*delta_t*delta_t * sin(yaw);
    v_p = v_p + v_aug*delta_t;

    yaw_p = yaw_p + 0.5*v_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + v_yawdd*delta_t;

    //write predicted sigma point into right column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }

  //predict state mean
  x_.fill(0.0);
  for(int i = 0; i < n_aug_sigma_; i++) {
    x_ += weights_(i) *  Xsig_pred_.col(i);
  }

  P_.fill(0.0);
  for(int i = 0; i < n_aug_sigma_; i++) {
    //predict state covariance matrix
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    //normalize angles
    while(x_diff(3) > M_PI)
      x_diff(3) -= 2 * M_PI;
    while(x_diff(3) < -M_PI)
      x_diff(3) += 2 * M_PI;

    P_ += weights_(i) * x_diff * x_diff.transpose();
  }
}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * Use lidar data to update the belief about the object's position.
   * Modify the state vector, x_, and covariance, P_.
   * also calculate the lidar NIS.
   */
  //set measurement dimension, lidar can measure px, py
  int n_z = 2;
  
  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  
  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  
  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  
  Zsig.fill(0.0);
  z_pred.fill(0.0);
  S.fill(0.0);
  
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    //transform sigma points into measurement space
    VectorXd state_vec = Xsig_pred_.col(i);
    double px = state_vec(0);
    double py = state_vec(1);
    
    Zsig.col(i) << px,
                   py;
    
    //calculate mean predicted measurement
    z_pred += weights_(i) * Zsig.col(i);
  }
  
  //calculate measurement covariance matrix S
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    VectorXd z_diff = Zsig.col(i) - z_pred;
    S += weights_(i) * z_diff * z_diff.transpose();
  }
  
  // Add R to S
  MatrixXd R = MatrixXd(2,2);
  R << std_laspx_*std_laspx_, 0,
       0, std_laspy_*std_laspy_;
  S += R;
  
  //vector for incoming radar measurement
  VectorXd z = VectorXd(n_z);
  
  double meas_px = meas_package.raw_measurements_(0);
  double meas_py = meas_package.raw_measurements_(1);
  
  z << meas_px,
       meas_py;
  
  //matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);
  
  //cross correlation matrix
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    
    //normalize angles
    while (x_diff(3) > M_PI) 
      x_diff(3) -= 2. * M_PI;
    while (x_diff(3) < -M_PI) 
      x_diff(3) += 2. * M_PI;
    
    
    VectorXd z_diff = Zsig.col(i) - z_pred;

    Tc += weights_(i) * x_diff * z_diff.transpose();

  }
  
  // difference
  VectorXd z_diff = z - z_pred;
  
  //calculate NIS
  NIS_laser_ = z_diff.transpose() * S.inverse() * z_diff;
  
  //calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();
  
  //update state mean and covariance matrix
  x_ += K*z_diff;
  P_ -= K*S*K.transpose();

}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * Use radar data to update the belief about the object's position. 
   * Modify the state vector, x_, and covariance, P_.
   * Also calculate the radar NIS.
   */
 //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 3;
  
  //matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  
  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  
  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  
  Zsig.fill(0.0);
  z_pred.fill(0.0);
  S.fill(0.0);
  double rho = 0;
  double phi = 0;
  double rho_d = 0;
  
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    //transform sigma points into measurement space
    VectorXd state_vec = Xsig_pred_.col(i);
    double px = state_vec(0);
    double py = state_vec(1);
    double v = state_vec(2);
    double yaw = state_vec(3);
    double yaw_d = state_vec(4);
    
    rho = sqrt(px*px+py*py);
    phi = atan2(py,px);
    rho_d = (px*cos(yaw)*v+py*sin(yaw)*v) / rho;
    
    Zsig.col(i) << rho,
                   phi,
                   rho_d;
    
    //mean predicted measurement
    z_pred += weights_(i) * Zsig.col(i);
  }
  
  //measurement covariance matrix S
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    VectorXd z_diff = Zsig.col(i) - z_pred;
    while (z_diff(1) > M_PI) 
      z_diff(1) -= 2. * M_PI;
    while (z_diff(1) < - M_PI) 
      z_diff(1) += 2. * M_PI;
    
    S += weights_(i) * z_diff * z_diff.transpose();
  }
  
  // Add R to S
  MatrixXd R = MatrixXd(3,3);
  R << std_radr_*std_radr_, 0, 0,
       0, std_radphi_*std_radphi_, 0,
       0, 0, std_radrd_*std_radrd_;
  S += R;
  
  //example vector for incoming radar measurement
  VectorXd z = VectorXd(n_z);
  
  double meas_rho = meas_package.raw_measurements_(0);
  double meas_phi = meas_package.raw_measurements_(1);
  double meas_rhod = meas_package.raw_measurements_(2);
  
  z << meas_rho,
       meas_phi,
       meas_rhod;
  
  //matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);
  
  //cross correlation matrix
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //normalize angles
    while (x_diff(3) > M_PI) 
      x_diff(3) -= 2. * M_PI;
    while (x_diff(3) < -M_PI) 
      x_diff(3) += 2. * M_PI;
    
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //normalize angles
    while (z_diff(1) > M_PI) 
      z_diff(1) -= 2. * M_PI;
    while (z_diff(1) < -M_PI) 
      z_diff(1) += 2. * M_PI;
  
    Tc += weights_(i) * x_diff * z_diff.transpose();
    
  }
  
  // difference
  VectorXd z_diff = z - z_pred;
  
  //normalize angles
  while (z_diff(1) > M_PI) 
    z_diff(1) -= 2. * M_PI;
  while (z_diff(1) < -M_PI) 
    z_diff(1) += 2. * M_PI;
  
  
  //calculate NIS
  NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;
  
  //calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();
  
  //update state mean and covariance matrix
  x_ += K*z_diff;
  P_ -= K*S*K.transpose();
}