#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {

    is_initialized_ = false;
    previous_timestamp_ = 0;

    // initializing matrices
    ekf_.x_ = VectorXd(4);

    ekf_.P_ = MatrixXd(4, 4);
    ekf_.P_ << 1, 0, 0,    0,
               0, 1, 0,    0,
               0, 0, 1000, 0,
               0, 0, 0,    1000;

    ekf_.F_ = MatrixXd(4, 4);

    H_laser_ = MatrixXd(2, 4);
    H_laser_ << 1, 0, 0, 0,
                0, 1, 0, 0;

    Hj_ = MatrixXd(3, 4);

    R_laser_ = MatrixXd(2, 2);
    R_laser_ << 0.0225, 0,
                0,      0.0225;
    R_radar_ = MatrixXd(3, 3);
    R_radar_ << 0.09, 0,      0,
                0,    0.0009, 0,
                0,    0,      0.09;

    ekf_.Q_ = MatrixXd(4, 4);

}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
    /*****************************************************************************
     *  Initialization
     ****************************************************************************/
    while (!is_initialized_){
        /**
        TODO:
          * Initialize the state ekf_.x_ with the first measurement.
          * Create the covariance matrix.
          * Remember: you'll need to convert radar from polar to cartesian coordinates.
        */
        // first measurement
        ekf_.x_ << 0,0,0,0;
        cout << "searching for non-zero measurement... " << endl;

        if ((measurement_pack.sensor_type_ == MeasurementPackage::RADAR) &&
            (measurement_pack.raw_measurements_[0]* measurement_pack.raw_measurements_[1] != 0)) {
            previous_timestamp_ = measurement_pack.timestamp_;

            double rho = measurement_pack.raw_measurements_[0];
            double phi = measurement_pack.raw_measurements_[1];

            ekf_.x_ << cos(phi) * rho, sin(phi) * rho, 1, 1;
            cout << "initial x (1st non-zero radar measurement):" << endl;
            cout << ekf_.x_ << endl;
            is_initialized_ = true;

        }
        else if ((measurement_pack.sensor_type_ == MeasurementPackage::LASER) &&
                (measurement_pack.raw_measurements_[0]* measurement_pack.raw_measurements_[1] != 0)){
            previous_timestamp_ = measurement_pack.timestamp_;
            cout << "initial time: " << previous_timestamp_ << endl;

            ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 1, 1;
            cout << "initial x (1st non-zero laser measurement):" << endl;
            cout << ekf_.x_ << endl;
            is_initialized_ = true;
        }
        cout << "-----------------------------" << endl;
        return;
    }

    /*****************************************************************************
     *  Prediction
     ****************************************************************************/

    /**
     TODO:
       * Update the state transition matrix F according to the new elapsed time.
        - Time is measured in seconds.
       * Update the process noise covariance matrix.
     */

    double dt = (double)(measurement_pack.timestamp_-previous_timestamp_)/1000000.0;
    cout << "delta t: " << dt << endl;
    previous_timestamp_ = measurement_pack.timestamp_;

    // update F matrix based on dt
    ekf_.F_ << 1, 0, dt, 0,
               0, 1, 0,  dt,
               0, 0, 1,  0,
               0, 0, 0,  1;

    // update Q matrix based on dt
    double dt_2 = dt * dt;
    double dt_3 = dt_2 * dt;
    double dt_4 = dt_3 * dt;

    int noise_ax = 9;
    int noise_ay = 9;

    ekf_.Q_ <<  dt_4/4*noise_ax,    0,               dt_3/2*noise_ax, 0,
                0,                  dt_4/4*noise_ay, 0,               dt_3/2*noise_ay,
                dt_3/2*noise_ax,    0,               dt_2*noise_ax,   0,
                0,                  dt_3/2*noise_ay, 0,               dt_2*noise_ay;

    ekf_.Predict();
    cout << "predicted x:" << endl;
    cout << ekf_.x_ << endl;

    /*****************************************************************************
     *  Update
     ****************************************************************************/

    /**
     TODO:
       * Use the sensor type to perform the update step.
       * Update the state and covariance matrices.
     */

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
        // Radar updates
        Hj_ = tools.CalculateJacobian(ekf_.x_);
        ekf_.H_ = Hj_;
        ekf_.R_ = R_radar_;
        cout << "radar raw measurements" << endl;
        cout << measurement_pack.raw_measurements_ << endl;
        ekf_.UpdateEKF(measurement_pack.raw_measurements_);
        cout << "radar updated " << endl;
    } else {
        // Laser updates
        ekf_.H_ = H_laser_;
        ekf_.R_ = R_laser_;
        cout << "lidar raw measurements" << endl;
        cout << measurement_pack.raw_measurements_ << endl;
        ekf_.Update(measurement_pack.raw_measurements_);
        cout << "lidar updated" << endl;
    }

    // print the output
    cout << "x_: " << endl;
    cout << ekf_.x_ << endl;
    cout << "P_: " << endl;
    cout << ekf_.P_ << endl;
    cout << "----------------------------------------------" << endl;
}
