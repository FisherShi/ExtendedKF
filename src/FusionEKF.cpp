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
    ekf_.F_ << 1, 0, 1, 0,
               0, 1, 0, 1,
               0, 0, 1, 0,
               0, 0, 0, 1;

    H_laser_ = MatrixXd(2, 4);
    H_laser_ << 1, 0, 0, 0,
                0, 1, 0, 0;

    Hj_ = MatrixXd(3, 4);
    Hj_ << 1, 1, 1, 1,
           1, 1, 1, 1,
           1, 1, 1, 1;

    R_laser_ = MatrixXd(2, 2);
    R_laser_ << 1, 1,
                1, 1;
    R_radar_ = MatrixXd(3, 3);
    R_radar_ << 1, 1, 1,
                1, 1, 1,
                1, 1, 1;

    ekf_.Q_ = MatrixXd(4, 4);
    ekf_.Q_ << 1, 1, 1, 1,
               1, 1, 1, 1,
               1, 1, 1, 1,
               1, 1, 1, 1;

}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
    /*****************************************************************************
     *  Initialization
     ****************************************************************************/
    if (!is_initialized_) {
        /**
        TODO:
          * Initialize the state ekf_.x_ with the first measurement.
          * Create the covariance matrix.
          * Remember: you'll need to convert radar from polar to cartesian coordinates.
        */
        // first measurement
        cout << "EKF: " << endl;
        ekf_.x_ << 1, 1, 1, 1;

        if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
            previous_timestamp_ = measurement_pack.timestamp_;

            float rho = measurement_pack.raw_measurements_[0];
            float phi = measurement_pack.raw_measurements_[1];
            ekf_.x_ << cos(phi)*rho, sin(phi)*rho, 1, 1;
            ekf_.Init(ekf_.x_,ekf_.P_,ekf_.F_,Hj_,R_radar_,ekf_.Q_);
            cout << "initial x radar:" << endl;
            cout << ekf_.x_ << endl;

        }
        else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
            previous_timestamp_ = measurement_pack.timestamp_;

            ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 1, 1;
            ekf_.Init(ekf_.x_,ekf_.P_,ekf_.F_,H_laser_,R_laser_,ekf_.Q_);
            cout << "initial x laser:" << endl;
            cout << ekf_.x_ << endl;
        }

        // done initializing, no need to predict or update
        is_initialized_ = true;
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
    ekf_.F_ << 1, 0, dt, 0,
               0, 1, 0,  dt,
               0, 0, 1,  0,
               0, 0, 0,  1;
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
        cout << "start radar update" << endl;
        ekf_.H_ = Hj_;
        cout << "radar raw measurements" << endl;
        cout << measurement_pack.raw_measurements_ << endl;
        ekf_.UpdateEKF(measurement_pack.raw_measurements_);
        cout << "radar updated " << endl;
    } else {
        // Laser updates
        cout << "start lidar update" << endl;
        ekf_.H_ = H_laser_;
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
    cout << "cycle done" << endl;
}
