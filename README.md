# SFND_Unscented_Kalman_Filter

<img src="media/ukf_highway_tracked.gif" width="700" height="400" />

In this project, an Unscented Kalman Filter is implemented to estimate the state of multiple cars on a highway using noisy lidar and radar measurements. 
Kalman filters work on 3 steps as below:
1. Initialization - Initialize state and covariance matrices
2. Prediction - Predict the location after a defined time interval
3. Update - Compare the predicted location with actual measurement and combine them to give updated location.

Standard Kalman Filter can predict only linear motion. Both Extended Kalman Filter and Unscented Kalman filter handles non-linear motion.
However their implementation is different. In Extended Kalman Filter, Jacobian Matrix is used to linearlize non-linear functions. Unscented Kalman Filter takes representative points from a Gaussian distribution instead of linearlizing non-linear functions.
In this project, Constant Turn Rate and Velocity Magnitude model (CTRV) is used for the Unscented Kalman Filter.
UKF proces chain performs the following steps for prediction and update.

Prediction:
1. Generate Sigma Points
2. Predict Sigma Points
3. Predict Mean and Covariance

Update:
1. Predict Measurement
2. Update State

The main program can be built and ran by doing the following from the project top directory.

1. mkdir build
2. cd build
3. cmake ..
4. make
5. ./ukf_highway

`main.cpp` is using `highway.h` to create a straight 3 lane highway environment with 3 traffic cars and the main ego car at the center. 
The viewer scene is centered around the ego car and the coordinate system is relative to the ego car as well. The ego car is green while the 
other traffic cars are blue. The traffic cars will be accelerating and altering their steering to change lanes. Each of the traffic car's has
it's own UKF object generated for it, and will update each indidual one during every time step. 

The red spheres above cars represent the (x,y) lidar detection and the purple lines show the radar measurements with the velocity magnitude along the detected angle. The Z axis is not taken into account for tracking, so you are only tracking along the X/Y axis.

---

## Other Important Dependencies
* cmake >= 3.5
  * All OSes: [click here for installation instructions](https://cmake.org/install/)
* make >= 4.1 (Linux, Mac), 3.81 (Windows)
  * Linux: make is installed by default on most Linux distros
  * Mac: [install Xcode command line tools to get make](https://developer.apple.com/xcode/features/)
  * Windows: [Click here for installation instructions](http://gnuwin32.sourceforge.net/packages/make.htm)
* gcc/g++ >= 5.4
  * Linux: gcc / g++ is installed by default on most Linux distros
  * Mac: same deal as make - [install Xcode command line tools](https://developer.apple.com/xcode/features/)
  * Windows: recommend using [MinGW](http://www.mingw.org/)
 * PCL 1.2

## Basic Build Instructions

1. Clone this repo.
2. Make a build directory: `mkdir build && cd build`
3. Compile: `cmake .. && make`
4. Run it: `./ukf_highway`

## Editor Settings

* indent using spaces
* set tab width to 2 spaces (keeps the matrices in source code aligned)

## Code Style

[Google's C++ style guide](https://google.github.io/styleguide/cppguide.html) 
