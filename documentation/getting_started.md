# Getting Started
This page explains how to download and install the software, run the MIT Controller, and how to create a controller of your own.

# Install dependencies

Packages:
```
sudo apt install mesa-common-dev freeglut3-dev coinor-libipopt-dev libblas-dev liblapack-dev gfortran liblapack-dev coinor-libipopt-dev cmake gcc build-essential libglib2.0-dev
```

Others:
- LCM 1.3.1 (it says Java 6, but you can use newer) (https://lcm-proj.github.io/)
- Qt 5.10.0 or newer (requires the gamepad library) (https://www.qt.io/download-qt-installer)
- Eigen (http://eigen.tuxfamily.org/)

NOTE: on Ubuntu 18.10 or 19.04, you may instead install Qt with
```
sudo apt install libqt5 libqt5gamepad5
```


# Download and build code

```
git clone https://github.com/mit-biomimetics/Cheetah-Software.git
cd Cheetah-Software
cd scripts # for now, you must actually go into this folder
./make_types.sh # you may see an error like `rm: cannot remove...` but this is okay
cd ..
mkdir build
cd build
cmake .. # there are still some warnings here
make -j
```


# Test
The code in the `common` folder has some tests.  From the build folder, these can be run with `common/test-common`.   There are two tests which commonly fail:

- OSQP - the solver itself is nondeterministic, so just run the test again and it should pass
- CASADI - the solver loads a library at runtime and sometimes has issues finding it.  It's fine if this fails as we aren't using this solver yet.

# Joystick
We use the Logitech F310 controller.  There's a switch in the back, which should be in the "X" position.  The controller needs to reconnected if you change the switch position.  Also, the LED on the front near the mode button should be off.
(https://www.amazon.com/Logitech-940-000110-Gamepad-F310/dp/B003VAHYQY)


# Simulation Example
The simulator default settings can be configured with `config/simulator-defaults.yaml` and `config/default-terrain.yaml`.  The default settings should be good for most uses, and the `default-terrain` file has commented out examples of how to add a mesh, box, and stairs.  The friction of the floor can be set from within the terrain file.

To launch the simulator, first plug in your joystick, then run `sim/sim` from within the `build` folder.    Select "Mini Cheetah" and "Simulator", then click "Start".  Somewhere in the output, you should see

```
[GameController] Found 1 joystick
```
The left panel allows you to change simulator settings.  The most useful setting is `simulation_speed`.  Defaults are loaded `simulator-defaults.yaml` when the simulator opens.  The settings file you load must have exactly the same set of parameters as is defined in the code.  You must recompile `sim` if you add or remove a parameter. You can save and load the parameters at any time, but note that the `use_spring_damper` setting will not take effect unless you restart the simulator.

The center panel allows you to change robot settings which are not specific to the controller.  Defaults are loaded `mini-cheetah-defatuls.yaml` when Mini Cheetah is selected and you click "Start".  If you stop and start the simulator, the file will be reloaded.  The settings file you load must have exactly the same set of parameters as is defined in the code.  You must recompile `sim` if you add or remove a parameter.  Currently most of these parameters do nothing and many will be removed. The most useful setting is `cheater_mode`, which sends the robot code the current position/orientation/velocity of the robot, and `controller_dt`, which changes how frequently the control code runs. (The `controller_dt` setting still needs to be tested).

The right panel allows you to change parameters which are specific to your controller, called "User Parameters".  If your control code does not use user parameters, it will be ignored.  If your code does use user parameters, then you must load user configuration parameters that match the parameters your code is expecting, or your controller will not start.


To start the robot control code, run `user/MIT_Controller/mit_ctrl m s`.  The `m` argument is for mini-cheetah, and the `s` indicates it should connect to the simulator. This uses shared memory to communicate with the simulator. The simulator should start running, and the robot should move to a ready position.  In the center column of the simulator window, set control mode to 10.  Once the robot has stopped moving, set control mode 1.  Then, set control mode to 4, and the robot will start trotting.

You can use the joysticks to drive the robot around.  You will see two robots - the gray one is the actual robot position from the simulation, and the red one is the estimate of the robot's position from our state estimator.  Turning on "cheater_mode" will make the estimated position equal to the actual position.  To adjust the simulation view, you can click and drag on the screen and scroll. Press and hold `t` to make the simulation run as fast as possible.  Press the spacebar to turn on free camera mode.  You can use the w,a,s,d,r,f keys to move the camera around, and click and drag to adjust the orientation.


# LCM
We use LCM (https://lcm-proj.github.io/) to connect the control interface to the actual mini cheetah hardware, and also as a debugging tool when running the simulator.  The `make_types.sh` script runs an LCM tool to generate C++ header files for the LCM data types.  When the simulator is running, you can run `scripts/launch_lcm_spy.sh` to open the LCM spy utility, which shows detailed information from the simulator and controller.  You can click on data streams to plot them, which is nice for debugging.  There is also a tool called `lcm-logger` which can save LCM data to a file.


# Writing a Robot Controller
To add your own robot controller, you should add a folder under `Cheetah-Software/user`, and add the folder to the `CMakeLists.txt` in `user`.  The `JPos_Controller` is an example of a very simple controller.  The `JPosUserParameters.h` file has an example of declaring two user parameters which can be adjusted from the simulator interface, but using user parameters is optional.  The `JPos_Controller.hpp` and `JPos_Controller.cpp` files are the actual controller, which should extend `RobotController`.  Notice that in the `JPos_Controller.hpp` file, the `getUserControlParameters` method retuns a pointer to the user parameters.  If you do not use user parameters, your `getUserControlParameters` should return `nullptr`.  Finally, your `main` function must be like the example main function in `main.cpp`.

The `runController` method of your controller will be called automatically, at 1 kHz.  Here, you have access to the following:

- `_quadruped` : contains constant parameters about the robot (link lengths, gear ratios, inertias...).  The `getHipLocation` function returns the location of the "hip" in the body coordinate system.  The x-axis points forward, y-axis to the left, and z-axis up.  The legs are ordered like this

```
FRONT
1 0  RIGHT
3 2
BACK
```

- `_model` : a dynamics model of the robot.  This can be used to compute forward kinematics, Jacobians, etc...
- `_legController`: Interface to the robot's legs. This data is syncronized with the hardware at around 700 Hz. There are multiple ways to control the legs, and the result from all the controllers are added together.
    - `commands[leg_id].tauFeedForward` : Leg torque (Nm, at the joint).  Order is ab/ad, hip, knee.
    - `commands[leg_id].forceFeedForward` : Force to apply at foot (N), in hip frame. (Same orientation as body frame, origin is the hip)
    - `commands[leg_id].qDes` : Desired joint position for joint PD controller (radians). Order is ab/ad, hip, knee.  `(0,0,0)` is leg pointing straight down.
    - `commands[leg_id].qdDes` : Desired joint velocity for joint PD controller (rad/sec).
    - `commands[leg_id].pDes, vDes` : Desired foot position/velocity for cartesian PD controller (meters, hip frame)
    - `commands[leg_id].kpCartesian, kdCartesian, kpJoint, kdJoint` : Gains for PD controllers (3x3 matrix).  Use the diagonal entries only.
    - `datas[leg_id].q` : Leg joint encoder (radians).  Order is ab/ad, hip, knee.  `(0,0,0)` is leg pointing straight down.
    - `datas[leg_id].qd` : Leg joint velocity (radians/sec).  Same order as `q`.
    - `datas[leg_id].p`  : Foot cartesian position, in hip frame. (Same orientation as body frame, origin is the hip)
    - `datas[leg_id].v`  : Foot cartesian velocity, in hip frame. 
    - `datas[leg_id].tau` : Estimate of motor torque from combination of all controllers
The joint PD control actually runs at 40 kHz on the motor controllers.
- `_stateEstimate, _stateEstimatorContainer` The result and interface for the provided state estimator.  If you provide the contact state of the robot (which feet are touching the ground), it will determine the robot's position/velocity in the world.
- `_driverCommand` : inputs from the game pad.
- `_controlParameters` : values from the center robot control parameters panel
- `_visualizationData` : interface to add debugging visualizations to the simulator window
- `_robotType` : If you are the mini Cheetah or Cheetah 3 robot.


If you would like to see more of how this works, look at the `robot` folder.  The `RobotRunner` class actually runs the control code, and connects it with either the `HardwareBridge` or `SimulationBridge`.  The code in the `rt` folder actually interacts with the hardware.

# User Parameters
User Parameters are settings which are specific to the controller you are running.  The list of user parameters and their values are defined in a yaml file.  On startup, the file `config/default-user.yaml` is loaded into the simulator. Currently you must manually make sure that your currently loaded list of user parameters (on the right of the simulator window) matches the user parameters you define in your controller.  If these do not match, the controller will see the mismatch and print an error like 

`parameter cmpc_gait wasn't found in parameter collection user-parameters`

In this case, you should edit `config/default-user.yaml` to have the same parameters as your controller's user parameters. After changing `default-user.yaml`, you must either restart the simulator, or use the "Load" button to reload the parameters.

If you do not want to use user parameters, you can simply put

`__collection-name__: user-parameters`

as the only thing in `config/default-user.yaml`.  Then, in your `RobotController` class, you should provide the method

```
  // indicate that there are no user control parameters for this controller
  virtual ControlParameters* getUserControlParameters() {
    return nullptr;
  }
```


# Errors
If the controller software encounters an error where it cannot determine what to do safely (such as incompatible control parameters), it will stop itself.  Currently it will crash with `terminate called after throwing an instance of ....` and should print the exception and an error message. 

If you see the controller crash with `CRASH: Caught 11` (segfault)  or `terminate called without active exception`, it is because the code has actually crashed.  If this crash is not caused by your controller code, this is a bug, and you should create a Github issue.

If the controller stops responding for more than 1 second (it has crashed, or is stuck in a loop), the simulator will go into an error state.  It will print `[ERROR] Timed out waiting for message from robot!  Did it crash?`. 

If the simulator itself throws an exception or crashes, this is a bug, and you should create a Github issue.

# Recovering from Errors
Currently, the most reliable way to recover from errors is:

- Kill the controller code
- Click "Stop" in the simulator
- Click "Start" in the simulator
- Restart the controller code


If you ever get the simulator in a state where you can't click "Stop" to reset the simulation (or the simulator crashes), it is a bug and you should open a Github issue.



# Unfinished
- Verify that the controller update rate can be changed in simulation/robot hardware
- Using the State Estimator (and disabling it if you don't want it)
- Using the Dynamics
- Final update rate to legs.
- Are other people going to use the wireless RC controller?  If so, add it to the robot controller.
- Safety checks
- Visualization data - not everything is implemented and they don't work on the robot


# Running on the robot differences
Running on the robot is very similar to running the simulator.  You will still have the gamepad, user parameters, and robot parameters.  In the simulation window, you will see only the state estimate of the robot, and cheater mode will not work.  Currently debugging visualizations don't work.




