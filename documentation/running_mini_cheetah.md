# Running Mini Cheetah

When viewed from the top, mini cheetah's two switches are in front.  The center switch is motor power and right switch is computer power.  Toggling the switch toward the back of the robot turns it on.

First turn on computer power.  It takes around 2 minutes for the computer to boot.
SSH into the computer with `ssh user@10.0.0.34`.  The password is the usual lab password.  You must set your computer's IP address to something like `10.0.0.xxx` for this to work.

Next, set up your computer for LCM.  In the scripts folder, run `./config_network_lcm.sh -I ` followed by the name of your network adapter you're connecting with.  If you want, you can edit the script to add a shortcut for your specific network adapter.

To build code which can run on the mini cheetah:

- `mkdir mc-build && cd mc-build`
- `cmake -DMINI_CHEETAH_BUILD=TRUE ..`
- `make -j`

Then, to copy the code to the robot, run `../scripts/send_to_mini_cheetah.sh user/yourController/your_controller`.  You will be prompted for the mini cheetah's password.

If you would like to open LCM spy, you can do this by running `../script/launch_lcm_spy.sh`.  If you receive an error about Java, try running `../scripts/launch_lcm_spy_jdk_fix.sh`, which has modified launch arguments to support newer versions of the JVM.


On the mini cheetah, you will find a folder in the home directory name `robot-software` and the date.  Inside this folder is a copy of the configure lcm script.  Run `./configure_network_lcm.sh mc-top` if you are connected with the ethernet cable on top of the robot.

To run the robot code, enter the build folder and run `(sudo) nohup mc_run.sh`.  The mini cheetah program will wait until the interface has been open and started in MiniCheetah, Robot mode.

The robot controller currently works with:

- State estimate (orientation only), access with `_stateEstimate`
- Leg control (torque feedforward, force feedforward, PD control cartesian, PD control joint), access through LegCommand.  Note that the leg command is zeroed for you on each iteration.
- Leg data (joint position/velocity, cartesian position/velocity in hip frame, jacobian in hip frame)
- Gamepad data
- Main Visualization (this is just a single cheetah model, the rest is still being implemented)
- Control parameters (set in the simulator gui)
- Configuration files (using the `THIS_COM`, or relative paths to the build folder)


It does not work with
- Full state estimator (position, velocity)
- Full visualization data
- Cheater state when running on the robot

Current LCM streams
- Raw spi data
- Raw vectornav data
- gamepad
- main visualization

Currently missing LCM streams
- Generic leg data
- State estimate
- Probably others I am forgetting

