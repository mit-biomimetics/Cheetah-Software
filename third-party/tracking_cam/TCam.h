#include <librealsense2/rs.hpp> 

struct Pose{
    // Wait for the next set of frames from the camera
    
    float x,y,z; // linear translation
    float rot_w, rot_x, rot_y,rot_z; // quaternions
    float vel_x, vel_y, vel_z; // linear velocity
    float acc_x, acc_y, acc_z; // linear acceleration
    float ang_vel_x, ang_vel_y, ang_vel_z; // angular acceleration
};

class TCam{
public:
    void init();
    void get_pose();
    Pose pose;


private:
    rs2::pipeline p;
};