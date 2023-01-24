#include "TCam.h"
#include <mutex>
#include <thread>

static std::mutex dataMutex;
void TCam::init(){
p.start();
}

void TCam::get_pose(){

    dataMutex.lock();
// Wait for the next set of frames from the camera
    auto frames = p.wait_for_frames();
    // Get a frame from the pose stream
    auto f = frames.first_or_default(RS2_STREAM_POSE);
    auto pose_data = f.as<rs2::pose_frame>().get_pose_data();
    pose.x = pose_data.translation.x;
    pose.y = pose_data.translation.y;
    pose.z = pose_data.translation.z;

    pose.rot_w = pose_data.rotation.w;
    pose.rot_x = pose_data.rotation.x;
    pose.rot_y= pose_data.rotation.y;
    pose.rot_z= pose_data.rotation.z;

    pose.vel_x = pose_data.velocity.x;
    pose.vel_y = pose_data.velocity.y;
    pose.vel_z = pose_data.velocity.z;

    pose.acc_x = pose_data.acceleration.x;
    pose.acc_y = pose_data.acceleration.y;
    pose.acc_z= pose_data.acceleration.z;

    pose.ang_vel_x = pose_data.angular_velocity.x;
    pose.ang_vel_y = pose_data.angular_velocity.y;
    pose.ang_vel_z = pose_data.angular_velocity.z;
    


    dataMutex.unlock();

}
