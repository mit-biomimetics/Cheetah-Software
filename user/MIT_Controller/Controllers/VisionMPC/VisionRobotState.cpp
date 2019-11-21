#include "VisionRobotState.h"
#include <eigen3/Eigen/Dense>
#include <iostream>
#include <math.h>

using std::cout;
using std::endl;

void VisionRobotState::set(flt* p_, flt* v_, flt* q_, flt* w_, flt* r_,flt yaw_)
{
    for(u8 i = 0; i < 3; i++)
    {
        this->p(i) = p_[i];
        this->v(i) = v_[i];
        this->w(i) = w_[i];
    }
    this->q.w() = q_[0];
    this->q.x() = q_[1];
    this->q.y() = q_[2];
    this->q.z() = q_[3];
    this->yaw = yaw_;

    //for(u8 i = 0; i < 12; i++)
    //    this->r_feet(i) = r[i];
    for(u8 rs = 0; rs < 3; rs++)
        for(u8 c = 0; c < 4; c++)
            this->r_feet(rs,c) = r_[rs*4 + c];

    R = this->q.toRotationMatrix();
    fpt yc = cos(yaw_);
    fpt ys = sin(yaw_);

    R_yaw <<  yc,  -ys,   0,
             ys,  yc,   0,
               0,   0,   1;

    Matrix<fpt,3,1> Id;
    Id << .07f, 0.26f, 0.242f;
    //Id << .07f, 0.71f, 0.712f;
    I_body.diagonal() = Id;

    //TODO: Consider normalizing quaternion??
}

void VisionRobotState::print()
{
   cout<<"Robot State:"<<endl<<"Position\n"<<p.transpose()
       <<"\nVelocity\n"<<v.transpose()<<"\nAngular Veloctiy\n"
       <<w.transpose()<<"\nRotation\n"<<R<<"\nYaw Rotation\n"
       <<R_yaw<<"\nFoot Locations\n"<<r_feet<<"\nInertia\n"<<I_body<<endl;
}



