#ifndef Vision_RobotState
#define Vision_RobotState

#include <eigen3/Eigen/Dense>
#include "../convexMPC/common_types.h"

using Eigen::Matrix;
using Eigen::Quaternionf;

class VisionRobotState
{
    public:
        void set(flt* p, flt* v, flt* q, flt* w, flt* r, flt yaw);
        //void compute_rotations();
        void print();
        Matrix<fpt,3,1> p,v,w;
        Matrix<fpt,3,4> r_feet;
        Matrix<fpt,3,3> R;
        Matrix<fpt,3,3> R_yaw;
        Matrix<fpt,3,3> I_body;
        Quaternionf q;
        fpt yaw;
        fpt m = 9;
    //private:
};
#endif
