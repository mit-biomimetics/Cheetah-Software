#ifndef vision_solver_mpc
#define vision_solver_mpc


#include <eigen3/Eigen/Dense>
#include "../convexMPC/common_types.h"
#include "VisionMPC_interface.h"
#include <iostream>
#include <stdio.h>

using Eigen::Matrix;
using Eigen::Quaternionf;
using Eigen::Quaterniond;


template <class T>
T vision_t_min(T a, T b) {
    if(a<b) return a;
    return b;
}

template <class T>
T vision_sq(T a) {
    return a*a;
}


void vision_solve_mpc(vision_mpc_update_data_t* update, vision_mpc_problem_setup* setup);
void vision_quat_to_rpy(Quaternionf q, Matrix<fpt,3,1>& rpy);
void vision_ct_ss_mats(Matrix<fpt,3,3> I_world, fpt m, Matrix<fpt,3,4> r_feet, Matrix<fpt,3,3> R_yaw, Matrix<fpt,13,13>& A, Matrix<fpt,13,12>& B);
void vision_resize_qp_mats(s16 horizon);
void vision_c2qp(Matrix<fpt,13,13> Ac, Matrix<fpt,13,12> Bc,fpt dt,s16 horizon);
mfp* vision_get_q_soln();
#endif
