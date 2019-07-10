#ifndef _solver_mpc
#define _solver_mpc


#include <eigen3/Eigen/Dense>
#include "common_types.h"
#include "convexMPC_interface.h"
#include <iostream>
#include <stdio.h>

using Eigen::Matrix;
using Eigen::Quaternionf;
using Eigen::Quaterniond;

template <class T>
void print_array(T* array, u16 rows, u16 cols)
{
    for(u16 r = 0; r < rows; r++)
    {
        for(u16 c = 0; c < cols; c++)
            std::cout<<(fpt)array[c+r*cols]<<" ";
        printf("\n");
    }
}

template <class T>
void print_named_array(const char* name, T* array, u16 rows, u16 cols)
{
    printf("%s:\n",name);
    print_array(array,rows,cols);
}

//print named variable
template <class T>
void pnv(const char* name, T v)
{
    printf("%s: ",name);
    std::cout<<v<<std::endl;
}

template <class T>
T t_min(T a, T b)
{
    if(a<b) return a;
    return b;
}

template <class T>
T sq(T a)
{
    return a*a;
}


void solve_mpc(update_data_t* update, problem_setup* setup);

void quat_to_rpy(Quaternionf q, Matrix<fpt,3,1>& rpy);
void ct_ss_mats(Matrix<fpt,3,3> I_world, fpt m, Matrix<fpt,3,4> r_feet, Matrix<fpt,3,3> R_yaw, Matrix<fpt,13,13>& A, Matrix<fpt,13,12>& B);
void resize_qp_mats(s16 horizon);
void c2qp(Matrix<fpt,13,13> Ac, Matrix<fpt,13,12> Bc,fpt dt,s16 horizon);
mfp* get_q_soln();
#endif
