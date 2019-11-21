#include "VisionMPC_interface.h"
#include "../convexMPC/common_types.h"
#include "VisionSolverMPC.h"
#include <eigen3/Eigen/Dense>
#include <pthread.h>
#include <stdio.h>
#include <string.h>

#define V_NUM_LEGS 4

vision_mpc_problem_setup v_problem_config;
vision_mpc_update_data_t v_update;

pthread_mutex_t vision_problem_cfg_mt;
pthread_mutex_t vision_mpc_update_mt;
pthread_t vision_mpc_solve_thread;

u8 v_first_run = 1;

void vision_initialize_mpc()
{
  //printf("Initializing MPC!\n");
  if(pthread_mutex_init(&vision_problem_cfg_mt,NULL)!=0)
    printf("[MPC ERROR] Failed to initialize problem configuration mutex.\n");

  if(pthread_mutex_init(&vision_mpc_update_mt,NULL)!=0)
    printf("[MPC ERROR] Failed to initialize update data mutex.\n");
}

void vision_setup_problem(double dt, int horizon, double mu, double f_max)
{
  if(v_first_run) {
    v_first_run = false;
    vision_initialize_mpc();
  }

  //pthread_mutex_lock(&problem_cfg_mt);

  v_problem_config.horizon = horizon;
  v_problem_config.f_max = f_max;
  v_problem_config.mu = mu;
  v_problem_config.dt = dt;

  //pthread_mutex_unlock(&problem_cfg_mt);
  vision_resize_qp_mats(horizon);
}

//inline to motivate gcc to unroll the loop in here.
inline void vision_mpf_to_flt(flt* dst, mfp* src, s32 n_items) {
  for(s32 i = 0; i < n_items; i++)
    *dst++ = *src++;
}

inline void vision_mint_to_u8(u8* dst, mint* src, s32 n_items) {
  for(s32 i = 0; i < n_items; i++)
    *dst++ = *src++;
}

int vision_has_solved = 0;

//safely copies problem data and starts the solver
void vision_update_problem_data(double* p, double* v, double* q, double* w, double* r, double yaw, double* weights, double* state_trajectory, double alpha, int* gait)
{
  vision_mpf_to_flt(v_update.p,p,3);
  vision_mpf_to_flt(v_update.v,v,3);
  vision_mpf_to_flt(v_update.q,q,4);
  vision_mpf_to_flt(v_update.w,w,3);
  vision_mpf_to_flt(v_update.r,r,12);
  v_update.yaw = yaw;
  vision_mpf_to_flt(v_update.weights,weights,12);
  //this is safe, the solver isn't running, and v_update_problem_data and setup_problem
  //are called from the same thread
  vision_mpf_to_flt(v_update.traj,state_trajectory,12*v_problem_config.horizon);
  v_update.alpha = alpha;
  vision_mint_to_u8(v_update.gait,gait,4*v_problem_config.horizon);

  vision_solve_mpc(&v_update, &v_problem_config);
  vision_has_solved = 1;
}

void vision_update_problem_data_floats(float* p, float* v, float* q, float* w,
                                float* r, float yaw, float* weights,
                                float* state_trajectory, float alpha, int* gait)
{
  v_update.alpha = alpha;
  v_update.yaw = yaw;
  vision_mint_to_u8(v_update.gait,gait,4*v_problem_config.horizon);
  memcpy((void*)v_update.p,(void*)p,sizeof(float)*3);
  memcpy((void*)v_update.v,(void*)v,sizeof(float)*3);
  memcpy((void*)v_update.q,(void*)q,sizeof(float)*4);
  memcpy((void*)v_update.w,(void*)w,sizeof(float)*3);
  memcpy((void*)v_update.r,(void*)r,sizeof(float)*12);
  memcpy((void*)v_update.weights,(void*)weights,sizeof(float)*12);
  memcpy((void*)v_update.traj,(void*)state_trajectory, sizeof(float) * 12 * v_problem_config.horizon);
  vision_solve_mpc(&v_update, &v_problem_config);
  vision_has_solved = 1;

}

double vision_get_solution(int index)
{
  if(!vision_has_solved) return 0.f;
  mfp* qs = vision_get_q_soln();
  return qs[index];
}
