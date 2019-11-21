#ifndef _Vision_mpc_interface
#define _Vision_mpc_interface

#define V_MAX_GAIT_SEGMENTS 36

#ifdef __cplusplus
#define EXTERNC extern "C"
#else
#define EXTERNC
#endif

struct vision_mpc_problem_setup
{
  float dt;
  float mu;
  float f_max;
  int horizon;
};

struct vision_mpc_update_data_t
{
  float p[3];
  float v[3];
  float q[4];
  float w[3];
  float r[12];
  float yaw;
  float weights[12];
  float traj[12*V_MAX_GAIT_SEGMENTS];
  float alpha;
  unsigned char gait[V_MAX_GAIT_SEGMENTS];
  unsigned char hack_pad[1000];
  int max_iterations;
  float x_drag;
};

EXTERNC void vision_setup_problem(double dt, int horizon, double mu, double f_max);
EXTERNC void vision_update_problem_data(double* p, double* v, double* q, double* w, double* r, double yaw, double* weights, double* state_trajectory, double alpha, int* gait);
EXTERNC double vision_get_solution(int index);
EXTERNC void vision_update_problem_data_floats(float* p, float* v, float* q, float* w,
                                        float* r, float yaw, float* weights,
                                        float* state_trajectory, float alpha, int* gait);

void update_x_drag(float x_drag);
#endif
