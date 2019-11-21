#ifndef JUMP_MPC_WRAPPER_H 
#define JUMP_MPC_WRAPPER_H 

#ifdef IPOPT_OPTION

#ifdef __cplusplus
extern "C" {
#endif

void JumpMPC_init();
void JumpMPC_Initialize();

void JumpMPC_SetCurrentState(double * current_state_in);
void JumpMPC_SetFootLocations(double * p_foot_0_in);
void JumpMPC_SetFirstRun();

// JumpMPC Settings
void JumpMPC_SetNumPredictions(double NUM_PREDICTIONS_in);
void JumpMPC_SetStateWeights(double * Q_in);
void JumpMPC_SetInputRegularization(double * R_in);
void JumpMPC_SetHeuristicGains(double * K_ref_in);
void JumpMPC_SetRobotParameters(double mass_in, double Inertia_in);
void JumpMPC_SetEnvironmentParameters(double gravity_in, double mu_in, double * z_g_in);
void JumpMPC_SetReady();

// IPOPT Solver
void JumpMPC_SetIPOPTOptions();
void JumpMPC_CreateSolveThread();
//void JumpMPC_SolvePrediction();

// Return values
double JumpMPC_GetSolution(int index);
double JumpMPC_GetSolutionTime(int index);
int JumpMPC_GetSolved();

#ifdef __cplusplus
}
#endif

#endif // IPOPT_OPTION

#endif
