#include "JumpNLPWrapper.h"

#ifdef IPOPT_OPTION
#include "JumpNLP.hpp"

#include "IpIpoptApplication.hpp"
#include "IpSolveStatistics.hpp"
#include <cstdlib>
#include <pthread.h>
#include <unistd.h>

#ifdef __cplusplus
extern "C" {
#endif

static SmartPtr<JumpNLP<double> > JumpMPCObj;
static SmartPtr<IpoptApplication> app;
static ApplicationReturnStatus status;
bool FIRST_RUN = false;
bool RUN_RPC = false;
bool INITIALIZED = false;

// For mutex stuff
static pthread_mutex_t solve_set_mutex = PTHREAD_MUTEX_INITIALIZER;


void JumpMPC_init()
{
    if (INITIALIZED == false)  {
        JumpMPCObj = new JumpNLP<double>();
        app = IpoptApplicationFactory();
        INITIALIZED = true;
        printf("[JumpMPC] Creating thread...\n");
        JumpMPC_CreateSolveThread();
    }
}


void JumpMPC_Initialize() {
    status = app->Initialize();
}


void JumpMPC_SetFirstRun() {
    pthread_mutex_lock(&solve_set_mutex);
    FIRST_RUN = true;
    pthread_mutex_unlock(&solve_set_mutex);
}


/*=========================== JumpMPC Settings ============================*/
void JumpMPC_SetCurrentState(double * current_state_in) {
    JumpMPCObj->SetCurrentState(current_state_in);
}

void JumpMPC_SetFootLocations(double * p_foot_0_in) {
    JumpMPCObj->SetFootLocations(p_foot_0_in);
}


/*=========================== JumpMPC Settings ============================*/
/**
 * Creates a new state error weighting matrix.
 *
 *  @param newQ
 *      the new matrix values for Q in an array.
 */
void JumpMPC_SetStateWeights(double * Q_in)
{
    JumpMPCObj->SetStateWeights(Q_in);
}

/**
 * Creates a new input regularization weighting matrix.
 *
 *  @param newR
 *      the new matrix values for R in an array.
 */
void JumpMPC_SetInputRegularization(double * R_in)
{
    JumpMPCObj->SetInputRegularization(R_in);
}

void JumpMPC_SetRobotParameters(double mass_in, double Inertia_in) {
    JumpMPCObj->SetRobotParameters(mass_in, Inertia_in);
}

void JumpMPC_SetEnvironmentParameters(double gravity_in, double mu_in, double * z_g_in) {
    JumpMPCObj->SetEnvironmentParameters(gravity_in, mu_in, z_g_in);
}


void JumpMPC_SetReady() {
    RUN_RPC = true;
}


/*============================ IPOPT Solver =============================*/

void JumpMPC_SetIPOPTOptions() {
    // Printing Options
    app->Options()->SetIntegerValue("file_print_level", 0);
    app->Options()->SetIntegerValue("print_level", 0);
    app->Options()->SetIntegerValue("print_frequency_iter", 20);
    app->Options()->SetStringValue("print_timing_statistics", "yes");


    app->Options()->SetNumericValue("bound_frac", 1e-6);
    app->Options()->SetNumericValue("bound_push", 1e-6);

    // Stopping conditions and tolerances
    app->Options()->SetIntegerValue("max_iter", 50);
    app->Options()->SetNumericValue("max_cpu_time",0.05);
    app->Options()->SetNumericValue("acceptable_tol", 0.00001);
    app->Options()->SetNumericValue("tol", 0.0001);

    app->Options()->SetStringValue("linear_solver", "ma27");
    //app->Options()->SetStringValue("mehrotra_algorithm","yes");
    app->Options()->SetStringValue("fixed_variable_treatment", "relax_bounds");
    //app->Options()->SetStringValue("honor_original_bounds","no")
    //app->Options()->SetStringValue("derivative_test", "first-order");
    //app->Options()->SetStringValue("derivative_test", "second-order");
    //app->Options()->SetNumericValue("derivative_test_tol",0.005);
    // Might actually be faster with hessian-approximation...
    //app->Options()->SetStringValue("hessian_approximation", "limited-memory");


    //app->Options()->SetStringValue("nlp_scaling_method","none");
    printf("[JumpMPC] Options Set\n");

}



static void * JumpMPC_SolvePrediction(void * arg) {
  (void)arg;
    // Microsecond delay
    unsigned int microseconds = 10000;

    // Solve for the rest of the time
    while (true) {
        if (FIRST_RUN == true) {
            RUN_RPC = true;
            if (RUN_RPC) {
              //printf("Run\n");
                // Run the optimization
                status = app->OptimizeTNLP(JumpMPCObj);
                RUN_RPC = false;
            }
        } else {
            // Do nothing
        }

        // Quick pause for stability
        usleep(microseconds);
    }
}

void JumpMPC_CreateSolveThread() {
    pthread_t prmpc_thread;
    int prmpc_thread_rc = pthread_create(&prmpc_thread, NULL, JumpMPC_SolvePrediction, NULL);
    if (prmpc_thread_rc) {
        printf("[JumpMPC] Failed to create JumpMPC solve thread!\n");
    } else {
        printf("[JumpMPC] Successfully created JumpMPC solve thread\n!");
    }
    //JumpMPC_SetMutex();
    //pthread_join(prmpc_thread, NULL);
}

double JumpMPC_GetSolution(int index) {
    return JumpMPCObj->GetSolution(index);
}


int JumpMPC_GetSolved() {
    return JumpMPCObj->GetSolved();
}

#ifdef __cplusplus
}
#endif

#endif // IPOPT_OPTION
